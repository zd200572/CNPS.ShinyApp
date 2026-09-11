# =====================================================================
# Reusable cycle results module (C / N / P / S).
# The module id IS the cycle key; results come from a shared
# reactiveValues filled by uploadServer.
# =====================================================================

cycleTabUI <- function(id) {
  ns <- NS(id)
  cfg <- CNPS_CONFIG[[id]]
  tagList(
    uiOutput(ns("summary_boxes")),
    navset_card_tab(
      id = ns("tabs"),
      nav_panel("基因丰度",
                card_body(
                  tags$p(class = "text-muted", "循环相关 KO 丰度与各过程（通路）聚合丰度。"),
                  h6("过程级丰度"), DT::DTOutput(ns("tbl_gene")),
                  downloadButton(ns("dl_gene"), "下载 CSV", class = "btn-sm btn-outline-primary mt-1"),
                  hr(), h6("KO 级丰度"), DT::DTOutput(ns("tbl_ko")),
                  downloadButton(ns("dl_ko"), "下载 CSV", class = "btn-sm btn-outline-primary mt-1")
                )),
      nav_panel("分组热图与差异检验",
                card_body(
                  tags$p(class = "text-muted",
                         "各组均值 log10 丰度热图；红点标记组间差异显著（p<0.05，ANOVA / Kruskal-Wallis）。"),
                  plotly::plotlyOutput(ns("plot_heat"), height = "480px"),
                  downloadButton(ns("dl_heat_png"), "下载 PNG", class = "btn-sm btn-outline-primary"),
                  hr(), h6("差异检验结果"), DT::DTOutput(ns("tbl_heat")),
                  downloadButton(ns("dl_heat_csv"), "下载 CSV", class = "btn-sm btn-outline-primary mt-1")
                )),
      nav_panel("倍数变化",
                card_body(
                  tags$p(class = "text-muted",
                         "以第一个分组为参照的相对倍数（log2 标度），全零过程自动省略。"),
                  plotly::plotlyOutput(ns("plot_fold"), height = "520px"),
                  downloadButton(ns("dl_fold_png"), "下载 PNG", class = "btn-sm btn-outline-primary"),
                  hr(), h6("倍数变化表"), DT::DTOutput(ns("tbl_fold")),
                  downloadButton(ns("dl_fold_csv"), "下载 CSV", class = "btn-sm btn-outline-primary mt-1")
                )),
      nav_panel("循环通路图",
                card_body(
                  tags$p(class = "text-muted",
                         "官方循环通路底图，各过程位置叠加该过程的组间倍数变化小热图（丰度全零的过程不叠加）。"),
                  imageOutput(ns("pathway_img"), height = "480px"),
                  downloadButton(ns("dl_pathway"), "下载 PDF", class = "btn-sm btn-outline-primary")
                )),
      nav_panel("宿主菌群",
                uiOutput(ns("host_selector")),
                card_body(
                  plotOutput(ns("plot_host"), height = "420px"),
                  downloadButton(ns("dl_host_png"), "下载 PNG", class = "btn-sm btn-outline-primary"),
                  hr(), h6("当前分类级相对丰度表"), DT::DTOutput(ns("tbl_host")),
                  downloadButton(ns("dl_host_csv"), "下载 CSV", class = "btn-sm btn-outline-primary mt-1")
                )),
      nav_panel("Beta 多样性",
                card_body(
                  layout_columns(col_widths = c(3, 3, 6),
                    selectInput(ns("ord_method"), "排序方法",
                                choices = c("PCoA", "PCA", "NMDS")),
                    selectInput(ns("ord_style"), "展示样式",
                                choices = c("点位" = "points", "95% 椭圆" = "ellipse",
                                            "样本标签" = "labels")),
                    NULL
                  ),
                  plotly::plotlyOutput(ns("plot_ord"), height = "440px"),
                  downloadButton(ns("dl_ord_png"), "下载 PNG", class = "btn-sm btn-outline-primary"),
                  hr(),
                  conditionalPanel(condition = sprintf("input['%s'] == 'PCoA'", ns("ord_method")),
                                   h6("Bray-Curtis 距离矩阵"), DT::DTOutput(ns("tbl_dist")),
                                   downloadButton(ns("dl_dist_csv"), "下载 CSV", class = "btn-sm btn-outline-primary mt-1"),
                                   h6("组间差异检验 (Adonis / ANOSIM / MRPP)"), DT::DTOutput(ns("tbl_diff"))),
                  conditionalPanel(condition = sprintf("input['%s'] == 'NMDS'", ns("ord_method")),
                                   h6("Stress"), verbatimTextOutput(ns("stress_out"))),
                  h6("样本坐标"), DT::DTOutput(ns("tbl_ord")),
                  downloadButton(ns("dl_ord_csv"), "下载 CSV", class = "btn-sm btn-outline-primary mt-1")
                ))
    )
  )
}

cycleTabServer <- function(id, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    cy <- id
    cfg <- CNPS_CONFIG[[cy]]

    res <- reactive({
      r <- results[[cy]]
      if (is.null(r)) NULL else r
    })

    # ---- summary value boxes ----
    output$summary_boxes <- renderUI({
      r <- req(res())
      n_detected <- sum(rowSums(as.matrix(r$gene_abun[, -1, drop = FALSE])) > 0)
      n_sig <- if (!is.null(r$heat_test) && "Normality.p.value" %in% names(r$heat_test))
        sum(r$heat_test$Normality.p.value < 0.05, na.rm = TRUE) else 0
      n_host <- length(r$host)
      layout_columns(col_widths = c(4, 4, 4),
        value_box("检出过程（丰度>0）", sprintf("%d / %d", n_detected, nrow(r$gene_abun)),
                  icon = icon("project-diagram"), theme = "primary"),
        value_box("差异显著过程 (p<0.05)", sprintf("%d", n_sig),
                  icon = icon("chart-line"), theme = "danger"),
        value_box("宿主菌群分析过程", sprintf("%d", n_host),
                  icon = icon("bacteria"), theme = "success")
      )
    })

    # ---- tab 1: abundance tables ----
    gene_dt <- reactive({
      r <- req(res())
      DT::datatable(r$gene_abun, rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE)) |>
        DT::formatRound(columns = 2:ncol(r$gene_abun), digits = 4)
    })
    output$tbl_gene <- DT::renderDT(gene_dt())
    output$dl_gene <- downloadHandler(
      filename = function() sprintf("%s_cycle_gene_abun.csv", cy),
      content = function(f) utils::write.csv(req(res())$gene_abun, f, row.names = FALSE))

    output$tbl_ko <- DT::renderDT({
      r <- req(res())
      DT::datatable(r$ko_abun, rownames = TRUE,
                    options = list(pageLength = 10, scrollX = TRUE)) |>
        DT::formatRound(columns = 1:ncol(r$ko_abun), digits = 4)
    })
    output$dl_ko <- downloadHandler(
      filename = function() sprintf("%s_cycle_ko_abun.csv", cy),
      content = function(f) utils::write.csv(req(res())$ko_abun, f, row.names = TRUE))

    # ---- tab 2: heatmap + diff test ----
    output$plot_heat <- plotly::renderPlotly({
      cnps_to_plotly(req(res())$heat_plot)
    })
    output$dl_heat_png <- downloadHandler(
      filename = function() sprintf("%s_cycle_gene_abun_group.png", cy),
      content = function(f) ggplot2::ggsave(f, req(res())$heat_plot,
        width = cfg$heatmap_w + 0.6 * res()$group_numb, height = cfg$heatmap_h))
    heat_dt <- reactive({
      r <- req(res())
      DT::datatable(r$heat_test, rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE))
    })
    output$tbl_heat <- DT::renderDT(heat_dt())
    output$dl_heat_csv <- downloadHandler(
      filename = function() sprintf("Diff_%s_gene_test.csv", cy),
      content = function(f) utils::write.csv(req(res())$heat_test, f, row.names = FALSE))

    # ---- tab 3: fold change ----
    output$plot_fold <- plotly::renderPlotly({
      cnps_to_plotly(req(res())$fold_plot)
    })
    output$dl_fold_png <- downloadHandler(
      filename = function() sprintf("%s_cycle_gene_fold_change.png", cy),
      content = function(f) ggplot2::ggsave(f, req(res())$fold_plot,
        width = cfg$fold_w + 0.6 * res()$group_numb, height = max(cfg$fold_h, 1.2)))
    output$tbl_fold <- DT::renderDT({
      r <- req(res())
      DT::datatable(r$fold_df, rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE)) |>
        DT::formatRound(columns = 2:ncol(r$fold_df), digits = 4)
    })
    output$dl_fold_csv <- downloadHandler(
      filename = function() sprintf("%s_gene_fold_change.csv", cy),
      content = function(f) utils::write.csv(req(res())$fold_df, f, row.names = FALSE))

    # ---- tab 4: pathway composite ----
    pathway_file <- reactive({
      r <- req(res())
      f <- tempfile(fileext = ".png")
      cnps_pathway_composite(r, f, device = "png")
      f
    })
    output$pathway_img <- renderImage({
      list(src = pathway_file(), contentType = "image/png", width = "100%")
    }, deleteFile = FALSE)
    output$dl_pathway <- downloadHandler(
      filename = function() sprintf("%s_cyc_fold_change.pdf", cy),
      content = function(f) cnps_pathway_composite(req(res()), f, device = "pdf"))

    # ---- tab 5: host taxonomy ----
    output$host_selector <- renderUI({
      r <- req(res())
      if (!length(r$host)) {
        return(card_body(tags$div(class = "alert alert-info",
          "样本中未检出可进行宿主分析的过程（丰度全零或该过程不提供宿主映射）。")))
      }
      codes <- names(r$host)
      lbls <- vapply(codes, function(x) {
        i <- which(cfg$processes$code == x)
        sprintf("%s (%s)", cfg$processes$cn[i], cfg$processes$label[i])
      }, "")
      tagList(
        card_body(
          layout_columns(col_widths = c(6, 6),
            selectInput(ns("host_proc"), "循环过程", choices = setNames(codes, lbls)),
            selectInput(ns("host_level"), "分类级",
                        choices = CNPS_TAX_LEVELS, selected = "Genus")
          )))
    })
    host_tbl <- reactive({
      r <- req(res())
      req(input$host_proc, input$host_level)
      j <- match(input$host_level, unname(CNPS_TAX_LEVELS))
      r$host[[input$host_proc]][[j]]
    })
    output$plot_host <- renderPlot({
      tbl <- req(host_tbl())
      i <- which(cfg$processes$code == input$host_proc)
      cnps_draw_host_plot(tbl, cfg$processes$title[i])
    })
    output$dl_host_png <- downloadHandler(
      filename = function() sprintf("%s_%s.png", input$host_proc, tolower(input$host_level)),
      content = function(f) {
        i <- which(cfg$processes$code == input$host_proc)
        cnps_save_host_png(req(host_tbl()), cfg$processes$title[i], f)
      })
    output$tbl_host <- DT::renderDT({
      tbl <- req(host_tbl())
      DT::datatable(round(as.matrix(tbl), 4), rownames = TRUE,
                    options = list(pageLength = 10, scrollX = TRUE))
    })
    output$dl_host_csv <- downloadHandler(
      filename = function() sprintf("%s_%s.csv", input$host_proc, tolower(input$host_level)),
      content = function(f) utils::write.csv(req(host_tbl()), f, row.names = TRUE))

    # ---- tab 6: beta diversity ----
    ord_data <- reactive({
      r <- req(res()); req(input$ord_method)
      b <- r$beta
      switch(input$ord_method,
             PCoA = if (!is.null(b$pcoa)) list(pd = b$pcoa[[3]], dist = b$pcoa[[1]],
                                               diff = b$pcoa[[2]]) else NULL,
             PCA  = if (!is.null(b$pca))  list(pd = b$pca[[1]]) else NULL,
             NMDS = if (!is.null(b$nmds)) list(pd = b$nmds[[2]], stress = b$nmds[[1]]) else NULL)
    })
    output$plot_ord <- plotly::renderPlotly({
      d <- req(ord_data())
      p <- cnps_ordination_plot(d$pd, input$ord_method, input$ord_style, paste0(cfg$en, " cycling"))
      cnps_to_plotly(p)
    })
    output$dl_ord_png <- downloadHandler(
      filename = function() sprintf("%s_%s_%s.png", cy, tolower(input$ord_method), input$ord_style),
      content = function(f) {
        d <- req(ord_data())
        p <- cnps_ordination_plot(d$pd, input$ord_method, input$ord_style, paste0(cfg$en, " cycling"))
        ggplot2::ggsave(f, p, width = 7.5, height = 5.4, dpi = 200)
      })
    output$tbl_dist <- DT::renderDT({
      d <- req(ord_data())
      m <- round(as.matrix(d$dist), 4)
      DT::datatable(m, rownames = TRUE, options = list(pageLength = 10, scrollX = TRUE))
    })
    output$dl_dist_csv <- downloadHandler(
      filename = function() sprintf("%s_distance_bray_curtis.csv", cy),
      content = function(f) utils::write.csv(round(as.matrix(req(ord_data())$dist), 6), f,
                                             row.names = TRUE))
    output$tbl_diff <- DT::renderDT({
      d <- req(ord_data())
      DT::datatable(as.data.frame(d$diff), rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE))
    })
    output$stress_out <- renderText({
      d <- req(ord_data())
      sprintf("stress = %.4g", as.numeric(d$stress[1]))
    })
    output$tbl_ord <- DT::renderDT({
      d <- req(ord_data())
      DT::datatable(d$pd, rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE)) |>
        DT::formatRound(columns = 2:3, digits = 4)
    })
    output$dl_ord_csv <- downloadHandler(
      filename = function() sprintf("%s_%s_coords.csv", cy, tolower(input$ord_method)),
      content = function(f) utils::write.csv(req(ord_data())$pd, f, row.names = FALSE))
  })
}
