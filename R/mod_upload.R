# =====================================================================
# Data upload / sample-data module + run control.
# Lightweight structural validation on upload (headers only); the heavy
# full read + preprocessing happens once on "开始分析".
# =====================================================================

uploadUI <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 380,
        card(
          card_header(class = "bg-primary text-white", "① 数据来源"),
          card_body(
            actionButton(ns("use_sample"), "使用内置示例数据",
                         icon = icon("magic"),
                         class = "btn btn-primary w-100 mb-3"),
          tags$label(class = "text-muted", "或上传自己的数据（Tab 分隔 txt）"),
          fileInput(ns("file_ko"), CNPS_FILE_SPECS$ko$cn, accept = c(".txt", ".tsv", ".csv"),
                    placeholder = "ko.txt"),
          tags$small(class = "text-muted", CNPS_FILE_SPECS$ko$hint),
          fileInput(ns("file_Gene"), CNPS_FILE_SPECS$Gene$cn, accept = c(".txt", ".tsv", ".csv"),
                    placeholder = "Gene.txt"),
          tags$small(class = "text-muted", CNPS_FILE_SPECS$Gene$hint),
          fileInput(ns("file_tax"), CNPS_FILE_SPECS$tax$cn, accept = c(".txt", ".tsv", ".csv"),
                    placeholder = "tax.txt"),
          tags$small(class = "text-muted", CNPS_FILE_SPECS$tax$hint),
          fileInput(ns("file_abundance"), CNPS_FILE_SPECS$abundance$cn,
                    accept = c(".txt", ".tsv", ".csv"), placeholder = "abundance.txt"),
          tags$small(class = "text-muted", CNPS_FILE_SPECS$abundance$hint),
          fileInput(ns("file_group"), CNPS_FILE_SPECS$group$cn, accept = c(".txt", ".tsv", ".csv"),
                    placeholder = "group.txt"),
          tags$small(class = "text-muted", CNPS_FILE_SPECS$group$hint)
        )
      ),
      card(
        card_header(class = "bg-primary text-white", "② 运行设置"),
        card_body(
          checkboxGroupInput(ns("cycles"), "选择要分析的循环：",
                             choices = c(碳循环 = "C", 氮循环 = "N", 磷循环 = "P", 硫循环 = "S"),
                             selected = c("C", "N", "P", "S"), inline = TRUE),
          actionButton(ns("run"), "开始分析",
                       icon = icon("play"),
                       class = "btn btn-success btn-lg w-100"),
          uiOutput(ns("run_state"))
        )
      )
    ),
    card(
      card_header(class = "bg-primary text-white", "③ 数据校验与预览"),
      card_body(
        uiOutput(ns("validation")),
        navset_card_tab(
          id = ns("preview_tabs"),
          nav_panel("KO 丰度表", DT::DTOutput(ns("prev_ko"))),
          nav_panel("Gene 注释", DT::DTOutput(ns("prev_Gene"))),
          nav_panel("Tax 注释", DT::DTOutput(ns("prev_tax"))),
          nav_panel("基因丰度表", DT::DTOutput(ns("prev_abundance"))),
          nav_panel("分组", DT::DTOutput(ns("prev_group")))
        )
      )
    ),
  )
}

uploadServer <- function(id, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # state: list(mode = "sample"/"upload", dat = preprocessed list, ready = bool)
    state <- reactiveValues(mode = NULL, dat = NULL, ready = FALSE, running = FALSE)

    upload_ids <- c("file_ko", "file_Gene", "file_tax", "file_abundance", "file_group")
    uploaded_paths <- reactive({
      setNames(lapply(upload_ids, function(x) {
        p <- input[[x]]$datapath; if (is.null(p)) "" else p
      }), c("ko", "Gene", "tax", "abundance", "group"))
    })

    # ---- sample data ----
    observeEvent(input$use_sample, {
      shinybusy::show_modal_spinner(spin = "folding-cube", text = "正在载入内置示例数据（约 400 MB，请稍候）……")
      dat <- tryCatch(cnps_load_sample_data(envir = new.env()),
                      error = function(e) e)
      shinybusy::remove_modal_spinner()
      if (inherits(dat, "error")) {
        showNotification(paste("示例数据载入失败：", conditionMessage(dat)), type = "error")
        return()
      }
      state$mode <- "sample"
      state$dat <- dat
      state$ready <- TRUE
      showNotification("示例数据已就绪，可直接点击「开始分析」", type = "message")
      for (x in upload_ids) try(shinyjs::disable(x), silent = TRUE)
    })

    # reset to upload mode when any file input changes (observeEvent handler
    # reads do NOT create reactive dependencies, unlike observe())
    observeEvent(uploaded_paths(), {
      if (identical(state$mode, "sample")) {
        state$mode <- "upload"; state$dat <- NULL; state$ready <- FALSE
        for (x in upload_ids) try(shinyjs::enable(x), silent = TRUE)
      }
    }, ignoreInit = TRUE)

    # ---- lightweight structural validation (headers only) ----
    light_check <- reactive({
      req(isTRUE(state$mode == "upload"))
      paths <- uploaded_paths()
      miss <- names(paths)[!nzchar(paths)]
      probs <- character(0)
      if (length(miss)) {
        probs <- c(probs, sprintf("尚未上传: %s",
                    paste(vapply(miss, function(n) CNPS_FILE_SPECS[[n]]$cn, ""), collapse = "、")))
      } else {
        heads <- lapply(paths, function(p) tryCatch(data.table::fread(p, nrows = 20, sep = "\t",
                          header = TRUE, check.names = FALSE, data.table = FALSE),
                          error = function(e) NULL))
        v <- cnps_validate_raw(heads)  # structural rules work on 20-row samples
        probs <- v$problems
        n_ok <- sum(vapply(heads, function(x) !is.null(x), logical(1)))
      }
      list(ok = length(probs) == 0, problems = probs,
           complete = !length(miss))
    })

    output$validation <- renderUI({
      if (isTRUE(state$mode == "sample")) {
        g <- dat_group_summary(state$dat)
        return(tags$div(class = "alert alert-success",
          tags$b("✓ 内置示例数据已就绪"),
          sprintf("（%d 个样本、%d 个分组）", g$n_sample, g$n_group)))
      }
      if (!isTRUE(state$mode == "upload")) {
        return(tags$div(class = "alert alert-secondary", "请上传 5 个数据文件，或点击「使用内置示例数据」快速体验。"))
      }
      ck <- light_check()
      if (!ck$complete) {
        return(tags$div(class = "alert alert-warning",
          tags$b("等待上传"), tags$ul(lapply(ck$problems, tags$li))))
      }
      if (ck$ok) {
        return(tags$div(class = "alert alert-success",
          tags$b("✓ 数据结构校验通过"), "可点击「开始分析」运行完整流程。"))
      }
      tags$div(class = "alert alert-danger",
        tags$b("校验未通过"), tags$ul(lapply(ck$problems, tags$li)))
    })

    output$run_state <- renderUI({
      if (isTRUE(state$running)) return(tags$div(class = "text-muted mt-2", icon("spinner", class = "fa-spin"), " 分析进行中……"))
      if (isTRUE(state$ready)) return(tags$div(class = "text-success mt-2", icon("circle-check"), " 数据就绪"))
      NULL
    })

    # ---- preview tables ----
    preview_data <- reactive({
      if (isTRUE(state$mode == "sample")) return(state$dat$raw_head)
      if (!isTRUE(state$mode == "upload")) return(NULL)
      paths <- uploaded_paths()
      if (!all(nzchar(paths))) return(NULL)
      lapply(paths, function(p) tryCatch(
        data.table::fread(p, nrows = 50, sep = "\t", header = TRUE, check.names = FALSE,
                          data.table = FALSE), error = function(e) NULL))
    })
    dat_group_summary <- function(dat) {
      list(n_sample = dat$Sample_numb, n_group = dat$Group_numb)
    }
    dn <- function(x) DT::datatable(x, options = list(pageLength = 5, scrollX = TRUE, dom = "tp"),
                                    rownames = FALSE)
    output$prev_ko <- DT::renderDT({ d <- preview_data(); if (!is.null(d$ko)) dn(d$ko) })
    output$prev_Gene <- DT::renderDT({ d <- preview_data(); if (!is.null(d$Gene)) dn(d$Gene) })
    output$prev_tax <- DT::renderDT({ d <- preview_data(); if (!is.null(d$tax)) dn(d$tax) })
    output$prev_abundance <- DT::renderDT({ d <- preview_data(); if (!is.null(d$abundance)) dn(d$abundance) })
    output$prev_group <- DT::renderDT({ d <- preview_data(); if (!is.null(d$group)) dn(d$group) })

    # ---- run ----
    observeEvent(input$run, {
      if (isTRUE(state$running)) return()
      if (!isTRUE(state$ready) && state$mode != "upload") {
        showNotification("请先载入示例数据或上传完整数据文件", type = "warning"); return()
      }
      cycles <- input$cycles
      if (!length(cycles)) { showNotification("请至少选择一个循环", type = "warning"); return() }
      state$running <- TRUE
      shinybusy::show_modal_spinner(spin = "folding-cube",
                                    text = "正在读取并预处理数据……首次分析约需数分钟")

      res_all <- list()
      ok_cycles <- character(0)
      tryCatch({
        if (isTRUE(state$mode == "sample")) {
          dat <- state$dat
        } else {
          paths <- uploaded_paths()
          if (!all(nzchar(paths))) stop("数据文件不完整")
          raw <- lapply(paths, cnps_read_table)
          v <- cnps_validate_raw(raw)
          if (!v$ok) stop(paste(v$problems, collapse = "；"))
          dat <- cnps_preprocess(raw)
        }
        export_root <- file.path(tempdir(), "CNPS_export")
        if (dir.exists(export_root)) unlink(export_root, recursive = TRUE)
        for (cy in cycles) {
          r1 <- tryCatch({
            res_all[[cy]] <- cnps_run_cycle(cy, dat, progress = function(m) {
              shinybusy::update_modal_spinner(text = m)
            })
            cnps_write_results(res_all[[cy]], export_root)
            cnps_write_plots(res_all[[cy]], export_root)
            TRUE
          }, error = function(e) {
            showNotification(sprintf("%s 分析失败：%s", CNPS_CONFIG[[cy]]$name,
                                     conditionMessage(e)), type = "error", duration = NULL)
            FALSE
          })
          if (isTRUE(r1)) ok_cycles <- c(ok_cycles, cy)
        }
      }, error = function(e) {
        showNotification(paste("分析失败：", conditionMessage(e)), type = "error", duration = NULL)
      })
      shinybusy::remove_modal_spinner()
      state$running <- FALSE
      if (length(ok_cycles)) {
        for (cy in ok_cycles) results[[cy]] <- res_all[[cy]]
        results$export_dir <- export_root
        results$ready <- TRUE
        showNotification(sprintf("分析完成：%s", paste(vapply(ok_cycles, function(c) CNPS_CONFIG[[c]]$name, ""), collapse = "、")),
                          type = "message")
        if (length(ok_cycles) < length(cycles)) {
          miss <- setdiff(cycles, ok_cycles)
          showNotification(sprintf("未完成：%s（详见日志）", paste(miss, collapse = "、")), type = "warning")
        }
      }
    })

    state
  })
}
