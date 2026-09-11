# =====================================================================
# CNPS 元素循环分析平台 - Shiny app entry point
# =====================================================================

library(shiny)
library(bslib)
library(DT)
library(plotly)
library(shinyjs)
library(shinybusy)
library(CNPS.cycle)

# app.R is sourced with cwd = app directory (shiny::runApp behavior)
for (f in list.files("R", full.names = TRUE, pattern = "\\.R$")) source(f, encoding = "UTF-8")

# ---- export tab ------------------------------------------------------
exportUI <- function(id) {
  ns <- NS(id)
  card(
    card_header(class = "bg-primary text-white", "结果打包导出"),
    card_body(
      tags$p("将本次分析的全部结果按官方 Results 目录结构打包为 zip："),
      tags$ul(
        tags$li("Gene/Abundance — KO 级与过程级丰度表"),
        tags$li("Gene/Heatmap — 分组热图 PDF 与差异检验表"),
        tags$li("Gene/Cycle image — 倍数变化图与循环通路合成图 PDF"),
        tags$li("Host_relative_Group — 宿主菌群各级相对丰度表与图"),
        tags$li("Beta diversity — 距离矩阵、差异检验、PCoA/PCA/NMDS 图表")
      ),
      uiOutput(ns("export_state")),
      downloadButton(ns("dl_zip"), "下载全部结果 (zip)",
                     class = "btn btn-primary btn-lg")
    )
  )
}
exportServer <- function(id, results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$export_state <- renderUI({
      if (isTRUE(results$ready) && !is.null(results$export_dir) && dir.exists(results$export_dir))
        tags$div(class = "alert alert-success", icon("check-circle"),
                 "分析结果已就绪，可打包下载。")
      else
        tags$div(class = "alert alert-secondary", "尚未运行分析——请先在「数据与运行」页完成分析。")
    })
    output$dl_zip <- downloadHandler(
      filename = function() sprintf("CNPS_Results_%s.zip", format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(f) {
        req(isTRUE(results$ready), results$export_dir, dir.exists(results$export_dir))
        owd <- setwd(results$export_dir); on.exit(setwd(owd), add = TRUE)
        zip::zip(zipfile = f, files = list.files("."), recurse = TRUE)
      },
      contentType = "application/zip")
  })
}

# ---- home tab --------------------------------------------------------
homeUI <- function() {
  shinyjs::useShinyjs()
  layout_columns(col_widths = c(4, 4, 4),
    card(card_body(
      h3(icon("info-circle"), "这是什么"),
      tags$p("基于", tags$b("CNPS.cycle"), " R 包的碳(C)/氮(N)/磷(P)/硫(S) 元素循环宏基因组分析平台。"),
      tags$p("上传 KEGG/NR 注释与丰度表后，一键完成：过程丰度聚合、分组差异检验、倍数变化、循环通路图、宿主菌群组成与 Beta 多样性分析。")
    )),
    card(card_body(
      h3(icon("list-ul"), "使用步骤"),
      tags$ol(
        tags$li("「数据与运行」页：点击", tags$b("使用内置示例数据"), "快速体验，或上传自己的 5 个数据文件"),
        tags$li("勾选要分析的循环，点击", tags$b("开始分析")),
        tags$li("在各循环页浏览交互图表，逐图下载或到「结果导出」页打包全部结果")
      )
    )),
    card(card_body(
      h3(icon("flask"), "示例数据"),
      tags$p("内置示例为 9 个样本、3 个分组（PC1/PC17/PC30）的宏基因组 KEGG/NR 注释数据。"),
      tags$p(class = "text-muted", "分析核心与官方自动化脚本完全一致，结果已通过与官方基准的逐值对拍验证。")
    ))
  )
}

# ---- theme & ui ------------------------------------------------------
theme <- bs_theme(
  version = 5, bootswatch = "flatly",
  primary = "#2c6e8f", success = "#3c9d6e"
)

ui <- page_navbar(
  title = tagList(icon("dna"), " CNPS 元素循环分析平台"),
  theme = theme,
  id = "main_nav",
  nav_panel("首页", homeUI()),
  nav_panel("数据与运行", uploadUI("upload")),
  nav_panel("碳循环", cycleTabUI("C")),
  nav_panel("氮循环", cycleTabUI("N")),
  nav_panel("磷循环", cycleTabUI("P")),
  nav_panel("硫循环", cycleTabUI("S")),
  nav_panel("结果导出", exportUI("export"))
)

server <- function(input, output, session) {
  results <- shiny::reactiveValues(C = NULL, N = NULL, P = NULL, S = NULL,
                                   ready = FALSE, export_dir = NULL)
  uploadServer("upload", results)
  cycleTabServer("C", results)
  cycleTabServer("N", results)
  cycleTabServer("P", results)
  cycleTabServer("S", results)
  exportServer("export", results)
}

shinyApp(ui, server)
