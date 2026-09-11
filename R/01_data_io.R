# =====================================================================
# Data input, validation and preprocessing.
# Reproduces the official "Data Preprocessing" chunk of
# SampleData_AutomatedExecutionScript.Rmd, with the paste -> paste0
# tax fix and robust column handling.
# =====================================================================

# ---- reading ---------------------------------------------------------

CNPS_FILE_SPECS <- list(
  ko        = list(fn = "ko.txt",        cn = "KEGG 丰度表 (KO)",
                   hint = "第1列 KO 编号，中间为各样本丰度，最后一列可为 Description"),
  Gene      = list(fn = "Gene.txt",      cn = "KEGG 基因注释 (Gene)",
                   hint = "至少含 GeneID、Entry 两列"),
  tax       = list(fn = "tax.txt",       cn = "NR 物种注释 (tax)",
                   hint = "两列：GeneID、完整 taxonomy 字符串"),
  abundance = list(fn = "abundance.txt", cn = "基因丰度表 (abundance)",
                   hint = "第1列 GeneID，其余为各样本丰度"),
  group     = list(fn = "group.txt",     cn = "分组文件 (group)",
                   hint = "两列：样本ID、分组")
)

# Read one table; returns data.frame or NULL on failure.
cnps_read_table <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  dt <- tryCatch(
    data.table::fread(path, sep = "\t", header = TRUE, encoding = "UTF-8",
                      check.names = FALSE, data.table = FALSE,
                      showProgress = FALSE, na.strings = c("", "NA")),
    error = function(e) NULL
  )
  if (is.null(dt) && requireNamespace("readr", quietly = TRUE)) {
    dt <- tryCatch(as.data.frame(readr::read_tsv(path, show_col_types = FALSE)),
                   error = function(e) NULL)
  }
  dt
}

# ---- validation ------------------------------------------------------

# Returns a list(ok = logical, problems = character vector)
cnps_validate_raw <- function(raw) {
  probs <- character(0)
  need <- c("ko", "Gene", "tax", "abundance", "group")
  for (n in need) {
    if (is.null(raw[[n]])) probs <- c(probs, sprintf("缺少 %s 表", CNPS_FILE_SPECS[[n]]$cn))
  }
  if (length(probs)) return(list(ok = FALSE, problems = probs))

  ko <- raw$ko; Gene <- raw$Gene; tax <- raw$tax
  ab <- raw$abundance; group <- raw$group

  if (ncol(ko) < 3) probs <- c(probs, "ko 表至少需要 KO 编号列 + 2 个样本列")
  if (ncol(group) != 2) probs <- c(probs, "group 表必须恰好两列（样本ID、分组）")
  if (ncol(tax) < 2) probs <- c(probs, "tax 表至少需要两列（GeneID、taxonomy）")
  if (ncol(Gene) < 2) probs <- c(probs, "Gene 表至少需要两列（GeneID、Entry）")
  if (ncol(ab) < 3) probs <- c(probs, "abundance 表至少需要 GeneID 列 + 2 个样本列")

  if (length(probs)) return(list(ok = FALSE, problems = probs))

  # sample ids consistency: group IDs must be sample columns of ko / abundance
  gids <- as.character(group[[1]])
  ko_cols <- as.character(colnames(ko)[-1])
  ab_cols <- as.character(colnames(ab)[-1])
  # ko 最后一列若是 Description 这类文本列，样本匹配时先剔除
  if (!is.numeric(ko[[ncol(ko)]]) && ncol(ko) > 2) ko_cols <- ko_cols[-length(ko_cols)]
  miss_ko <- setdiff(gids, ko_cols)
  miss_ab <- setdiff(gids, ab_cols)
  if (length(miss_ko)) probs <- c(probs, sprintf("group 中样本不在 ko 表样本列中: %s",
                                                 paste(utils::head(miss_ko, 5), collapse = ", ")))
  if (length(miss_ab)) probs <- c(probs, sprintf("group 中样本不在 abundance 表样本列中: %s",
                                                 paste(utils::head(miss_ab, 5), collapse = ", ")))
  if (anyDuplicated(gids)) probs <- c(probs, "group 表存在重复样本ID")

  # numeric check for ko / abundance sample columns
  ko_body <- if (!is.numeric(ko[[ncol(ko)]]) && ncol(ko) > 2) ko[, -c(1, ncol(ko)), drop = FALSE] else ko[, -1, drop = FALSE]
  if (any(vapply(ko_body, function(x) !is.numeric(x) && any(is.na(suppressWarnings(as.numeric(as.character(x))))), logical(1))))
    probs <- c(probs, "ko 表样本列存在无法转为数值的内容")
  ab_body <- ab[, -1, drop = FALSE]
  if (any(vapply(ab_body, function(x) !is.numeric(x) && any(is.na(suppressWarnings(as.numeric(as.character(x))))), logical(1))))
    probs <- c(probs, "abundance 表样本列存在无法转为数值的内容")

  list(ok = length(probs) == 0, problems = probs)
}

# ---- preprocessing (faithful to the official script) ------------------

cnps_preprocess <- function(raw) {
  ko <- raw$ko
  # 最后一列为非数值的 Description 列则剔除
  if (!is.numeric(ko[[ncol(ko)]]) && ncol(ko) > 2) ko <- ko[, -ncol(ko), drop = FALSE]
  rownames(ko) <- as.character(ko[, 1])
  ko <- ko[, -1, drop = FALSE]
  for (i in seq_len(ncol(ko))) ko[[i]] <- as.numeric(ko[[i]])

  group <- raw$group[, 1:2]
  colnames(group) <- c("ID", "Group")
  group$ID <- as.character(group$ID)
  group$Group <- factor(as.character(group$Group),
                        levels = unique(as.character(group$Group)))
  Group_numb  <- length(unique(group$Group))
  Sample_numb <- length(unique(group$ID))

  Gene <- raw$Gene[!duplicated(raw$Gene[, 1:2]), 1:2, drop = FALSE]
  colnames(Gene) <- c("GeneID", "Entry")
  Gene <- Gene %>%
    dplyr::group_by(Entry) %>%
    dplyr::mutate(index = dplyr::row_number()) %>%
    tidyr::pivot_wider(names_from = Entry, values_from = GeneID) %>%
    dplyr::select(-index) %>%
    as.data.frame()

  tax <- raw$tax[, 1:2]
  colnames(tax) <- c("V1", "V2")
  tax$V1 <- as.character(tax$V1)
  tax$V2 <- paste0("k__", gsub(".*;k__", "", as.character(tax$V2)))

  abundance <- raw$abundance
  rownames(abundance) <- as.character(abundance[, 1])
  abundance <- abundance[, -1, drop = FALSE]
  for (i in seq_len(ncol(abundance))) abundance[[i]] <- as.numeric(abundance[[i]])
  abundance$V1 <- rownames(abundance)
  abundance <- abundance[, c("V1", setdiff(colnames(abundance), "V1")), drop = FALSE]

  # 只保留 group 中出现的样本列，顺序与 group 一致
  gids <- group$ID
  ko <- ko[, gids, drop = FALSE]
  abundance <- abundance[, c("V1", gids), drop = FALSE]

  list(ko = ko, Gene = Gene, tax = tax, abundance = abundance,
       group = group, Group_numb = Group_numb, Sample_numb = Sample_numb)
}

# ---- sample data (bundled with the package) ---------------------------

cnps_load_sample_data <- function(envir = globalenv()) {
  for (n in c("ko", "Gene", "tax", "abundance", "group")) {
    if (!exists(n, envir = envir, inherits = FALSE)) {
      data(list = n, package = "CNPS.cycle", envir = envir)
    }
  }
  raw <- list(ko = envir$ko, Gene = envir$Gene, tax = envir$tax,
              abundance = envir$abundance, group = envir$group)
  # 包内置数据即为已预处理形态的直接来源，但 ko/Gene/tax/abundance 与
  # 官方脚本开头一致：ko 带 Entry+Description 列，Gene 长表，tax 完整字符串
  out <- cnps_preprocess(raw)
  out$raw_head <- lapply(raw, function(x) utils::head(as.data.frame(x), 50))
  out
}
