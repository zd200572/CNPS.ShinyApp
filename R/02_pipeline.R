# =====================================================================
# Core pipeline adapter.
#
# IMPORTANT (verified against package source + runtime experiment):
# Ccyc/Ncyc/Pcyc/Scyc.abundance and all 42 host-taxonomy functions
# IGNORE their formal arguments and resolve the free names
# ko / Gene / tax / abundance from the calling scope chain, which for
# package functions terminates at the global environment. The user's
# preprocessed tables must therefore be assigned to .GlobalEnv under
# those exact names before calling them (exactly what the official Rmd
# does), and must be removed afterwards. `group` alone is a real formal.
# abun.heatmap.g / fold.change / pcoa.arg / pca.arg / nmds.arg honor
# their arguments, but the ordination functions need a global
# `cbbPalette` which the package does not define.
# =====================================================================

CNPS_GLOBAL_NAMES <- c("ko", "Gene", "tax", "abundance", "group", "cbbPalette")

cnps_bind_globals <- function(dat) {
  assign("ko", dat$ko, envir = .GlobalEnv)
  assign("Gene", dat$Gene, envir = .GlobalEnv)
  assign("tax", dat$tax, envir = .GlobalEnv)
  assign("abundance", dat$abundance, envir = .GlobalEnv)
  assign("group", dat$group, envir = .GlobalEnv)
  assign("cbbPalette", CNPS_PALETTE, envir = .GlobalEnv)
  invisible(TRUE)
}

cnps_unbind_globals <- function() {
  rm(list = intersect(CNPS_GLOBAL_NAMES, ls(envir = .GlobalEnv)), envir = .GlobalEnv)
  invisible(TRUE)
}

# Package namespace handle (functions resolved from there, not the search path)
cnps_ns <- function() asNamespace("CNPS.cycle")

# ---- one full cycle ----------------------------------------------------
# cycle: "C" / "N" / "P" / "S"
# dat:   preprocessed list from cnps_preprocess()
# progress: optional function(msg) for progress updates
# Returns a structured list with all tables and plots.
cnps_run_cycle <- function(cycle, dat, progress = function(msg) NULL) {
  stopifnot(cycle %in% CNPS_CYCLE_KEYS)
  cfg <- CNPS_CONFIG[[cycle]]
  ns <- cnps_ns()

  cnps_bind_globals(dat)
  on.exit(cnps_unbind_globals(), add = TRUE)

  abun_fun <- get(cfg$abun_fun, envir = ns)
  host_flag_fun <- get(cfg$host_fun, envir = ns)

  # 1) KO-level abundance table (subset of global ko, as the official script)
  progress(sprintf("%s：提取循环相关 KO 丰度", cfg$name))
  kos <- CNPS_CYCLE_KOS[[cycle]]
  ko_abun <- dat$ko[rownames(dat$ko) %in% kos, , drop = FALSE]

  # 2) process-level abundance aggregation (reads global `ko`)
  progress(sprintf("%s：聚合 %d 个循环过程丰度", cfg$name, nrow(cfg$processes)))
  gene_abun <- abun_fun()
  gene_abun <- as.data.frame(gene_abun)

  # 3) group heatmap + significance tests (real formals)
  progress(sprintf("%s：分组热图与差异检验", cfg$name))
  hm <- abun.heatmap.g(gene_abun, dat$group, dat$Group_numb)

  # 4) fold change (real formals); needs >= 2 groups
  progress(sprintf("%s：倍数变化分析", cfg$name))
  fc <- fold.change(gene_abun, dat$group)

  # 5) host-taxonomy flags, then per-process host analysis
  progress(sprintf("%s：检测宿主分析可用性", cfg$name))
  flags <- host_flag_fun(dat$Gene)
  host <- list(); host_skipped <- character(0)
  for (i in seq_len(nrow(cfg$processes))) {
    pr <- cfg$processes[i, ]
    has_abun <- sum(as.numeric(gene_abun[i, 2:ncol(gene_abun)])) > 0
    flag_ok <- identical(unname(flags[i]), 1)
    if (!has_abun || !pr$host || !flag_ok) {
      if (has_abun && pr$host && !flag_ok) host_skipped <- c(host_skipped, pr$code)
      next
    }
    progress(sprintf("%s：宿主菌群分析 - %s", cfg$name, pr$label))
    host_fn <- get(pr$code, envir = ns)
    res <- tryCatch(
      host_fn(dat$Gene, dat$tax, dat$abundance, dat$group),
      error = function(e) {
        warning(sprintf("宿主分析 %s 失败: %s", pr$code, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res)) host[[pr$code]] <- res
  }

  # 6) beta diversity on the process-level matrix (samples x processes)
  progress(sprintf("%s：Beta 多样性 (PCoA / PCA / NMDS)", cfg$name))
  sub <- gene_abun[, 2:ncol(gene_abun), drop = FALSE]
  beta <- list()
  beta$pcoa <- tryCatch(pcoa.arg(sub, dat$group), error = function(e) NULL)
  beta$pca  <- tryCatch(pca.arg(sub, dat$group),  error = function(e) NULL)
  set.seed(123)  # metaMDS random starts -> reproducible runs
  beta$nmds <- tryCatch(nmds.arg(sub, dat$group), error = function(e) NULL)

  list(cycle = cycle, cfg = cfg, group_numb = dat$Group_numb,
       ko_abun = ko_abun, gene_abun = gene_abun,
       heat_test = hm[[1]], heat_plot = hm[[2]],
       fold_df = fc[[1]], fold_heat = fc[[2]], fold_plot = fc[[3]],
       host_flags = flags, host = host, host_skipped = host_skipped,
       beta = beta)
}

# ---- export the full official-style Results tree ------------------------
# root: output directory (created recursively)
cnps_write_results <- function(res, root) {
  cycle <- res$cycle
  cfg <- res$cfg
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, cycle_name_dir(cycle), "Gene", "Abundance"),
             recursive = TRUE, showWarnings = FALSE)
  base <- file.path(root, cycle_name_dir(cycle))

  utils::write.table(res$ko_abun, file.path(base, "Gene", "Abundance",
                    sprintf("%s_cycle_ko_abun.txt", cycle)), sep = "\t")
  utils::write.table(res$gene_abun, file.path(base, "Gene", "Abundance",
                     sprintf("%s_cycle_gene_abun.txt", cycle)), sep = "\t", row.names = FALSE)

  dir.create(file.path(base, "Gene", "Heatmap"), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(res$heat_test, file.path(base, "Gene", "Heatmap",
                     sprintf("Diff_%s_gene_test.txt", cycle)), sep = "\t", row.names = FALSE)

  dir.create(file.path(base, "Gene", "Cycle image"), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(res$fold_df, file.path(base, "Gene", "Cycle image",
                     "Gene_fold_change.txt"), sep = "\t", row.names = FALSE)

  dir.create(file.path(base, "Host_relative_Group"), recursive = TRUE, showWarnings = FALSE)
  for (code in names(res$host)) {
    d <- file.path(base, "Host_relative_Group", code)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    for (j in seq_along(res$host[[code]])) {
      utils::write.csv(res$host[[code]][[j]],
                       file.path(d, sprintf("%s_%s.csv", code,
                                            c("phylum","class","order","family","genus","species")[j])),
                       row.names = FALSE)
    }
  }
  invisible(base)
}

cycle_name_dir <- function(cycle) c(C = "Carbon", N = "Nitrogen",
                                    P = "Phosphorus", S = "Sulfur")[cycle]
