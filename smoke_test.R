# Smoke test for the CNPS Shiny app pipeline.
#
# Runs all four cycles through the app's data-preprocessing + pipeline
# adapter and, when the official Rmd baseline outputs are available,
# compares them value-by-value.
#
# Environment overrides (all optional):
#   CNPS_APP_DIR    - app directory           (default: this script's dir)
#   CNPS_BASELINE   - official Results tree   (default: D:/share/Projects/CNPS/run/Results)
#   CNPS_SAMPLE_TXT - dir with ko/Gene/tax/abundance/group .txt files;
#                     when absent, the package's built-in sample data is used.
#
# Run (UTF-8 native locale keeps the Chinese messages readable on Windows):
#   set LC_CTYPE=.UTF-8
#   Rscript smoke_test.R

options(warn = 1)
APP_DIR <- Sys.getenv("CNPS_APP_DIR", getwd())
BASELINE <- Sys.getenv("CNPS_BASELINE", "D:/share/Projects/CNPS/run/Results")
DATA_DIR <- Sys.getenv("CNPS_SAMPLE_TXT", "D:/share/Projects/CNPS/data")

for (f in list.files(file.path(APP_DIR, "R"), full.names = TRUE)) source(f, encoding = "UTF-8")
library(CNPS.cycle)

HAS_BASELINE <- dir.exists(BASELINE)
failures <- character(0)
check_num <- function(name, got, want, tol = 1e-6) {
  g <- suppressWarnings(as.numeric(got)); w <- suppressWarnings(as.numeric(want))
  if (length(g) != length(w) || any(is.na(g)) || any(is.na(w)) || max(abs(g - w)) > tol) {
    cat("  MISMATCH", name, ": len", length(g), "vs", length(w),
        " maxdiff", if (length(g) == length(w) && !any(is.na(g)) && !any(is.na(w))) max(abs(g - w)) else NA, "\n")
    failures <<- c(failures, name)
  } else cat("  ok", name, "\n")
}

cat("== 1. load input data ==\n")
txt_files <- file.path(DATA_DIR, c("ko.txt", "Gene.txt", "tax.txt", "abundance.txt", "group.txt"))
if (all(file.exists(txt_files))) {
  cat("  using uploaded-style txt files from", DATA_DIR, "\n")
  raw <- list(
    ko        = cnps_read_table(txt_files[1]),
    Gene      = cnps_read_table(txt_files[2]),
    tax       = cnps_read_table(txt_files[3]),
    abundance = cnps_read_table(txt_files[4]),
    group     = cnps_read_table(txt_files[5])
  )
  for (n in names(raw)) cat(sprintf("  %s: %d x %d\n", n, nrow(raw[[n]]), ncol(raw[[n]])))
  v <- cnps_validate_raw(raw)
  cat("  validate ok:", v$ok, "\n")
  if (!v$ok) { cat(paste("  PROBLEM:", v$problems, collapse = "\n"), "\n"); stop("validation failed") }
  dat <- cnps_preprocess(raw)
} else {
  cat("  txt files not found - using the package's built-in sample data\n")
  dat <- cnps_load_sample_data(envir = new.env())
}
cat(sprintf("  preprocessed: ko %d x %d | gene_wide %d x %d | tax %d x %d | abundance %d x %d | groups %d\n",
            nrow(dat$ko), ncol(dat$ko), nrow(dat$Gene), ncol(dat$Gene),
            nrow(dat$tax), ncol(dat$tax), nrow(dat$abundance), ncol(dat$abundance),
            dat$Group_numb))
cat("  group levels:", paste(levels(dat$group$Group), collapse = ","), "\n")

results <- list()
for (cy in CNPS_CYCLE_KEYS) {
  cat(sprintf("== 2. run cycle %s ==\n", cy))
  t0 <- Sys.time()
  results[[cy]] <- cnps_run_cycle(cy, dat, progress = function(m) cat("   .", m, "\n"))
  cat(sprintf("   done in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

if (HAS_BASELINE) {
  cat("== 3. compare against official baseline ==\n")
  cmp_gene_abun <- function(cy) {
    p <- file.path(BASELINE, cycle_name_dir(cy), "Gene", "Abundance",
                   sprintf("%s_cycle_gene_abun.txt", cy))
    b <- read.delim(p, check.names = FALSE)
    r <- results[[cy]]$gene_abun
    check_num(paste0(cy, " gene_abun values"),
              as.matrix(r[, -1]), as.matrix(b[, -1]))
    if (!identical(as.character(r[[1]]), as.character(b[[1]])))
      { cat("  LABEL MISMATCH", cy, "\n"); failures <<- c(failures, paste0(cy, " labels")) }
    else cat("  ok", cy, "labels\n")
    if (!identical(colnames(r)[-1], colnames(b)[-1]))
      { cat("  SAMPLE COL MISMATCH", cy, "\n"); failures <<- c(failures, paste0(cy, " samples")) }
    else cat("  ok", cy, "sample cols\n")
  }
  cmp_ko_abun <- function(cy) {
    p <- file.path(BASELINE, cycle_name_dir(cy), "Gene", "Abundance",
                   sprintf("%s_cycle_ko_abun.txt", cy))
    b <- read.delim(p, check.names = FALSE, row.names = 1)
    r <- results[[cy]]$ko_abun
    common <- intersect(rownames(b), rownames(r))
    check_num(paste0(cy, " ko_abun"), as.matrix(r[common, ]), as.matrix(b[common, ]))
  }
  cmp_heat <- function(cy) {
    p <- file.path(BASELINE, cycle_name_dir(cy), "Gene", "Heatmap",
                   sprintf("Diff_%s_gene_test.txt", cy))
    b <- read.delim(p, check.names = FALSE)
    r <- results[[cy]]$heat_test
    m <- merge(r, b, by = "Type", suffixes = c(".r", ".b"))
    check_num(paste0(cy, " heat p-values"), m$Normality.p.value.r, m$Normality.p.value.b)
  }
  cmp_fold <- function(cy) {
    p <- file.path(BASELINE, cycle_name_dir(cy), "Gene", "Cycle image", "Gene_fold_change.txt")
    b <- read.delim(p, check.names = FALSE)
    r <- results[[cy]]$fold_df
    check_num(paste0(cy, " fold change"), as.matrix(r[, -1]), as.matrix(b[, -1]))
  }
  cmp_beta <- function(cy) {
    b_pcoa <- read.delim(file.path(BASELINE, cycle_name_dir(cy), "Beta diversity", "PCoA",
                                   sprintf("%s_pcoa.txt", cy)), check.names = FALSE)
    r_pcoa <- results[[cy]]$beta$pcoa[[3]]
    # compare numeric coordinate columns only (both sides end with a text Group col)
    check_num(paste0(cy, " pcoa scores"), as.matrix(r_pcoa[, 2:3]), as.matrix(b_pcoa[, 2:3]))
    b_dist <- read.delim(file.path(BASELINE, cycle_name_dir(cy), "Beta diversity", "Distance",
                                   "distance_bray_curtis.txt"), check.names = FALSE, row.names = 1)
    check_num(paste0(cy, " bray-curtis"), as.matrix(results[[cy]]$beta$pcoa[[1]]), as.matrix(b_dist))
    r_nmds <- results[[cy]]$beta$nmds[[2]]
    # metaMDS is stochastic (random starts): coordinates are NOT comparable
    # across sessions. Verify finiteness here; reproducibility comes from
    # set.seed() in the pipeline.
    ok_nmds <- all(is.finite(as.matrix(r_nmds[, 2:3])))
    cat(sprintf("  %s nmds finite: %s (stress=%.4g)\n", cy, ok_nmds,
                as.numeric(results[[cy]]$beta$nmds[[1]][1])))
    if (!ok_nmds) failures <<- c(failures, paste0(cy, " nmds not finite"))
  }
  cmp_host <- function(cy) {
    n_expected <- length(list.dirs(file.path(BASELINE, cycle_name_dir(cy), "Host_relative_Group"),
                                   recursive = FALSE))
    n_got <- length(results[[cy]]$host)
    cat(sprintf("  %s host processes: got %d, baseline dirs %d\n", cy, n_got, n_expected))
    if (n_got == 0) failures <<- c(failures, paste0(cy, " host empty"))
    for (code in names(results[[cy]]$host)) {
      p <- file.path(BASELINE, cycle_name_dir(cy), "Host_relative_Group", code,
                     sprintf("%s_genus.csv", code))
      if (file.exists(p)) {
        b <- read.csv(p, check.names = FALSE, row.names = 1)
        r <- results[[cy]]$host[[code]][[5]]
        rn <- intersect(rownames(b), rownames(r))
        check_num(paste0(cy, "/", code, " genus"), as.matrix(r[rn, ]), as.matrix(b[rn, ]))
        break
      }
    }
  }
  for (cy in CNPS_CYCLE_KEYS) {
    cmp_gene_abun(cy); cmp_ko_abun(cy); cmp_heat(cy); cmp_fold(cy); cmp_beta(cy); cmp_host(cy)
  }
} else {
  cat("== 3. baseline not available - structural sanity checks only ==\n")
  for (cy in CNPS_CYCLE_KEYS) {
    r <- results[[cy]]
    ok <- is.data.frame(r$gene_abun) && nrow(r$gene_abun) > 0 &&
          nrow(r$ko_abun) > 0 && !is.null(r$heat_test) &&
          !is.null(r$fold_df) && !is.null(r$beta$pcoa)
    cat(sprintf("  %s outputs present: %s (host processes: %d)\n", cy, ok, length(r$host)))
    if (!ok) failures <- c(failures, paste0(cy, " incomplete"))
  }
}

# globals must be cleaned up (inherits=FALSE: the package's own lazy data
# sits on the search path and would otherwise cause false positives)
leftover <- CNPS_GLOBAL_NAMES[vapply(CNPS_GLOBAL_NAMES, function(n)
  exists(n, envir = .GlobalEnv, inherits = FALSE), logical(1))]
cat("== 4. globals cleaned:", length(leftover) == 0,
    if (length(leftover)) paste("leftover:", paste(leftover, collapse = ",")) else "", "\n")
if (length(leftover)) failures <- c(failures, "globals not cleaned")

cat("\n====================================\n")
if (length(failures)) {
  cat("SMOKE TEST FAILURES:", length(failures), "\n", paste(unique(failures), collapse = "\n"), "\n")
  quit(status = 1)
} else {
  cat("SMOKE TEST PASSED\n")
}
