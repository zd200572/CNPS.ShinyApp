# =====================================================================
# Plot helpers: pathway composite figure, host-ratio drawing, and
# rebuilt ordination plots (package plotdata is verified by the smoke
# test; the app redraws them with correct cycle titles).
# =====================================================================

# ---- pathway composite: background PDF + fold-change mini heatmaps ----
# Faithful port of the official "Visualization of the Cycle Pathway
# Diagram" chunks; viewport coordinates come from the config table.
cnps_pathway_composite <- function(res, outfile, width = 13, height = 7,
                                   device = c("pdf", "png")) {
  device <- match.arg(device)
  cfg <- res$cfg
  abun <- res$gene_abun
  img <- system.file("data", cfg$bg_pdf, package = "CNPS.cycle")

  if (device == "pdf") {
    grDevices::pdf(outfile, width = width, height = height, onefile = TRUE)
  } else {
    grDevices::png(outfile, width = width * 150, height = height * 150, res = 150)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  gg <- ggplot2::ggplot()
  gg <- ggimage::ggbackground(gg, img)
  print(gg)
  for (i in seq_len(nrow(cfg$processes))) {
    if (sum(as.numeric(abun[i, 2:ncol(abun)])) <= 0) next
    mini <- res$fold_heat + ggplot2::ylim(cfg$processes$label[i])
    print(mini, vp = grid::viewport(width = 0.035 * res$group_numb,
                                    height = 0.05,
                                    x = cfg$processes$vx[i], y = cfg$processes$vy[i]))
  }
  invisible(outfile)
}

# ---- host taxonomy figure (pheatmap + grid title, device side effect) --
cnps_draw_host_plot <- function(host_tbl, title) {
  CNPS.cycle::host.ratio(host_tbl, title)
  invisible(NULL)
}

cnps_save_host_png <- function(host_tbl, title, outfile, width = 1600, height = 900) {
  grDevices::png(outfile, width = width, height = height, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  CNPS.cycle::host.ratio(host_tbl, title)
  invisible(outfile)
}

# ---- rebuilt ordination plots ----------------------------------------
cnps_ordination_plot <- function(plotdata, method = c("PCoA", "PCA", "NMDS"),
                                 style = c("points", "ellipse", "labels"),
                                 cycle_name = "") {
  method <- match.arg(method)
  style <- match.arg(style)
  pd <- as.data.frame(plotdata)
  xcol <- if (method == "NMDS") "NMDS1" else "PC1"
  ycol <- if (method == "NMDS") "NMDS2" else "PC2"
  xlab <- if (method == "NMDS") "NMDS1" else paste0(method, "1 (", x_percent(pd, 1), "%)")
  ylab <- if (method == "NMDS") "NMDS2" else paste0(method, "2 (", x_percent(pd, 2), "%)")

  p <- ggplot2::ggplot(pd, ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]],
                                        colour = Group, fill = Group)) +
    ggplot2::geom_point(size = 2.6, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = CNPS_PALETTE) +
    ggplot2::scale_fill_manual(values = CNPS_PALETTE) +
    theme_cnps() +
    ggplot2::labs(title = sprintf("%s - Function genes of %s", method, cycle_name),
                  x = xlab, y = ylab, colour = "Group", fill = "Group")
  if (style == "ellipse")
    p <- p + ggplot2::stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95,
                                   linetype = 2, linewidth = 0.4)
  if (style == "labels")
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = 3,
                                      show.legend = FALSE, seed = 1)
  p
}

# axis percent: rough eigenvalue share from coordinate spread, matching
# the floor() display style of the package plots
x_percent <- function(pd, k) {
  coords <- as.matrix(pd[, 2:3])
  if (ncol(coords) < 2 || nrow(coords) < 2) return(0)
  v <- tryCatch({
    cm <- scale(coords, center = TRUE, scale = FALSE)
    eig <- svd(cm)$d^2
    round(100 * eig[k] / sum(eig))
  }, error = function(e) 0)
  v
}

# ---- plotly conversions ----------------------------------------------
cnps_to_plotly <- function(p) {
  tryCatch(plotly::ggplotly(p, tooltip = c("all")), error = function(e) p)
}

# ---- shared ggplot theme ---------------------------------------------
theme_cnps <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}

# ---- write the full official-style Results tree (tables + figures) ----
cnps_write_plots <- function(res, root) {
  cy <- res$cycle; cfg <- res$cfg; G <- res$group_numb
  base <- file.path(root, cycle_name_dir(cy))

  hm_dir <- file.path(base, "Gene", "Heatmap")
  dir.create(hm_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(hm_dir, sprintf("%s_cycle_gene_abun_group.pdf", cy)),
                  res$heat_plot, width = cfg$heatmap_w + 0.6 * G, height = cfg$heatmap_h)
  utils::write.csv(res$heat_test, file.path(hm_dir, sprintf("Diff_%s_gene_test.csv", cy)),
                   row.names = FALSE)

  img_dir <- file.path(base, "Gene", "Cycle image")
  dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(img_dir, "Gene_fold_change.pdf"), res$fold_plot,
                  width = cfg$fold_w + 0.6 * G, height = cfg$fold_h)
  cnps_pathway_composite(res, file.path(img_dir, sprintf("%s_cyc_fold_change.pdf", cy)))

  host_dir <- file.path(base, "Host_relative_Group")
  dir.create(host_dir, recursive = TRUE, showWarnings = FALSE)
  for (code in names(res$host)) {
    d <- file.path(host_dir, code)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    for (j in seq_along(res$host[[code]])) {
      cnps_save_host_png(res$host[[code]][[j]], cfg$processes$title[cfg$processes$code == code],
                         file.path(d, sprintf("%s_%s.png", code,
                                              c("phylum","class","order","family","genus","species")[j])))
    }
  }
  # shared 0-100% legend (official style)
  if (length(res$host)) {
    df <- t(matrix(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1), nrow = 2))
    df1 <- matrix(c("0%","20%","40%","60%","80%","100%"), nrow = 1, ncol = 6)
    bk <- seq(0, 1, by = 0.01)
    try(pheatmap::pheatmap(df, fontsize = 30, cluster_rows = FALSE, fontface = "bold",
        cluster_cols = FALSE, cellwidth = 80, cellheight = 50, legend = FALSE,
        breaks = bk, show_rownames = FALSE, show_colnames = FALSE,
        color = grDevices::colorRampPalette(c("white", "Red"))(100),
        display_numbers = df1, number_color = "black", border_color = "black",
        filename = file.path(host_dir, "Legend_relative_Group.pdf"),
        width = 7, height = 2), silent = TRUE)
  }

  bd <- file.path(base, "Beta diversity")
  for (s in c("Distance", "PCA", "PCoA", "NMDS"))
    dir.create(file.path(bd, s), recursive = TRUE, showWarnings = FALSE)
  b <- res$beta
  if (!is.null(b$pcoa)) {
    utils::write.table(as.matrix(b$pcoa[[1]]),
           file.path(bd, "Distance", "distance_bray_curtis.txt"), sep = "\t")
    utils::write.table(b$pcoa[[2]], file.path(bd, "Distance", "diff_test.txt"), sep = "\t")
    utils::write.table(b$pcoa[[3]], file.path(bd, "PCoA", sprintf("%s_pcoa.txt", cy)), sep = "\t")
    for (st in c("group", "ellipse", "label")) {
      p <- cnps_ordination_plot(b$pcoa[[3]], "PCoA",
                                switch(st, group = "points", label = "labels", st), paste0(cfg$en, " cycling"))
      ggplot2::ggsave(file.path(bd, "PCoA", sprintf("%s_pcoa_%s.pdf", cy, st)), p,
                      width = 7.5, height = 5.4)
    }
  }
  if (!is.null(b$pca)) {
    utils::write.table(b$pca[[1]], file.path(bd, "PCA", sprintf("%s_pca.txt", cy)), sep = "\t")
    for (st in c("group", "ellipse", "label")) {
      p <- cnps_ordination_plot(b$pca[[1]], "PCA",
                                switch(st, group = "points", label = "labels", st), paste0(cfg$en, " cycling"))
      ggplot2::ggsave(file.path(bd, "PCA", sprintf("%s_pca_%s.pdf", cy, st)), p,
                      width = 7.5, height = 5.4)
    }
  }
  if (!is.null(b$nmds)) {
    utils::write.table(data.frame(x = b$nmds[[1]]),
                       file.path(bd, "NMDS", sprintf("%s_stress.txt", cy)), sep = "\t")
    utils::write.table(b$nmds[[2]], file.path(bd, "NMDS", sprintf("%s_nmds.txt", cy)), sep = "\t")
    for (st in c("group", "ellipse", "label")) {
      p <- cnps_ordination_plot(b$nmds[[2]], "NMDS",
                                switch(st, group = "points", label = "labels", st), paste0(cfg$en, " cycling"))
      ggplot2::ggsave(file.path(bd, "NMDS", sprintf("%s_nmds_%s.pdf", cy, st)), p,
                      width = 7.5, height = 5.4)
    }
  }
  invisible(base)
}
