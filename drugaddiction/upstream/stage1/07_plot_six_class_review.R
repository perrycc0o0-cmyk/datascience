#!/usr/bin/env Rscript
# Six-class display review after user-requested nearest-class matching.
# The original broad-class labels and matching uncertainty remain in source data.
# This renderer does not alter matrices, embeddings or matching assignments.

stage1_sixclass_plot_main <- function(root, config = file.path(root, "config", "processing_config.json")) {
  root <- normalizePath(root, mustWork = TRUE)
  cfg <- jsonlite::read_json(config, simplifyVector = FALSE)
  naming_path <- file.path(root, "config", "sample_naming.json")
  if (!file.exists(naming_path)) stop("Canonical sample-naming manifest is required.")
  naming <- jsonlite::read_json(naming_path, simplifyVector = FALSE)
  if (!is.null(cfg$stage_gate) && cfg$stage_gate != "global_celltype_review")
    stop("This renderer is restricted to the first user-review gate.")
  libs <- unlist(cfg$runtime$r_libraries, use.names = FALSE)
  if (length(libs)) .libPaths(unique(c(libs, .libPaths())))
  if (!is.null(cfg$runtime$fontconfig_file)) Sys.setenv(FONTCONFIG_FILE = cfg$runtime$fontconfig_file)
  req <- c("ggplot2", "patchwork", "ggrepel", "Matrix", "jsonlite")
  absent <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing R plotting dependencies: ", paste(absent, collapse = ", "))
  if (!isTRUE(capabilities("cairo"))) stop("Editable-text PDF requires Cairo support.")
  fc <- Sys.which("fc-match")
  if (!nzchar(fc)) stop("fc-match is required to verify Arial rather than substitute a font.")
  font <- system2(fc, c("-f", shQuote("%{family}\t%{file}"), "Arial"), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(font, "status")) || !length(font) ||
      !grepl("(^|[,;])\\s*Arial\\s*([,;]|$)", strsplit(font[[1L]], "\t", fixed = TRUE)[[1L]][1L], ignore.case = TRUE))
    stop("Fontconfig did not resolve Arial: ", paste(font, collapse = " "))
  if (requireNamespace("showtext", quietly = TRUE)) showtext::showtext_auto(FALSE)
  data_dir <- file.path(root, "05_six_class_review", "source_data")
  payload <- file.path(data_dir, "global_marker_review_payload.rds")
  if (!file.exists(file.path(root, "03_global", "complete.json")) || !file.exists(payload))
    stop("Global stage completion and marker payload are required before plotting.")
  cache <- readRDS(payload)
  d <- cache$coordinates; x <- cache$expression; st <- cache$marker_stats
  proposals <- cache$proposals; panels <- cache$marker_panels
  if (!identical(colnames(x), d$cell_id) || anyDuplicated(d$cell_id)) stop("Marker payload/coordinate alignment is invalid.")
  six_classes <- c("Excitatory", "Inhibitory", "Oligodendrocytes", "OPC", "Astrocytes", "Microglia")
  if (!setequal(unique(as.character(d$proposed_major_class)), six_classes) ||
      !all(as.character(proposals$proposed_major_class) %in% six_classes))
    stop("The matched payload must contain exactly the six requested display classes.")
  xy_before <- d[, c("cell_id", "UMAP_1", "UMAP_2")]
  out <- file.path(root, "05_six_class_review")
  if (file.exists(file.path(out, "complete.json"))) stop("Review figures already complete; use an explicit distinct plot revision directory.")
  for (z in c("01_UMAP", "02_DotPlot", "03_marker_UMAP", "source_data", "provenance"))
    dir.create(file.path(out, z), recursive = TRUE, showWarnings = FALSE)
  preset <- list(base = 8, axis = 7.5, title = 9, legend = 6.5, tag = 10,
    label_mm = 2.6, vector_point_mm = .09, raster_pointsize = 1.2, raster_pixels = 1200L,
    dot_radius = 4.0, dpi = 600L)
  palette <- c(Excitatory = "#3B6FBF", Inhibitory = "#E85D5D", Oligodendrocytes = "#B5895A",
    OPC = "#9A6FD1", Astrocytes = "#4DB894", Microglia = "#F0986D")
  cluster_ids <- unique(as.character(d$cluster))
  cluster_ids <- cluster_ids[order(as.numeric(cluster_ids))]
  cluster_colors <- stats::setNames(grDevices::colorRampPalette(c("#1E3A8A", "#DC2626", "#059669",
    "#7C2D12", "#6D28D9", "#0F766E", "#BE185D", "#374151", "#D97706", "#0369A1", "#4F46E5", "#A16207"))(length(cluster_ids)), cluster_ids)
  theme_review <- function(umap = FALSE) {
    p <- ggplot2::theme_classic(base_family = "Arial", base_size = preset$base) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = preset$title, hjust = .5),
        axis.text = ggplot2::element_text(size = preset$axis, color = "black"),
        axis.title = ggplot2::element_text(size = preset$axis),
        axis.line = ggplot2::element_line(linewidth = .3, color = "black"),
        axis.ticks = ggplot2::element_line(linewidth = .3, color = "black"),
        legend.text = ggplot2::element_text(size = preset$legend),
        legend.title = ggplot2::element_text(size = preset$legend, face = "bold"),
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.key = ggplot2::element_rect(fill = "white", color = NA),
        plot.margin = ggplot2::margin(5, 5, 5, 5), plot.tag = ggplot2::element_blank())
    if (umap) p <- p + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
    p
  }
  manifest <- list()
  save_plot <- function(p, stem, width_mm, height_mm) {
    stem <- file.path(out, stem)
    grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
      family = "Arial", bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(paste0(stem, ".png"), width = width_mm, height = height_mm, units = "mm",
        res = preset$dpi, background = "white")
    } else grDevices::png(paste0(stem, ".png"), width = width_mm, height = height_mm, units = "mm",
      res = preset$dpi, type = "cairo", family = "Arial", bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
    manifest[[length(manifest) + 1L]] <<- data.frame(stem = stem, width_mm = width_mm,
      height_mm = height_mm, png_dpi = preset$dpi, pdf_editable_text = TRUE)
  }
  points <- function() {
    if (requireNamespace("scattermore", quietly = TRUE))
      scattermore::geom_scattermore(pointsize = preset$raster_pointsize,
        pixels = c(preset$raster_pixels, preset$raster_pixels))
    else ggplot2::geom_point(size = preset$vector_point_mm, stroke = 0)
  }
  arrows <- function(p) {
    xr <- range(d$UMAP_1); yr <- range(d$UMAP_2); dx <- diff(xr); dy <- diff(yr)
    x0 <- xr[1L] - .05 * dx; y0 <- yr[1L] - .05 * dy
    p + ggplot2::annotate("segment", x = x0, xend = x0 + .15 * dx, y = y0, yend = y0,
      linewidth = .35, color = "black", arrow = grid::arrow(length = grid::unit(1.7, "mm"), type = "closed")) +
      ggplot2::annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + .15 * dy,
      linewidth = .35, color = "black", arrow = grid::arrow(length = grid::unit(1.7, "mm"), type = "closed")) +
      ggplot2::annotate("text", x = x0 + .03 * dx, y = y0 - .03 * dy, label = "UMAP-1",
      hjust = 0, vjust = 1, family = "Arial", size = 2) +
      ggplot2::annotate("text", x = x0 - .03 * dx, y = y0 + .12 * dy, label = "UMAP-2",
      angle = 90, hjust = 1, vjust = 0, family = "Arial", size = 2) + ggplot2::coord_fixed(clip = "off")
  }
  umap <- function(group, colors, title, label = TRUE) {
    dd <- d; dd$group <- factor(as.character(dd[[group]]), levels = names(colors))
    if (anyNA(dd$group)) stop("Unmapped group color in UMAP: ", group)
    dd <- dd[order((seq_len(nrow(dd)) * 104729) %% 1000003), , drop = FALSE]
    p <- ggplot2::ggplot(dd, ggplot2::aes(UMAP_1, UMAP_2, color = group)) + points() +
      ggplot2::scale_color_manual(values = colors, drop = FALSE) + ggplot2::labs(title = title, color = NULL) + theme_review(TRUE)
    if (label) {
      lab <- stats::aggregate(dd[, c("UMAP_1", "UMAP_2")], list(group = dd$group), stats::median)
      p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = group), family = "Arial",
        color = "black", size = preset$label_mm, seed = 42L, max.overlaps = Inf,
        min.segment.length = 0, segment.size = .2, show.legend = FALSE) + ggplot2::theme(legend.position = "none")
    } else p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1.7)))
    arrows(p)
  }
  dotplot <- function(tab, genes, levels, title) {
    dd <- tab[tab$gene %in% genes, , drop = FALSE]
    dd$gene <- factor(dd$gene, levels = genes); dd$cluster <- factor(dd$cluster, levels = levels)
    dd$pct_display <- ifelse(dd$pct_expressed < 5, NA_real_, dd$pct_expressed)
    ggplot2::ggplot(dd, ggplot2::aes(gene, cluster, size = pct_display, color = avg_expression_scaled)) +
      ggplot2::geom_point(na.rm = TRUE) +
      ggplot2::scale_color_gradientn(colors = c("grey80", "#DC143C"), limits = c(-2.5, 2.5),
        name = "Average Expression") +
      ggplot2::scale_radius(range = c(0, preset$dot_radius), limits = c(0, 100),
        breaks = c(25, 50, 75), name = "Percent Expressed") +
      ggplot2::labs(title = title, x = NULL, y = NULL) + theme_review() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = .5, face = "italic"),
        axis.text.y = ggplot2::element_text(face = "bold"), legend.position = "bottom", legend.box = "horizontal") +
      ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top", direction = "horizontal",
        barwidth = grid::unit(28, "mm"), barheight = grid::unit(2.2, "mm")),
        size = ggplot2::guide_legend(title.position = "top", nrow = 1))
  }
  combine <- function(plots, ncol) patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(tag_levels = "A", theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = preset$tag))) &
    ggplot2::theme(plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = preset$tag))
  used_palette <- palette[names(palette) %in% d$proposed_major_class]
  display_names <- c(Excitatory = "Excitatory_Neurons", Inhibitory = "Inhibitory_Neurons",
    Oligodendrocytes = "Oligodendrocytes", OPC = "OPC", Astrocytes = "Astrocytes",
    Microglia = "Microglia")
  d$major_display <- unname(display_names[as.character(d$proposed_major_class)])
  display_colors <- stats::setNames(unname(used_palette), display_names[names(used_palette)])
  pA <- umap("cluster", cluster_colors, "Unsupervised clusters")
  pB <- umap("major_display", display_colors, "Nearest-class grouping")
  save_plot(pA, "01_UMAP/A_global_clusters", 110, 100)
  save_plot(pB, "01_UMAP/B_proposed_cell_types", 110, 100)
  sample_ids <- vapply(naming$sample_map, `[[`, "", "sample_id")
  if (!setequal(unique(as.character(d$sample_id)), sample_ids)) stop("Plot sample IDs are not canonical paper sample names.")
  sample_colors <- stats::setNames(c("#377EB8", "#66A9D3", "#A6CEE3", "#984EA3", "#C08FD0", "#DCC6E7", "#E41A1C", "#F4777F", "#FBB4AE"), sample_ids)
  pE <- umap("sample_id", sample_colors, "Library", FALSE)
  pF <- umap("condition_label", c(WD = "#888888", Ext = "#4BBAD3", `Ret-Ext` = "#EE4B35"), "Condition", FALSE)
  save_plot(pE, "01_UMAP/E_sample_diagnostic", 140, 110)
  save_plot(pF, "01_UMAP/F_condition_diagnostic", 125, 110)
  genes <- unique(unlist(panels, use.names = FALSE)); genes <- genes[genes %in% rownames(x)]
  labels <- stats::setNames(paste0(proposals$cluster, " | ", proposals$proposed_major_class), proposals$cluster)
  cluster_stats <- st; cluster_stats$cluster <- labels[as.character(cluster_stats$cluster)]
  cluster_levels <- rev(unname(labels[cluster_ids]))
  pC <- dotplot(cluster_stats, genes, cluster_levels, "Canonical markers across new clusters")
  save_plot(pC, "02_DotPlot/C_canonical_markers_by_cluster", 183, max(100, 35 + length(cluster_ids) * 3.6))
  # Recompute gene-wise color scaling within this displayed grouping, as in Seurat DotPlot.
  source(file.path(root, "scripts", "02_global_celltype.R"), local = TRUE)
  group_levels <- names(used_palette)
  group_stats <- do.call(rbind, lapply(rownames(x), function(g) stage1_gene_stats(
    as.numeric(x[g, ]), g, as.character(d$proposed_major_class), group_levels)))
  display_stats <- group_stats
  display_stats$cluster <- unname(display_names[as.character(display_stats$cluster)])
  pD <- dotplot(display_stats, genes, unname(display_names[group_levels]), "Markers by nearest-class grouping")
  save_plot(pD, "02_DotPlot/D_canonical_markers_by_celltype", 183, 80)
  audit_dot_points <- function(plot, plot_id, filename) {
    built <- ggplot2::ggplot_build(plot)
    source <- plot$data
    actual <- built$data[[1L]]
    if (nrow(actual) != nrow(source)) stop("DotPlot audit cannot align source rows and point layer.")
    mapped_x <- as.integer(source$gene); mapped_y <- as.integer(source$cluster)
    if (!identical(as.numeric(actual$x), as.numeric(mapped_x)) ||
        !identical(as.numeric(actual$y), as.numeric(mapped_y))) stop("DotPlot point layer order does not match source keys.")
    get_actual <- function(nm) if (nm %in% names(actual)) actual[[nm]] else rep(NA, nrow(actual))
    audit <- data.frame(plot_id = plot_id, layer_index = 1L, geom = "GeomPoint",
      source_key = paste(as.character(source$gene), as.character(source$cluster), sep = "::"),
      gene = as.character(source$gene), group = as.character(source$cluster),
      pct_source = source$pct_expressed, scaled_source = source$avg_expression_scaled,
      pct_display_source = source$pct_display,
      plotted_x = get_actual("x"), plotted_y = get_actual("y"),
      plotted_size = get_actual("size"), plotted_colour = get_actual("colour"),
      plotted_alpha = get_actual("alpha"), plotted_shape = get_actual("shape"),
      plotted_stroke = get_actual("stroke"), stringsAsFactors = FALSE)
    write.csv(audit, file.path(out, "source_data", filename), row.names = FALSE)
  }
  audit_dot_points(pC, "C", "C_actual_DotPlot_point_layer.csv")
  audit_dot_points(pD, "D", "D_actual_DotPlot_point_layer.csv")
  extra <- setdiff(rownames(x), genes)
  if (length(extra)) {
    pExtra <- dotplot(cluster_stats, extra, cluster_levels, "Supplemental marker review")
    save_plot(pExtra, "02_DotPlot/G_supplemental_marker_review", 210, max(100, 35 + length(cluster_ids) * 3.6))
  }
  save_plot(combine(list(pA, pB, pD), 2), "Celltype_marker_review_combined", 220, 205)
  save_plot(combine(list(pA, pB), 2), "01_UMAP/Combined_global_UMAPs", 220, 100)
  save_plot(combine(list(pE, pF), 2), "01_UMAP/Combined_sample_condition_diagnostics", 240, 110)
  feature_genes <- intersect(c("Slc17a7", "Satb2", "Gad1", "Gad2", "Mog", "Mbp",
    "Cacng4", "Pdgfra", "Slco1c1", "Gfap", "P2ry12", "Tmem119"), rownames(x))
  features <- list()
  for (i in seq_along(feature_genes)) {
    g <- feature_genes[[i]]; dd <- d; dd$value <- as.numeric(x[g, ])
    dd <- dd[order(dd$value), , drop = FALSE]
    p <- ggplot2::ggplot(dd, ggplot2::aes(UMAP_1, UMAP_2, color = value)) + points() +
      ggplot2::scale_color_gradientn(colors = c("grey60", "purple"), name = "Expression") +
      ggplot2::labs(title = g) + theme_review(TRUE) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "italic"), legend.position = "bottom") +
      ggplot2::guides(color = ggplot2::guide_colorbar(direction = "horizontal", title.position = "top",
        barwidth = grid::unit(18, "mm"), barheight = grid::unit(2, "mm"))) + ggplot2::coord_fixed()
    features[[i]] <- p
    save_plot(p, paste0("03_marker_UMAP/", LETTERS[[i]], "_", g), 89, 85)
  }
  if (length(features)) save_plot(combine(features, 4), "03_marker_UMAP/Combined_canonical_marker_UMAP", 240, ceiling(length(features) / 4) * 70)
  write.csv(cluster_stats, file.path(out, "source_data", "DotPlot_canonical_and_review_by_cluster.csv"), row.names = FALSE)
  write.csv(group_stats, file.path(out, "source_data", "DotPlot_by_proposed_celltype.csv"), row.names = FALSE)
  write.csv(do.call(rbind, manifest), file.path(out, "provenance", "figure_export_manifest.csv"), row.names = FALSE)
  if (!identical(xy_before, d[, c("cell_id", "UMAP_1", "UMAP_2")])) stop("Plotting changed source coordinates.")
  status <- list(status = "complete", completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stage_gate = "global_celltype_review", n_cells = nrow(d), n_clusters = length(cluster_ids),
    n_proposed_classes = length(used_palette), annotation_approved = FALSE,
    waiting_for = "User marker review of six-class nearest-match display before any Exc/Inh reclustering",
    subclustering_started = FALSE, downstream_started = FALSE,
    font = font, backend = "R", single_panels_untagged = TRUE, combined_labels = "uppercase",
    marker_values_and_coordinates_changed = FALSE,
    DotPlot_encoding = "percent detected (below 5% not plotted); radius proportional to detection; gene-wise z of log1p(mean linear normalized expression), clipped only at display limits -2.5 to 2.5",
    annotation = "User-requested nearest-class matching into six display classes; pending marker review; source labels and matching uncertainty retained, not confirmed biological lineage")
  jsonlite::write_json(status, file.path(out, "complete.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)
  cat("Six-class matched display ready. STOP: await user marker confirmation.\n")
}

if (sys.nframe() == 0L && Sys.getenv("STAGE1_SOURCE_ONLY") != "1") {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) < 2L || a[[1L]] != "--root") stop("Usage: Rscript 07_plot_six_class_review.R --root NEW_ROOT [--config CONFIG]")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required to read configuration.")
  cp <- if ("--config" %in% a) a[[match("--config", a) + 1L]] else file.path(a[[2L]], "config", "processing_config.json")
  stage1_sixclass_plot_main(a[[2L]], cp)
}
