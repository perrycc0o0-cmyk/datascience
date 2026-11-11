#!/usr/bin/env Rscript
# Render saved neuronal-subcluster annotation review only; never fit or relabel.

nr_flat <- function(x) as.character(unlist(x, use.names = FALSE))
nr_named <- function(x) unlist(x, use.names = TRUE)
nr_columns <- function(x, needed, label) {
  absent <- setdiff(needed, names(x))
  if (length(absent)) stop(label, " lacks fields: ", paste(absent, collapse = ", "))
}
nr_sort_ids <- function(x) {
  x <- unique(as.character(x))
  if (all(grepl("^[0-9]+$", x))) x[order(as.numeric(x))] else sort(x)
}
nr_warning_present <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  z <- trimws(tolower(as.character(x)))
  !is.na(z) & nzchar(z) & !z %in% c("false", "0", "none", "no", "not_flagged")
}

neuron_review_plot_main <- function(root, runtime_config = NULL, output_name = "04_marker_review") {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
  if (!is.null(runtime_config)) {
    rt <- jsonlite::read_json(runtime_config, simplifyVector = FALSE)
    if (!is.null(rt$runtime$r_libraries)) .libPaths(unique(c(nr_flat(rt$runtime$r_libraries), .libPaths())))
    if (!is.null(rt$runtime$fontconfig_file)) Sys.setenv(FONTCONFIG_FILE = rt$runtime$fontconfig_file)
  }
  needed <- c("ggplot2", "ggrepel", "patchwork", "ragg", "Matrix", "scales", "jsonlite")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing R plotting packages: ", paste(absent, collapse = ", "))
  if (!isTRUE(capabilities("cairo"))) stop("Cairo is required for editable PDF text")
  fc <- Sys.which("fc-match")
  if (!nzchar(fc)) stop("fc-match is required to verify Arial")
  font <- system2(fc, c("-f", shQuote("%{family}\t%{file}"), "Arial"), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(font, "status")) || !length(font) ||
      !grepl("(^|[,;])\\s*Arial\\s*([,;]|$)", strsplit(font[[1L]], "\t", fixed = TRUE)[[1L]][1L], ignore.case = TRUE))
    stop("Arial was not resolved: ", paste(font, collapse = " "))
  if (requireNamespace("showtext", quietly = TRUE)) showtext::showtext_auto(FALSE)
  config <- file.path(root, "config", "subtype_markers_and_style.json")
  cfg <- jsonlite::read_json(config, simplifyVector = FALSE)
  style <- cfg$style
  if (!identical(style$backend, "R")) stop("This renderer is R-only")
  if (length(output_name) != 1L || !grepl("^[A-Za-z0-9_]+$", output_name)) stop("Invalid output directory")
  out <- file.path(root, output_name)
  existing <- if (dir.exists(out)) list.files(out, recursive = TRUE, full.names = TRUE, pattern = "\\.(png|pdf)$") else character()
  if (file.exists(file.path(out, "complete.json")) || length(existing))
    stop("Review exports already exist; select an explicit new output directory")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  nums <- function(x) as.numeric(unlist(x, use.names = FALSE))
  excluded <- nr_flat(cfg$display_excluded_genes)
  fs <- style$font_sizes
  palette <- nr_flat(style$colors)
  knots <- nums(style$color_knots)
  color_limits <- nums(style$color_limits)
  palette_values <- (knots - color_limits[1L]) / diff(color_limits)
  size_limits <- nums(style$size_limits)
  cutoff <- as.numeric(style$hide_below_detection_percent)
  raster <- requireNamespace("scattermore", quietly = TRUE)
  base_theme <- function(umap = FALSE) {
    p <- ggplot2::theme_classic(base_family = "Arial", base_size = as.numeric(fs$base)) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = as.numeric(fs$title), hjust = .5),
        plot.subtitle = ggplot2::element_text(size = as.numeric(fs$legend), hjust = .5),
        axis.text = ggplot2::element_text(size = as.numeric(fs$axis), color = "black"),
        axis.title = ggplot2::element_text(size = as.numeric(fs$axis)),
        axis.line = ggplot2::element_line(linewidth = .3, color = "black"),
        axis.ticks = ggplot2::element_line(linewidth = .3, color = "black"),
        legend.text = ggplot2::element_text(size = as.numeric(fs$legend)),
        legend.title = ggplot2::element_text(size = as.numeric(fs$legend), face = "bold"),
        panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.key = ggplot2::element_rect(fill = "white", color = NA),
        plot.margin = ggplot2::margin(5, 5, 5, 5), plot.tag = ggplot2::element_blank())
    if (umap) p <- p + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
    p
  }
  point_layer <- function() {
    if (raster) scattermore::geom_scattermore(pointsize = as.numeric(style$umap_raster_pointsize),
      pixels = rep(as.integer(style$umap_raster_pixels), 2L))
    else ggplot2::geom_point(size = as.numeric(style$umap_vector_point_mm), stroke = 0)
  }
  tag_theme <- ggplot2::theme(plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = as.numeric(fs$tag)))
  combine <- function(plots, ncol) {
    patchwork::wrap_plots(plots, ncol = ncol) +
      patchwork::plot_annotation(tag_levels = "A", theme = tag_theme) & tag_theme
  }
  all_status <- list()
  for (lineage in c("Exc", "Inh")) {
    spec <- cfg$lineages[[lineage]]
    input <- file.path(root, "03_annotation", paste0(lineage, "_marker_review_payload.rds"))
    if (!file.exists(input)) stop("Missing annotation review payload: ", lineage)
    protected <- c(input, config)
    hashes_before <- tools::md5sum(protected)
    cache <- readRDS(input)
    nr_columns(cache, c("expression", "coordinates", "marker_stats_by_cluster", "marker_stats_by_subtype",
      "proposals", "display_panels", "subtype_order", "colors"), paste(lineage, "payload"))
    d <- cache$coordinates
    nr_columns(d, c("cell_id", "sample_id", "condition_label", "subcluster", "subtype_proposed", "UMAP_1", "UMAP_2",
      "original_global_cluster", "initial_major_class", "forced_display_only", "vascular_review_flag", "lineage_identity_warning"),
      paste(lineage, "coordinates"))
    if (!identical(colnames(cache$expression), as.character(d$cell_id)) || anyDuplicated(d$cell_id) ||
        any(!is.finite(d$UMAP_1)) || any(!is.finite(d$UMAP_2))) stop("Invalid marker/coordinate alignment")
    if (!inherits(cache$expression, "sparseMatrix")) stop("Expected sparse normalized marker expression")
    if (any(!is.finite(cache$expression@x)) || any(cache$expression@x < 0)) stop("Invalid normalized marker values")
    original_d <- d
    xy_before <- d[, c("cell_id", "UMAP_1", "UMAP_2")]
    raw_ids <- nr_sort_ids(d$subcluster)
    prop <- cache$proposals
    nr_columns(prop, c("subcluster", "subtype_proposed", "n_cells", "annotation_status", "review_notes"), "proposals")
    if (anyDuplicated(prop$subcluster) || !setequal(as.character(prop$subcluster), raw_ids)) stop("Proposal/cluster mismatch")
    expected_labels <- as.character(prop$subtype_proposed[match(as.character(d$subcluster), as.character(prop$subcluster))])
    if (!identical(expected_labels, as.character(d$subtype_proposed))) stop("Cell/proposal labels disagree")
    desired_order <- nr_flat(spec$order)
    declared_order <- nr_flat(cache$subtype_order)
    known_present <- desired_order[desired_order %in% unique(d$subtype_proposed)]
    if (!identical(declared_order[declared_order %in% desired_order], desired_order[desired_order %in% declared_order]))
      stop("Payload changes approved subtype order")
    subtype_order <- unique(c(known_present, declared_order[declared_order %in% unique(d$subtype_proposed)]))
    subtype_order <- c(subtype_order, setdiff(unique(as.character(d$subtype_proposed)), subtype_order))
    colors <- nr_named(spec$colors)
    colors <- c(colors, nr_named(cfg$review_colors)[!names(nr_named(cfg$review_colors)) %in% names(colors)])
    payload_colors <- nr_named(cache$colors)
    shared_colors <- intersect(names(payload_colors), names(colors))
    if (length(shared_colors) && any(toupper(payload_colors[shared_colors]) != toupper(colors[shared_colors])))
      stop("Payload changes original subtype palette")
    if (length(setdiff(subtype_order, names(colors)))) stop("No declared color for subtype: ", paste(setdiff(subtype_order, names(colors)), collapse = ","))
    colors <- colors[subtype_order]
    expected_panels <- lapply(spec$display_panels, nr_flat)
    supplied_panels <- lapply(cache$display_panels, nr_flat)
    if (!identical(expected_panels, supplied_panels)) stop("Payload changes the selected B/A display panel")
    requested_display <- unique(nr_flat(spec$display_panels))
    requested_context <- unique(c(nr_flat(spec$context_markers),
      nr_flat(lapply(spec$annotation_panels, function(v) v$genes)), nr_flat(cfg$lineage_review_markers)))
    gene_names <- rownames(cache$expression)
    genes <- intersect(setdiff(requested_display, excluded), gene_names)
    context <- intersect(setdiff(requested_context, excluded), gene_names)
    if (!length(genes) || !length(context)) stop("No available display or context genes")
    dest <- file.path(out, lineage)
    for (p in c("01_UMAP", "02_DotPlot", "source_data", "provenance")) dir.create(file.path(dest, p), recursive = TRUE, showWarnings = FALSE)
    write.csv(data.frame(gene = unique(c(requested_display, requested_context)),
      available = unique(c(requested_display, requested_context)) %in% gene_names,
      excluded_from_display_only = unique(c(requested_display, requested_context)) %in% excluded),
      file.path(dest, "source_data", "marker_availability.csv"), row.names = FALSE)
    con <- gzfile(file.path(dest, "source_data", "coordinates_and_review_flags.csv.gz"), "wt")
    write.csv(original_d, con, row.names = FALSE); close(con)
    write.csv(prop, file.path(dest, "source_data", "annotation_proposals.csv"), row.names = FALSE)
    manifest <- list()
    save_plot <- function(p, stem, width, height) {
      path <- file.path(dest, stem)
      grDevices::cairo_pdf(paste0(path, ".pdf"), width = width / 25.4, height = height / 25.4, family = "Arial", bg = "white")
      tryCatch(print(p), finally = grDevices::dev.off())
      ragg::agg_png(paste0(path, ".png"), width = width, height = height, units = "mm", res = as.integer(style$png_dpi), background = "white")
      tryCatch(print(p), finally = grDevices::dev.off())
      manifest[[length(manifest) + 1L]] <<- data.frame(stem = stem, width_mm = width, height_mm = height,
        png_dpi = as.integer(style$png_dpi), pdf_editable_text = TRUE, stringsAsFactors = FALSE)
    }
    add_arrows <- function(p) {
      xr <- range(d$UMAP_1); yr <- range(d$UMAP_2); dx <- diff(xr); dy <- diff(yr)
      if (dx <= 0 || dy <= 0) stop("UMAP has a zero-range axis")
      x0 <- xr[1L] - .05 * dx; y0 <- yr[1L] - .05 * dy
      p + ggplot2::annotate("segment", x = x0, xend = x0 + .15 * dx, y = y0, yend = y0,
        linewidth = .35, color = "black", arrow = grid::arrow(length = grid::unit(1.7, "mm"), type = "closed")) +
        ggplot2::annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + .15 * dy,
          linewidth = .35, color = "black", arrow = grid::arrow(length = grid::unit(1.7, "mm"), type = "closed")) +
        ggplot2::annotate("text", x = x0 + .03 * dx, y = y0 - .03 * dy, label = "UMAP-1", hjust = 0, vjust = 1,
          family = "Arial", size = 2) +
        ggplot2::annotate("text", x = x0 - .03 * dx, y = y0 + .12 * dy, label = "UMAP-2", angle = 90, hjust = 1,
          vjust = 0, family = "Arial", size = 2) + ggplot2::coord_fixed(clip = "off")
    }
    umap <- function(group, palette, title = NULL, label = TRUE) {
      dd <- d; dd$plot_group <- factor(as.character(dd[[group]]), levels = names(palette))
      if (anyNA(dd$plot_group)) stop("Unmapped UMAP group: ", group)
      dd <- dd[order((seq_len(nrow(dd)) * 104729) %% 1000003), , drop = FALSE]
      p <- ggplot2::ggplot(dd, ggplot2::aes(UMAP_1, UMAP_2, color = plot_group)) + point_layer() +
        ggplot2::scale_color_manual(values = palette, drop = FALSE) + ggplot2::labs(title = title, color = NULL) + base_theme(TRUE)
      if (label) {
        lab <- stats::aggregate(dd[, c("UMAP_1", "UMAP_2")], list(plot_group = dd$plot_group), stats::median)
        p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = plot_group), color = "black", family = "Arial",
          size = as.numeric(style$umap_label_mm), seed = 42L, max.overlaps = Inf, min.segment.length = 0,
          segment.size = .2, show.legend = FALSE) + ggplot2::theme(legend.position = "none")
      } else p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1.7)))
      add_arrows(p)
    }
    ramp <- nr_flat(style$numeric_cluster_ramp)
    raw_colors <- stats::setNames(grDevices::colorRampPalette(ramp)(length(raw_ids)), raw_ids)
    pA <- umap("subcluster", raw_colors, paste(spec$name, "subclusters"))
    pB <- umap("subtype_proposed", colors, NULL)
    save_plot(pA, "01_UMAP/A_numeric_clusters", 110, 100)
    save_plot(pB, "01_UMAP/B_proposed_subtypes", 110, 100)
    sample_order <- c(paste0("WD", 1:3), paste0("Ext", 1:3), paste0("Ret-Ext", 1:3))
    if (!all(as.character(d$sample_id) %in% sample_order) ||
        !all(as.character(d$condition_label) %in% c("WD", "Ext", "Ret-Ext"))) stop("Canonical sample/condition names are invalid")
    sample_colors <- stats::setNames(c("#377EB8", "#66A9D3", "#A6CEE3", "#984EA3", "#C08FD0", "#DCC6E7", "#E41A1C", "#F4777F", "#FBB4AE"), sample_order)
    sample_colors <- sample_colors[sample_order %in% d$sample_id]
    condition_colors <- c(WD = "#888888", Ext = "#4BBAD3", `Ret-Ext` = "#EE4B35")
    condition_colors <- condition_colors[names(condition_colors) %in% d$condition_label]
    pE <- umap("sample_id", sample_colors, "Library", FALSE)
    pF <- umap("condition_label", condition_colors, "Condition", FALSE)
    parent_ids <- nr_sort_ids(d$original_global_cluster)
    parent_colors <- stats::setNames(grDevices::colorRampPalette(ramp)(length(parent_ids)), parent_ids)
    pG <- umap("original_global_cluster", parent_colors, "Parent global clusters")
    if (lineage == "Inh" && "26" %in% parent_ids)
      pG <- pG + ggplot2::labs(subtitle = "Global cluster 26 retained for lineage-identity review")
    save_plot(pE, "01_UMAP/E_sample", 140, 110)
    save_plot(pF, "01_UMAP/F_condition", 125, 110)
    save_plot(pG, "01_UMAP/G_parent_global_cluster", 110, 100)
    stats_c <- cache$marker_stats_by_cluster
    stats_d <- cache$marker_stats_by_subtype
    for (tt in list(stats_c, stats_d)) {
      nr_columns(tt, c("gene", "cluster", "n_cells", "pct_expressed", "avg_expression_scaled"), "marker statistics")
      if (any(!is.finite(tt$pct_expressed)) || any(tt$pct_expressed < 0 | tt$pct_expressed > 100) ||
          any(!is.finite(tt$avg_expression_scaled)) || any(tt$avg_expression_scaled < color_limits[1L] - 1e-10 | tt$avg_expression_scaled > color_limits[2L] + 1e-10))
        stop("Invalid saved marker statistic or color scale")
    }
    if (!setequal(as.character(unique(stats_c$cluster)), raw_ids) ||
        !setequal(as.character(unique(stats_d$cluster)), unique(as.character(d$subtype_proposed)))) stop("Marker statistics have different groups")
    raw_labels <- stats::setNames(paste0(prop$subcluster, " | ", prop$subtype_proposed), as.character(prop$subcluster))
    stats_c$source_group <- as.character(stats_c$cluster)
    stats_c$cluster <- unname(raw_labels[stats_c$source_group])
    stats_d$source_group <- as.character(stats_d$cluster)
    dotplot <- function(tab, gene_order, group_order, title, max_size, id) {
      dd <- tab[tab$gene %in% gene_order, , drop = FALSE]
      dd$gene <- factor(dd$gene, levels = gene_order); dd$cluster <- factor(dd$cluster, levels = group_order)
      if (anyNA(dd$gene) || anyNA(dd$cluster) || any(as.character(dd$gene) %in% excluded)) stop("DotPlot gene/group mismatch")
      dd$pct_display <- ifelse(dd$pct_expressed < cutoff, NA_real_, dd$pct_expressed)
      p <- ggplot2::ggplot(dd, ggplot2::aes(gene, cluster, size = pct_display, color = avg_expression_scaled)) +
        ggplot2::geom_point(shape = 16, stroke = 0, na.rm = TRUE) +
        ggplot2::scale_color_gradientn(colors = palette, values = palette_values, limits = color_limits,
          breaks = c(-2, -1, 0, 1, 2), name = "Average Expression") +
        ggplot2::scale_size_area(max_size = max_size, limits = size_limits, breaks = nums(style$size_breaks), name = "Percent Expressed") +
        ggplot2::labs(title = title, x = NULL, y = NULL) + base_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = .5, face = "italic"),
          axis.text.y = ggplot2::element_text(face = "bold"), legend.position = "bottom", legend.box = "horizontal") +
        ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top", direction = "horizontal",
          barwidth = grid::unit(28, "mm"), barheight = grid::unit(2.2, "mm")), size = ggplot2::guide_legend(title.position = "top", nrow = 1))
      actual <- ggplot2::ggplot_build(p)$data[[1L]]
      if (nrow(actual) != nrow(dd) ||
          !isTRUE(all.equal(as.numeric(actual$x), as.numeric(dd$gene))) ||
          !isTRUE(all.equal(as.numeric(actual$y), as.numeric(dd$cluster))) ||
          !isTRUE(all.equal(actual$size, max_size * sqrt(dd$pct_display / 100), tolerance = 1e-10)))
        stop("Actual point placement or area differs from the saved source")
      expected_color <- scales::gradient_n_pal(palette, values = palette_values)((dd$avg_expression_scaled - color_limits[1L]) / diff(color_limits))
      if (!identical(toupper(actual$colour), toupper(expected_color))) stop("Actual color mapping differs from true saved Z")
      source_points <- data.frame(plot_id = id, source_key = paste(as.character(dd$gene), dd$source_group, sep = "::"),
        gene = as.character(dd$gene), group = dd$source_group, group_label = as.character(dd$cluster),
        pct_expressed = dd$pct_expressed, avg_expression_scaled = dd$avg_expression_scaled,
        pct_display = dd$pct_display, plotted_x = as.numeric(actual$x), plotted_y = as.numeric(actual$y),
        plotted_size = actual$size, plotted_colour = actual$colour, plotted_shape = actual$shape,
        plotted_stroke = actual$stroke, stringsAsFactors = FALSE)
      write.csv(source_points, file.path(dest, "source_data", paste0(id, "_actual_points.csv")), row.names = FALSE)
      p
    }
    cluster_height <- max(110, 45 + length(raw_ids) * 5)
    pC <- dotplot(stats_c, genes, rev(unname(raw_labels[raw_ids])), "Selected markers across fresh subclusters",
      as.numeric(style$max_cluster_point_mm), "C")
    pD <- dotplot(stats_d, genes, rev(subtype_order), "Candidate subtype markers",
      as.numeric(style$max_group_point_mm), "D")
    pH <- dotplot(stats_c, context, rev(unname(raw_labels[raw_ids])), "Subtype context and lineage review",
      as.numeric(style$max_cluster_point_mm), "H")
    save_plot(pC, "02_DotPlot/C_by_cluster", 183, cluster_height)
    save_plot(pD, "02_DotPlot/D_by_subtype", 183, if (lineage == "Exc") 100 else 115)
    save_plot(pH, "02_DotPlot/H_context", max(240, 40 + length(context) * 4.5), cluster_height)
    save_plot(combine(list(pA, pB), 2), "01_UMAP/Combined_numeric_and_subtypes", 220, 100)
    save_plot(combine(list(pE, pF), 2), "01_UMAP/Combined_sample_condition", 240, 110)
    write.csv(cache$marker_stats_by_cluster, file.path(dest, "source_data", "saved_marker_stats_by_cluster.csv"), row.names = FALSE)
    write.csv(cache$marker_stats_by_subtype, file.path(dest, "source_data", "saved_marker_stats_by_subtype.csv"), row.names = FALSE)
    write.csv(do.call(rbind, manifest), file.path(dest, "provenance", "figure_export_manifest.csv"), row.names = FALSE)
    if (!identical(original_d, cache$coordinates) || !identical(xy_before, d[, c("cell_id", "UMAP_1", "UMAP_2")]) ||
        !identical(hashes_before, tools::md5sum(protected))) stop("Plotting changed source values, metadata or coordinates")
    status <- list(status = "complete", completed_UTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      lineage = lineage, n_cells = nrow(d), n_subclusters = length(raw_ids), candidate_subtypes = subtype_order,
      annotation_pending = TRUE, annotation_approved = FALSE, plotting_only = TRUE,
      source_values_labels_coordinates_unchanged = TRUE, source_md5 = as.list(hashes_before),
      font = font, backend = "R", png_dpi = as.integer(style$png_dpi), single_panels_untagged = TRUE, combined_labels = "uppercase",
      displayed_markers = genes, displayed_context = context, display_excluded_genes = excluded,
      missing_display_genes = setdiff(setdiff(requested_display, excluded), gene_names),
      missing_context_genes = setdiff(setdiff(requested_context, excluded), gene_names),
      point_size = "Point area proportional to unchanged detection percentage; below 5 percent not displayed",
      color = "Saved gene-wise scaled expression, true limits -2.5 to 2.5; common continuous palette",
      subclustering_recomputed_by_renderer = FALSE, UMAP_recomputed_by_renderer = FALSE,
      parent_global_cluster26_retained = lineage == "Inh" && "26" %in% parent_ids,
      retained_identity_warning_nuclei = sum(nr_warning_present(d$lineage_identity_warning)),
      downstream_started = FALSE, stop_before_downstream = TRUE,
      waiting_for = "User review of fresh Exc/Inh UMAP and marker DotPlots before downstream analysis")
    jsonlite::write_json(status, file.path(dest, "complete.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
    all_status[[lineage]] <- status
    cat(lineage, "marker review plots completed; stop before downstream.\n")
  }
  final <- list(status = "awaiting_user_subtype_marker_review", completed_UTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    lineages = all_status, annotation_approved = FALSE, downstream_started = FALSE, stop_before_downstream = TRUE,
    session_info = capture.output(utils::sessionInfo()))
  jsonlite::write_json(final, file.path(out, "complete.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
  invisible(final)
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) < 2L || a[1L] != "--root")
    stop("Usage: Rscript 03_plot_neuron_review.R --root RUN_ROOT [--runtime-config JSON] [--output-name REVISION]")
  opt <- function(flag, default = NULL) {
    j <- match(flag, a)
    if (is.na(j)) return(default)
    if (j == length(a)) stop("Missing value for ", flag)
    a[j + 1L]
  }
  neuron_review_plot_main(a[2L], runtime_config = opt("--runtime-config"), output_name = opt("--output-name", "04_marker_review"))
}
