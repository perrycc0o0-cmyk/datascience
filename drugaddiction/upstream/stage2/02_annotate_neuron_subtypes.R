#!/usr/bin/env Rscript
# Marker-supported proposals for the newly clustered neuronal objects.
# Reads existing RNA/data; does not refilter, normalize, integrate or recluster.

nr_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
nr_log <- function(...) { cat(nr_now(), ..., "\n"); flush.console() }
nr_chars <- function(x) as.character(unlist(x, use.names = FALSE))
nr_levels <- function(x) {
  v <- unique(as.character(x)); n <- suppressWarnings(as.numeric(v))
  v[order(is.na(n), n, v)]
}
nr_file_meta <- function(path) {
  z <- file.info(path)
  if (nrow(z) != 1L || !is.finite(z$size) || z$isdir) stop("Missing input file: ", path)
  list(path = normalizePath(path, mustWork = TRUE), bytes = unname(z$size),
       mtime_numeric = unname(as.numeric(z$mtime)), mtime = as.character(z$mtime))
}
nr_write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, null = "null", na = "null", dataframe = "rows", digits = 16)
nr_atomic_rds <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite an existing review payload: ", path)
  tmp <- tempfile(".incomplete-", tmpdir = dirname(path))
  saveRDS(x, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) stop("Cannot finalize review payload: ", tmp)
}
nr_require_columns <- function(x, columns, description) {
  absent <- setdiff(columns, names(x))
  if (length(absent)) stop(description, " missing columns: ", paste(absent, collapse = ", "))
}
nr_gene_stats <- function(value, gene, group, levels) {
  group <- as.character(group); levels <- as.character(levels)
  if (length(value) != length(group) || anyNA(group) || anyDuplicated(levels) ||
      !setequal(unique(group), levels) || any(!is.finite(value)) || any(value < 0))
    stop("Invalid normalized marker values or grouping: ", gene)
  n <- as.integer(table(factor(group, levels = levels)))
  ss <- rowsum(cbind(detected = as.numeric(value > 0), log_value = value,
                    linear_value = expm1(value)), group, reorder = FALSE)
  ss <- ss[match(levels, rownames(ss)), , drop = FALSE]
  if (any(!is.finite(ss))) stop("Nonfinite marker sums: ", gene)
  av <- ss[, "linear_value"] / n
  lv <- log1p(av); s <- stats::sd(lv)
  z <- if (is.finite(s) && s > 0) (lv - mean(lv)) / s else rep(0, length(lv))
  data.frame(gene = gene, cluster = levels, n_cells = n,
    n_cells_expressed = ss[, "detected"], pct_expressed = 100 * ss[, "detected"] / n,
    avg_normalized_data = ss[, "log_value"] / n, avg_linear_normalized = av,
    log1p_avg_linear_normalized = lv, gene_z_unclipped = z,
    avg_expression_scaled = pmax(-2.5, pmin(2.5, z)), row.names = NULL,
    stringsAsFactors = FALSE)
}
nr_all_gene_stats <- function(expression, group, levels) {
  if (!nrow(expression)) stop("No available marker genes for review")
  do.call(rbind, lapply(rownames(expression), function(g)
    nr_gene_stats(as.numeric(expression[g, ]), g, group, levels)))
}
nr_panel_spec <- function(x, order, lineage) {
  if (!is.list(x) || !length(x) || is.null(names(x)) || anyDuplicated(names(x)) ||
      !setequal(names(x), order)) stop("Invalid annotation_panels for ", lineage)
  x <- x[order]
  out <- lapply(names(x), function(label) {
    p <- x[[label]]
    if (!is.list(p) || is.null(p$genes) || is.null(p$use_for_identity_signature))
      stop("Each annotation panel needs genes and use_for_identity_signature: ", label)
    genes <- unique(nr_chars(p$genes))
    if (!length(genes) || anyNA(genes) || any(!nzchar(genes))) stop("Empty panel: ", label)
    enabled <- p$use_for_identity_signature
    if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled))
      stop("Identity-signature switch must be a scalar Boolean: ", label)
    list(genes = genes, use_for_identity_signature = enabled)
  })
  names(out) <- names(x)
  if (lineage == "Exc" && !all(vapply(out, `[[`, logical(1), "use_for_identity_signature")))
    stop("All six excitatory panels must participate as identity signatures")
  if (lineage == "Inh" && (!"GABA_Amb" %in% names(out) ||
      isTRUE(out$GABA_Amb$use_for_identity_signature) ||
      sum(vapply(out, `[[`, logical(1), "use_for_identity_signature")) != 7L))
    stop("Inhibitory annotation requires seven known signatures and disabled GABA_Amb")
  out
}
nr_panel_evidence <- function(st, clusters, panels) {
  available_genes <- unique(as.character(st$gene))
  records <- list()
  for (cl in clusters) {
    one <- st[as.character(st$cluster) == cl, , drop = FALSE]
    if (anyDuplicated(one$gene)) stop("Duplicate gene/cluster statistic")
    for (label in names(panels)) {
      genes <- panels[[label]]$genes
      available <- genes[genes %in% available_genes]
      missing <- setdiff(genes, available)
      k <- match(available, one$gene)
      if (anyNA(k)) stop("Incomplete gene-by-cluster statistic grid")
      z <- one$gene_z_unclipped[k]; pct <- one$pct_expressed[k]
      n_expected <- length(genes); n_available <- length(available)
      required_available <- min(2L, n_expected)
      required_detected <- min(2L, n_available)
      score <- if (n_available) mean(z) else NA_real_
      enough_available <- n_available >= required_available
      support <- n_available > 0L && sum(pct >= 5) >= required_detected && any(pct >= 10)
      qualifies <- isTRUE(panels[[label]]$use_for_identity_signature) &&
        enough_available && support && is.finite(score) && score > 0
      records[[length(records) + 1L]] <- data.frame(subcluster = cl, panel = label,
        use_for_identity_signature = panels[[label]]$use_for_identity_signature,
        n_expected_markers = n_expected, n_available_markers = n_available,
        required_available_markers = required_available,
        marker_genes = paste(genes, collapse = ";"), available_genes = paste(available, collapse = ";"),
        missing_genes = paste(missing, collapse = ";"), mean_gene_z = score,
        n_markers_detected_ge5pct = sum(pct >= 5),
        n_markers_detected_ge10pct = sum(pct >= 10),
        required_detected_ge5pct = required_detected,
        enough_available_markers = enough_available, detection_support = support,
        positive_mean_gene_z = is.finite(score) && score > 0,
        qualifies_for_identity_assignment = qualifies,
        stringsAsFactors = FALSE, row.names = NULL)
    }
  }
  do.call(rbind, records)
}
nr_propose <- function(st, clusters, panels, lineage, coordinates) {
  evidence <- nr_panel_evidence(st, clusters, panels)
  get <- function(cl, gene) {
    j <- as.character(st$cluster) == cl & st$gene == gene
    if (sum(j) != 1L) return(NA_real_)
    as.numeric(st$pct_expressed[j])
  }
  decisions <- vector("list", length(clusters))
  for (i in seq_along(clusters)) {
    cl <- clusters[[i]]
    ev <- evidence[evidence$subcluster == cl, , drop = FALSE]
    known <- ev[ev$use_for_identity_signature, , drop = FALSE]
    ranked <- known[is.finite(known$mean_gene_z), , drop = FALSE]
    ranked <- ranked[order(-ranked$mean_gene_z, match(ranked$panel, names(panels))), , drop = FALSE]
    qualified <- known[known$qualifies_for_identity_assignment, , drop = FALSE]
    qualified <- qualified[order(-qualified$mean_gene_z, match(qualified$panel, names(panels))), , drop = FALSE]
    candidate <- if (nrow(qualified)) qualified$panel[[1L]] else NA_character_
    runner <- if (nrow(qualified) >= 2L) qualified$panel[[2L]] else NA_character_
    candidate_score <- if (nrow(qualified)) qualified$mean_gene_z[[1L]] else NA_real_
    runner_score <- if (nrow(qualified) >= 2L) qualified$mean_gene_z[[2L]] else NA_real_
    margin <- if (nrow(qualified) >= 2L) candidate_score - runner_score else NA_real_
    low_margin <- is.finite(margin) && margin < .15
    pct <- vapply(c("Gad1", "Gad2", "Slc17a7", "Slc17a6"), function(g) get(cl, g), numeric(1))
    gate_evaluable <- all(is.finite(pct))
    missing_lineage <- names(pct)[!is.finite(pct)]
    inhibitory_dominant <- gate_evaluable && pct[["Gad1"]] >= 10 && pct[["Gad2"]] >= 10 &&
      pct[["Slc17a7"]] < 10 && pct[["Slc17a6"]] < 10
    amb_evaluable <- all(is.finite(pct[c("Gad1", "Gad2", "Slc17a6")]))
    amb_support <- is.finite(pct[["Slc17a6"]]) && pct[["Slc17a6"]] >= 10
    amb_without_gad <- amb_evaluable && amb_support && pct[["Gad1"]] < 10 && pct[["Gad2"]] < 10
    proposed <- if (is.na(candidate)) "Unknown" else candidate
    warning <- FALSE; warning_reason <- character(); amb_reason <- NA_character_
    notes <- c(if (!nrow(qualified)) "No known subtype panel passes the fixed marker-support gates",
      if (low_margin) "Top two qualifying panels differ by <0.15 mean gene-Z; review required",
      if (length(missing_lineage)) paste("Lineage gate not fully evaluable; missing:", paste(missing_lineage, collapse = ", ")))
    if (!gate_evaluable) {
      warning <- TRUE; warning_reason <- c(warning_reason, "lineage_gate_incomplete")
    }
    if (lineage == "Exc" && inhibitory_dominant) {
      proposed <- "Unknown"; warning <- TRUE
      warning_reason <- c(warning_reason, "Gad1_Gad2_supported_without_Slc17a7_Slc17a6_support")
      notes <- c(notes, "Both Gad1 and Gad2 are detected in >=10%, with both Slc17a7 and Slc17a6 <10%; excitatory subtype not assigned")
    }
    if (lineage == "Inh" && amb_without_gad) {
      proposed <- "GABA_Amb"; warning <- TRUE; amb_reason <- "lineage_conflict"
      warning_reason <- c(warning_reason, "Slc17a6_supported_without_Gad1_Gad2_support")
      notes <- c(notes, "Slc17a6 >=10% with Gad1 and Gad2 <10%; GABA_Amb is a review container, not evidence of inhibitory identity")
    } else if (lineage == "Inh" && !nrow(qualified) && amb_support) {
      proposed <- "GABA_Amb"; warning <- TRUE; amb_reason <- "unresolved_container"
      warning_reason <- c(warning_reason, "known_subtypes_unresolved_Slc17a6_supported")
      notes <- c(notes, "No qualifying known inhibitory subtype; Slc17a6 supports the legacy ambiguous container only")
    }
    cells <- as.character(coordinates$subcluster) == cl
    n26 <- sum(as.character(coordinates$original_global_cluster[cells]) == "26")
    n11 <- sum(as.character(coordinates$original_global_cluster[cells]) == "11")
    n_forced <- sum(coordinates$forced_display_only[cells])
    n_vascular <- sum(coordinates$vascular_review_flag[cells])
    if (n_forced) notes <- c(notes, paste0(n_forced, " nuclei retain a nearest-six-class parent assignment flag"))
    if (n26) notes <- c(notes, paste0(n26, " nuclei originate from global cluster 26; lineage remains subject to marker review"))
    if (n_vascular) notes <- c(notes, paste0(n_vascular, " nuclei retain inherited vascular-review flags"))
    decisions[[i]] <- data.frame(subcluster = cl, subtype_proposed = proposed,
      n_cells = sum(cells), annotation_status = "pending_user_marker_review",
      candidate = candidate, runner_up = runner, candidate_mean_gene_z = candidate_score,
      runner_up_mean_gene_z = runner_score, qualified_score_margin = margin,
      low_score_margin = low_margin, n_qualifying_panels = nrow(qualified),
      qualifying_panels = paste(qualified$panel, collapse = ";"),
      highest_scoring_known_panel = if (nrow(ranked)) ranked$panel[[1L]] else NA_character_,
      highest_known_panel_score = if (nrow(ranked)) ranked$mean_gene_z[[1L]] else NA_real_,
      lineage_identity_warning = warning,
      lineage_warning_reason = paste(unique(warning_reason), collapse = ";"),
      lineage_gate_evaluable = gate_evaluable, missing_lineage_markers = paste(missing_lineage, collapse = ";"),
      Gad1_pct = unname(pct[["Gad1"]]), Gad2_pct = unname(pct[["Gad2"]]),
      Slc17a7_pct = unname(pct[["Slc17a7"]]), Slc17a6_pct = unname(pct[["Slc17a6"]]),
      Amb_marker_supported = amb_support, Amb_container_reason = amb_reason,
      n_original_global_26 = n26, n_original_global_11 = n11,
      n_inherited_forced_display_only = n_forced, n_inherited_vascular_review = n_vascular,
      annotation_is_probability = FALSE, review_notes = paste(notes, collapse = "; "),
      stringsAsFactors = FALSE, row.names = NULL)
  }
  list(proposals = do.call(rbind, decisions), panel_evidence = evidence)
}
nr_gene_availability <- function(spec, requested, available, lineage_genes) {
  records <- list()
  add <- function(role, panel, genes, identity = FALSE) {
    genes <- unique(nr_chars(genes))
    if (!length(genes)) return(invisible(NULL))
    records[[length(records) + 1L]] <<- data.frame(role = role, panel = panel, gene = genes,
      available = genes %in% available, use_for_identity_signature = identity,
      stringsAsFactors = FALSE)
  }
  for (label in names(spec$annotation_panels)) add("annotation", label,
    spec$annotation_panels[[label]]$genes, spec$annotation_panels[[label]]$use_for_identity_signature)
  for (label in names(spec$display_panels)) add("display", label, spec$display_panels[[label]])
  for (label in names(spec$context_markers)) add("context", label, spec$context_markers[[label]])
  add("lineage_review", "transmitter_and_neuron_context", lineage_genes)
  out <- do.call(rbind, records)
  if (!setequal(unique(out$gene), requested)) stop("Availability roles do not cover all requested genes")
  out
}
nr_get_layer <- function(obj, assay, layer, features = NULL, cells = NULL) {
  a <- obj[[assay]]
  if (is.null(a)) stop("Missing assay: ", assay)
  if (inherits(a, "Assay5")) {
    found <- grep(paste0("^", layer, "($|\\.)"), SeuratObject::Layers(a, search = NA), value = TRUE)
    if (length(found) != 1L || found[[1L]] != layer) stop("Expected one joined ", assay, "/", layer, " layer")
    x <- SeuratObject::LayerData(a, layer = layer, features = features, cells = cells)
  } else {
    x <- SeuratObject::GetAssayData(a, slot = layer)
    if (!is.null(features)) x <- x[features, , drop = FALSE]
    if (!is.null(cells)) x <- x[, cells, drop = FALSE]
  }
  if (!inherits(x, "sparseMatrix") || any(!is.finite(x@x)) || any(x@x < 0))
    stop("Expected finite nonnegative sparse matrix: ", assay, "/", layer)
  if (layer == "counts" && any(x@x != round(x@x))) stop("Counts must remain integer UMIs")
  x
}
nr_write_csv <- function(x, path) {
  if (grepl("\\.gz$", path)) {
    con <- gzfile(path, "wt"); on.exit(close(con))
    utils::write.csv(x, con, row.names = FALSE, na = "NA")
  } else utils::write.csv(x, path, row.names = FALSE, na = "NA")
}

run_neuron_annotation <- function(root, config = file.path(root, "config", "stage2_config.json"),
                                  marker_config = file.path(root, "config", "subtype_markers_and_style.json")) {
  root <- normalizePath(root, mustWork = TRUE)
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required to read configuration")
  cfg <- jsonlite::read_json(config, simplifyVector = FALSE)
  if (!identical(cfg$stage_gate, "neuron_subtype_review")) stop("Incorrect stage gate")
  libs <- nr_chars(cfg$runtime$r_libraries)
  if (length(libs)) .libPaths(unique(c(libs, .libPaths())))
  needed <- c("Seurat", "SeuratObject", "Matrix", "jsonlite")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing existing runtime packages: ", paste(missing, collapse = ", "))
  marker_cfg <- jsonlite::read_json(marker_config, simplifyVector = FALSE)
  if (!all(c("Exc", "Inh") %in% names(marker_cfg$lineages))) stop("Marker configuration needs Exc and Inh")
  if (length(cfg$lineages) != 2L || !setequal(vapply(cfg$lineages, `[[`, "", "id"), c("Exc", "Inh")))
    stop("Stage 2 must contain exactly the two requested lineages")
  recluster_done <- file.path(root, "reclustering_complete.json")
  if (!file.exists(recluster_done) ||
      !identical(jsonlite::read_json(recluster_done)$status, "complete"))
    stop("Both new neuronal reclustering stages must be complete before marker review")
  out <- file.path(root, "03_annotation")
  if (dir.exists(out)) stop("Refusing to overwrite an existing annotation directory")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  status <- list(status = "running", started = nr_now(), stage = "input_validation",
    annotation_approved = FALSE, awaiting_user_marker_confirmation = TRUE,
    condition_DE_started = FALSE, downstream_started = FALSE)
  update <- function(stage, ...) {
    status$stage <<- stage; status$updated <<- nr_now()
    extra <- list(...); for (nm in names(extra)) status[[nm]] <<- extra[[nm]]
    nr_write_json(status, file.path(out, "status.json"))
  }
  small_inputs <- c(config, marker_config, recluster_done,
    file.path(root, "scripts", "02_annotate_neuron_subtypes.R"))
  small_inputs <- small_inputs[file.exists(small_inputs)]
  small_before <- tools::md5sum(small_inputs)
  all_records <- list()
  tryCatch({
    for (lineage_cfg in cfg$lineages) {
      lineage <- lineage_cfg$id; spec <- marker_cfg$lineages[[lineage]]
      lineage_dir <- file.path(root, lineage_cfg$output_directory)
      input <- file.path(lineage_dir, "neuron_reclustered.rds")
      input_meta <- file.path(lineage_dir, "source_data", "cell_annotations_UMAP.csv.gz")
      before_file <- nr_file_meta(input)
      if (!file.exists(input_meta)) stop("Missing original neuron coordinate/source manifest: ", input_meta)
      lineage_done <- file.path(lineage_dir, "complete.json")
      if (!file.exists(lineage_done) ||
          !identical(jsonlite::read_json(lineage_done)$status, "complete")) stop("Incomplete neuronal lineage: ", lineage)
      small_lineage <- c(input_meta, lineage_done,
        file.path(lineage_dir, "source_data", "parameters.json"))
      small_lineage <- small_lineage[file.exists(small_lineage)]
      small_lineage_before <- tools::md5sum(small_lineage)
      update(paste0(lineage, "_read_existing_normalized_markers"))
      nr_log("Read new neuronal object for marker-only annotation:", lineage)
      obj <- readRDS(input)
      if (!inherits(obj, "Seurat")) stop("Expected Seurat neuronal input")
      cells <- colnames(obj); md <- obj[[]]
      nr_require_columns(md, c("sample_id", "condition_label", "subcluster", "original_global_cluster",
        "initial_major_class", "forced_display_only", "vascular_review_flag"), "Neuronal metadata")
      if (length(cells) != as.integer(lineage_cfg$expected_nuclei) || anyDuplicated(cells) ||
          !identical(rownames(md), cells) || anyNA(md$subcluster) || anyNA(md$original_global_cluster) ||
          anyNA(md$forced_display_only) || anyNA(md$vascular_review_flag) ||
          !is.logical(md$forced_display_only) || !is.logical(md$vascular_review_flag))
        stop("Unexpected neuron membership/order or missing inherited review flags")
      if (anyNA(md$sample_id) || anyNA(md$condition_label) ||
          any(!grepl("^(WD|Ext|Ret-Ext)[1-3]$", md$sample_id))) stop("Invalid canonical sample labels")
      before_cluster <- as.character(md$subcluster)
      before_idents <- as.character(SeuratObject::Idents(obj))
      before_xy <- SeuratObject::Embeddings(obj, "umap.neuron")
      if (!identical(rownames(before_xy), cells) || ncol(before_xy) != 2L || any(!is.finite(before_xy)))
        stop("Invalid existing neuronal UMAP")
      before_counts <- nr_get_layer(obj, "RNA", "counts")
      before_observed <- nr_get_layer(obj, "observed", "counts")
      if (!identical(dimnames(before_counts), dimnames(before_observed)) ||
          !identical(colnames(before_counts), cells)) stop("Corrected and observed counts are not aligned")
      source <- utils::read.csv(input_meta, stringsAsFactors = FALSE, check.names = FALSE)
      nr_require_columns(source, c("cell_id", "subcluster", "UMAP_1", "UMAP_2"), "Original neuron coordinate manifest")
      if (!identical(as.character(source$cell_id), cells) ||
          !identical(as.character(source$subcluster), before_cluster) ||
          !isTRUE(all.equal(as.numeric(source$UMAP_1), as.numeric(before_xy[, 1L]), tolerance = 1e-12)) ||
          !isTRUE(all.equal(as.numeric(source$UMAP_2), as.numeric(before_xy[, 2L]), tolerance = 1e-12)))
        stop("Source coordinates and new neuronal object disagree")
      subtype_order <- nr_chars(spec$order)
      if (!length(subtype_order) || anyDuplicated(subtype_order) || "Unknown" %in% subtype_order)
        stop("Subtype order must contain the unique named manuscript subtypes, without Unknown")
      panels <- nr_panel_spec(spec$annotation_panels, subtype_order, lineage)
      spec$annotation_panels <- panels
      if (is.null(spec$display_panels) || is.null(names(spec$display_panels))) stop("Missing named display panels")
      display_panels <- lapply(spec$display_panels, nr_chars)
      if (!all(names(display_panels) %in% subtype_order)) stop("Display panels contain an unknown subtype")
      display_excluded <- nr_chars(marker_cfg$display_excluded_genes)
      if (any(unlist(display_panels, use.names = FALSE) %in% display_excluded))
        stop("A display-excluded gene remains in the requested DotPlot panel")
      context <- lapply(spec$context_markers, nr_chars)
      lineage_genes <- unique(c(nr_chars(marker_cfg$lineage_review_markers),
                                "Slc17a6", "Gad1", "Gad2", "Slc17a7"))
      requested <- unique(c(unlist(lapply(panels, `[[`, "genes"), use.names = FALSE),
        unlist(display_panels, use.names = FALSE), unlist(context, use.names = FALSE), lineage_genes))
      available <- requested[requested %in% rownames(obj[["RNA"]])]
      if (!length(available)) stop("No configured neuronal marker gene is available")
      # Keep only the small marker matrix; never densify the complete RNA matrix.
      expression <- nr_get_layer(obj, "RNA", "data", features = available, cells = cells)
      expression <- expression[available, cells, drop = FALSE]
      if (!inherits(expression, "sparseMatrix") || !identical(colnames(expression), cells) ||
          any(!is.finite(expression@x)) || any(expression@x < 0)) stop("Invalid sparse marker-only RNA/data")
      clusters <- nr_levels(before_cluster)
      st_cluster <- nr_all_gene_stats(expression, before_cluster, clusters)
      coordinates <- md
      coordinates$cell_id <- cells
      coordinates$UMAP_1 <- as.numeric(before_xy[, 1L]); coordinates$UMAP_2 <- as.numeric(before_xy[, 2L])
      coordinates$subcluster <- before_cluster
      coordinates <- coordinates[, c("cell_id", setdiff(names(coordinates), "cell_id")), drop = FALSE]
      rownames(coordinates) <- NULL
      update(paste0(lineage, "_marker_support_proposals"))
      ann <- nr_propose(st_cluster, clusters, panels, lineage, coordinates)
      map <- match(coordinates$subcluster, ann$proposals$subcluster)
      if (anyNA(map)) stop("A nucleus lacks a cluster-level proposal")
      coordinates$subtype_proposed <- ann$proposals$subtype_proposed[map]
      coordinates$lineage_identity_warning <- ann$proposals$lineage_identity_warning[map]
      coordinates$lineage_warning_reason <- ann$proposals$lineage_warning_reason[map]
      coordinates$subtype_annotation_status <- ann$proposals$annotation_status[map]
      coordinates$subtype_low_score_margin <- ann$proposals$low_score_margin[map]
      coordinates$Amb_container_reason <- ann$proposals$Amb_container_reason[map]
      if (nrow(coordinates) != length(cells) || anyNA(coordinates$subtype_proposed)) stop("Incomplete nucleus mapping")
      subtype_order <- c(subtype_order, "Unknown")
      present_subtypes <- subtype_order[subtype_order %in% coordinates$subtype_proposed]
      st_subtype <- nr_all_gene_stats(expression, coordinates$subtype_proposed, present_subtypes)
      availability <- nr_gene_availability(spec, requested, available, lineage_genes)
      colors <- nr_chars(spec$colors); names(colors) <- names(spec$colors)
      if (anyNA(colors) || !all(setdiff(subtype_order, "Unknown") %in% names(colors))) stop("Missing named subtype colors")
      if (!"Unknown" %in% names(colors)) colors <- c(colors, Unknown = "#D8D8D8")
      colors <- colors[subtype_order]
      source_counts <- as.data.frame(table(subcluster = factor(coordinates$subcluster, levels = clusters),
        original_global_cluster = factor(as.character(coordinates$original_global_cluster),
          levels = nr_levels(coordinates$original_global_cluster))), stringsAsFactors = FALSE)
      source_counts <- source_counts[source_counts$Freq > 0, , drop = FALSE]
      names(source_counts)[names(source_counts) == "Freq"] <- "n_cells"
      sample_counts <- as.data.frame(table(sample_id = coordinates$sample_id,
        subcluster = factor(coordinates$subcluster, levels = clusters)), stringsAsFactors = FALSE)
      names(sample_counts)[names(sample_counts) == "Freq"] <- "n_cells"
      subtype_counts <- as.data.frame(table(subtype = factor(coordinates$subtype_proposed,
        levels = subtype_order)), stringsAsFactors = FALSE)
      names(subtype_counts)[names(subtype_counts) == "Freq"] <- "n_cells"
      source26 <- as.character(md$original_global_cluster) == "26"
      if (lineage == "Inh" && !any(source26)) stop("Previously included global cluster 26 is absent; do not silently exclude it")
      provenance <- list(lineage = lineage, input_file = before_file,
        small_input_md5 = as.list(small_lineage_before), config_md5 = as.list(small_before),
        n_cells = ncol(obj), n_clusters = length(clusters), n_available_review_genes = length(available),
        expected_nuclei = lineage_cfg$expected_nuclei,
        n_original_global_26 = sum(source26),
        n_original_global_11 = sum(as.character(md$original_global_cluster) == "11"),
        n_inherited_forced_display_only = sum(md$forced_display_only),
        n_inherited_vascular_review = sum(md$vascular_review_flag),
        score_definition = "Mean of available curated marker gene-Z scores; gene-Z is calculated across new subclusters from log1p(mean(expm1(RNA/data))); unclipped Z for annotation, clipped +/-2.5 only for display",
        missing_gene_policy = "Missing genes are recorded, never imputed as zero; a fully missing panel has NA score",
        qualification = list(available_markers = "at least min(2, expected marker count)",
          detection_ge5pct = "at least min(2, available marker count)", detection_ge10pct = "at least one marker",
          mean_gene_z = "strictly greater than zero", low_margin = "difference between top two qualifying known panels <0.15"),
        Amb_definition = "Disabled identity signature; Slc17a6>=10% with both Gad1/Gad2<10% is a lineage-conflict container; otherwise Slc17a6>=10% permits an unresolved container only if no known panel qualifies",
        annotation_backend = "Rule-based cluster-level marker evidence; not atlas label transfer or calibrated probabilities",
        scientific_cautions = spec$scientific_cautions, annotation_approved = FALSE,
        awaiting_user_marker_confirmation = TRUE, cell_membership_unchanged = TRUE,
        counts_unchanged = TRUE, UMAP_unchanged = TRUE, cluster_assignments_unchanged = TRUE,
        inherited_flags_preserved = TRUE, normalization_recomputed = FALSE, latent_space_recomputed = FALSE,
        cells_removed = 0L, cells_reassigned_to_other_lineage = 0L,
        downstream_started = FALSE, condition_DE_started = FALSE, completed = nr_now())
      # This script never mutates the Seurat object. Verify the core data and source flags directly.
      if (!identical(cells, colnames(obj)) || !identical(before_cluster, as.character(obj$subcluster)) ||
          !identical(before_idents, as.character(SeuratObject::Idents(obj))) ||
          !identical(before_xy, SeuratObject::Embeddings(obj, "umap.neuron")) ||
          !identical(before_counts, nr_get_layer(obj, "RNA", "counts")) ||
          !identical(before_observed, nr_get_layer(obj, "observed", "counts")) ||
          !identical(coordinates$forced_display_only, md$forced_display_only) ||
          !identical(coordinates$vascular_review_flag, md$vascular_review_flag) ||
          !identical(coordinates$original_global_cluster, md$original_global_cluster) ||
          !identical(coordinates$initial_major_class, md$initial_major_class) ||
          !identical(coordinates$cell_id, cells) ||
          !identical(coordinates$UMAP_1, as.numeric(before_xy[, 1L])) ||
          !identical(coordinates$UMAP_2, as.numeric(before_xy[, 2L])))
        stop("Marker annotation changed protected cells, counts, flags, clusters or coordinates")
      if (!identical(before_file, nr_file_meta(input)) ||
          !identical(small_lineage_before, tools::md5sum(small_lineage))) stop("Neuronal inputs changed during review")
      payload <- list(expression = expression, coordinates = coordinates,
        marker_stats_by_cluster = st_cluster, marker_stats_by_subtype = st_subtype,
        proposals = ann$proposals, panel_evidence = ann$panel_evidence,
        marker_availability = availability, display_panels = display_panels,
        annotation_panels = panels, context_markers = context, subtype_order = subtype_order,
        colors = colors, lineage = lineage, provenance = provenance)
      nr_atomic_rds(payload, file.path(out, paste0(lineage, "_marker_review_payload.rds")))
      nr_write_csv(st_cluster, file.path(out, paste0(lineage, "_marker_stats_by_cluster.csv")))
      nr_write_csv(st_subtype, file.path(out, paste0(lineage, "_marker_stats_by_subtype.csv")))
      nr_write_csv(availability, file.path(out, paste0(lineage, "_marker_availability.csv")))
      nr_write_csv(availability[!availability$available, , drop = FALSE], file.path(out, paste0(lineage, "_missing_markers.csv")))
      nr_write_csv(ann$panel_evidence, file.path(out, paste0(lineage, "_panel_scores_and_support.csv")))
      nr_write_csv(ann$proposals, file.path(out, paste0(lineage, "_subtype_proposals.csv")))
      nr_write_csv(source_counts, file.path(out, paste0(lineage, "_global_sources_by_subcluster.csv")))
      nr_write_csv(sample_counts, file.path(out, paste0(lineage, "_sample_by_subcluster_counts.csv")))
      nr_write_csv(subtype_counts, file.path(out, paste0(lineage, "_subtype_counts.csv")))
      nr_write_csv(coordinates, file.path(out, paste0(lineage, "_cell_annotations_UMAP.csv.gz")))
      metadata_mapping <- list(cell_id = cells, subcluster = before_cluster,
        subtype_proposed = coordinates$subtype_proposed, metadata = coordinates,
        proposals = ann$proposals, provenance = provenance)
      nr_atomic_rds(metadata_mapping, file.path(out, paste0(lineage, "_subtype_metadata_mapping.rds")))
      nr_write_json(provenance, file.path(out, paste0(lineage, "_complete.json")))
      all_records[[lineage]] <- provenance
      nr_log("Marker review prepared:", lineage, ncol(obj), "nuclei; pending human confirmation")
      rm(obj, expression, before_counts, before_observed, before_xy, source, coordinates,
        payload, metadata_mapping, md, st_cluster, st_subtype, ann); gc(FALSE)
    }
    if (!identical(small_before, tools::md5sum(small_inputs))) stop("Configuration/script changed during review")
    status$status <- "complete"; status$stage <- "neuron_subtype_marker_review_ready"
    status$completed <- nr_now(); status$lineages <- all_records
    status$package_versions <- stats::setNames(lapply(needed, function(p) as.character(utils::packageVersion(p))), needed)
    status$large_Seurat_object_written <- FALSE; status$atlas_label_transfer <- FALSE
    nr_write_json(status, file.path(out, "status.json"))
    nr_write_json(status, file.path(out, "complete.json"))
    writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
    nr_log("Stop: review subtype UMAP/DotPlot before any IEG, condition DE, enrichment or other downstream analysis")
    invisible(status)
  }, error = function(e) {
    status$status <- "failed"; status$error <- conditionMessage(e); status$failed_at <- nr_now()
    nr_write_json(status, file.path(out, "status.json")); stop(e)
  })
}

neuron_annotation_unit_tests <- function() {
  tests <- 0L
  check <- function(x) { stopifnot(isTRUE(x)); tests <<- tests + 1L }
  # Pure synthetic helper tests: no biological matrices or Seurat functions are read/run.
  s <- nr_gene_stats(log1p(c(0, 2, 4, 0)), "g", c("a", "a", "b", "b"), c("a", "b"))
  check(identical(s$n_cells, c(2L, 2L)) && identical(s$pct_expressed, c(50, 50)))
  check(isTRUE(all.equal(s$avg_linear_normalized, c(1, 2))) &&
    isTRUE(all.equal(s$log1p_avg_linear_normalized, log1p(c(1, 2)))))
  check(isTRUE(all.equal(s$gene_z_unclipped, as.numeric(scale(log1p(c(1, 2)))))))
  panels <- list(A = list(genes = c("a1", "a2"), use_for_identity_signature = TRUE),
    B = list(genes = c("b1", "b2"), use_for_identity_signature = TRUE),
    Missing = list(genes = c("missing1", "missing2"), use_for_identity_signature = TRUE))
  st <- data.frame(gene = c("a1", "a2", "b1", "b2", "Gad1", "Gad2", "Slc17a7", "Slc17a6"),
    cluster = "0", gene_z_unclipped = c(1, 1, .9, .9, 0, 0, 0, 0),
    pct_expressed = c(20, 6, 20, 6, 20, 20, 20, 0))
  d <- data.frame(subcluster = "0", original_global_cluster = "26", forced_display_only = TRUE,
    vascular_review_flag = FALSE)
  e <- nr_panel_evidence(st, "0", panels)
  check(is.na(e$mean_gene_z[e$panel == "Missing"]) && !e$qualifies_for_identity_assignment[e$panel == "Missing"])
  a <- nr_propose(st, "0", panels, "Exc", d)
  check(a$proposals$subtype_proposed == "A" && a$proposals$runner_up == "B" && a$proposals$low_score_margin)
  check(a$proposals$n_original_global_26 == 1L && a$proposals$n_inherited_forced_display_only == 1L)
  st$pct_expressed[st$gene == "Slc17a7"] <- 0
  a <- nr_propose(st, "0", panels, "Exc", d)
  check(a$proposals$subtype_proposed == "Unknown" && a$proposals$lineage_identity_warning)
  st$pct_expressed[st$gene %in% c("Gad1", "Gad2")] <- 0
  st$pct_expressed[st$gene == "Slc17a6"] <- 20
  panels$GABA_Amb <- list(genes = "Slc17a6", use_for_identity_signature = FALSE)
  a <- nr_propose(st, "0", panels, "Inh", d)
  check(a$proposals$subtype_proposed == "GABA_Amb" && a$proposals$Amb_container_reason == "lineage_conflict")
  check(!a$panel_evidence$qualifies_for_identity_assignment[a$panel_evidence$panel == "GABA_Amb"])
  st$pct_expressed[st$gene %in% c("Gad1", "Gad2")] <- 20
  st$gene_z_unclipped[st$gene %in% c("a1", "a2", "b1", "b2")] <- -1
  a <- nr_propose(st, "0", panels, "Inh", d)
  check(a$proposals$subtype_proposed == "GABA_Amb" && a$proposals$Amb_container_reason == "unresolved_container")
  st$pct_expressed[st$gene == "Slc17a6"] <- 0
  a <- nr_propose(st, "0", panels, "Inh", d)
  check(a$proposals$subtype_proposed == "Unknown")
  a <- nr_propose(st[st$gene != "Gad1", ], "0", panels, "Exc", d)
  check(!a$proposals$lineage_gate_evaluable && a$proposals$lineage_identity_warning && is.na(a$proposals$Gad1_pct))
  one <- list(Single = list(genes = "a1", use_for_identity_signature = TRUE))
  st$gene_z_unclipped[st$gene == "a1"] <- 1
  e <- nr_panel_evidence(st, "0", one)
  check(e$qualifies_for_identity_assignment && e$required_available_markers == 1L && e$required_detected_ge5pct == 1L)
  incomplete <- list(Two = list(genes = c("a1", "unavailable"), use_for_identity_signature = TRUE))
  e <- nr_panel_evidence(st, "0", incomplete)
  check(!e$qualifies_for_identity_assignment && is.finite(e$mean_gene_z) && e$n_available_markers == 1L)
  cat("PURE_SYNTHETIC_ANNOTATION_HELPER_TESTS_OK ", tests, "\n", sep = "")
  invisible(tests)
}

if (sys.nframe() == 0L && Sys.getenv("NEURON_ANNOTATION_SOURCE_ONLY") != "1") {
  args <- commandArgs(trailingOnly = TRUE)
  if (identical(args, "--self-test")) {
    neuron_annotation_unit_tests()
  } else {
    if (length(args) < 2L || args[[1L]] != "--root")
      stop("Usage: Rscript --vanilla 02_annotate_neuron_subtypes.R --root RUN_ROOT [--config CONFIG] [--marker-config MARKER_CONFIG]; or --self-test")
    take <- function(flag, default) {
      j <- match(flag, args)
      if (is.na(j)) return(default)
      if (j == length(args)) stop("Missing path after ", flag)
      args[[j + 1L]]
    }
    run_neuron_annotation(args[[2L]], take("--config", file.path(args[[2L]], "config", "stage2_config.json")),
      take("--marker-config", file.path(args[[2L]], "config", "subtype_markers_and_style.json")))
  }
}
