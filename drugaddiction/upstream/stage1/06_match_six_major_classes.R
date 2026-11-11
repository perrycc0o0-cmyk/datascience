#!/usr/bin/env Rscript
# Assign previously unresolved global clusters to their nearest of six reference classes.
# This changes annotation metadata only; it does not normalize, integrate or recluster.

sixclass_main <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  cfg <- jsonlite::read_json(file.path(root, "config", "processing_config.json"))
  .libPaths(unique(c(unlist(cfg$runtime$r_libraries), .libPaths())))
  needed <- c("Seurat", "SeuratObject", "Matrix", "jsonlite")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing dependencies: ", paste(absent, collapse = ", "))
  if (!identical(cfg$stage_gate, "global_celltype_review")) stop("Unexpected review gate")
  classes <- c("Excitatory", "Inhibitory", "Oligodendrocytes", "OPC", "Astrocytes", "Microglia")
  source_dir <- file.path(root, "03_global", "source_data")
  out <- file.path(root, "05_six_class_review")
  if (file.exists(file.path(out, "assignment_complete.json")) ||
      file.exists(file.path(out, "global_celltype_six_class.rds")))
    stop("Six-class annotation already exists; do not overwrite it")
  dir.create(file.path(out, "source_data"), recursive = TRUE, showWarnings = FALSE)
  status_file <- file.path(out, "assignment_status.json")
  update <- function(stage, status = "running", ...) {
    jsonlite::write_json(c(list(status = status, stage = stage,
      updated_UTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      subclustering_started = FALSE, downstream_started = FALSE), list(...)),
      status_file, auto_unbox = TRUE, pretty = TRUE, digits = 16)
    cat(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), stage, "\n")
  }
  tryCatch({
    update("loading_completed_global_object")
    if (!file.exists(file.path(root, "03_global", "complete.json"))) stop("Original global stage is incomplete")
    obj <- readRDS(file.path(root, "03_global", "global_celltype_review.rds"))
    cache <- readRDS(file.path(source_dir, "global_marker_review_payload.rds"))
    d <- cache$coordinates
    proposals <- cache$proposals
    if (!identical(colnames(obj), d$cell_id)) stop("Object/payload cell order differs")
    if (!identical(as.character(obj$proposed_major_class), as.character(d$proposed_major_class)))
      stop("Original annotations differ between object and payload")
    cluster <- as.character(d$cluster)
    initial <- as.character(d$proposed_major_class)
    samples <- unique(as.character(d$sample_id))
    if (length(samples) != 9L || any(!grepl("^(WD|Ext|Ret-Ext)[1-3]$", samples)))
      stop("Canonical nine-sample naming is not preserved")
    unknown <- as.character(proposals$cluster[!proposals$proposed_major_class %in% classes])
    if (!length(unknown)) stop("There are no unresolved clusters to match")
    if (!all(unique(initial) %in% c(classes, "Mixed", "Unresolved"))) stop("Unexpected initial annotation")
    before_xy <- SeuratObject::Embeddings(obj, "umap.global")
    before_cells <- colnames(obj)
    before_clusters <- as.character(SeuratObject::Idents(obj))
    before_counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
    before_count_sum <- sum(before_counts)
    before_count_nnz <- length(before_counts@x)

    update("library_balanced_reference_profiles")
    features <- unique(c(SeuratObject::VariableFeatures(obj), unlist(cache$marker_panels, use.names = FALSE)))
    fmap <- read.csv(file.path(source_dir, "gene_feature_union_map.csv"), stringsAsFactors = FALSE)
    mt_features <- unique(fmap$seurat_feature[fmap$gene_id %in% unlist(cfg$qc$mitochondrial_gene_ids)])
    technical <- features %in% mt_features | grepl("^(Rpl|Rps|mt[-.])", features, ignore.case = TRUE)
    features <- features[!technical]
    data <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
    features <- intersect(features, rownames(data))
    if (length(features) < 100L) stop("Insufficient annotation features")
    x <- data[features, , drop = FALSE]
    if (!identical(colnames(x), d$cell_id)) stop("Expression/cell order differs")
    # LogNormalize values are inverted before averaging, then log1p-transformed.
    x@x <- expm1(x@x)
    if (any(!is.finite(x@x))) stop("Nonfinite expression in matching input")
    key <- paste(d$sample_id, cluster, sep = "::")
    group_keys <- unique(key)
    group_index <- match(key, group_keys)
    design <- Matrix::sparseMatrix(i = seq_along(group_index), j = group_index,
      x = 1, dims = c(length(group_index), length(group_keys)))
    totals <- as.matrix(x %*% design)
    group_n <- as.numeric(Matrix::colSums(design))
    first <- match(group_keys, key)
    group_meta <- data.frame(key = group_keys, sample_id = d$sample_id[first],
      cluster = cluster[first], initial_class = initial[first], n_cells = group_n,
      stringsAsFactors = FALSE)
    rownames(totals) <- features
    colnames(totals) <- group_keys
    # Reference classes contain only the originally unambiguous clusters.
    ref_linear <- list(); ref_n <- list()
    for (sid in samples) {
      mat <- matrix(NA_real_, nrow = length(features), ncol = length(classes),
        dimnames = list(features, classes))
      nn <- stats::setNames(numeric(length(classes)), classes)
      for (cl in classes) {
        j <- which(group_meta$sample_id == sid & group_meta$initial_class == cl)
        nn[cl] <- sum(group_meta$n_cells[j])
        if (nn[cl] >= 20L) mat[, cl] <- rowSums(totals[, j, drop = FALSE]) / nn[cl]
      }
      ref_linear[[sid]] <- mat; ref_n[[sid]] <- nn
    }
    refs <- do.call(cbind, ref_linear)
    eligible_features <- rowSums(is.finite(refs)) > 0L &
      apply(log1p(refs), 1L, function(v) stats::sd(v[is.finite(v)]) > 1e-8)
    if (sum(eligible_features) < 100L) stop("Insufficient variable reference genes")
    features <- features[eligible_features]
    totals <- totals[eligible_features, , drop = FALSE]
    ref_linear <- lapply(ref_linear, function(v) v[eligible_features, , drop = FALSE])
    write.csv(data.frame(gene = features), file.path(out, "source_data", "nearest_class_feature_genes.csv"), row.names = FALSE)
    write.csv(group_meta, file.path(out, "source_data", "cluster_by_library_counts.csv"), row.names = FALSE)

    update("nearest_class_spearman_matching")
    per_library <- list(); decisions <- list()
    for (q in unknown) {
      per_q <- list()
      for (sid in samples) {
        j <- which(group_meta$sample_id == sid & group_meta$cluster == q)
        if (length(j) != 1L || group_meta$n_cells[j] < 10L || any(ref_n[[sid]] < 20L)) next
        query <- log1p(totals[, j] / group_meta$n_cells[j])
        score <- as.numeric(stats::cor(query, log1p(ref_linear[[sid]]), method = "spearman"))
        if (length(score) != length(classes) || any(!is.finite(score))) stop("Invalid per-library correlation")
        per_q[[sid]] <- data.frame(cluster = q, sample_id = sid, query_n = group_meta$n_cells[j],
          reference_class = classes, spearman_rho = score, stringsAsFactors = FALSE)
      }
      if (length(per_q) < 2L) stop("Fewer than two eligible libraries for unresolved cluster ", q)
      tab <- do.call(rbind, per_q)
      means <- vapply(classes, function(cl) mean(tab$spearman_rho[tab$reference_class == cl]), numeric(1))
      ranking <- order(-means, match(names(means), classes))
      best <- classes[ranking[1L]]; second <- classes[ranking[2L]]
      votes <- vapply(per_q, function(z) z$reference_class[which.max(z$spearman_rho)], "")
      pr <- proposals[as.character(proposals$cluster) == q, , drop = FALSE]
      decisions[[q]] <- data.frame(cluster = q, n_cells = pr$n_cells,
        initial_class = as.character(pr$proposed_major_class), assigned_class = best,
        runner_up = second, mean_spearman = means[best], runner_up_mean_spearman = means[second],
        score_margin = means[best] - means[second], eligible_libraries = length(per_q),
        library_winner_agreement = mean(votes == best),
        canonical_marker_candidate = as.character(pr$dominant_canonical_candidate),
        matches_canonical_candidate = best == pr$dominant_canonical_candidate,
        vascular_review_flag = pr$vascular_review_flag,
        original_review_notes = pr$review_notes,
        interpretation = "Nearest of six requested classes; not independent proof of cell identity",
        stringsAsFactors = FALSE)
      per_library[[q]] <- tab
    }
    matching <- do.call(rbind, decisions); rownames(matching) <- NULL
    score_table <- do.call(rbind, per_library); rownames(score_table) <- NULL
    write.csv(matching, file.path(out, "source_data", "nearest_class_cluster_assignments.csv"), row.names = FALSE)
    write.csv(score_table, file.path(out, "source_data", "nearest_class_scores_by_library.csv"), row.names = FALSE)
    map <- stats::setNames(as.character(proposals$proposed_major_class), as.character(proposals$cluster))
    map[matching$cluster] <- matching$assigned_class
    assigned <- unname(map[cluster])
    if (anyNA(assigned) || !setequal(unique(assigned), classes)) stop("Six-class assignment is incomplete")
    cache$coordinates$initial_major_class <- initial
    cache$coordinates$proposed_major_class <- assigned
    cache$coordinates$annotation_origin <- ifelse(cluster %in% unknown, "nearest_class_match", "initial_marker_supported")
    cache$coordinates$forced_display_only <- cluster %in% unknown
    cache$coordinates$vascular_review_flag <- proposals$vascular_review_flag[match(cluster, as.character(proposals$cluster))]
    cache$proposals$initial_major_class <- as.character(cache$proposals$proposed_major_class)
    cache$proposals$proposed_major_class <- unname(map[as.character(cache$proposals$cluster)])
    cache$proposals$annotation_status <- "six_class_nearest_match_pending_user_review"
    if (!identical(cache$coordinates[, c("cell_id", "UMAP_1", "UMAP_2")],
      d[, c("cell_id", "UMAP_1", "UMAP_2")])) stop("Annotation changed coordinates")
    saveRDS(cache, file.path(out, "source_data", "global_marker_review_payload.rds"))
    con <- gzfile(file.path(out, "source_data", "global_cell_annotations_UMAP.csv.gz"), "wt")
    write.csv(cache$coordinates, con, row.names = FALSE); close(con)
    write.csv(cache$proposals, file.path(out, "source_data", "cluster_annotation_proposals.csv"), row.names = FALSE)
    write.csv(as.data.frame(table(factor(assigned, levels = classes))),
      file.path(out, "source_data", "six_class_nucleus_counts.csv"), row.names = FALSE)

    update("saving_six_class_metadata_object")
    obj$initial_major_class <- initial
    obj$proposed_major_class <- assigned
    obj$major_class_six <- assigned
    obj$annotation_origin <- cache$coordinates$annotation_origin
    obj$forced_display_only <- cache$coordinates$forced_display_only
    obj$vascular_review_flag <- cache$coordinates$vascular_review_flag
    obj@misc$six_class_assignment <- list(method = "Per-library Spearman nearest-centroid matching; equal library weights",
      initial_reference_classes_only = TRUE, features = features,
      minimum_query_nuclei_per_library = 10L, minimum_reference_nuclei_per_library_class = 20L,
      cluster_assignments = matching, user_requested_six_class_assignment = TRUE,
      marker_confirmation_pending = TRUE, source_global_object = file.path(root, "03_global", "global_celltype_review.rds"))
    after_counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
    if (!identical(before_cells, colnames(obj)) || !identical(before_xy, SeuratObject::Embeddings(obj, "umap.global")) ||
        !identical(before_clusters, as.character(SeuratObject::Idents(obj))) ||
        !identical(before_counts, after_counts)) stop("Metadata assignment changed cells, counts, clusters or UMAP")
    saveRDS(obj, file.path(out, "global_celltype_six_class.rds"))
    report <- list(status = "complete", updated_UTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      n_cells = nrow(d), n_clusters = nrow(proposals), reassigned_clusters = matching$cluster,
      reassigned_nuclei = sum(matching$n_cells), n_classes = 6L, counts_unchanged = TRUE,
      count_sum = before_count_sum, count_nonzeros = before_count_nnz,
      cell_names_unchanged = TRUE, clusters_unchanged = TRUE, UMAP_unchanged = TRUE,
      method = "Spearman correlation of log1p mean linear-normalized expression against original clear-class profiles within each library; arithmetic mean correlations across eligible libraries",
      feature_selection = "Original 2000 HVGs plus canonical broad-class markers, excluding mitochondrial/ribosomal genes and invariant reference genes",
      minimum_query_nuclei_per_library = 10L, minimum_reference_nuclei_per_library_class = 20L,
      annotation_is_probability = FALSE, vascular_flags_preserved = TRUE,
      subclustering_started = FALSE, downstream_started = FALSE,
      script_md5 = unname(tools::md5sum(file.path(root, "scripts", "06_match_six_major_classes.R"))))
    jsonlite::write_json(report, file.path(out, "assignment_complete.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
    update("six_class_annotation_ready", "complete", reassigned_clusters = matching$cluster)
    print(matching[, c("cluster", "initial_class", "assigned_class", "mean_spearman", "score_margin", "library_winner_agreement", "vascular_review_flag")], row.names = FALSE)
    invisible(report)
  }, error = function(e) {
    update("failed", "failed", error = conditionMessage(e)); stop(e)
  })
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) != 2L || a[[1L]] != "--root") stop("Usage: Rscript 06_match_six_major_classes.R --root RUN_ROOT")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
  sixclass_main(a[[2L]])
}
