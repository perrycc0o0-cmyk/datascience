# 06 match six major classes

# 1. Packages
library(jsonlite, include.only = character())
library(Seurat, include.only = character())
library(SeuratObject, include.only = character())
library(Matrix, include.only = character())

# 2. Paths
input_dir <- "/path/to/matrices"
work_dir <- "/path/to/work"
result_dir <- file.path(work_dir, "stage1")
figure_dir <- result_dir
fontconfig_file <- NULL
script_path <- "/path/to/06_match_six_major_classes.R"

# 3. Samples and conditions
sample_cols <- c("WD1", "WD2", "WD3", "Ext1", "Ext2", "Ext3", "Ret-Ext1", "Ret-Ext2", "Ret-Ext3")
metadata <- data.frame(
  pool_id = sample_cols,
  source_sample_id = c("Con1", "Con2", "Con3", "Ext1", "Ext2", "Ext3", "Ret_ext1", "Ret_ext2", "Ret_ext3"),
  model = rep("Rat IL", 9L),
  condition = factor(rep(c("WD", "Ext", "Ret-Ext"), each = 3L), levels = c("WD", "Ext", "Ret-Ext")),
  rats_per_pool = rep(5L, 9L), stringsAsFactors = FALSE)

# 4. Parameters and marker definitions
processing_settings <- list(
  `schema_version` = "1.0",
  `run_id` = "20261101_Manuscript_matched_stage1_01",
  `stage_gate` = "global_celltype_review",
  `analysis_seed` = 42L,
  `threads` = 4L,
  `samples` = list(
    list(
      `sample_id` = "Con1",
      `condition_label` = "WD"
    ),
    list(
      `sample_id` = "Con2",
      `condition_label` = "WD"
    ),
    list(
      `sample_id` = "Con3",
      `condition_label` = "WD"
    ),
    list(
      `sample_id` = "Ext1",
      `condition_label` = "Ext"
    ),
    list(
      `sample_id` = "Ext2",
      `condition_label` = "Ext"
    ),
    list(
      `sample_id` = "Ext3",
      `condition_label` = "Ext"
    ),
    list(
      `sample_id` = "Ret_ext1",
      `condition_label` = "Ret-Ext"
    ),
    list(
      `sample_id` = "Ret_ext2",
      `condition_label` = "Ret-Ext"
    ),
    list(
      `sample_id` = "Ret_ext3",
      `condition_label` = "Ret-Ext"
    )
  ),
  `runtime` = list(
    `threads` = 4L,
    `blas_threads` = 1L,
    `r_libraries` = list(),
    `fontconfig_file` = NULL,
    `font_family` = "Arial"
  ),
  `initial_eligibility` = list(
    `method` = "native Seurat CreateSeuratObject on original called counts",
    `min_cells_per_gene` = 3L,
    `min_features_per_cell` = 300L,
    `preserve_all_called_barcode_audit` = TRUE
  ),
  `qc` = list(
    `metric_assay` = "observed original UMI counts",
    `metric_gene_scope` = "genes retained by native initial min.cells=3/min.features=300 qualification",
    `nFeature_lower_exclusive` = 700L,
    `nFeature_upper_exclusive` = 5000L,
    `nCount_lower_exclusive` = 500L,
    `nCount_upper_exclusive` = 30000L,
    `percent_mt_upper_exclusive` = 5L,
    `log10GenesPerUMI_lower_exclusive` = 0.75,
    `complexity_formula` = "log10(nFeature) / log10(nCount)",
    `nonfinite_metrics_fail` = TRUE,
    `MAD_filters` = FALSE,
    `mitochondrial_gene_ids` = list(
      "ENSRNOG00000030644",
      "ENSRNOG00000031033",
      "ENSRNOG00000034234",
      "ENSRNOG00000030371",
      "ENSRNOG00000033299",
      "ENSRNOG00000031979",
      "ENSRNOG00000030700",
      "ENSRNOG00000033615",
      "ENSRNOG00000031053",
      "ENSRNOG00000029707",
      "ENSRNOG00000029971",
      "ENSRNOG00000029042",
      "ENSRNOG00000031766"
    ),
    `mitochondrial_symbols` = list(
      "ND1",
      "ND2",
      "COX1",
      "COX2",
      "ATP8",
      "ATP6",
      "COX3",
      "ND3",
      "ND4L",
      "ND4",
      "ND5",
      "ND6",
      "CYTB"
    ),
    `ribo_symbol_pattern` = "^(Rpl|Rps)",
    `ribo_filter` = FALSE,
    `corrected_metrics_role` = "audit only; cannot rescue observed-QC failures",
    `full_GEX_metrics_role` = "diagnostic only; not used for the manuscript thresholds"
  ),
  `ambient` = list(
    `method` = "SoupX",
    `input` = "fresh original raw/filtered H5, all called barcodes and input genes",
    `precluster_genes_min` = 200L,
    `precluster_umi_min` = 300L,
    `precluster_n_hvg` = 2000L,
    `precluster_n_pcs` = 30L,
    `precluster_neighbor_dims` = 20L,
    `precluster_resolution` = 0.6,
    `round_to_integer` = TRUE,
    `rho_review_max` = 0.3,
    `force_accept` = FALSE,
    `reuse_previous_analysis_counts` = FALSE
  ),
  `doublet` = list(
    `method` = "scDblFinder",
    `per_capture` = TRUE,
    `input` = "corrected counts of initial-qualified genes and observed-QC-passing nuclei",
    `clusters` = FALSE,
    `dbr_per_1000` = 0.004,
    `dbr_denominator` = "original Cell Ranger-called barcode count per capture before initial eligibility and QC",
    `dbr_sd` = "package default",
    `xgb_nthreads` = 4L,
    `BPPARAM` = "SerialParam",
    `chemistry` = "Chromium GEM-X Single Cell 3 prime v4",
    `vendor_doublet_calls_role` = "not combined with exclusions"
  ),
  `global` = list(
    `n_hvg` = 2000L,
    `n_pcs` = 50L,
    `pc_max` = 30L,
    `regress_vars` = list(
      "nCount_RNA",
      "percent.mt"
    ),
    `harmony_group_by` = "sample_id",
    `regress_condition` = FALSE,
    `neighbors_k` = 20L,
    `resolutions` = list(
      0.1,
      0.2,
      0.3,
      0.4,
      0.5,
      0.6,
      0.7,
      0.8,
      0.9,
      1.0
    ),
    `primary_resolution` = 1.0,
    `umap_n_neighbors` = 30L,
    `umap_min_dist` = 0.3,
    `umap_metric` = "cosine",
    `seed` = 42L,
    `threads` = 4L,
    `percent_mt_metadata_source` = "percent_mt_observed in initial-qualified gene set",
    `annotation_policy` = "fresh clusters; major-cell-type marker review before any lineage subclustering"
  ),
  `output_interface` = list(
    `sample_directory` = "02_qc_doublet/<sample_id>",
    `singlet_object` = "02_qc_doublet/<sample_id>/singlet_clean.rds",
    `all_called_QC` = "02_qc_doublet/<sample_id>/all_called_barcodes_QC.csv.gz",
    `aggregate_summary` = "02_qc_doublet/attrition_summary.csv",
    `primary_assay` = "RNA: SoupX-corrected integer UMI counts",
    `audit_assay` = "observed: matching original UMI counts",
    `metadata_join_key` = "sample_id__10x-barcode",
    `no_old_embeddings_or_labels` = TRUE
  ),
  `input_root` = NULL
)
processing_settings$input_root <- input_dir
processing_settings$runtime["fontconfig_file"] <- list(fontconfig_file)
processing_settings$samples <- lapply(seq_len(nrow(metadata)), function(i) list(
  sample_id = metadata$source_sample_id[i], condition_label = as.character(metadata$condition[i])))
sample_naming <- list(
  `schema_version` = "1.0",
  `requested_on` = "2026-11-01",
  `source_files_read_only` = TRUE,
  `naming_only_no_count_changes` = TRUE,
  `cell_id_format` = "<sample_id>__<original_10x_barcode>",
  `canonical_conditions` = list(
    "WD",
    "Ext",
    "Ret-Ext"
  ),
  `sample_map` = list(
    list(
      `source_sample_id` = "Con1",
      `sample_id` = "WD1",
      `condition_label` = "WD"
    ),
    list(
      `source_sample_id` = "Con2",
      `sample_id` = "WD2",
      `condition_label` = "WD"
    ),
    list(
      `source_sample_id` = "Con3",
      `sample_id` = "WD3",
      `condition_label` = "WD"
    ),
    list(
      `source_sample_id` = "Ext1",
      `sample_id` = "Ext1",
      `condition_label` = "Ext"
    ),
    list(
      `source_sample_id` = "Ext2",
      `sample_id` = "Ext2",
      `condition_label` = "Ext"
    ),
    list(
      `source_sample_id` = "Ext3",
      `sample_id` = "Ext3",
      `condition_label` = "Ext"
    ),
    list(
      `source_sample_id` = "Ret_ext1",
      `sample_id` = "Ret-Ext1",
      `condition_label` = "Ret-Ext"
    ),
    list(
      `source_sample_id` = "Ret_ext2",
      `sample_id` = "Ret-Ext2",
      `condition_label` = "Ret-Ext"
    ),
    list(
      `source_sample_id` = "Ret_ext3",
      `sample_id` = "Ret-Ext3",
      `condition_label` = "Ret-Ext"
    )
  )
)
sample_naming$sample_map <- lapply(seq_len(nrow(metadata)), function(i) list(
  source_sample_id = metadata$source_sample_id[i], sample_id = metadata$pool_id[i],
  condition_label = as.character(metadata$condition[i])))

# 5. Helper functions
save_settings <- function(value, filename) {
  path <- file.path(result_dir, "provenance", filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    previous <- jsonlite::read_json(path, simplifyVector = FALSE)
    if (!isTRUE(all.equal(previous, value, tolerance = 0))) stop("Existing settings differ: ", path)
  } else {
    jsonlite::write_json(value, path, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = 17)
  }
  path
}

# 6. Analysis and output
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
config_snapshot <- save_settings(processing_settings, "processing_settings.json")
naming_snapshot <- save_settings(sample_naming, "sample_naming.json")
# Assign previously unresolved global clusters to their nearest of six reference classes.
# This changes annotation metadata only; it does not normalize, integrate or recluster.

sixclass_main <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  cfg <- processing_settings
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
      script_md5 = unname(tools::md5sum(script_path)))
    jsonlite::write_json(report, file.path(out, "assignment_complete.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
    update("six_class_annotation_ready", "complete", reassigned_clusters = matching$cluster)
    print(matching[, c("cluster", "initial_class", "assigned_class", "mean_spearman", "score_margin", "library_winner_agreement", "vascular_review_flag")], row.names = FALSE)
    invisible(report)
  }, error = function(e) {
    update("failed", "failed", error = conditionMessage(e)); stop(e)
  })
}

sixclass_main(result_dir)
