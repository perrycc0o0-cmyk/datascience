#!/usr/bin/env Rscript
# Independently recompute Exc/Inh latent spaces from the approved stage-1 cells.
# No new QC, cell exclusion, subtype assignment or condition comparison is performed.

s2_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
s2_json <- function(value, path) jsonlite::write_json(value, path, pretty = TRUE,
  auto_unbox = TRUE, null = "null", na = "null", dataframe = "rows", digits = 16)
s2_atomic_rds <- function(value, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  tmp <- tempfile(".incomplete-", tmpdir = dirname(path))
  saveRDS(value, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) stop("Atomic rename failed; temporary object retained: ", tmp)
}
s2_levels <- function(value) {
  ids <- unique(as.character(value)); nums <- suppressWarnings(as.numeric(ids))
  ids[order(is.na(nums), nums, ids)]
}
s2_counts <- function(obj, assay) {
  a <- obj[[assay]]
  if (is.null(a)) stop("Missing assay: ", assay)
  layers <- SeuratObject::Layers(a, search = NA)
  cc <- grep("^counts($|\\.)", layers, value = TRUE)
  if (length(cc) != 1L) stop("Expected one counts layer in ", assay)
  x <- SeuratObject::LayerData(a, layer = cc)
  if (!inherits(x, "sparseMatrix") || any(!is.finite(x@x)) || any(x@x < 0) ||
      any(x@x != round(x@x))) stop("Counts must be sparse nonnegative integer UMIs: ", assay)
  x
}
s2_copy_feature_metadata <- function(from, to) {
  if (!identical(rownames(from), rownames(to))) stop("Feature order differs during metadata transfer")
  fm <- from[[]]
  keep <- intersect(c("gene_id", "original_symbol", "seurat_feature", "duplicate_symbol",
    "passes_initial_min_cells", "source_sample_id"), names(fm))
  if (!length(keep)) stop("Original gene identity metadata is missing")
  fm <- fm[, keep, drop = FALSE]; rownames(fm) <- rownames(to)
  SeuratObject::AddMetaData(to, metadata = fm)
}
s2_write_gz_csv <- function(tab, path) {
  con <- gzfile(path, "wt"); on.exit(close(con))
  write.csv(tab, con, row.names = FALSE)
}

run_stage2_recluster <- function(root, config = file.path(root, "config", "stage2_config.json")) {
  root <- normalizePath(root, mustWork = TRUE)
  cfg <- jsonlite::read_json(config, simplifyVector = FALSE)
  libs <- unlist(cfg$runtime$r_libraries, use.names = FALSE)
  .libPaths(unique(c(libs, .libPaths())))
  required <- c("Seurat", "SeuratObject", "Matrix", "harmony", "uwot", "jsonlite", "future")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required dependencies: ", paste(missing, collapse = ", "))
  if (!identical(cfg$run_id, "20261101_Manuscript_matched_stage2_01") ||
      !identical(cfg$stage_gate, "neuron_subtype_review")) stop("Unexpected run identity or review gate")
  if (!isTRUE(cfg$keep_previously_forced_clusters) || !isTRUE(cfg$do_not_refilter))
    stop("The approved cells must not be refiltered or silently excluded")
  processing <- cfg$processing
  expected <- list(normalization = "LogNormalize", scale_factor = 10000, hvg_method = "vst", n_hvg = 2000,
    pca_components_max = 100, harmony_pca_dimensions_max = 30, harmony_group_by = "sample_id",
    harmony_max_iter = 10, harmony_project_dim = FALSE, neighbor_harmony_dimensions_max = 19,
    clustering_algorithm = 1, resolution = .5, umap_engine = "uwot::umap", umap_metric = "cosine",
    umap_repulsion_strength = 1, umap_sgd_threads = 1)
  for (nm in names(expected)) if (!isTRUE(all.equal(processing[[nm]], expected[[nm]], check.attributes = FALSE)))
    stop("Manuscript subclustering parameter mismatch: ", nm)
  if (length(unlist(processing$scale_vars_to_regress))) stop("Subset Methods specify scaling without regression")
  marker_cfg <- cfg$annotation_marker_discovery
  marker_expected <- list(test = "wilcox", min_pct = .25, logfc_threshold = .25,
    max_cells_per_ident = 1000, only_pos = FALSE, return_thresh = 1)
  for (nm in names(marker_expected)) if (!isTRUE(all.equal(marker_cfg[[nm]], marker_expected[[nm]], check.attributes = FALSE)))
    stop("Annotation-marker setting mismatch: ", nm)
  threads <- as.integer(cfg$runtime$threads); seed <- as.integer(cfg$analysis_seed)
  if (threads != 4L || seed != 42L || as.integer(cfg$runtime$blas_threads) != 1L) stop("Unexpected reproducibility budget")
  Sys.setenv(OMP_NUM_THREADS = threads, OMP_THREAD_LIMIT = threads, OPENBLAS_NUM_THREADS = 1,
    MKL_NUM_THREADS = 1, RCPP_PARALLEL_NUM_THREADS = threads, CUDA_VISIBLE_DEVICES = "")
  future::plan(future::sequential)
  options(Seurat.object.assay.version = "v5")
  for (z in c("logs", "provenance")) dir.create(file.path(root, z), recursive = TRUE, showWarnings = FALSE)
  status_path <- file.path(root, "logs", "01_recluster_neurons.status.json")
  done_path <- file.path(root, "reclustering_complete.json")
  if (file.exists(done_path)) stop("Stage-2 reclustering already completed; do not overwrite")
  existing <- vapply(cfg$lineages, function(z) file.exists(file.path(root, z$output_directory,
    "neuron_reclustered.rds")) || file.exists(file.path(root, z$output_directory, "complete.json")), logical(1))
  if (any(existing)) stop("Existing lineage outputs require explicit recovery, not an automatic retry")
  state <- list(status = "running", stage = "preflight", lineage = NULL, started_UTC = s2_now(),
    updated_UTC = s2_now(), completed_lineages = character(), stage_gate = cfg$stage_gate,
    subtype_annotation_started = FALSE, downstream_started = FALSE)
  lineage_path <- NULL
  update <- function(stage, lineage = NULL) {
    state$stage <<- stage; state$lineage <<- lineage; state$updated_UTC <<- s2_now()
    s2_json(state, status_path)
    if (!is.null(lineage_path)) s2_json(state, file.path(lineage_path, "status.json"))
    cat(s2_now(), if (is.null(lineage)) "ALL" else lineage, stage, "\n"); flush.console()
  }
  tryCatch({
    update("loading_approved_global_object_once")
    input <- normalizePath(cfg$input_rds, mustWork = TRUE)
    parent_complete <- file.path(cfg$parent_stage1_root, "05_six_class_review", "assignment_complete.json")
    if (!file.exists(parent_complete)) stop("Approved six-class source is not complete")
    input_before <- file.info(input)[, c("size", "mtime")]
    parent <- readRDS(input)
    parent_md <- parent[[]]
    necessary <- c("sample_id", "condition_label", "source_sample_id", "source_cell_id", "global_cluster",
      "major_class_six", "initial_major_class", "forced_display_only", "vascular_review_flag", "percent.mt")
    if (!all(necessary %in% names(parent_md))) stop("Required source metadata or flags missing")
    if (!identical(rownames(parent_md), colnames(parent)) || anyDuplicated(colnames(parent))) stop("Invalid parent cell IDs")
    expected_samples <- c(paste0("WD", 1:3), paste0("Ext", 1:3), paste0("Ret-Ext", 1:3))
    if (!setequal(unique(as.character(parent_md$sample_id)), expected_samples) ||
        any(!grepl("^(WD|Ext|Ret-Ext)[1-3]__", colnames(parent)))) stop("Canonical sample/cell names are not preserved")
    if (!all(as.character(parent_md$condition_label) == sub("[1-3]$", "", as.character(parent_md$sample_id))))
      stop("Condition labels disagree with manuscript sample naming")
    parent_rna <- s2_counts(parent, "RNA")
    parent_observed <- s2_counts(parent, "observed")
    if (!identical(dimnames(parent_rna), dimnames(parent_observed))) stop("Observed/corrected count alignment differs")
    runtime <- list(R = R.version.string, package_versions = stats::setNames(lapply(required,
      function(p) as.character(utils::packageVersion(p))), required), library_paths = .libPaths(),
      future_plan = "sequential", threads = threads, BLAS_threads = 1L, seed = seed)
    s2_json(runtime, file.path(root, "provenance", "runtime.json"))
    s2_json(list(input_rds = input, input_bytes = input_before$size, input_mtime = as.character(input_before$mtime),
      parent_nuclei = ncol(parent), config_md5 = unname(tools::md5sum(config)),
      script_md5 = unname(tools::md5sum(file.path(root, "scripts", "01_recluster_neurons.R")))),
      file.path(root, "provenance", "input_provenance.json"))
    summaries <- list()
    for (spec in cfg$lineages) {
      id <- spec$id
      if (!id %in% c("Exc", "Inh")) stop("Only Exc and Inh are authorized")
      expected_spec <- if (id == "Exc") list(source_class = "Excitatory", expected_nuclei = 62991,
        neighbors_k = 20, umap_n_neighbors = 10, umap_min_dist = .8, umap_spread = 6) else
        list(source_class = "Inhibitory", expected_nuclei = 21693, neighbors_k = 30,
          umap_n_neighbors = 30, umap_min_dist = .3, umap_spread = 1)
      for (nm in names(expected_spec)) if (!isTRUE(all.equal(spec[[nm]], expected_spec[[nm]], check.attributes = FALSE)))
        stop("Lineage Methods parameter mismatch: ", id, "/", nm)
      if (basename(spec$output_directory) != spec$output_directory) stop("Invalid lineage output directory")
      lineage_path <- file.path(root, spec$output_directory)
      for (z in c(lineage_path, file.path(lineage_path, "source_data"), file.path(lineage_path, "checkpoints")))
        dir.create(z, recursive = TRUE, showWarnings = FALSE)
      update("extracting_approved_lineage_counts", id)
      cells <- colnames(parent)[as.character(parent_md[[cfg$source_class_column]]) == spec$source_class]
      if (length(cells) != as.integer(spec$expected_nuclei)) stop("Approved lineage cell count mismatch: ", id)
      md <- parent_md[cells, , drop = FALSE]
      md$original_global_cluster <- as.character(md$global_cluster)
      md$original_global_seurat_clusters <- as.character(md$seurat_clusters)
      md$seurat_clusters <- NULL
      md$lineage <- id
      corrected <- parent_rna[, cells, drop = FALSE]
      observed <- parent_observed[, cells, drop = FALSE]
      fresh_assay <- SeuratObject::CreateAssay5Object(counts = corrected, min.cells = 0, min.features = 0)
      obj <- SeuratObject::CreateSeuratObject(counts = fresh_assay, project = paste0("stage2_", id),
        min.cells = 0, min.features = 0, meta.data = md)
      obj[["RNA"]] <- s2_copy_feature_metadata(parent[["RNA"]], obj[["RNA"]])
      obj[["observed"]] <- SeuratObject::CreateAssay5Object(counts = observed, min.cells = 0, min.features = 0)
      obj[["observed"]] <- s2_copy_feature_metadata(parent[["observed"]], obj[["observed"]])
      SeuratObject::DefaultAssay(obj) <- "RNA"
      obj$orig.ident <- as.character(md$sample_id)
      if (length(obj@reductions) || length(obj@graphs) || length(obj@neighbors) ||
          length(SeuratObject::VariableFeatures(obj))) stop("A fresh subset inherited global latent structures")
      if (!identical(colnames(obj), cells) || !identical(s2_counts(obj, "RNA"), corrected) ||
          !identical(s2_counts(obj, "observed"), observed)) stop("Fresh subset changed counts or cell order")
      for (nm in necessary) if (!identical(as.character(obj[[]][[nm]]), as.character(md[[nm]])))
        stop("Source metadata changed: ", nm)
      input_audit <- data.frame(cell_id = cells, sample_id = md$sample_id, condition_label = md$condition_label,
        original_global_cluster = md$original_global_cluster, initial_major_class = md$initial_major_class,
        forced_display_only = md$forced_display_only, vascular_review_flag = md$vascular_review_flag,
        corrected_UMI = as.numeric(Matrix::colSums(corrected)),
        observed_UMI = as.numeric(Matrix::colSums(observed)), stringsAsFactors = FALSE)
      s2_write_gz_csv(input_audit, file.path(lineage_path, "source_data", "input_cells_audit.csv.gz"))
      write.csv(as.data.frame(table(sample_id = factor(md$sample_id, levels = expected_samples))),
        file.path(lineage_path, "source_data", "input_nuclei_by_sample.csv"), row.names = FALSE)
      update("fresh_normalization_and_2000_HVG", id)
      obj <- Seurat::NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
      obj <- Seurat::FindVariableFeatures(obj, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      hvg <- SeuratObject::VariableFeatures(obj, assay = "RNA")
      if (length(hvg) != 2000L) stop("Expected exactly 2000 subset HVGs")
      update("scale_without_regression_and_PCA100", id)
      obj <- Seurat::ScaleData(obj, assay = "RNA", features = hvg, vars.to.regress = NULL,
        block.size = 500, verbose = TRUE)
      npcs <- as.integer(min(100L, length(hvg) - 1L, ncol(obj) - 1L))
      obj <- Seurat::RunPCA(obj, assay = "RNA", features = hvg, npcs = npcs,
        seed.use = seed, reduction.name = "pca", reduction.key = "PC_", verbose = TRUE)
      pc <- SeuratObject::Embeddings(obj[["pca"]])
      if (ncol(pc) < 30L || !identical(rownames(pc), cells) || any(!is.finite(pc))) stop("Invalid fresh PCA")
      harmony_dims <- seq_len(min(30L, ncol(pc)))
      update("sample_only_Harmony_first30_PCs", id)
      obj <- harmony::RunHarmony(obj, group.by.vars = "sample_id", reduction.use = "pca",
        dims.use = harmony_dims, reduction.save = "harmony.sample", project.dim = FALSE,
        ncores = threads, max_iter = 10L, plot_convergence = FALSE, verbose = TRUE)
      h <- SeuratObject::Embeddings(obj[["harmony.sample"]])
      if (!identical(rownames(h), cells) || ncol(h) != length(harmony_dims) || any(!is.finite(h)))
        stop("Invalid Harmony embedding; no PCA fallback")
      graph_dims <- seq_len(min(19L, ncol(h)))
      if (length(graph_dims) != 19L) stop("Nineteen Harmony dimensions are required")
      update("neighbors_and_resolution_point5", id)
      graph_prefix <- tolower(id)
      graphs <- paste0(graph_prefix, c("_nn", "_snn"))
      obj <- Seurat::FindNeighbors(obj, reduction = "harmony.sample", dims = graph_dims,
        k.param = as.integer(spec$neighbors_k), graph.name = graphs, verbose = TRUE)
      obj <- Seurat::FindClusters(obj, graph.name = graphs[2L], resolution = .5,
        algorithm = 1L, cluster.name = "subcluster", random.seed = seed, verbose = TRUE)
      obj$subcluster <- factor(as.character(obj$subcluster), levels = s2_levels(obj$subcluster))
      obj$seurat_clusters <- obj$subcluster
      SeuratObject::Idents(obj) <- "subcluster"
      update("fresh_lineage_UMAP", id)
      umap_args <- list(n_neighbors = as.integer(spec$umap_n_neighbors), n_components = 2L, metric = "cosine",
        min_dist = spec$umap_min_dist, spread = spec$umap_spread, set_op_mix_ratio = 1,
        local_connectivity = 1, repulsion_strength = 1, negative_sample_rate = 5,
        learning_rate = 1, n_epochs = NULL, init = "spectral", n_threads = threads,
        n_sgd_threads = 1L, fast_sgd = FALSE, ret_model = FALSE, verbose = TRUE)
      set.seed(seed)
      xy <- do.call(uwot::umap, c(list(X = h[, graph_dims, drop = FALSE]), umap_args))
      if (!is.matrix(xy) || nrow(xy) != length(cells) || ncol(xy) != 2L || any(!is.finite(xy))) stop("Invalid fresh UMAP")
      rownames(xy) <- cells; colnames(xy) <- c("UMAP_1", "UMAP_2")
      obj[["umap.neuron"]] <- SeuratObject::CreateDimReducObject(embeddings = xy,
        key = "UMAP_", assay = "RNA", misc = list(engine = "uwot::umap", input_reduction = "harmony.sample",
          input_dims = graph_dims, parameters = umap_args, seed = seed))
      s2_atomic_rds(list(pca = obj[["pca"]], harmony = obj[["harmony.sample"]], umap = obj[["umap.neuron"]],
        graphs = obj@graphs), file.path(lineage_path, "checkpoints", "latent_and_graphs.rds"))
      current_md <- obj[[]]
      coord <- data.frame(cell_id = cells, sample_id = current_md$sample_id,
        condition_label = current_md$condition_label, subcluster = as.character(obj$subcluster),
        original_global_cluster = current_md$original_global_cluster,
        initial_major_class = current_md$initial_major_class, approved_major_class = current_md$major_class_six,
        forced_display_only = current_md$forced_display_only, vascular_review_flag = current_md$vascular_review_flag,
        source_sample_id = current_md$source_sample_id, source_cell_id = current_md$source_cell_id,
        UMAP_1 = xy[, 1L], UMAP_2 = xy[, 2L], stringsAsFactors = FALSE)
      s2_write_gz_csv(coord, file.path(lineage_path, "source_data", "cell_annotations_UMAP.csv.gz"))
      write.csv(data.frame(gene = hvg), file.path(lineage_path, "source_data", "HVG_2000.csv"), row.names = FALSE)
      sd <- SeuratObject::Stdev(obj[["pca"]]); pct <- 100 * sd^2 / sum(sd^2)
      write.csv(data.frame(PC = seq_along(sd), stdev = sd, variance = sd^2,
        percent_of_computed_PCs = pct, cumulative_percent = cumsum(pct),
        input_to_Harmony = seq_along(sd) %in% harmony_dims,
        used_Harmony_dimension_count = length(graph_dims)),
        file.path(lineage_path, "source_data", "PCA_variance_profile.csv"), row.names = FALSE)
      write.csv(as.data.frame(table(subcluster = obj$subcluster)),
        file.path(lineage_path, "source_data", "cluster_counts.csv"), row.names = FALSE)
      write.csv(obj[["RNA"]][[]], file.path(lineage_path, "source_data", "gene_feature_map.csv"), row.names = FALSE)
      update("annotation_marker_discovery_not_condition_DE", id)
      set.seed(seed)
      markers <- Seurat::FindAllMarkers(obj, assay = "RNA", slot = "data", test.use = "wilcox",
        only.pos = FALSE, min.pct = .25, logfc.threshold = .25, max.cells.per.ident = 1000,
        random.seed = seed, return.thresh = 1, verbose = TRUE)
      if (!is.data.frame(markers) || !all(c("gene", "cluster") %in% names(markers))) stop("Marker discovery did not return its expected table")
      write.csv(markers, file.path(lineage_path, "source_data", "global_cluster_markers.csv"), row.names = FALSE)
      if (!identical(colnames(obj), cells) || !identical(s2_counts(obj, "RNA"), corrected) ||
          !identical(s2_counts(obj, "observed"), observed)) stop("Reprocessing changed counts or membership")
      for (nm in necessary) if (!identical(as.character(obj[[]][[nm]]), as.character(md[[nm]])))
        stop("Reprocessing changed retained source metadata: ", nm)
      pars <- list(lineage = id, n_cells = ncol(obj), n_genes = nrow(obj),
        source_class_column = cfg$source_class_column, source_class = spec$source_class,
        source_global_object = input, new_QC_or_cell_exclusions = FALSE, previously_forced_cells_retained = sum(md$forced_display_only),
        NormalizeData_scale_factor = 10000, HVG_method = "vst", HVGs = length(hvg), regression = character(),
        PCA_computed = ncol(pc), Harmony_input_PCs = harmony_dims, Harmony_group = "sample_id",
        Harmony_max_iter = 10L, Harmony_project_dim = FALSE, graph_Harmony_dims = graph_dims,
        neighbors_k = spec$neighbors_k, clustering_algorithm = "Louvain, Seurat algorithm=1", resolution = .5,
        UMAP = list(engine = "uwot::umap", input_Harmony_dims = graph_dims, parameters = umap_args),
        annotation_marker_discovery = marker_cfg, runtime = runtime, seed = seed,
        count_sum_corrected = sum(corrected), count_sum_observed = sum(observed),
        source_counts_and_metadata_unchanged = TRUE, subtype_annotation_started = FALSE, downstream_started = FALSE)
      obj@misc$manuscript_stage2 <- pars
      s2_json(pars, file.path(lineage_path, "source_data", "parameters.json"))
      update("saving_atomic_lineage_object_without_compression", id)
      s2_atomic_rds(obj, file.path(lineage_path, "neuron_reclustered.rds"))
      writeLines(capture.output(utils::sessionInfo()), file.path(lineage_path, "sessionInfo.txt"))
      one <- list(status = "complete", lineage = id, completed_UTC = s2_now(), n_cells = ncol(obj),
        n_clusters = length(levels(obj$subcluster)), n_HVGs = length(hvg), n_PCs = ncol(pc),
        Harmony_dimensions = ncol(h), graph_dimensions = length(graph_dims),
        counts_and_metadata_unchanged = TRUE, subtype_annotation_started = FALSE, downstream_started = FALSE,
        awaiting = "Subtype-marker annotation and UMAP/DotPlot user review; no downstream analysis")
      s2_json(one, file.path(lineage_path, "complete.json"))
      summaries[[id]] <- one
      state$completed_lineages <- c(state$completed_lineages, id)
      update("lineage_reclustering_complete", id)
      rm(obj, corrected, observed, fresh_assay, input_audit, coord, markers, pc, h, xy, current_md, md); gc(FALSE)
      lineage_path <- NULL
    }
    if (!identical(input_before, file.info(input)[, c("size", "mtime")])) stop("Parent input file changed during run")
    state$status <- "complete"; state$stage <- "independent_neuron_reclustering_ready"
    state$lineage <- NULL; state$updated_UTC <- s2_now(); state$lineages <- summaries
    state$awaiting <- "Subtype annotation/plots, followed by the user's second review gate"
    s2_json(state, status_path); s2_json(state, done_path)
    cat(s2_now(), "Reclustering complete; no subtype or downstream analyses were run.\n")
    invisible(state)
  }, error = function(e) {
    state$status <- "failed"; state$error <- conditionMessage(e); state$updated_UTC <- s2_now()
    s2_json(state, status_path)
    if (!is.null(lineage_path)) s2_json(state, file.path(lineage_path, "status.json"))
    stop(e)
  })
}

if (sys.nframe() == 0L && Sys.getenv("STAGE2_SOURCE_ONLY") != "1") {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) != 2L || a[1L] != "--root") stop("Usage: Rscript 01_recluster_neurons.R --root NEW_STAGE2_ROOT")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required to read configuration")
  run_stage2_recluster(a[2L])
}
