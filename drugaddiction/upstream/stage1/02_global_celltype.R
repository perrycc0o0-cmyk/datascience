#!/usr/bin/env Rscript
# Stage 1: fresh global latent space and provisional broad-class marker evidence.
# This script deliberately contains no neuronal subclustering or condition DE.

stage1_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
stage1_log <- function(...) { cat(stage1_now(), ..., "\n"); flush.console() }
stage1_write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, null = "null", na = "null", dataframe = "rows", digits = 16)
stage1_atomic_rds <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite completed object: ", path)
  tmp <- tempfile(".incomplete-", tmpdir = dirname(path))
  saveRDS(x, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) stop("Atomic rename failed; temporary object retained: ", tmp)
}
stage1_levels <- function(x) {
  v <- unique(as.character(x)); n <- suppressWarnings(as.numeric(v))
  v[order(is.na(n), n, v)]
}
stage1_panel <- function() list(
  Excitatory = c("Slc17a7", "Satb2", "Slc30a3"),
  Inhibitory = c("Gad1", "Gad2", "Dlx1"),
  Oligodendrocytes = c("Mog", "Mbp", "St18"),
  OPC = c("Cacng4", "Pdgfra", "Vcan"),
  Astrocytes = c("Slco1c1", "Gfap", "Aqp4"),
  Microglia = c("P2ry12", "Tmem119", "Cx3cr1"))
stage1_review_genes <- function() c("Snap25", "Rbfox3", "Plp1", "Sox10", "Cnp", "Mag",
  "Aldh1l1", "Cldn5", "Pecam1", "Kdr", "Flt1", "Cdh5", "Rgs5", "Pdgfrb",
  "Kcnj8", "Abcc9", "Acta2", "Col1a1", "Col1a2", "Dcn", "Lum", "Ptprc",
  "Csf1r", "Mrc1", "Cd163", "Lyz2")

stage1_get_counts <- function(obj, assay) {
  a <- obj[[assay]]
  if (is.null(a)) stop("Missing assay: ", assay)
  if (inherits(a, "Assay5")) {
    layers <- SeuratObject::Layers(a, search = NA)
    cc <- grep("^counts($|\\.)", layers, value = TRUE)
    if (length(cc) != 1L) stop("Each input must have one counts layer: ", assay)
    x <- SeuratObject::LayerData(a, layer = cc[[1L]])
  } else x <- SeuratObject::GetAssayData(a, slot = "counts")
  if (!inherits(x, "sparseMatrix") || any(!is.finite(x@x)) || any(x@x < 0) ||
      any(x@x != round(x@x))) stop("Invalid sparse integer counts in ", assay)
  x
}

stage1_gene_stats <- function(value, gene, group, levels) {
  if (length(value) != length(group) || any(!is.finite(value)) || any(value < 0))
    stop("Invalid normalized marker values: ", gene)
  n <- as.integer(table(factor(group, levels = levels)))
  ss <- rowsum(cbind(detected = as.numeric(value > 0), log_value = value,
                    linear_value = expm1(value)), group, reorder = FALSE)
  ss <- ss[match(levels, rownames(ss)), , drop = FALSE]
  av <- ss[, "linear_value"] / n
  lv <- log1p(av); s <- stats::sd(lv)
  z <- if (is.finite(s) && s > 0) (lv - mean(lv)) / s else rep(0, length(lv))
  data.frame(gene = gene, cluster = levels, n_cells = n,
    n_cells_expressed = ss[, "detected"], pct_expressed = 100 * ss[, "detected"] / n,
    avg_normalized_data = ss[, "log_value"] / n, avg_linear_normalized = av,
    log1p_avg_linear_normalized = lv, gene_z_unclipped = z,
    avg_expression_scaled = pmax(-2.5, pmin(2.5, z)), row.names = NULL)
}

stage1_annotation <- function(stats, clusters) {
  panels <- stage1_panel()
  get <- function(cluster, gene, column) {
    i <- stats$cluster == cluster & stats$gene == gene
    if (sum(i) != 1L) return(NA_real_)
    stats[[column]][i]
  }
  decisions <- evidence <- vector("list", length(clusters))
  for (i in seq_along(clusters)) {
    cl <- clusters[[i]]
    pp <- vapply(names(panels), function(label) {
      mean(vapply(panels[[label]], function(g) get(cl, g, "gene_z_unclipped"), numeric(1)), na.rm = TRUE)
    }, numeric(1))
    detected <- function(g) vapply(g, function(z) get(cl, z, "pct_expressed"), numeric(1))
    enough <- function(g, cutoff, n = 2L) sum(detected(g) >= cutoff, na.rm = TRUE) >= n
    pass <- c(
      Excitatory = isTRUE(detected("Slc17a7") >= 10) &&
        enough(c("Satb2", "Slc30a3"), 5, 1L),
      Inhibitory = enough(c("Gad1", "Gad2"), 10),
      Oligodendrocytes = isTRUE(detected("Mbp") >= 10) &&
        enough(c("Mog", "St18"), 5, 1L),
      OPC = isTRUE(detected("Pdgfra") >= 5) &&
        enough(c("Cacng4", "Vcan"), 10, 1L),
      Astrocytes = enough(c("Slco1c1", "Gfap", "Aqp4"), 5),
      Microglia = enough(c("P2ry12", "Tmem119", "Cx3cr1"), 2))
    available <- vapply(panels, function(g) sum(g %in% stats$gene), integer(1))
    pass <- pass & available >= 2L & is.finite(pp) & pp > 0
    candidates <- names(pass)[pass]
    candidates <- candidates[order(pp[candidates], decreasing = TRUE)]
    dominant <- if (length(candidates)) candidates[[1L]] else "Unresolved"
    margin <- if (length(candidates) >= 2L) pp[candidates[1L]] - pp[candidates[2L]] else NA_real_
    mixed <- length(candidates) >= 2L && is.finite(margin) && margin <= 0.25
    endothelial <- enough(c("Cldn5", "Pecam1", "Kdr", "Flt1", "Cdh5"), 15)
    mural <- enough(c("Rgs5", "Pdgfrb", "Kcnj8", "Abcc9", "Acta2"), 15)
    fibroblast <- enough(c("Col1a1", "Col1a2", "Dcn", "Lum"), 15)
    vascular_flag <- endothelial || mural || fibroblast
    myeloid_flag <- isTRUE(detected("Ptprc") >= 20) &&
      enough(c("Csf1r", "Mrc1", "Cd163", "Lyz2"), 10) && !isTRUE(pass["Microglia"])
    proposed <- if (mixed) "Mixed" else dominant
    if (vascular_flag || myeloid_flag) proposed <- "Unresolved"
    cross_neuronal <- enough(c("Gad1", "Gad2"), 10) &&
      isTRUE(detected("Slc17a7") >= 10)
    if (cross_neuronal) proposed <- "Mixed"
    review <- c(if (mixed) "Competing canonical panels",
      if (vascular_flag) "Vascular-associated marker evidence; not forced into six classes",
      if (myeloid_flag) "Other myeloid marker evidence; not forced into Microglia",
      if (cross_neuronal) "Excitatory and inhibitory marker cross-signal",
      if (dominant == "Unresolved") "Insufficient canonical marker evidence")
    evidence[[i]] <- data.frame(cluster = cl, panel = names(panels),
      n_markers_available = available, marker_genes = vapply(panels, paste, "", collapse = ";"),
      mean_gene_z = pp, detection_gate_passed = pass, row.names = NULL)
    decisions[[i]] <- data.frame(cluster = cl,
      n_cells = stats$n_cells[match(cl, stats$cluster)],
      proposed_major_class = proposed, dominant_canonical_candidate = dominant,
      qualifying_canonical_panels = paste(candidates, collapse = ";"),
      top_panel_score = if (dominant == "Unresolved") NA_real_ else unname(pp[dominant]),
      top_two_score_margin = unname(margin), cross_neuronal_markers = cross_neuronal,
      vascular_review_flag = vascular_flag, other_myeloid_review_flag = myeloid_flag,
      annotation_status = "pending_user_marker_review", annotation_is_probability = FALSE,
      review_notes = paste(review, collapse = "; "), row.names = NULL)
  }
  list(decisions = do.call(rbind, decisions), evidence = do.call(rbind, evidence))
}

stage1_select_pcs <- function(sd, cap = 30L) {
  variance <- sd^2
  pct <- 100 * variance / sum(variance)
  cum <- cumsum(pct)
  eligible <- which(cum > 90 & pct < 5)
  drops <- which(head(pct, -1L) - tail(pct, -1L) > 0.1)
  co1 <- if (length(eligible)) eligible[[1L]] else Inf
  co2 <- if (length(drops)) max(drops) + 1L else Inf
  k <- as.integer(min(co1, co2, cap, length(sd)))
  if (!is.finite(k) || k < 2L) stop("Variance-profile PC rule returned too few dimensions.")
  list(npcs = k, criterion_cumulative90_pct5 = co1, criterion_last_drop_over_point1 = co2,
    cap = cap, definition = "100 * stdev^2 / sum(stdev^2) across the 50 computed PCs; min(first cumulative>90% and PC<5%, last adjacent percent drop>0.1 plus one, 30)",
    table = data.frame(PC = seq_along(sd), stdev = sd, variance = variance,
      percent_of_computed_PCs = pct, cumulative_percent = cum, used = seq_along(sd) <= k))
}

run_stage1_global <- function(root, config = file.path(root, "config", "processing_config.json")) {
  root <- normalizePath(root, mustWork = TRUE)
  cfg <- jsonlite::read_json(config, simplifyVector = FALSE)
  naming_path <- file.path(root, "config", "sample_naming.json")
  if (!file.exists(naming_path)) stop("Canonical sample-naming manifest is required.")
  naming <- jsonlite::read_json(naming_path, simplifyVector = FALSE)
  sample_map <- naming$sample_map
  source_ids <- vapply(sample_map, `[[`, "", "source_sample_id")
  paper_ids <- vapply(sample_map, `[[`, "", "sample_id")
  if (length(source_ids) != 9L || anyDuplicated(source_ids) || anyDuplicated(paper_ids))
    stop("Canonical sample map must contain nine unique source and paper sample IDs.")
  names(paper_ids) <- source_ids
  setting <- function(name, default) if (is.null(cfg[[name]])) default else cfg[[name]]
  seed <- as.integer(setting("analysis_seed", 42)); threads <- as.integer(setting("threads", 4))
  if (!is.null(cfg$stage_gate) && cfg$stage_gate != "global_celltype_review")
    stop("This executable is restricted to the global-celltype review gate.")
  if (threads < 1L || threads > 8L) stop("Stage 1 thread budget must be 1 to 8.")
  libs <- unlist(cfg$runtime$r_libraries, use.names = FALSE)
  if (length(libs)) .libPaths(unique(c(libs, .libPaths())))
  expected_global <- list(n_hvg = 2000, n_pcs = 50, pc_max = 30, neighbors_k = 20,
    primary_resolution = 1, umap_n_neighbors = 30, umap_min_dist = .3,
    umap_metric = "cosine", harmony_group_by = "sample_id", regress_condition = FALSE)
  for (nm in names(expected_global)) {
    if (is.null(cfg$global[[nm]]) || !isTRUE(all.equal(cfg$global[[nm]], expected_global[[nm]], check.attributes = FALSE)))
      stop("Manuscript global configuration mismatch: ", nm)
  }
  if (!isTRUE(all.equal(as.numeric(unlist(cfg$global$resolutions)), seq(.1, 1, .1), tolerance = 1e-10)))
    stop("Expected clustering resolutions from 0.1 to 1.0.")
  if (!identical(as.character(unlist(cfg$global$regress_vars)), c("nCount_RNA", "percent.mt")))
    stop("Expected nCount_RNA and raw percent.mt regression only.")
  Sys.setenv(OMP_NUM_THREADS = as.character(threads), OMP_THREAD_LIMIT = as.character(threads),
    OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = as.character(threads), CUDA_VISIBLE_DEVICES = "")
  required <- c("Seurat", "SeuratObject", "Matrix", "harmony", "future", "uwot", "jsonlite")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing runtime packages: ", paste(missing, collapse = ", "))
  future::plan(future::sequential)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1L)
  options(future.globals.maxSize = 20 * 1024^3)
  out <- file.path(root, "03_global")
  for (d in c(out, file.path(out, "source_data"), file.path(out, "checkpoints"), file.path(root, "provenance")))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  done <- file.path(out, "complete.json")
  if (file.exists(done) || file.exists(file.path(out, "global_celltype_review.rds")))
    stop("Stage 1 already has completed outputs; use a distinct run directory.")
  status <- list(status = "running", started = stage1_now(), stage = "input_validation",
    no_previous_embeddings_or_cluster_maps = TRUE, downstream_started = FALSE)
  status_path <- file.path(out, "status.json")
  update <- function(stage) { status$stage <<- stage; status$updated <<- stage1_now(); stage1_write_json(status, status_path) }
  update("input_validation")
  tryCatch({
    samples <- cfg$samples
    if (is.null(samples) || length(samples) != 9L) stop("Exactly nine sample records are required.")
    sample_ids <- vapply(samples, `[[`, "", "sample_id")
    if (anyDuplicated(sample_ids)) stop("Duplicated sample identities.")
    if (!setequal(sample_ids, source_ids)) stop("Input samples do not match canonical naming map.")
    objs <- vector("list", length(samples)); names(objs) <- sample_ids
    maps <- audits <- vector("list", length(samples))
    rename_audit <- rename_cells <- vector("list", length(samples))
    for (i in seq_along(samples)) {
      s <- samples[[i]]; sid <- s$sample_id
      p <- file.path(root, "02_qc_doublet", sid, "singlet_clean.rds")
      completion <- file.path(root, "02_qc_doublet", sid, "complete.json")
      if (!file.exists(p) || !file.exists(completion)) stop("Incomplete library: ", sid)
      cc <- jsonlite::read_json(completion)
      if (!identical(cc$status, "complete")) stop("Library completion status is not complete: ", sid)
      stage1_log("Loading new singlets:", sid)
      obj <- readRDS(p)
      if (!inherits(obj, "Seurat") || !ncol(obj)) stop("Invalid singlet object: ", sid)
      md <- obj[[]]
      if (!all(c("sample_id", "condition_label", "percent.mt") %in% names(md)))
        stop("Required sample and raw mitochondrial QC metadata missing: ", sid)
      if (!all(as.character(md$sample_id) == sid) ||
          !all(as.character(md$condition_label) == s$condition_label) || anyNA(md$percent.mt))
        stop("Sample or condition or mitochondrial metadata mismatch: ", sid)
      if ("passes_QC" %in% names(md) && !all(md$passes_QC)) stop("Nonpassing cells in singlet input.")
      x <- stage1_get_counts(obj, "RNA"); observed <- stage1_get_counts(obj, "observed")
      if (!identical(dimnames(x), dimnames(observed))) stop("RNA/observed dimnames disagree: ", sid)
      paper_id <- unname(paper_ids[[sid]])
      name_record <- sample_map[[match(sid, source_ids)]]
      if (!identical(name_record$condition_label, s$condition_label)) stop("Canonical condition mismatch: ", sid)
      source_cells <- colnames(obj)
      barcodes <- if ("barcode_10x" %in% names(md)) as.character(md$barcode_10x) else
        sub(paste0("^", sid, "__"), "", source_cells)
      if (anyNA(barcodes) || any(!nzchar(barcodes)) || anyDuplicated(barcodes)) stop("Invalid original barcodes.")
      canonical_cells <- paste(paper_id, barcodes, sep = "__")
      obj <- SeuratObject::RenameCells(obj, new.names = canonical_cells)
      obj$source_cell_id <- source_cells; obj$source_sample_id <- sid
      obj$cell_id <- canonical_cells; obj$sample_id <- paper_id; obj$orig.ident <- paper_id
      renamed_x <- stage1_get_counts(obj, "RNA"); renamed_obs <- stage1_get_counts(obj, "observed")
      unchanged <- identical(x@x, renamed_x@x) && identical(x@i, renamed_x@i) && identical(x@p, renamed_x@p) &&
        identical(observed@x, renamed_obs@x) && identical(observed@i, renamed_obs@i) && identical(observed@p, renamed_obs@p)
      if (!unchanged || !identical(colnames(renamed_x), canonical_cells) ||
          !identical(dimnames(renamed_x), dimnames(renamed_obs))) stop("Naming changed count values/order or assay alignment.")
      rename_cells[[i]] <- data.frame(source_cell_id = source_cells, cell_id = canonical_cells,
        source_sample_id = sid, sample_id = paper_id, barcode_10x = barcodes, stringsAsFactors = FALSE)
      rename_audit[[i]] <- data.frame(source_sample_id = sid, sample_id = paper_id, n_cells = ncol(x),
        corrected_sum_before = sum(x), corrected_sum_after = sum(renamed_x),
        observed_sum_before = sum(observed), observed_sum_after = sum(renamed_obs),
        corrected_nnz_before = length(x@x), corrected_nnz_after = length(renamed_x@x),
        observed_nnz_before = length(observed@x), observed_nnz_after = length(renamed_obs@x),
        sparse_slots_identical = unchanged)
      rm(renamed_x, renamed_obs)
      fm <- obj[["RNA"]][[]]
      if (nrow(fm) != nrow(x)) stop("Feature map length mismatch: ", sid)
      if (!"seurat_feature" %in% names(fm)) fm$seurat_feature <- rownames(x)
      if (!identical(as.character(fm$seurat_feature), rownames(x))) stop("Feature map order mismatch.")
      fm$source_sample_id <- sid
      maps[[i]] <- fm
      # No derived embedding, graph or inherited cluster label is retained.
      obj@reductions <- list(); obj@graphs <- list(); obj@neighbors <- list()
      for (nm in names(md)[grepl("snn_res|cluster|subtype|major_class|celltype|annotation", names(md), ignore.case = TRUE)])
        obj[[nm]] <- NULL
      SeuratObject::DefaultAssay(obj) <- "RNA"
      for (assay in c("RNA", "observed")) {
        if (inherits(obj[[assay]], "Assay5")) {
          layers <- SeuratObject::Layers(obj[[assay]], search = NA)
          for (layer in layers[!grepl("^counts($|\\.)", layers)])
            SeuratObject::LayerData(obj[[assay]], layer = layer) <- NULL
        }
      }
      obj$nCount_RNA <- as.numeric(Matrix::colSums(x))
      obj$nFeature_RNA <- as.numeric(Matrix::colSums(x > 0))
      audits[[i]] <- data.frame(sample_id = paper_id, source_sample_id = sid, condition_label = s$condition_label,
        n_cells = ncol(x), n_features = nrow(x), corrected_total = sum(x), observed_total = sum(observed),
        input_file = p, input_bytes = file.info(p)$size, input_mtime = as.character(file.info(p)$mtime))
      objs[[i]] <- obj
      rm(obj, x, observed); gc(FALSE)
    }
    audit <- do.call(rbind, audits)
    all_cells <- unlist(lapply(objs, colnames), use.names = FALSE)
    if (anyDuplicated(all_cells)) stop("Sample-prefixed cell IDs must be globally unique.")
    fmap <- do.call(rbind, maps)
    if ("gene_id" %in% names(fmap)) {
      by_feature <- split(as.character(fmap$gene_id), fmap$seurat_feature)
      if (any(vapply(by_feature, function(v) length(unique(v)) > 1L, logical(1))))
        stop("A feature maps to different gene IDs across libraries.")
    }
    fmap <- fmap[!duplicated(fmap$seurat_feature), , drop = FALSE]
    update("merge_normalize_hvg")
    obj <- merge(objs[[1L]], y = objs[-1L], merge.data = FALSE, merge.dr = FALSE,
                 project = "manuscript_global_stage1")
    rm(objs, maps); gc(FALSE)
    for (assay in c("RNA", "observed")) {
      if (inherits(obj[[assay]], "Assay5")) obj[[assay]] <- SeuratObject::JoinLayers(obj[[assay]])
      features <- rownames(obj[[assay]])
      idx <- match(features, fmap$seurat_feature)
      if (anyNA(idx)) stop("Merged features missing from explicit feature map.")
      fmd <- fmap[idx, , drop = FALSE]; rownames(fmd) <- NULL
      obj[[assay]]@meta.data <- fmd
    }
    if (ncol(obj) != nrow(obj[[]]) || !setequal(colnames(obj), all_cells)) stop("Cell membership changed during merge.")
    obj <- Seurat::NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    hvg <- SeuratObject::VariableFeatures(obj, assay = "RNA")
    if (length(hvg) != 2000L) stop("Expected exactly 2,000 selected HVGs.")
    update("scale_regress_and_pca")
    obj <- Seurat::ScaleData(obj, assay = "RNA", features = hvg,
      vars.to.regress = c("nCount_RNA", "percent.mt"), block.size = 500, verbose = TRUE)
    obj <- Seurat::RunPCA(obj, assay = "RNA", features = hvg, npcs = 50,
      seed.use = seed, reduction.name = "pca", reduction.key = "PC_", verbose = TRUE)
    pcs <- stage1_select_pcs(SeuratObject::Stdev(obj[["pca"]]), 30L)
    dims <- seq_len(pcs$npcs)
    update("sample_only_harmony")
    obj <- harmony::RunHarmony(obj, group.by.vars = "sample_id", reduction.use = "pca",
      dims.use = dims, reduction.save = "harmony.sample", project.dim = FALSE,
      ncores = threads, max_iter = 10L, plot_convergence = FALSE, verbose = TRUE)
    h <- SeuratObject::Embeddings(obj[["harmony.sample"]])
    if (!identical(rownames(h), colnames(obj)) || ncol(h) != length(dims) || any(!is.finite(h)))
      stop("Invalid Harmony output; no PCA fallback.")
    update("graph_clusters_and_umap")
    obj <- Seurat::FindNeighbors(obj, reduction = "harmony.sample", dims = dims, k.param = 20,
      graph.name = c("global_nn", "global_snn"), verbose = TRUE)
    resolutions <- seq(.1, 1, .1)
    for (res in resolutions) {
      nm <- paste0("global_cluster_res_", sprintf("%.1f", res))
      obj <- Seurat::FindClusters(obj, graph.name = "global_snn", resolution = res,
        algorithm = 1L, cluster.name = nm, random.seed = seed, verbose = TRUE)
    }
    obj$global_cluster <- factor(as.character(obj$global_cluster_res_1.0),
                                 levels = stage1_levels(obj$global_cluster_res_1.0))
    obj$seurat_clusters <- obj$global_cluster
    SeuratObject::Idents(obj) <- "global_cluster"
    # Native uwot exposes explicit neighbor-search and SGD thread budgets.
    # The parameters below match the Seurat UMAP defaults except the explicitly
    # requested manuscript neighbor count, minimum distance and cosine metric.
    umap_args <- list(n_neighbors = 30L, n_components = 2L, metric = "cosine",
      min_dist = .3, spread = 1, set_op_mix_ratio = 1, local_connectivity = 1,
      repulsion_strength = 1, negative_sample_rate = 5, learning_rate = 1,
      n_epochs = NULL, init = "spectral", n_threads = threads, n_sgd_threads = 1L,
      fast_sgd = FALSE, ret_model = FALSE, verbose = TRUE)
    set.seed(seed)
    embedding <- do.call(uwot::umap, c(list(X = h), umap_args))
    if (!is.matrix(embedding) || nrow(embedding) != ncol(obj) || ncol(embedding) != 2L ||
        any(!is.finite(embedding))) stop("Native uwot returned an invalid embedding.")
    rownames(embedding) <- rownames(h); colnames(embedding) <- c("UMAP_1", "UMAP_2")
    obj[["umap.global"]] <- SeuratObject::CreateDimReducObject(embeddings = embedding,
      key = "UMAP_", assay = "RNA", misc = list(engine = "uwot::umap",
        input_reduction = "harmony.sample", input_dims = dims, parameters = umap_args, seed = seed))
    xy <- SeuratObject::Embeddings(obj[["umap.global"]])
    if (!identical(rownames(xy), colnames(obj)) || any(!is.finite(xy))) stop("Invalid fresh UMAP.")
    stage1_atomic_rds(list(pca = obj[["pca"]], harmony = obj[["harmony.sample"]],
      umap = obj[["umap.global"]], graphs = obj@graphs), file.path(out, "checkpoints", "fresh_latent_and_graphs.rds"))
    update("canonical_marker_evidence")
    panels <- stage1_panel(); genes <- unique(c(unlist(panels, use.names = FALSE), stage1_review_genes()))
    available <- intersect(genes, rownames(obj[["RNA"]]))
    availability <- data.frame(gene = genes, available = genes %in% available,
      role = ifelse(genes %in% unlist(panels), "manuscript_major_class", "supplemental_artifact_review"))
    normalized <- SeuratObject::LayerData(obj[["RNA"]], layer = "data", features = available)
    if (!identical(colnames(normalized), colnames(obj))) normalized <- normalized[, colnames(obj), drop = FALSE]
    cl <- as.character(obj$global_cluster); levels <- stage1_levels(cl)
    st <- do.call(rbind, lapply(available, function(g) stage1_gene_stats(as.numeric(normalized[g, ]), g, cl, levels)))
    ann <- stage1_annotation(st, levels)
    md <- obj[[]]
    obj$proposed_major_class <- ann$decisions$proposed_major_class[match(cl, ann$decisions$cluster)]
    obj$annotation_status <- "pending_user_marker_review"
    d <- data.frame(cell_id = colnames(obj), source_cell_id = obj$source_cell_id,
      sample_id = obj$sample_id, source_sample_id = obj$source_sample_id, condition_label = obj$condition_label,
      cluster = cl, proposed_major_class = obj$proposed_major_class,
      UMAP_1 = xy[, 1L], UMAP_2 = xy[, 2L], stringsAsFactors = FALSE, row.names = NULL)
    write.csv(audit, file.path(out, "source_data", "input_singlets_audit.csv"), row.names = FALSE)
    write.csv(do.call(rbind, rename_audit), file.path(out, "source_data", "sample_naming_count_invariance.csv"), row.names = FALSE)
    naming_con <- gzfile(file.path(out, "source_data", "cell_sample_naming_map.csv.gz"), "wt")
    write.csv(do.call(rbind, rename_cells), naming_con, row.names = FALSE); close(naming_con)
    write.csv(availability, file.path(out, "source_data", "marker_availability.csv"), row.names = FALSE)
    write.csv(st, file.path(out, "source_data", "marker_by_cluster.csv"), row.names = FALSE)
    write.csv(ann$decisions, file.path(out, "source_data", "cluster_annotation_proposals.csv"), row.names = FALSE)
    write.csv(ann$evidence, file.path(out, "source_data", "canonical_panel_evidence.csv"), row.names = FALSE)
    write.csv(pcs$table, file.path(out, "source_data", "PCA_variance_profile.csv"), row.names = FALSE)
    write.csv(data.frame(feature = hvg), file.path(out, "source_data", "HVG_2000.csv"), row.names = FALSE)
    write.csv(fmap, file.path(out, "source_data", "gene_feature_union_map.csv"), row.names = FALSE)
    write.csv(do.call(rbind, lapply(resolutions, function(res) data.frame(resolution = res,
      n_clusters = length(unique(obj[[]][[paste0("global_cluster_res_", sprintf("%.1f", res))]])),
      primary = res == 1))), file.path(out, "source_data", "resolution_summary.csv"), row.names = FALSE)
    con <- gzfile(file.path(out, "source_data", "global_cell_annotations_UMAP.csv.gz"), "wt")
    write.csv(d, con, row.names = FALSE); close(con)
    stage1_atomic_rds(list(expression = normalized, coordinates = d, marker_panels = panels,
      marker_stats = st, proposals = ann$decisions), file.path(out, "source_data", "global_marker_review_payload.rds"))
    if (is.null(cfg$run_cluster_marker_discovery) || isTRUE(cfg$run_cluster_marker_discovery)) {
      update("annotation_marker_discovery")
      set.seed(seed)
      cluster_markers <- Seurat::FindAllMarkers(obj, assay = "RNA", slot = "data", test.use = "wilcox",
        only.pos = FALSE, min.pct = .25, logfc.threshold = .25, max.cells.per.ident = 1000,
        random.seed = seed, return.thresh = 1, verbose = TRUE)
      write.csv(cluster_markers, file.path(out, "source_data", "global_cluster_markers.csv"), row.names = FALSE)
    }
    obj@misc$manuscript_stage1 <- list(config_md5 = unname(tools::md5sum(config)), parameters = list(
      NormalizeData_scale_factor = 10000, HVG_method = "vst", HVGs = 2000L,
      regression = c("nCount_RNA corrected UMIs", "percent.mt observed initial-eligible genes"), PCA_computed = 50L,
      PC_selection = pcs[names(pcs) != "table"], Harmony_group = "sample_id", graph_k = 20L,
      clustering_algorithm = "Louvain, Seurat algorithm=1", resolutions = resolutions, primary_resolution = 1,
      UMAP = list(engine = "uwot::umap", neighbors = 30L, min_dist = .3, metric = "cosine",
        seed = seed, thread_limits = c(neighbor_search = threads, SGD = 1L), parameters = umap_args)),
      annotations = "Provisional cluster-level marker review; not atlas probabilities or final user approval",
      stop_after = "global_celltype_review", subclustering_started = FALSE, downstream_started = FALSE)
    obj@misc$manuscript_stage1$script_md5 <- unname(tools::md5sum(file.path(root, "scripts", "02_global_celltype.R")))
    obj@misc$manuscript_stage1$sample_naming_md5 <- unname(tools::md5sum(naming_path))
    obj@misc$manuscript_stage1$sample_naming <- sample_map
    obj@misc$manuscript_stage1$runtime <- list(library_paths = .libPaths(),
      package_versions = stats::setNames(lapply(required, function(p) as.character(utils::packageVersion(p))), required),
      seed = seed, threads = threads, BLAS_threads = 1L)
    stage1_write_json(obj@misc$manuscript_stage1, file.path(out, "source_data", "global_parameters.json"))
    update("saving_global_review_object")
    final <- file.path(out, "global_celltype_review.rds")
    stage1_atomic_rds(obj, final)
    writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
    status$status <- "complete"; status$stage <- "global_celltype_review_ready"
    status$completed <- stage1_now(); status$n_cells <- ncol(obj); status$n_clusters <- length(levels)
    status$n_pcs <- pcs$npcs; status$output_rds <- final; status$awaiting_user_marker_confirmation <- TRUE
    status$subclustering_started <- FALSE; status$downstream_started <- FALSE
    stage1_write_json(status, status_path); stage1_write_json(status, done)
    stage1_log("Stage 1 complete. Stop for user marker confirmation; no subclustering or downstream analysis.")
    invisible(final)
  }, error = function(e) {
    status$status <- "failed"; status$error <- conditionMessage(e); status$failed_at <- stage1_now()
    stage1_write_json(status, status_path)
    stop(e)
  })
}

if (sys.nframe() == 0L && Sys.getenv("STAGE1_SOURCE_ONLY") != "1") {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) < 2L || a[[1L]] != "--root") stop("Usage: Rscript 02_global_celltype.R --root NEW_ROOT [--config CONFIG]")
  cp <- if ("--config" %in% a) a[[match("--config", a) + 1L]] else file.path(a[[2L]], "config", "processing_config.json")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required to read configuration.")
  run_stage1_global(a[[2L]], cp)
}
