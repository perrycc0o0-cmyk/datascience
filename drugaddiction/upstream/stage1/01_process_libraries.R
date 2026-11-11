#!/usr/bin/env Rscript
# Original per-library H5 input through fixed-QC singlets only.
# Usage: Rscript --vanilla 01_process_libraries.R --root RUN [--samples Con1,Ext1] [--resume]

args <- commandArgs(trailingOnly = TRUE)
arg <- function(key, default = NULL) {
  i <- match(key, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], "--")) stop("Missing value for ", key)
  args[i + 1L]
}
root_arg <- arg("--root")
if (is.null(root_arg)) stop("--root is required")
root <- normalizePath(root_arg, mustWork = TRUE)
resume <- "--resume" %in% args
cfgfile <- file.path(root, "config", "processing_config.json")
if (!file.exists(cfgfile)) stop("Missing processing_config.json")
suppressPackageStartupMessages(library(jsonlite))
cfg <- read_json(cfgfile, simplifyVector = FALSE)
flat <- function(x) as.character(unlist(x, use.names = FALSE))
.libPaths(unique(c(flat(cfg$runtime$r_libraries), .libPaths())))
Sys.setenv(OMP_NUM_THREADS = as.character(cfg$threads), OMP_THREAD_LIMIT = as.character(cfg$threads),
           OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", CUDA_VISIBLE_DEVICES = "",
           RCPP_PARALLEL_NUM_THREADS = as.character(cfg$threads))
if (!is.null(cfg$runtime$fontconfig_file)) Sys.setenv(FONTCONFIG_FILE = cfg$runtime$fontconfig_file)
options(stringsAsFactors = FALSE, warn = 1, Seurat.object.assay.version = "v5")
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(hdf5r)
  library(SoupX); library(scDblFinder); library(SingleCellExperiment); library(ggplot2)
})
future::plan(future::sequential)
xgboost::xgb.set.config(nthread = as.integer(cfg$doublet$xgb_nthreads))
stopifnot(xgboost::xgb.get.config()$nthread == cfg$doublet$xgb_nthreads)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1L)
stopifnot(isFALSE(cfg$qc$MAD_filters), isFALSE(cfg$ambient$force_accept),
          isTRUE(cfg$ambient$round_to_integer), isTRUE(cfg$doublet$per_capture),
          identical(cfg$stage_gate, "global_celltype_review"),
          identical(cfg$qc$metric_assay, "observed original UMI counts"))
if (!requireNamespace("ragg", quietly = TRUE)) stop("ragg is required for QC PNGs")
font_match <- system2("fc-match", c("-f", shQuote("%{family}"), "Arial"), stdout = TRUE)
if (!length(font_match) || !grepl("^Arial([, ]|$)", font_match[1])) stop("Arial font resolution failed")

for (d in c("01_ambient", "02_qc_doublet", "logs", "provenance"))
  dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
scriptfile <- file.path(root, "scripts", "01_process_libraries.R")
stopifnot(file.exists(scriptfile))
cfgmd5 <- unname(tools::md5sum(cfgfile)); scriptmd5 <- unname(tools::md5sum(scriptfile))
now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
logmsg <- function(...) { cat(now(), ..., "\n"); flush.console() }
writej <- function(x, path) write_json(x, path, auto_unbox = TRUE, pretty = TRUE,
                                      null = "null", na = "null", digits = 17)
atomic_rds <- function(x, path) {
  if (file.exists(path) || file.exists(paste0(path, ".incomplete")))
    stop("Refusing to overwrite RDS/checkpoint: ", path)
  tmp <- paste0(path, ".incomplete"); saveRDS(x, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) stop("Checkpoint rename failed: ", path)
}
writegz <- function(x, path) {
  con <- gzfile(path, "wt"); on.exit(close(con)); write.csv(x, con, row.names = FALSE, na = "NA")
}
status <- function(state, sample = NULL, stage = NULL, details = NULL) {
  writej(list(status = state, updated_UTC = now(), sample = sample, stage = stage,
              details = details, pid = Sys.getpid(), config_md5 = cfgmd5, script_md5 = scriptmd5),
         file.path(root, "logs", "01_process_libraries.status.json"))
}

read_h5 <- function(path) {
  if (!file.exists(path)) stop("Missing original H5: ", path)
  h <- H5File$new(path, "r"); on.exit(h$close_all())
  ids <- as.character(h[["matrix/features/id"]][])
  symbols <- as.character(h[["matrix/features/name"]][])
  bc <- as.character(h[["matrix/barcodes"]][]); shape <- as.integer(h[["matrix/shape"]][])
  if (h$exists("matrix/features/feature_type")) {
    ft <- as.character(h[["matrix/features/feature_type"]][])
    if (any(ft != "Gene Expression")) stop("Non-GEX features require explicit modality selection")
  }
  stopifnot(length(shape) == 2L, shape[1] == length(ids), shape[2] == length(bc),
            length(symbols) == length(ids), !anyDuplicated(ids), !anyDuplicated(bc))
  counts <- sparseMatrix(i = as.integer(h[["matrix/indices"]][]) + 1L,
                         p = as.integer(h[["matrix/indptr"]][]),
                         x = as.numeric(h[["matrix/data"]][]), dims = shape,
                         dimnames = list(ids, bc), repr = "C")
  stopifnot(all(is.finite(counts@x)), all(counts@x >= 0), all(counts@x == round(counts@x)))
  list(counts = counts, ids = ids, symbols = symbols)
}
qcstats <- function(m, mt_ids, ribo_ids) {
  n <- as.numeric(Matrix::colSums(m)); g <- as.numeric(Matrix::colSums(m > 0))
  fraction <- function(ids) {
    ids <- intersect(ids, rownames(m))
    numerator <- if (length(ids)) as.numeric(Matrix::colSums(m[ids, , drop = FALSE])) else rep(0, ncol(m))
    ifelse(n > 0, 100 * numerator / n, NA_real_)
  }
  complexity <- rep(NA_real_, length(n)); ok <- n > 1 & g > 0
  complexity[ok] <- log10(g[ok]) / log10(n[ok])
  data.frame(nCount = n, nFeature = g, percent_mt = fraction(mt_ids),
             percent_ribo = fraction(ribo_ids), log10GenesPerUMI = complexity,
             row.names = colnames(m))
}
prefixed <- function(x, suffix) { names(x) <- paste0(names(x), suffix); x }

plot_qc <- function(q, sample, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  th <- theme_classic(base_size = 12, base_family = cfg$runtime$font_family) +
    theme(plot.title = element_text(face = "bold"), plot.margin = margin(8, 12, 8, 8))
  export <- function(p, stem) {
    ggsave(file.path(outdir, paste0(stem, ".png")), p, width = 7.5, height = 5.4,
           dpi = 300, device = ragg::agg_png, bg = "white")
    ggsave(file.path(outdir, paste0(stem, ".pdf")), p, width = 7.5, height = 5.4,
           device = grDevices::cairo_pdf, family = cfg$runtime$font_family, bg = "white")
  }
  rules <- list(
    list("nFeature_observed", "Detected genes", c(cfg$qc$nFeature_lower_exclusive, cfg$qc$nFeature_upper_exclusive), "A_observed_genes"),
    list("nCount_observed", "Original UMIs", c(cfg$qc$nCount_lower_exclusive, cfg$qc$nCount_upper_exclusive), "B_observed_UMI"),
    list("percent_mt_observed", "Original mitochondrial UMIs (%)", cfg$qc$percent_mt_upper_exclusive, "C_observed_mt"),
    list("log10GenesPerUMI_observed", "log10(genes) / log10(UMIs)", cfg$qc$log10GenesPerUMI_lower_exclusive, "D_observed_complexity"))
  for (r in rules) {
    z <- data.frame(value = q[[r[[1]]]]); z <- z[is.finite(z$value), , drop = FALSE]
    p <- ggplot(z, aes(value)) + geom_histogram(bins = 80, fill = "#3B6FBF", colour = "white", linewidth = .1) +
      geom_vline(xintercept = r[[3]], linetype = "dashed", colour = "#DC143C") +
      labs(title = paste(sample, r[[2]]), subtitle = "Observed counts; initial-qualified gene set", x = r[[2]], y = "Called barcodes") + th
    if (r[[1]] %in% c("nFeature_observed", "nCount_observed")) p <- p + scale_x_log10()
    export(p, r[[4]])
  }
  z <- q[q$nCount_observed > 0 & q$nFeature_observed > 0, , drop = FALSE]
  p <- ggplot(z, aes(nCount_observed, nFeature_observed, colour = final_call)) +
    geom_point(size = .35, alpha = .5) + scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c(singlet = "#3B6FBF", doublet = "#E85D5D", QC_failed = "#BDBDBD")) +
    labs(title = paste(sample, "fixed QC and doublet calls"), x = "Original UMIs", y = "Original detected genes", colour = NULL) + th
  export(p, "E_QC_doublet_calls")
}

samples <- cfg$samples
requested <- arg("--samples")
if (!is.null(requested)) {
  requested <- strsplit(requested, ",", fixed = TRUE)[[1L]]
  if (any(!requested %in% vapply(samples, `[[`, "", "sample_id"))) stop("Unknown requested sample")
  samples <- Filter(function(s) s$sample_id %in% requested, samples)
}
stopifnot(length(samples) > 0L)
mt_ids <- flat(cfg$qc$mitochondrial_gene_ids); mt_symbols <- flat(cfg$qc$mitochondrial_symbols)

run_one <- function(s) {
  sample <- s$sample_id
  seed <- as.integer(cfg$analysis_seed + match(sample, vapply(cfg$samples, `[[`, "", "sample_id")))
  ad <- file.path(root, "01_ambient", sample); qd <- file.path(root, "02_qc_doublet", sample)
  dir.create(ad, recursive = TRUE, showWarnings = FALSE); dir.create(qd, recursive = TRUE, showWarnings = FALSE)
  fp <- file.path(cfg$input_root, sample, "filtered_feature_bc_matrix.h5")
  rp <- file.path(cfg$input_root, sample, "raw_feature_bc_matrix.h5")
  if (!all(file.exists(c(fp, rp)))) stop("Missing original library H5s: ", sample)
  info_before <- file.info(c(fp, rp))[, c("size", "mtime"), drop = FALSE]
  input_md5 <- as.list(tools::md5sum(c(fp, rp)))
  identity <- list(config_md5 = cfgmd5, script_md5 = scriptmd5, input_md5 = input_md5)
  check_identity <- function(x) {
    stopifnot(identical(x$config_md5, cfgmd5), identical(x$script_md5, scriptmd5), identical(x$input_md5, input_md5))
  }
  done <- file.path(qd, "complete.json")
  if (file.exists(done)) {
    if (!resume) stop("Completed sample exists; explicit --resume is required: ", sample)
    d <- read_json(done); check_identity(d)
    stopifnot(file.exists(file.path(qd, "singlet_clean.rds")))
    logmsg(sample, "same-identity completed checkpoint retained"); return(invisible(NULL))
  }
  manifest <- data.frame(sample_id = sample, condition_label = s$condition_label,
    input_kind = c("filtered_called", "raw_droplets"), path = c(fp, rp), bytes = info_before$size,
    mtime_UTC = format(info_before$mtime, tz = "UTC", usetz = TRUE), md5 = unlist(input_md5, use.names = FALSE))
  write.csv(manifest, file.path(qd, "input_manifest.csv"), row.names = FALSE)
  set.seed(seed); status("running", sample, "read_original"); logmsg(sample, "reading original filtered H5")
  obs <- read_h5(fp); n0 <- ncol(obs$counts)
  stopifnot(all(mt_ids %in% obs$ids), setequal(obs$ids[obs$symbols %in% mt_symbols], mt_ids))
  symbols <- obs$symbols; empty <- is.na(symbols) | !nzchar(symbols); symbols[empty] <- obs$ids[empty]
  featuremap <- data.frame(gene_id = obs$ids, original_symbol = obs$symbols,
    seurat_feature = make.unique(gsub("_", "-", symbols)),
    duplicate_symbol = duplicated(symbols) | duplicated(symbols, fromLast = TRUE))
  ribo_ids <- obs$ids[grepl(cfg$qc$ribo_symbol_pattern, symbols)]
  full_raw <- qcstats(obs$counts, mt_ids, ribo_ids)

  # Native Seurat determines initial qualification; no fixed-QC value is computed on corrected counts.
  initial <- CreateSeuratObject(counts = obs$counts, project = sample,
    min.cells = as.integer(cfg$initial_eligibility$min_cells_per_gene),
    min.features = as.integer(cfg$initial_eligibility$min_features_per_cell))
  eligible_genes <- rownames(initial); eligible_cells <- colnames(initial)
  stopifnot(length(eligible_genes) > 2000L, length(eligible_cells) > 1000L,
    all(eligible_genes %in% obs$ids), all(eligible_cells %in% colnames(obs$counts)))
  featuremap$passes_initial_min_cells <- featuremap$gene_id %in% eligible_genes
  featuremap$n_detecting_called_barcodes <- as.numeric(Matrix::rowSums(obs$counts > 0))
  featuremap$n_detecting_initial_cells <- as.numeric(Matrix::rowSums(obs$counts[, eligible_cells, drop = FALSE] > 0))
  write.csv(featuremap, file.path(ad, "feature_mapping.csv"), row.names = FALSE)
  mt_eligible <- intersect(mt_ids, eligible_genes)
  if (!length(mt_eligible)) stop("No reference mitochondrial genes survive initial qualification")
  write.csv(featuremap[match(mt_ids, featuremap$gene_id), , drop = FALSE],
            file.path(qd, "mitochondrial_gene_scope_audit.csv"), row.names = FALSE)
  rawstats <- qcstats(obs$counts[eligible_genes, , drop = FALSE], mt_ids, ribo_ids)
  stopifnot(all(rawstats[eligible_cells, "nCount"] == initial$nCount_RNA),
            all(rawstats[eligible_cells, "nFeature"] == initial$nFeature_RNA))
  rm(initial); gc(FALSE)

  ambdone <- file.path(ad, "complete.json"); acpath <- file.path(ad, "corrected_counts_gene_id.rds")
  if (file.exists(ambdone)) {
    if (!resume) stop("Ambient checkpoint exists; explicit --resume is required: ", sample)
    am <- read_json(ambdone); check_identity(am)
    corrected <- readRDS(acpath); rho <- am$estimated_rho
    stopifnot(identical(dimnames(corrected), dimnames(obs$counts)))
  } else {
    status("running", sample, "SoupX_precluster")
    use <- full_raw$nFeature >= cfg$ambient$precluster_genes_min & full_raw$nCount >= cfg$ambient$precluster_umi_min
    stopifnot(sum(use) > 1000L)
    pre <- CreateSeuratObject(obs$counts[, use, drop = FALSE], min.cells = 0, min.features = 0, project = sample)
    pre <- NormalizeData(pre, verbose = FALSE)
    pre <- FindVariableFeatures(pre, nfeatures = as.integer(cfg$ambient$precluster_n_hvg), verbose = FALSE)
    VariableFeatures(pre) <- setdiff(VariableFeatures(pre), mt_ids)
    pre <- ScaleData(pre, features = VariableFeatures(pre), verbose = FALSE)
    pre <- RunPCA(pre, npcs = as.integer(cfg$ambient$precluster_n_pcs), features = VariableFeatures(pre), seed.use = seed, verbose = FALSE)
    pre <- FindNeighbors(pre, dims = seq_len(cfg$ambient$precluster_neighbor_dims), verbose = FALSE)
    pre <- FindClusters(pre, resolution = cfg$ambient$precluster_resolution, algorithm = 1, random.seed = seed, verbose = FALSE)
    cl <- setNames(rep("low_quality", n0), colnames(obs$counts)); cl[colnames(pre)] <- as.character(Idents(pre))
    write.csv(data.frame(barcode = names(cl), SoupX_precluster = unname(cl)), file.path(ad, "preliminary_clusters.csv"), row.names = FALSE)
    rm(pre); gc(FALSE)
    status("running", sample, "SoupX_autoEstCont"); logmsg(sample, "fresh raw-droplet SoupX estimation")
    raw <- read_h5(rp)
    stopifnot(identical(raw$ids, obs$ids), identical(raw$symbols, obs$symbols), all(colnames(obs$counts) %in% colnames(raw$counts)))
    sc <- SoupChannel(tod = raw$counts, toc = obs$counts); sc <- setClusters(sc, cl)
    rm(raw); gc(FALSE)
    sc <- autoEstCont(sc, doPlot = FALSE, forceAccept = FALSE, verbose = TRUE)
    rho <- unique(sc$metaData$rho)
    stopifnot(length(rho) == 1L, is.finite(rho), rho > 0)
    writej(c(identity, list(sample_id = sample, estimated_rho = rho, updated_UTC = now())), file.path(ad, "rho_estimate.json"))
    if (rho > cfg$ambient$rho_review_max) stop("SoupX rho exceeds review boundary; inspect before accepting: ", rho)
    write.csv(sc$soupProfile, file.path(ad, "soup_profile.csv"), row.names = TRUE)
    if (!is.null(sc$fit)) atomic_rds(sc$fit, file.path(ad, "autoEstCont_fit.rds"))
    status("running", sample, "SoupX_adjustCounts"); set.seed(seed)
    corrected <- adjustCounts(sc, roundToInt = TRUE, verbose = 1)
    stopifnot(identical(dimnames(corrected), dimnames(obs$counts)), all(is.finite(corrected@x)),
      all(corrected@x >= 0), all(corrected@x == round(corrected@x)))
    atomic_rds(corrected, acpath)
    write.csv(data.frame(featuremap, observed_UMI = as.numeric(Matrix::rowSums(obs$counts)),
      corrected_UMI = as.numeric(Matrix::rowSums(corrected)),
      removed_UMI = as.numeric(Matrix::rowSums(obs$counts) - Matrix::rowSums(corrected))),
      file.path(ad, "gene_UMI_correction.csv"), row.names = FALSE)
    writej(c(identity, list(status = "complete", sample_id = sample, estimated_rho = rho,
      method = "SoupX", version = as.character(packageVersion("SoupX")), completed_UTC = now())), ambdone)
    rm(sc); gc(FALSE)
  }
  status("running", sample, "observed_fixed_QC")
  cs <- qcstats(corrected[eligible_genes, , drop = FALSE], mt_ids, ribo_ids)
  full_corrected <- qcstats(corrected, mt_ids, ribo_ids)
  qc <- data.frame(cell_id = paste(sample, colnames(obs$counts), sep = "__"), sample_id = sample,
    condition_label = s$condition_label, barcode_10x = colnames(obs$counts),
    passes_initial_eligibility = colnames(obs$counts) %in% eligible_cells,
    prefixed(rawstats, "_observed"), prefixed(cs, "_corrected"),
    prefixed(full_raw, "_fullGEX_observed"), prefixed(full_corrected, "_fullGEX_corrected"),
    estimated_SoupX_rho = rho, row.names = NULL)
  qc$removed_UMI <- qc$nCount_observed - qc$nCount_corrected
  qc$removed_UMI_fraction <- ifelse(qc$nCount_observed > 0, qc$removed_UMI / qc$nCount_observed, NA_real_)
  qc$fail_nonfinite_metric <- !is.finite(qc$nFeature_observed) | !is.finite(qc$nCount_observed) |
    !is.finite(qc$percent_mt_observed) | !is.finite(qc$log10GenesPerUMI_observed)
  qc$fail_low_genes <- qc$nFeature_observed <= cfg$qc$nFeature_lower_exclusive
  qc$fail_high_genes <- qc$nFeature_observed >= cfg$qc$nFeature_upper_exclusive
  qc$fail_low_UMI <- qc$nCount_observed <= cfg$qc$nCount_lower_exclusive
  qc$fail_high_UMI <- qc$nCount_observed >= cfg$qc$nCount_upper_exclusive
  qc$fail_high_observed_mt <- !is.finite(qc$percent_mt_observed) | qc$percent_mt_observed >= cfg$qc$percent_mt_upper_exclusive
  qc$fail_low_complexity <- !is.finite(qc$log10GenesPerUMI_observed) |
    qc$log10GenesPerUMI_observed <= cfg$qc$log10GenesPerUMI_lower_exclusive
  fail_cols <- grep("^fail_", names(qc), value = TRUE)
  stopifnot(!anyNA(qc[fail_cols]))
  qc$passes_core_QC <- qc$passes_initial_eligibility & !apply(qc[fail_cols], 1L, any)
  if (sum(qc$passes_core_QC) <= 1000L) stop("Too few nuclei pass fixed QC; inspect rather than change thresholds")
  qc$scDblFinder_score <- NA_real_; qc$scDblFinder_class <- NA_character_
  qc$doublet_partition <- ifelse(qc$passes_core_QC, sample, NA_character_)
  writegz(qc, file.path(qd, "all_called_barcodes_QC_before_doublet.csv.gz"))
  writej(list(primary_metric_assay = cfg$qc$metric_assay, primary_gene_scope = cfg$qc$metric_gene_scope,
    applied = cfg$qc, n_initial_genes = length(eligible_genes), n_initial_cells = length(eligible_cells),
    n_reference_mt_genes = length(mt_ids), n_initial_mt_genes = length(mt_eligible),
    mt_ids_not_initial_qualified = setdiff(mt_ids, mt_eligible), no_MAD_filter = TRUE,
    no_corrected_metric_filter = TRUE), file.path(qd, "applied_thresholds.json"))
  dp <- file.path(qd, "doublet_calls.rds")
  if (file.exists(dp)) {
    if (!resume) stop("Doublet checkpoint exists; explicit --resume is required: ", sample)
    dc <- readRDS(dp); check_identity(dc)
    stopifnot(identical(dc$barcodes, qc$barcode_10x[qc$passes_core_QC]))
  } else {
    status("running", sample, "scDblFinder"); logmsg(sample, "scDblFinder on", sum(qc$passes_core_QC), "fixed-QC nuclei")
    sce <- SingleCellExperiment(assays = list(counts = corrected[eligible_genes, qc$passes_core_QC, drop = FALSE]))
    expected_dbr <- cfg$doublet$dbr_per_1000 * n0 / 1000
    stopifnot(expected_dbr > 0, expected_dbr < .2)
    set.seed(seed)
    sce <- scDblFinder(sce, clusters = FALSE, dbr = expected_dbr, dbr.per1k = cfg$doublet$dbr_per_1000,
      BPPARAM = BiocParallel::SerialParam(), verbose = TRUE)
    stopifnot(identical(colnames(sce), qc$barcode_10x[qc$passes_core_QC]))
    dc <- c(identity, list(barcodes = colnames(sce), score = as.numeric(sce$scDblFinder.score),
      class = as.character(sce$scDblFinder.class), expected_dbr = expected_dbr,
      original_called_denominator = n0, BPPARAM = "SerialParam", seed = seed,
      xgboost_global_config = xgboost::xgb.get.config(), completed_UTC = now()))
    stopifnot(all(is.finite(dc$score)), all(dc$class %in% c("singlet", "doublet")))
    atomic_rds(dc, dp); rm(sce); gc(FALSE)
  }
  qc$scDblFinder_score[qc$passes_core_QC] <- dc$score
  qc$scDblFinder_class[qc$passes_core_QC] <- dc$class
  qc$passes_QC <- qc$passes_core_QC & !is.na(qc$scDblFinder_class) & qc$scDblFinder_class == "singlet"
  qc$final_call <- ifelse(!qc$passes_core_QC, "QC_failed", qc$scDblFinder_class)
  stopifnot(!anyNA(qc$passes_QC), sum(qc$passes_QC) > 0L,
    sum(!qc$passes_core_QC) + sum(dc$class == "doublet") + sum(qc$passes_QC) == n0)
  qc$percent.mt <- qc$percent_mt_observed
  qc$percent.ribo <- qc$percent_ribo_observed
  qc$log10GenesPerUMI <- qc$log10GenesPerUMI_observed
  writegz(qc, file.path(qd, "all_called_barcodes_QC.csv.gz"))
  plot_qc(qc, sample, file.path(qd, "QC_after"))
  status("running", sample, "save_singlets")
  keep <- which(qc$passes_QC)
  cor <- Matrix::drop0(corrected[eligible_genes, keep, drop = FALSE])
  rawkeep <- Matrix::drop0(obs$counts[eligible_genes, keep, drop = FALSE])
  fmap <- featuremap[match(eligible_genes, featuremap$gene_id), , drop = FALSE]
  dimnames(cor) <- list(fmap$seurat_feature, qc$cell_id[keep]); dimnames(rawkeep) <- dimnames(cor)
  meta <- qc[keep, , drop = FALSE]; rownames(meta) <- meta$cell_id
  obj <- CreateSeuratObject(counts = cor, min.cells = 0, min.features = 0, project = sample, meta.data = meta)
  obj[["observed"]] <- CreateAssay5Object(counts = rawkeep)
  # Assay5 uses positional feature-metadata row names; order is checked against explicit feature names.
  rownames(fmap) <- NULL
  obj[["RNA"]]@meta.data <- fmap; obj[["observed"]]@meta.data <- fmap
  stopifnot(isTRUE(validObject(obj)), identical(rownames(obj[["RNA"]]), fmap$seurat_feature),
    identical(rownames(obj[["RNA"]]), rownames(obj[["observed"]])), identical(colnames(obj), rownames(obj[[]])))
  obj@misc$processing <- c(identity, list(sample_id = sample, seed = seed, original_H5 = c(fp, rp),
    primary_assay = cfg$output_interface$primary_assay, audit_assay = cfg$output_interface$audit_assay,
    QC_metric_gene_scope = cfg$qc$metric_gene_scope, QC_metric_assay = cfg$qc$metric_assay,
    initial_min_cells = cfg$initial_eligibility$min_cells_per_gene,
    initial_min_features = cfg$initial_eligibility$min_features_per_cell,
    no_old_embeddings_or_labels = TRUE, no_normalization_or_integration_yet = TRUE))
  atomic_rds(obj, file.path(qd, "singlet_clean.rds"))
  summary <- data.frame(sample_id = sample, condition_label = s$condition_label, CellRanger_called = n0,
    initial_eligible_cells = length(eligible_cells), initial_eligible_genes = length(eligible_genes),
    initial_ineligible_cells = sum(!qc$passes_initial_eligibility), SoupX_rho = rho,
    min_genes_exclusive = cfg$qc$nFeature_lower_exclusive, max_genes_exclusive = cfg$qc$nFeature_upper_exclusive,
    min_UMI_exclusive = cfg$qc$nCount_lower_exclusive, max_UMI_exclusive = cfg$qc$nCount_upper_exclusive,
    max_observed_mt_exclusive = cfg$qc$percent_mt_upper_exclusive,
    min_observed_complexity_exclusive = cfg$qc$log10GenesPerUMI_lower_exclusive,
    core_QC_failed = sum(!qc$passes_core_QC), core_QC_passed = sum(qc$passes_core_QC),
    scDblFinder_doublets = sum(dc$class == "doublet"), final_singlets = sum(qc$passes_QC),
    expected_dbr = dc$expected_dbr, n_reference_mt_genes = length(mt_ids), n_initial_mt_genes = length(mt_eligible))
  for (nm in fail_cols) summary[[paste0(nm, "_nonexclusive")]] <- sum(qc[[nm]])
  write.csv(summary, file.path(qd, "sample_summary.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(qd, "sessionInfo.txt"))
  stopifnot(identical(info_before, file.info(c(fp, rp))[, c("size", "mtime"), drop = FALSE]))
  writej(c(identity, list(status = "complete", sample_id = sample, completed_UTC = now(), summary = summary,
    no_old_embeddings_or_labels = TRUE, no_old_results_overwritten = TRUE)), done)
  logmsg(sample, "COMPLETE:", sum(qc$passes_QC), "singlets;", sum(dc$class == "doublet"), "predicted doublets")
  rm(obj, obs, corrected, cor, rawkeep); gc(FALSE)
}

tryCatch({
  status("running", stage = "original_libraries")
  for (s in samples) run_one(s)
  tabs <- lapply(cfg$samples, function(s) {
    p <- file.path(root, "02_qc_doublet", s$sample_id, "sample_summary.csv")
    if (file.exists(p)) read.csv(p, check.names = FALSE) else NULL
  })
  tabs <- Filter(Negate(is.null), tabs)
  write.csv(do.call(rbind, tabs), file.path(root, "02_qc_doublet", "attrition_summary.csv"), row.names = FALSE)
  status("complete", stage = "requested_original_libraries_completed")
  writej(list(status = "complete", config_md5 = cfgmd5, script_md5 = scriptmd5,
    requested_samples = vapply(samples, `[[`, "", "sample_id"),
    all_nine_libraries_complete = all(vapply(cfg$samples, function(s) file.exists(file.path(root, "02_qc_doublet", s$sample_id, "complete.json")), logical(1))),
    completed_UTC = now(), stop_before_global = TRUE, no_subclustering_or_downstream_analysis = TRUE),
    file.path(root, "provenance", "01_libraries_complete.json"))
  logmsg("Requested captures completed; global/lineage/downstream analyses were not executed")
}, error = function(e) {
  status("failed", stage = "inspect_log", details = conditionMessage(e))
  logmsg("FAILED:", conditionMessage(e)); quit(status = 1)
})
