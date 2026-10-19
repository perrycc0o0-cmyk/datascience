suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SoupX)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(ggplot2)
})
source('R/common.R')
source('R/preprocessing.R')

ctx <- pipeline_context()
cfg <- ctx$cfg
out <- new_stage(ctx, '01_preprocessing')
samples <- read_samples(ctx)
require_version('SoupX', cfg$ambient$version)
require_version('scDblFinder', cfg$doublet$version)
summaries <- vector('list', nrow(samples))
objects <- character(nrow(samples))

for (i in seq_len(nrow(samples))) {
  s <- samples[i, ]
  message('Processing ', s$sample)
  filtered <- read_gene_expression_h5(s$filtered_h5)
  raw <- read_gene_expression_h5(s$raw_h5)
  original_called <- ncol(filtered)
  if (!all(rownames(filtered) %in% rownames(raw)) || !all(colnames(filtered) %in% colnames(raw))) {
    stop('Raw and filtered matrices do not have matching features/barcodes: ', s$sample)
  }
  raw <- raw[rownames(filtered), , drop = FALSE]
  obj <- CreateSeuratObject(filtered, project = s$sample,
                            min.cells = cfg$qc$initial_min_cells,
                            min.features = cfg$qc$initial_min_features)
  if (ncol(obj) < 10L) stop('Too few initial nuclei: ', s$sample)
  observed <- rna_matrix(obj)
  mt <- grep(cfg$qc$mt_pattern, rownames(observed), value = TRUE)
  if (!length(mt)) stop('No features match mitochondrial pattern ', cfg$qc$mt_pattern,
                        ' in ', s$sample, '; check the reference gene symbols')
  ribo <- grep(cfg$qc$ribo_pattern, rownames(observed), value = TRUE)
  obj$percent.mt <- as.numeric(Matrix::colSums(observed[mt, , drop = FALSE]) /
                                Matrix::colSums(observed) * 100)
  obj$percent.ribo <- if (length(ribo)) as.numeric(Matrix::colSums(observed[ribo, , drop = FALSE]) /
                                                   Matrix::colSums(observed) * 100) else 0
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  qc <- data.frame(barcode = colnames(obj), sample = s$sample, condition = s$condition,
                   nFeature_observed = obj$nFeature_RNA, nCount_observed = obj$nCount_RNA,
                   percent.mt = obj$percent.mt, percent.ribo = obj$percent.ribo,
                   log10GenesPerUMI = obj$log10GenesPerUMI, stringsAsFactors = FALSE)
  qc$passes_genes <- qc$nFeature_observed > cfg$qc$genes_min & qc$nFeature_observed < cfg$qc$genes_max
  qc$passes_UMI <- qc$nCount_observed > cfg$qc$umi_min & qc$nCount_observed < cfg$qc$umi_max
  qc$passes_mt <- qc$percent.mt < cfg$qc$mt_max
  qc$passes_complexity <- is.finite(qc$log10GenesPerUMI) & qc$log10GenesPerUMI > cfg$qc$complexity_min
  qc$passes_QC <- with(qc, passes_genes & passes_UMI & passes_mt & passes_complexity)

  pre <- NormalizeData(obj, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
  pre <- FindVariableFeatures(pre, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
  pre <- ScaleData(pre, features = VariableFeatures(pre), verbose = FALSE)
  npcs <- min(30L, available_pcs(pre, 30L))
  pre <- RunPCA(pre, features = VariableFeatures(pre), npcs = npcs,
                seed.use = cfg$seed, verbose = FALSE)
  pre <- FindNeighbors(pre, reduction = 'pca', dims = seq_len(npcs), verbose = FALSE)
  pre <- FindClusters(pre, resolution = cfg$ambient$precluster_resolution,
                      random.seed = cfg$seed, verbose = FALSE)
  soup <- SoupChannel(tod = raw, toc = filtered[, colnames(obj), drop = FALSE])
  soup <- setClusters(soup, stats::setNames(as.character(Idents(pre)), colnames(pre)))
  soup <- autoEstCont(soup, forceAccept = cfg$ambient$force_accept)
  corrected <- adjustCounts(soup, roundToInt = cfg$ambient$round_to_integer)
  corrected <- corrected[rownames(observed), colnames(observed), drop = FALSE]
  integer_counts(corrected)
  qc$ambient_fraction <- as.numeric(soup$metaData[qc$barcode, 'rho'])
  save_object(soup$fit, file.path(out, 'ambient_models', paste0(s$sample, '.rds')))
  rm(pre, raw, filtered, soup)
  invisible(gc())

  kept <- qc$barcode[qc$passes_QC]
  if (length(kept) < 10L) stop('Too few nuclei pass fixed QC: ', s$sample)
  sce <- SingleCellExperiment(assays = list(counts = corrected[, kept, drop = FALSE]))
  expected_dbr <- cfg$doublet$dbr_per_1000 * original_called / 1000
  if (!is.finite(expected_dbr) || expected_dbr <= 0 || expected_dbr >= 1) stop('Invalid expected doublet rate')
  set.seed(as.integer(cfg$seed))
  sce <- scDblFinder(sce, clusters = cfg$doublet$clusters, dbr = expected_dbr,
                     verbose = TRUE)
  d <- as.data.frame(SummarizedExperiment::colData(sce))
  qc$scDblFinder.score <- d[match(qc$barcode, rownames(d)), 'scDblFinder.score']
  qc$scDblFinder.class <- as.character(d[match(qc$barcode, rownames(d)), 'scDblFinder.class'])
  qc$passes_doublet <- !is.na(qc$scDblFinder.class) & qc$scDblFinder.class == 'singlet'
  qc$retained <- qc$passes_QC & qc$passes_doublet
  final_cells <- qc$barcode[qc$retained]
  if (!length(final_cells)) stop('No singlets retained: ', s$sample)
  md <- qc[match(final_cells, qc$barcode), , drop = FALSE]
  rownames(md) <- final_cells
  final <- CreateSeuratObject(corrected[, final_cells, drop = FALSE], project = s$sample,
                              min.cells = 0, min.features = 0, meta.data = md)
  final[['observed']] <- CreateAssayObject(counts = observed[, final_cells, drop = FALSE])
  final <- RenameCells(final, new.names = paste(s$sample, final_cells, sep = '_'))
  final$orig.ident <- s$sample
  final$batch <- s$sample
  final$sample <- s$sample
  final$condition <- s$condition
  objects[i] <- file.path(out, 'libraries', paste0(s$sample, '.rds'))
  save_object(final, objects[i])
  write_table(qc, file.path(out, 'QC', paste0(s$sample, '_nuclei.csv')))
  summaries[[i]] <- data.frame(sample = s$sample, condition = s$condition,
    original_called = original_called, after_initial_filter = ncol(observed),
    initial_removed = original_called - ncol(observed),
    fixed_QC_removed = sum(!qc$passes_QC), after_fixed_QC = sum(qc$passes_QC),
    doublets_removed = sum(qc$scDblFinder.class == 'doublet', na.rm = TRUE),
    retained_singlets = length(final_cells), expected_doublet_rate = expected_dbr,
    median_ambient_fraction = stats::median(qc$ambient_fraction))
  long <- rbind(data.frame(value = qc$nFeature_observed, metric = 'Detected genes'),
                data.frame(value = qc$nCount_observed, metric = 'UMIs'),
                data.frame(value = qc$percent.mt, metric = 'Mitochondrial (%)'),
                data.frame(value = qc$log10GenesPerUMI, metric = 'log10 genes / log10 UMIs'))
  lines <- data.frame(metric = c('Detected genes', 'Detected genes', 'UMIs', 'UMIs',
                                 'Mitochondrial (%)', 'log10 genes / log10 UMIs'),
                       threshold = c(cfg$qc$genes_min, cfg$qc$genes_max, cfg$qc$umi_min,
                                     cfg$qc$umi_max, cfg$qc$mt_max, cfg$qc$complexity_min))
  p <- ggplot(long, aes(value)) + geom_histogram(bins = 60, fill = '#7B9CB9', colour = 'white') +
    geom_vline(data = lines, aes(xintercept = threshold), linetype = 'dashed', colour = '#C34A42') +
    facet_wrap(~metric, scales = 'free', ncol = 2) + annotation_theme(ctx) +
    labs(title = s$sample, x = NULL, y = 'Nuclei')
  plot_export(p, file.path(out, 'QC', paste0(s$sample, '_thresholds')), 10, 7, cfg$plot$dpi)
  rm(obj, observed, corrected, sce, final, qc, md, d, long)
  invisible(gc())
}

write_table(do.call(rbind, summaries), file.path(out, 'library_summary.csv'))
singlets <- lapply(objects, readRDS)
merged <- merge(singlets[[1L]], y = singlets[-1L], merge.data = FALSE)
merged <- join_rna_layers(merged)
save_object(merged, file.path(out, 'merged.rds'))
write_table(data.frame(sample = samples$sample, object = objects), file.path(out, 'objects.csv'))
write_session(file.path(out, 'sessionInfo.txt'))
