suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
})
source('R/common.R')

ctx <- pipeline_context()
cfg <- ctx$cfg
out <- new_stage(ctx, '05_DESeq2')
samples <- read_samples(ctx)
conditions <- c('WD', 'Ext', 'Ret-Ext')
contrasts <- list(Ext_vs_WD = c('Ext', 'WD'), RetExt_vs_WD = c('Ret-Ext', 'WD'),
                  RetExt_vs_Ext = c('Ret-Ext', 'Ext'))
threshold <- log2(as.numeric(cfg$de$fold_change))
de_summaries <- list()
filter_summaries <- list()

for (lineage in c('Excitatory', 'Inhibitory')) {
  obj <- readRDS(stage_path(ctx, '04_ieg', paste0(lineage, '.rds')))
  require_neuron_metadata(obj)
  obj <- join_rna_layers(obj)
  meta <- obj[[]]
  if (!'IEGhigh' %in% names(meta) || anyNA(meta$IEGhigh)) stop('Missing IEG-high classification')
  high_cells <- rownames(meta)[meta$IEGhigh]
  if (!length(high_cells)) stop('No IEG-high nuclei in ', lineage)
  meta <- meta[high_cells, , drop = FALSE]
  sample_index <- match(as.character(meta$sample), samples$sample)
  if (anyNA(sample_index) || any(as.character(meta$condition) != samples$condition[sample_index])) {
    stop('Sample and condition metadata do not match config/samples.csv')
  }
  n_nuclei <- tabulate(sample_index, nbins = nrow(samples))
  if (any(n_nuclei == 0L)) stop('Every biological sample must contribute IEG-high nuclei: ', lineage)
  counts <- rna_matrix(obj, 'counts')[, high_cells, drop = FALSE]
  integer_counts(counts)
  membership <- Matrix::sparseMatrix(i = seq_along(sample_index), j = sample_index, x = 1,
                                     dims = c(length(high_cells), nrow(samples)))
  pseudobulk <- as.matrix(counts %*% membership)
  rownames(pseudobulk) <- rownames(counts)
  colnames(pseudobulk) <- samples$sample
  integer_counts(pseudobulk)
  if (any(colSums(pseudobulk) == 0)) stop('Empty pseudobulk library in ', lineage)
  if (any(pseudobulk > .Machine$integer.max)) stop('Pseudobulk count exceeds R integer capacity')
  storage.mode(pseudobulk) <- 'integer'
  sample_meta <- data.frame(sample = samples$sample,
                            condition = factor(samples$condition, levels = conditions),
                            n_nuclei = n_nuclei, row.names = samples$sample)
  keep <- rowSums(pseudobulk) >= as.numeric(cfg$de$min_total_counts) &
    rowSums(pseudobulk > 0) >= as.integer(cfg$de$min_samples)
  if (!any(keep)) stop('No genes pass the configured count filter in ', lineage)
  dds <- DESeqDataSetFromMatrix(countData = pseudobulk[keep, , drop = FALSE],
                                colData = sample_meta, design = stats::as.formula(cfg$de$design))
  dds <- DESeq(dds, test = 'Wald', fitType = 'parametric', sfType = 'ratio',
               minReplicatesForReplace = Inf, parallel = FALSE, quiet = TRUE)
  count_table <- data.frame(gene = rownames(pseudobulk), pseudobulk, check.names = FALSE)
  write_table(count_table, file.path(out, 'counts', paste0(lineage, '_pseudobulk_counts.csv')))
  write_table(sample_meta, file.path(out, 'counts', paste0(lineage, '_sample_metadata.csv')))
  save_object(dds, file.path(out, 'objects', paste0(lineage, '_DESeq2.rds')))
  for (contrast_name in names(contrasts)) {
    comparison <- contrasts[[contrast_name]]
    result <- results(dds, contrast = c('condition', comparison), alpha = as.numeric(cfg$de$padj),
                      pAdjustMethod = 'BH', independentFiltering = TRUE)
    result_table <- as.data.frame(result)
    full <- data.frame(gene = rownames(pseudobulk), tested = keep,
                       total_count = rowSums(pseudobulk),
                       n_samples_expressed = rowSums(pseudobulk > 0), stringsAsFactors = FALSE)
    matched <- match(full$gene, rownames(result_table))
    for (column in names(result_table)) full[[column]] <- result_table[[column]][matched]
    full$neuron_type <- lineage
    full$contrast <- contrast_name
    full$numerator <- comparison[1L]
    full$denominator <- comparison[2L]
    full$significant <- !is.na(full$padj) & full$padj <= as.numeric(cfg$de$padj) &
      !is.na(full$log2FoldChange) & abs(full$log2FoldChange) >= threshold
    full$direction <- ifelse(full$significant,
                              ifelse(full$log2FoldChange > 0, 'Up', 'Down'), 'Not significant')
    table_dir <- file.path(out, 'tables', lineage)
    write_table(full, file.path(table_dir, paste0(contrast_name, '_all_genes.csv')))
    write_table(full[full$significant, , drop = FALSE],
                file.path(table_dir, paste0(contrast_name, '_significant.csv')))
    de_summaries[[paste(lineage, contrast_name)]] <- data.frame(
      neuron_type = lineage, contrast = contrast_name, n_genes = nrow(full),
      n_count_filter = sum(keep), n_with_pvalue = sum(!is.na(full$pvalue)),
      n_with_padj = sum(!is.na(full$padj)), n_up = sum(full$direction == 'Up'),
      n_down = sum(full$direction == 'Down'), padj_cutoff = as.numeric(cfg$de$padj),
      absolute_log2FC_cutoff = threshold
    )
    filtering <- S4Vectors::metadata(result)
    filter_summaries[[paste(lineage, contrast_name)]] <- data.frame(
      neuron_type = lineage, contrast = contrast_name,
      min_total_counts = as.numeric(cfg$de$min_total_counts),
      min_samples_expressed = as.integer(cfg$de$min_samples),
      independent_filtering = TRUE,
      independent_filter_threshold = if (length(filtering$filterThreshold)) as.numeric(filtering$filterThreshold) else NA_real_,
      independent_filter_quantile = if (length(filtering$filterTheta)) as.numeric(filtering$filterTheta) else NA_real_,
      p_adjust_method = 'BH', design = cfg$de$design
    )
  }
}
write_table(do.call(rbind, de_summaries), file.path(out, 'DE_summary.csv'))
write_table(do.call(rbind, filter_summaries), file.path(out, 'filtering_parameters.csv'))
write_session(file.path(out, 'sessionInfo.txt'))
