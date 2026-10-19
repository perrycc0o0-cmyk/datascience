suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
source('R/common.R')

ctx <- pipeline_context()
cfg <- ctx$cfg
out <- new_stage(ctx, '04_ieg')
samples <- read_samples(ctx)
lineages <- c('Excitatory', 'Inhibitory')
ieg_genes <- as.character(cfg$ieg$genes)
cutoffs <- list()
metadata <- list()

for (lineage in lineages) {
  obj <- readRDS(neuron_object_path(ctx, lineage))
  require_neuron_metadata(obj)
  DefaultAssay(obj) <- 'RNA'
  obj <- join_rna_layers(obj)
  integer_counts(rna_matrix(obj, 'counts'))
  obj <- NormalizeData(obj, assay = 'RNA', normalization.method = 'LogNormalize',
                       scale.factor = as.numeric(cfg$global$scale_factor), verbose = FALSE)
  available_ieg <- intersect(ieg_genes, rownames(obj[['RNA']]))
  write_table(data.frame(neuron_type = lineage, gene = ieg_genes,
                         available = ieg_genes %in% available_ieg),
              file.path(out, paste0(lineage, '_IEG_availability.csv')))
  if (!length(available_ieg)) stop('None of the prespecified IEG genes is available in ', lineage)
  meta <- obj[[]]
  sample_index <- match(as.character(meta$sample), samples$sample)
  if (anyNA(sample_index) || any(as.character(meta$condition) != samples$condition[sample_index])) {
    stop('Sample and condition metadata do not match config/samples.csv')
  }
  expression <- rna_matrix(obj, 'data')
  ieg_sum <- as.numeric(Matrix::colSums(expression[available_ieg, , drop = FALSE]))
  cutoff <- as.numeric(stats::quantile(ieg_sum, probs = as.numeric(cfg$ieg$quantile),
                                       type = as.integer(cfg$ieg$quantile_type), names = FALSE))
  obj$neuron_type <- lineage
  obj$ieg_sum <- ieg_sum
  obj$IEGhigh <- ieg_sum >= cutoff
  obj$IEGclass <- ifelse(obj$IEGhigh, 'IEGhigh', 'IEGlow')
  m <- obj[[]]
  metadata[[lineage]] <- data.frame(
    Cell = colnames(obj), sample = as.character(m$sample), condition = as.character(m$condition),
    neuron_type = lineage, subtype = as.character(m$subtype), ieg_sum = m$ieg_sum,
    IEGhigh = m$IEGhigh, IEGclass = m$IEGclass, stringsAsFactors = FALSE
  )
  cutoffs[[lineage]] <- data.frame(
    neuron_type = lineage, quantile = as.numeric(cfg$ieg$quantile),
    quantile_type = as.integer(cfg$ieg$quantile_type), cutoff = cutoff,
    n_IEG_genes = length(available_ieg), IEG_genes = paste(available_ieg, collapse = ';'),
    n_nuclei = ncol(obj), n_IEGhigh = sum(obj$IEGhigh),
    fraction_IEGhigh = mean(obj$IEGhigh), stringsAsFactors = FALSE
  )
  save_object(obj, file.path(out, paste0(lineage, '.rds')))
}

metadata <- do.call(rbind, metadata)
rownames(metadata) <- NULL
write_table(metadata, file.path(out, 'nucleus_metadata.csv'))
write_table(do.call(rbind, cutoffs), file.path(out, 'IEG_cutoffs.csv'))
summary <- as.data.frame(with(metadata, table(neuron_type, condition, IEGclass)), stringsAsFactors = FALSE)
names(summary)[ncol(summary)] <- 'n_nuclei'
write_table(summary, file.path(out, 'IEG_counts.csv'))
write_session(file.path(out, 'sessionInfo.txt'))
