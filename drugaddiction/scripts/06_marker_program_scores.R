suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
source('R/common.R')

ctx <- pipeline_context()
cfg <- ctx$cfg
out <- new_stage(ctx, '06_scores')
lineages <- c('Excitatory', 'Inhibitory')
objects <- lapply(lineages, function(lineage) {
  obj <- readRDS(stage_path(ctx, '04_ieg', paste0(lineage, '.rds')))
  require_neuron_metadata(obj)
  if (!all(c('IEGhigh', 'IEGclass', 'neuron_type') %in% names(obj[[]]))) {
    stop('Missing IEG-high annotation in ', lineage)
  }
  obj$original_cell <- colnames(obj)
  DefaultAssay(obj) <- 'RNA'
  join_rna_layers(obj)
})
neuron <- merge(objects[[1L]], y = objects[[2L]], add.cell.ids = c('Exc', 'Inh'),
                 merge.data = TRUE, project = 'Neuron')
rm(objects)
neuron <- join_rna_layers(neuron)
neuron <- NormalizeData(neuron, assay = 'RNA', normalization.method = 'LogNormalize',
                        scale.factor = as.numeric(cfg$global$scale_factor), verbose = FALSE)
signatures <- list()
signature_summary <- list()

for (condition in c('Ext', 'Ret-Ext')) {
  cells <- colnames(neuron)[as.character(neuron$condition) == condition]
  local <- subset(neuron, cells = cells)
  local <- join_rna_layers(local)
  Idents(local) <- factor(as.character(local$IEGclass), levels = c('IEGhigh', 'IEGlow'))
  group_counts <- table(Idents(local))
  if (any(group_counts < 3L)) stop('At least three nuclei per IEG class are required in ', condition)
  markers <- FindMarkers(local, assay = 'RNA', ident.1 = 'IEGhigh', ident.2 = 'IEGlow',
                          test.use = 'wilcox', slot = 'data', min.pct = as.numeric(cfg$programs$min_pct),
                          logfc.threshold = as.numeric(cfg$programs$logfc_threshold),
                          only.pos = FALSE, base = 2, random.seed = as.integer(cfg$seed),
                          verbose = FALSE)
  if (!all(c('p_val_adj', 'avg_log2FC') %in% names(markers))) {
    stop('FindMarkers must return Bonferroni-adjusted P values and average log2 fold changes')
  }
  markers$gene <- rownames(markers)
  markers$condition <- condition
  marker_order <- order(markers$p_val_adj, markers$gene, na.last = TRUE)
  markers <- markers[marker_order, , drop = FALSE]
  rownames(markers) <- NULL
  candidate <- markers[!is.na(markers$p_val_adj) & markers$p_val_adj < 0.05 &
                         !is.na(markers$avg_log2FC) & markers$avg_log2FC > 0 &
                         !markers$gene %in% as.character(cfg$ieg$genes), , drop = FALSE]
  selected <- head(as.character(candidate$gene), as.integer(cfg$programs$top_genes))
  if (!length(selected)) stop('No positive non-IEG genes are available for the ', condition, ' program')
  program <- if (condition == 'Ext') 'Ext_IEGhigh_UP' else 'RetExt_IEGhigh_UP'
  signatures[[program]] <- selected
  slug <- if (condition == 'Ext') 'Ext' else 'RetExt'
  write_table(markers, file.path(out, 'signature_DE', paste0(slug, '_IEGhigh_vs_IEGlow_all.csv')))
  write_table(candidate, file.path(out, 'signature_DE', paste0(slug, '_IEGhigh_vs_IEGlow_UP.csv')))
  write_table(data.frame(rank = seq_along(selected), gene = selected, program = program),
              file.path(out, 'gene_sets', paste0(program, '.csv')))
  dir.create(file.path(out, 'gene_sets'), recursive = TRUE, showWarnings = FALSE)
  writeLines(selected, file.path(out, 'gene_sets', paste0('geneset_', slug, '_IEGhigh_vs_low_UP.txt')))
  signature_summary[[program]] <- data.frame(
    program = program, source_condition = condition, n_IEGhigh = as.integer(group_counts['IEGhigh']),
    n_IEGlow = as.integer(group_counts['IEGlow']), n_positive_candidates = nrow(candidate),
    n_selected = length(selected), requested_top_genes = as.integer(cfg$programs$top_genes),
    discovery_neuron_types = 'Excitatory + Inhibitory', p_adjust_method = 'Bonferroni',
    min_pct = as.numeric(cfg$programs$min_pct),
    logfc_threshold = as.numeric(cfg$programs$logfc_threshold), stringsAsFactors = FALSE
  )
}

high <- subset(neuron, cells = colnames(neuron)[neuron$IEGhigh])
scopes <- list(All = colnames(high),
                Excitatory = colnames(high)[high$neuron_type == 'Excitatory'],
                Inhibitory = colnames(high)[high$neuron_type == 'Inhibitory'])
all_scores <- list()
statistics <- list()
for (scope in names(scopes)) {
  if (!length(scopes[[scope]])) stop('No IEG-high nuclei for scope: ', scope)
  score_obj <- subset(high, cells = scopes[[scope]])
  score_obj <- join_rna_layers(score_obj)
  DefaultAssay(score_obj) <- 'RNA'
  for (program in names(signatures)) {
    available <- intersect(signatures[[program]], rownames(score_obj[['RNA']]))
    if (length(available) != length(signatures[[program]])) stop('Program genes are missing from the scoring object')
    score_obj <- AddModuleScore(score_obj, assay = 'RNA', features = list(available),
                                pool = rownames(score_obj[['RNA']]),
                                nbin = as.integer(cfg$programs$nbin),
                                ctrl = as.integer(cfg$programs$control_genes),
                                name = 'marker_program', seed = as.integer(cfg$seed),
                                search = FALSE, slot = 'data')
    m <- score_obj[[]]
    score_data <- data.frame(
      Cell = rownames(m), original_cell = as.character(m$original_cell),
      sample = as.character(m$sample), condition = as.character(m$condition),
      neuron_type = as.character(m$neuron_type), subtype = as.character(m$subtype),
      ieg_sum = m$ieg_sum, IEGhigh = m$IEGhigh, scope = scope, program = program,
      score = m$marker_program1, stringsAsFactors = FALSE
    )
    all_scores[[paste(scope, program)]] <- score_data
    comparison <- score_data[score_data$condition %in% c('Ext', 'Ret-Ext'), , drop = FALSE]
    ext <- comparison$score[comparison$condition == 'Ext']
    ret <- comparison$score[comparison$condition == 'Ret-Ext']
    if (!length(ext) || !length(ret)) stop('Ext and Ret-Ext nuclei are required for nuclear score comparison')
    result <- stats::wilcox.test(ext, ret, alternative = 'two.sided', paired = FALSE,
                                 exact = FALSE, correct = TRUE)
    statistics[[paste(scope, program)]] <- data.frame(
      scope = scope, program = program, comparison = 'Ext vs Ret-Ext',
      n_Ext_nuclei = length(ext), n_RetExt_nuclei = length(ret),
      median_Ext = stats::median(ext), median_RetExt = stats::median(ret),
      Wilcoxon_W = unname(result$statistic), pvalue = result$p.value,
      test = 'Two-sided nuclear Wilcoxon rank-sum', stringsAsFactors = FALSE
    )
    prism_n <- max(length(ext), length(ret))
    length(ext) <- prism_n
    length(ret) <- prism_n
    write_table(data.frame(Ext = ext, Ret_Ext = ret),
                file.path(out, 'Prism', paste0(scope, '_', program, '.csv')))
  }
}
write_table(do.call(rbind, signature_summary), file.path(out, 'signature_summary.csv'))
write_table(do.call(rbind, all_scores), file.path(out, 'scores_long.csv'))
write_table(do.call(rbind, statistics), file.path(out, 'score_statistics.csv'))
write_session(file.path(out, 'sessionInfo.txt'))
