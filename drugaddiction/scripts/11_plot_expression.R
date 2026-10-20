source('R/common.R')
source('R/plotting.R')
suppressPackageStartupMessages({library(ggrepel); library(Matrix)})
ctx <- pipeline_context()
dest <- new_stage(ctx, '11_expression_figures')
comparisons <- expand.grid(neuron_type = c('Excitatory', 'Inhibitory'),
                           contrast = c('RetExt_vs_Ext', 'Ext_vs_WD', 'RetExt_vs_WD'),
                           stringsAsFactors = FALSE)
for (i in seq_len(nrow(comparisons))) {
  z <- comparisons[i, ]
  d <- read_csv_table(stage_path(ctx, '05_DESeq2', 'tables', z$neuron_type,
                                paste0(z$contrast, '_all_genes.csv')))
  title <- paste(z$neuron_type, gsub('_', ' ', z$contrast))
  p <- volcano_plot(d, ctx, title)
  plot_export(p, file.path(dest, 'volcano', paste0(z$neuron_type, '_', z$contrast)),
              ctx$cfg$plot$volcano_width, ctx$cfg$plot$volcano_height, ctx$cfg$plot$dpi)
  e <- collect_enrichment(stage_path(ctx, '07_enrichment', 'tables', z$neuron_type, z$contrast))
  ep <- enrichment_plot(e, ctx, title)
  path <- paste0('enrichment/', z$neuron_type, '_', z$contrast)
  plot_export(ep$plot, file.path(dest, path), 11, max(5, 2 + .32 * nrow(ep$data)), ctx$cfg$plot$dpi)
  write_table(ep$data, file.path(dest, paste0(path, '_plot_data.csv')))
}

scores <- read_csv_table(stage_path(ctx, '06_scores', 'scores_long.csv'))
if (!all(c('scope', 'program', 'condition', 'score') %in% names(scores))) {
  stop('Score table requires scope, program, condition and score')
}
score_statistics <- read_csv_table(stage_path(ctx, '06_scores', 'score_statistics.csv'))
required_statistics <- c('scope', 'program', 'pvalue', 'padj', 'p_adjust_method',
                          'planned_tests', 'comparison_family')
if (!all(required_statistics %in% names(score_statistics))) {
  stop('Score statistics require BH-adjusted results from step 06')
}
if (anyNA(score_statistics[, c('scope', 'program', 'p_adjust_method', 'planned_tests',
                               'comparison_family')]) ||
    anyDuplicated(score_statistics[, c('scope', 'program')]) ||
    any(score_statistics$p_adjust_method != 'BH') ||
    any(score_statistics$planned_tests != 2L)) {
  stop('Score statistics must contain unique comparisons with BH adjustment across two programs')
}
expected_programs <- c('Ext_IEGhigh_UP', 'RetExt_IEGhigh_UP')
for (scope in unique(scores$scope)) {
  family <- score_statistics[score_statistics$scope == scope, , drop = FALSE]
  if (nrow(family) != length(expected_programs) ||
      !setequal(family$program, expected_programs) ||
      length(unique(family$comparison_family)) != 1L) {
    stop('Incomplete marker-program correction family for scope: ', scope)
  }
}
if (any(!is.na(score_statistics$padj) &
        (!is.finite(score_statistics$padj) | score_statistics$padj < 0 |
           score_statistics$padj > 1))) {
  stop('Invalid BH-adjusted score P values')
}
write_table(score_statistics, file.path(dest, 'scores', 'score_statistics.csv'))
for (scope in unique(scores$scope)) for (program in unique(scores$program)) {
  d <- scores[scores$scope == scope & scores$program == program &
                scores$condition %in% c('Ext', 'Ret-Ext'), , drop = FALSE]
  if (!nrow(d)) next
  d$condition <- factor(d$condition, levels = c('Ext', 'Ret-Ext'))
  p <- ggplot(d, aes(condition, score, fill = condition)) +
    geom_violin(trim = FALSE, linewidth = .35, alpha = .8) +
    geom_boxplot(width = .12, outlier.shape = NA, fill = 'white', linewidth = .35) +
    scale_fill_manual(values = condition_colors) +
    labs(x = NULL, y = 'Marker-program score', title = paste(scope, program)) +
    paper_theme(ctx) + theme(legend.position = 'none')
  statistic <- score_statistics[score_statistics$scope == scope &
                                   score_statistics$program == program, , drop = FALSE]
  if (nrow(statistic) != 1L) stop('A unique score comparison is required for ', scope, ' / ', program)
  label <- if (is.na(statistic$padj)) 'NA' else stars(statistic$padj)
  span <- diff(range(d$score, finite = TRUE)); if (!is.finite(span) || span == 0) span <- 1
  p <- p + annotate('text', x = 1.5, y = max(d$score, na.rm = TRUE) + .1 * span,
                    label = label, family = ctx$cfg$plot$font, size = 5)
  base <- paste0('scores/', scope, '_', program)
  plot_export(p, file.path(dest, base), 5, 6, ctx$cfg$plot$dpi)
}

genes <- read_csv_table(ctx$cfg$paths$heatmap_genes)
de <- read_csv_table(stage_path(ctx, '05_DESeq2', 'tables', 'Excitatory', 'RetExt_vs_Ext_all_genes.csv'))
sig <- de$gene[!is.na(de$padj) & de$padj <= ctx$cfg$de$padj]
genes$passes_DE <- genes$gene %in% sig
obj <- readRDS(stage_path(ctx, '04_ieg', 'Excitatory.rds'))
obj <- join_rna_layers(obj)
if (!all(c('condition', 'IEGhigh') %in% names(obj[[]]))) stop('Missing IEG metadata')
genes$available <- genes$gene %in% rownames(obj)
write_table(genes, file.path(dest, 'heatmap/gene_membership.csv'))
genes <- genes[genes$passes_DE & genes$available, , drop = FALSE]
genes <- genes[order(genes$gene_order), , drop = FALSE]
if (nrow(genes)) {
  m <- obj[[]]
  expr <- rna_matrix(obj, 'data')
  means <- vapply(c('Ext', 'Ret-Ext'), function(condition) {
    cells <- rownames(m)[m$IEGhigh & m$condition == condition]
    if (!length(cells)) stop('No IEG-high nuclei in ', condition)
    Matrix::rowMeans(expr[genes$gene, cells, drop = FALSE])
  }, numeric(nrow(genes)))
  rownames(means) <- genes$gene
  sds <- apply(means, 1, sd)
  z <- sweep(means, 1, rowMeans(means), '-')
  good <- is.finite(sds) & sds > 0
  z[good, ] <- sweep(z[good, , drop = FALSE], 1, sds[good], '/')
  z[!good, ] <- 0
  long <- do.call(rbind, lapply(seq_len(ncol(z)), function(j) {
    data.frame(gene = genes$gene, module = genes$module, condition = colnames(z)[j],
               mean_log_expression = means[, j], gene_Z = z[, j], stringsAsFactors = FALSE)
  }))
  long$gene <- factor(long$gene, levels = rev(genes$gene))
  long$condition <- factor(long$condition, levels = c('Ext', 'Ret-Ext'))
  long$module <- factor(long$module, levels = unique(genes$module))
  p <- ggplot(long, aes(condition, gene, fill = gene_Z)) + geom_tile() +
    facet_grid(module ~ ., scales = 'free_y', space = 'free_y') +
    scale_fill_gradient2(low = '#2166AC', mid = 'white', high = '#B2182B', midpoint = 0, name = 'Gene Z') +
    labs(x = NULL, y = NULL) + paper_theme(ctx) +
    theme(axis.text.y = element_text(face = 'italic'), strip.text.y = element_text(angle = 0))
  write_table(long, file.path(dest, 'heatmap/module_expression.csv'))
} else p <- empty_panel('No selected genes meet the adjusted-P criterion', ctx)
plot_export(p, file.path(dest, 'heatmap/module_heatmap'), 9, max(5, nrow(genes) * .22), ctx$cfg$plot$dpi)
write_session(file.path(dest, 'sessionInfo.txt'))
