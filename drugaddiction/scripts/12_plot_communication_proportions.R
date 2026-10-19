source('R/common.R')
source('R/plotting.R')
suppressPackageStartupMessages({library(Matrix); library(igraph)})
ctx <- pipeline_context()
dest <- new_stage(ctx, '12_communication_figures')
counts <- read_csv_table(stage_path(ctx, '08_GLMM', 'pool_counts.csv'))
tests <- read_csv_table(stage_path(ctx, '08_GLMM', 'model_statistics.csv'))

for (nt in c('Excitatory', 'Inhibitory')) {
  obj <- readRDS(stage_path(ctx, '04_ieg', paste0(nt, '.rds')))
  m <- obj[[]]
  m <- m[m$IEGhigh & m$condition %in% c('Ext', 'Ret-Ext'), , drop = FALSE]
  comp <- as.data.frame(table(condition = m$condition, subtype = m$subtype), stringsAsFactors = FALSE)
  names(comp)[3] <- 'n_IEGhigh'
  comp$percent <- ave(comp$n_IEGhigh, comp$condition, FUN = function(x) 100 * x / sum(x))
  colors <- if (nt == 'Excitatory') exc_colors else inh_colors
  p <- ggplot(comp, aes(condition, percent, fill = subtype)) + geom_col(width = .65) +
    scale_fill_manual(values = colors) + labs(x = NULL, y = 'IEG-high composition (%)', fill = NULL) + paper_theme(ctx)
  base <- paste0('proportions/', nt, '_IEGhigh_composition')
  plot_export(p, file.path(dest, base), 7, 6, ctx$cfg$plot$dpi)
  write_table(comp, file.path(dest, paste0(base, '_data.csv')))
  analysis <- if (nt == 'Excitatory') 'excitatory_IEGhigh_composition' else 'inhibitory_IEGhigh_fraction'
  d <- counts[counts$analysis == analysis & !is.na(counts$percent), , drop = FALSE]
  if ('excluded' %in% names(d)) d <- d[!d$excluded, , drop = FALSE]
  d$condition <- factor(d$condition, levels = c('Ext', 'Ret-Ext'))
  t <- tests[tests$analysis == analysis, , drop = FALSE]
  t$label <- stars(t$padj)
  t$label[is.na(t$label)] <- 'NA'
  height <- aggregate(percent ~ subtype, d, max)
  names(height)[2] <- 'top'
  t <- merge(t, height, by = 'subtype', all.x = TRUE, sort = FALSE)
  span <- max(d$percent, na.rm = TRUE); if (!is.finite(span) || span == 0) span <- 1
  t$y <- t$top + .13 * span
  p <- ggplot(d, aes(condition, percent, color = condition)) +
    stat_summary(fun.min = min, fun.max = max, geom = 'linerange', linewidth = .55) +
    geom_point(position = position_jitter(width = .065, height = 0, seed = ctx$cfg$seed), size = 2.7) +
    facet_wrap(~ subtype, ncol = 3, scales = 'free_y') +
    scale_color_manual(values = condition_colors) +
    geom_text(data = t, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE,
              family = ctx$cfg$plot$font, size = 4.5) +
    labs(x = NULL, y = if (nt == 'Excitatory') 'Within-IEG-high composition (%)' else 'Within-subtype IEG-high (%)') +
    paper_theme(ctx) + theme(legend.position = 'none')
  base <- paste0('proportions/', nt, '_GLMM')
  plot_export(p, file.path(dest, base), 10, if (nt == 'Excitatory') 7 else 9, ctx$cfg$plot$dpi)
  rm(obj)
}

outgoing <- read_csv_table(stage_path(ctx, '09_CellChat', 'outgoing_activity_self_excluded.csv'))
comp <- read_csv_table(stage_path(ctx, '09_CellChat', 'RetExt_vs_Ext_outgoing_comparison.csv'))
if ('outgoing_weight_percent_change' %in% names(comp)) {
  d <- comp[comp$scope %in% c('All', 'Excitatory') & comp$subtype %in% c('GLU_L6_CLA', 'GLU_L5_6'), , drop = FALSE]
  d$scope <- factor(d$scope, levels = c('All', 'Excitatory'), labels = c('All neuronal targets', 'Excitatory targets'))
  p <- ggplot(d, aes(subtype, outgoing_weight_percent_change, fill = subtype)) +
    geom_col(width = .65, na.rm = TRUE) + facet_wrap(~ scope) + scale_fill_manual(values = exc_colors) +
    geom_hline(yintercept = 0, linewidth = .3) + labs(x = NULL, y = 'Ret-Ext vs Ext outgoing change (%)') +
    paper_theme(ctx) + theme(legend.position = 'none')
  plot_export(p, file.path(dest, 'communication/outgoing_percent_change'), 9, 5, ctx$cfg$plot$dpi)
}
d <- outgoing[outgoing$condition %in% c('Ext', 'Ret-Ext') & outgoing$scope %in% c('All', 'Excitatory') &
                outgoing$subtype %in% names(exc_colors), , drop = FALSE]
p <- ggplot(d, aes(subtype, outgoing_weight, fill = condition)) +
  geom_col(position = position_dodge(.8), width = .7) + facet_wrap(~ scope, scales = 'free_y') +
  scale_fill_manual(values = condition_colors) + labs(x = NULL, y = 'Inferred outgoing strength', fill = NULL) +
  paper_theme(ctx) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_export(p, file.path(dest, 'communication/outgoing_strength'), 12, 6, ctx$cfg$plot$dpi)

exc_ext <- read_csv_table(stage_path(ctx, '09_CellChat', 'Excitatory', 'Ext', 'edges.csv'))
exc_ret <- read_csv_table(stage_path(ctx, '09_CellChat', 'Excitatory', 'Ret-Ext', 'edges.csv'))
edges <- merge(exc_ext[, c('source', 'target', 'interaction_weight')],
               exc_ret[, c('source', 'target', 'interaction_weight')], by = c('source', 'target'),
               all = TRUE, suffixes = c('_Ext', '_RetExt'))
edges$interaction_weight_Ext[is.na(edges$interaction_weight_Ext)] <- 0
edges$interaction_weight_RetExt[is.na(edges$interaction_weight_RetExt)] <- 0
edges$difference <- edges$interaction_weight_RetExt - edges$interaction_weight_Ext
edges <- edges[edges$source != edges$target & edges$difference != 0, , drop = FALSE]
write_table(edges, file.path(dest, 'communication/differential_network_edges.csv'))
network_base <- file.path(dest, 'communication/differential_network')
dir.create(dirname(network_base), recursive = TRUE, showWarnings = FALSE)
network_draw <- function() {
  if (!nrow(edges)) { plot.new(); text(.5, .5, 'No differential edges'); return(invisible(NULL)) }
  vertices <- data.frame(name = names(exc_colors))
  graph <- igraph::graph_from_data_frame(edges[, c('source', 'target')], directed = TRUE, vertices = vertices)
  ew <- abs(edges$difference)
  width <- .4 + 4 * ew / max(ew)
  plot(graph, layout = igraph::layout_in_circle(graph), vertex.color = exc_colors[igraph::V(graph)$name],
       vertex.size = 22, vertex.label.cex = .9, vertex.label.family = ctx$cfg$plot$font,
       vertex.label.color = 'black', edge.color = ifelse(edges$difference > 0, '#D73027', '#4575B4'),
       edge.width = width, edge.arrow.size = .35, margin = .3)
}
grDevices::png(paste0(network_base, '.png'), width = 8, height = 8, units = 'in',
               res = ctx$cfg$plot$dpi, type = 'cairo', family = ctx$cfg$plot$font)
network_draw(); grDevices::dev.off()
grDevices::cairo_pdf(paste0(network_base, '.pdf'), width = 8, height = 8, family = ctx$cfg$plot$font)
network_draw(); grDevices::dev.off()

pathway_diff <- list()
sum_pathways <- function(x) {
  if (!nrow(x)) return(data.frame(source = character(), pathway = character(), probability = numeric()))
  aggregate(probability ~ source + pathway, x, sum)
}
for (scope in c('All', 'Excitatory')) {
  ext <- read_csv_table(stage_path(ctx, '09_CellChat', scope, 'Ext', 'outgoing_pathways.csv'))
  ret <- read_csv_table(stage_path(ctx, '09_CellChat', scope, 'Ret-Ext', 'outgoing_pathways.csv'))
  ext <- ext[ext$source != ext$target, , drop = FALSE]
  ret <- ret[ret$source != ret$target, , drop = FALSE]
  e <- sum_pathways(ext)
  r <- sum_pathways(ret)
  x <- merge(e, r, by = c('source', 'pathway'), all = TRUE, suffixes = c('_Ext', '_RetExt'))
  x$probability_Ext[is.na(x$probability_Ext)] <- 0
  x$probability_RetExt[is.na(x$probability_RetExt)] <- 0
  x$difference <- x$probability_RetExt - x$probability_Ext
  x$scope <- scope
  pathway_diff[[scope]] <- x
}
pd <- do.call(rbind, pathway_diff)
write_table(pd, file.path(dest, 'communication/pathway_differences.csv'))
for (direction in c('RetExt_higher', 'Ext_higher')) {
  if (direction == 'RetExt_higher') {
    z <- pd[pd$source == 'GLU_L6_CLA' & pd$difference > 0, , drop = FALSE]
    z$value <- z$difference
  } else {
    z <- pd[pd$source %in% c('GLU_L5', 'GLU_L5_6') & pd$difference < 0, , drop = FALSE]
    z$value <- -z$difference
  }
  z <- z[order(-z$value), , drop = FALSE]
  z <- utils::head(z, 20L)
  if (nrow(z)) {
    z$key <- make.unique(paste(z$source, z$pathway, z$scope, sep = ' | '))
    z$key <- factor(z$key, levels = rev(z$key))
    p <- ggplot(z, aes(value, key, fill = scope)) + geom_col(width = .7) +
      labs(x = if (direction == 'RetExt_higher') 'Ret-Ext minus Ext' else 'Ext minus Ret-Ext', y = NULL, fill = NULL) + paper_theme(ctx)
  } else p <- empty_panel('No pathways in the specified direction', ctx)
  plot_export(p, file.path(dest, paste0('communication/', direction, '_outgoing_pathways')), 11, max(5, .3 * nrow(z) + 2), ctx$cfg$plot$dpi)
}
z <- pd[pd$source %in% c('GLU_L6_CLA', 'GLU_L5', 'GLU_L5_6'), , drop = FALSE]
if (nrow(z)) {
  paths <- aggregate(abs(difference) ~ pathway, z, sum)
  keep <- utils::head(paths$pathway[order(-paths[[2]])], 20L)
  z <- z[z$pathway %in% keep, , drop = FALSE]
  z$scaled_difference <- 1000 * z$difference
  p <- ggplot(z, aes(pathway, source, fill = scaled_difference)) + geom_tile() +
    geom_text(aes(label = sprintf('%.2f', scaled_difference)), family = ctx$cfg$plot$font, size = 3) +
    facet_wrap(~ scope) + scale_fill_gradient2(low = '#4575B4', mid = 'white', high = '#D73027', midpoint = 0) +
    labs(x = NULL, y = NULL, fill = expression(10^3*Delta)) + paper_theme(ctx) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
} else p <- empty_panel('No outgoing pathway differences', ctx)
plot_export(p, file.path(dest, 'communication/pathway_matrix'), 13, 6, ctx$cfg$plot$dpi)
write_session(file.path(dest, 'sessionInfo.txt'))
