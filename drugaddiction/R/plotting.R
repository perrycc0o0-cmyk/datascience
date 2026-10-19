suppressPackageStartupMessages({library(ggplot2); library(Seurat)})

major_colors <- c(Excitatory = '#3B6FBF', Inhibitory = '#E85D5D',
                  Oligodendrocytes = '#B5895A', OPC = '#9A6FD1',
                  Astrocytes = '#4DB894', Microglia = '#F0986D')
condition_colors <- c(WD = '#9E9E9E', Ext = '#49B8D1', 'Ret-Ext' = '#EF4936')
exc_colors <- c(GLU_L6_CLA = '#8DD3C7', GLU_L5_6 = '#FFFFB3',
                GLU_L2_L3 = '#BEBADA', GLU_L5 = '#FB8072',
                GLU_L4_5 = '#80B1D3', GLU_ADGM = '#FDB462')
inh_colors <- c(GABA_Amb = '#66C2A5', GABA_Etv1 = '#FC8D62',
                GABA_Lamp5 = '#8DA0CB', GABA_MSN = '#E78AC3',
                GABA_Sst = '#A6D854', GABA_Hs3st2 = '#FFD92F',
                GABA_Pvalb = '#E5C494', GABA_Vip = '#FF6347')

paper_theme <- function(ctx) {
  theme_classic(base_size = ctx$cfg$plot$base_size, base_family = ctx$cfg$plot$font) +
    theme(axis.text = element_text(color = 'black'),
          plot.title = element_text(face = 'bold', hjust = 0),
          legend.title = element_text(face = 'plain'))
}

stars <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 1e-4, '****', ifelse(p < 1e-3, '***',
                ifelse(p < 1e-2, '**', ifelse(p < .05, '*', 'ns')))))
}

empty_panel <- function(message, ctx) {
  ggplot() + annotate('text', x = 0, y = 0, label = message,
                     family = ctx$cfg$plot$font, size = 4) + theme_void()
}

object_umap <- function(obj, group, colors, ctx, title = NULL) {
  if (!'umap' %in% names(obj@reductions)) stop('UMAP reduction is absent')
  Seurat::DimPlot(obj, reduction = 'umap', group.by = group,
                  cols = colors, pt.size = .15, raster = FALSE) +
    paper_theme(ctx) + labs(title = title, x = 'UMAP 1', y = 'UMAP 2')
}

marker_dotplot <- function(obj, markers, group, ctx) {
  features <- unique(unlist(markers, use.names = FALSE))
  features <- intersect(features, rownames(obj))
  if (!length(features)) return(empty_panel('No marker genes available', ctx))
  Seurat::DotPlot(obj, features = features, group.by = group,
                  cols = c('#DDEFF4', '#EF4936'), dot.scale = 8,
                  scale = TRUE, col.min = -2, col.max = 2) +
    paper_theme(ctx) + labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, face = 'italic'))
}

read_csv_table <- function(path) {
  if (!file.exists(path)) stop('Missing input table: ', path)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

volcano_plot <- function(d, ctx, title) {
  need <- c('gene', 'log2FoldChange', 'pvalue', 'padj')
  if (!all(need %in% names(d))) stop('DE table requires: ', paste(need, collapse = ', '))
  d <- d[is.finite(d$log2FoldChange) & is.finite(d$padj) & d$padj >= 0, , drop = FALSE]
  if (!nrow(d)) return(empty_panel('No adjusted-P values available', ctx))
  d$minus_log10_padj <- -log10(pmax(d$padj, .Machine$double.xmin))
  cutoff <- log2(ctx$cfg$de$fold_change)
  d$significant <- d$padj <= ctx$cfg$de$padj & abs(d$log2FoldChange) >= cutoff
  label <- d[d$significant, , drop = FALSE]
  label <- label[order(label$padj, -abs(label$log2FoldChange), label$gene), , drop = FALSE]
  label <- utils::head(label, 20L)
  p <- ggplot(d, aes(log2FoldChange, minus_log10_padj)) +
    geom_point(aes(color = log2FoldChange, size = minus_log10_padj), alpha = .72) +
    scale_color_gradientn(colors = c('#3288BD', '#66C2A5', '#FFFFBF', '#F46D43', '#9E0142'),
                          values = scales::rescale(c(-2, -1, 0, 1, 2)), limits = c(-2, 2),
                          oob = scales::squish, name = expression(log[2]*FC)) +
    scale_size_continuous(range = c(.7, 3.2), guide = 'none') +
    geom_hline(yintercept = -log10(ctx$cfg$de$padj), linetype = 'dashed', color = 'black', linewidth = .35) +
    geom_vline(xintercept = c(-cutoff, cutoff), linetype = 'dashed', color = 'grey50', linewidth = .35) +
    labs(title = title, x = expression(log[2]*' fold change'),
         y = expression(-log[10]*' BH-adjusted P')) + paper_theme(ctx)
  if (nrow(label)) p <- p + ggrepel::geom_text_repel(data = label, aes(label = gene),
          family = ctx$cfg$plot$font, fontface = 'italic', size = 4.5,
          max.overlaps = Inf, box.padding = .45, point.padding = .2, seed = ctx$cfg$seed)
  p
}

collect_enrichment <- function(directory) {
  files <- list.files(directory, pattern = '\\.csv$', full.names = TRUE)
  files <- files[!grepl('_significant\\.csv$', basename(files))]
  rows <- lapply(files, function(f) {
    x <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    if (!nrow(x) || !all(c('Description', 'p.adjust', 'Count') %in% names(x))) return(NULL)
    if (!'direction' %in% names(x)) x$direction <- if (grepl('^Up_', basename(f))) 'Up' else 'Down'
    if (!'database' %in% names(x)) x$database <- if (grepl('KEGG', basename(f))) 'KEGG' else 'GO'
    data.frame(term = x$Description, padj = x$p.adjust, Count = x$Count,
               direction = x$direction, database = x$database, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame(term = character(), padj = numeric(), Count = integer(),
                                      direction = character(), database = character()))
  do.call(rbind, rows)
}

enrichment_plot <- function(d, ctx, title) {
  d <- d[is.finite(d$padj) & d$padj <= ctx$cfg$ora$padj, , drop = FALSE]
  patterns <- utils::read.csv('config/enrichment_terms.csv', stringsAsFactors = FALSE)$pattern
  if (nrow(d)) d <- d[grepl(paste(patterns, collapse = '|'), d$term, ignore.case = TRUE), , drop = FALSE]
  if (!nrow(d)) return(list(plot = empty_panel('No significant selected pathways', ctx), data = d))
  d <- d[order(d$padj, -d$Count, d$term), , drop = FALSE]
  d <- do.call(rbind, lapply(split(d, d$direction), function(x) utils::head(x, 12L)))
  d$key <- make.unique(paste(d$term, d$direction, sep = ' | '))
  d$log_padj <- -log10(pmax(d$padj, .Machine$double.xmin))
  d$signed_log_padj <- ifelse(d$direction == 'Up', d$log_padj, -d$log_padj)
  d$key <- factor(d$key, levels = d$key[order(d$signed_log_padj)])
  labels <- setNames(d$term, as.character(d$key))
  p <- ggplot(d, aes(signed_log_padj, key)) +
    geom_segment(aes(x = 0, xend = signed_log_padj, yend = key), linewidth = 5, color = '#DDE5EA') +
    geom_point(aes(size = Count, color = padj, shape = database)) +
    scale_y_discrete(labels = labels) + scale_color_gradient(low = '#D73027', high = '#4575B4', name = 'BH-adjusted P') +
    scale_size_continuous(range = c(3, 8)) +
    labs(title = title, x = expression('Direction × '*-log[10]*' BH-adjusted P'), y = NULL) + paper_theme(ctx)
  list(plot = p, data = d)
}
