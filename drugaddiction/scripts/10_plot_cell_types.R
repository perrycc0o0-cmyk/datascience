source('R/common.R')
source('R/plotting.R')
ctx <- pipeline_context()
dest <- new_stage(ctx, '10_celltype_figures')
markers <- jsonlite::read_json('config/markers.json', simplifyVector = FALSE)
obj <- readRDS(stage_path(ctx, '02_annotation', 'annotated.rds'))
obj <- join_rna_layers(obj)
DefaultAssay(obj) <- 'RNA'

p <- object_umap(obj, 'cell_type', major_colors, ctx)
plot_export(p, file.path(dest, 'cell_types/major_cell_UMAP'), 8, 7, ctx$cfg$plot$dpi)
p <- marker_dotplot(obj, markers$major, 'cell_type', ctx)
plot_export(p, file.path(dest, 'cell_types/major_marker_dotplot'), 12, 5, ctx$cfg$plot$dpi)
m <- obj[[]]
composition <- as.data.frame(table(condition = m$condition, cell_type = m$cell_type), stringsAsFactors = FALSE)
names(composition)[3] <- 'n_nuclei'
composition$percent <- ave(composition$n_nuclei, composition$condition,
                          FUN = function(x) 100 * x / sum(x))
write_table(composition, file.path(dest, 'major_composition.csv'))
p <- ggplot(subset(composition, condition %in% c('Ext', 'Ret-Ext')),
            aes(condition, percent, fill = cell_type)) +
  geom_col(width = .65) + scale_fill_manual(values = major_colors) +
  labs(x = NULL, y = 'Nuclei (%)', fill = NULL) + paper_theme(ctx)
plot_export(p, file.path(dest, 'cell_types/major_cell_composition'), 6, 6, ctx$cfg$plot$dpi)

metrics <- c('nFeature_RNA', 'percent.ribo', 'nCount_RNA')
for (metric in metrics) {
  observed <- switch(metric, nFeature_RNA = 'nFeature_observed',
                      nCount_RNA = 'nCount_observed', percent.ribo = 'percent.ribo_observed')
  use <- if (observed %in% names(m)) observed else metric
  if (!use %in% names(m)) stop('QC metric not available: ', metric)
  d <- data.frame(sample = m$sample, value = m[[use]])
  p <- ggplot(d, aes(sample, value, fill = sample)) +
    geom_violin(scale = 'width', linewidth = .25, show.legend = FALSE) +
    labs(x = NULL, y = metric) + paper_theme(ctx) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot_export(p, file.path(dest, paste0('QC/', metric, '_by_sample')), 10, 5, ctx$cfg$plot$dpi)
}
major_features <- c('Slc17a7', 'Satb2', 'Gad1', 'Gad2', 'Mog', 'Mbp',
                    'Cacng4', 'Pdgfra', 'Slco1c1', 'Gfap', 'P2ry12', 'Tmem119')
for (gene in intersect(major_features, rownames(obj))) {
  p <- FeaturePlot(obj, features = gene, reduction = 'umap',
                   cols = c('#EEEEEE', '#753B8E'), raster = FALSE) + paper_theme(ctx)
  plot_export(p, file.path(dest, paste0('cell_types/marker_', gene)), 6, 5.5, ctx$cfg$plot$dpi)
}
rm(obj)
for (nt in c('Excitatory', 'Inhibitory')) {
  obj <- readRDS(neuron_object_path(ctx, nt))
  DefaultAssay(obj) <- 'RNA'
  cluster_col <- if ('neuron_cluster' %in% names(obj[[]])) 'neuron_cluster' else 'seurat_clusters'
  colors <- if (nt == 'Excitatory') exc_colors else inh_colors
  p <- object_umap(obj, cluster_col, NULL, ctx)
  plot_export(p, file.path(dest, paste0('neurons/', nt, '_clusters')), 7, 6, ctx$cfg$plot$dpi)
  p <- object_umap(obj, 'subtype', colors, ctx)
  plot_export(p, file.path(dest, paste0('neurons/', nt, '_subtypes')), 8, 6, ctx$cfg$plot$dpi)
  p <- marker_dotplot(obj, markers[[nt]], 'subtype', ctx)
  plot_export(p, file.path(dest, paste0('neurons/', nt, '_markers')), 10, 5, ctx$cfg$plot$dpi)
  for (gene in intersect(unique(unlist(markers[[nt]], use.names = FALSE)), rownames(obj))) {
    p <- FeaturePlot(obj, features = gene, reduction = 'umap',
                     cols = c('#EEEEEE', '#753B8E'), raster = FALSE) + paper_theme(ctx)
    plot_export(p, file.path(dest, paste0('Marker_UMAP/', nt, '_', gene)), 6, 5.5, ctx$cfg$plot$dpi)
  }
  emb <- as.data.frame(Embeddings(obj, 'umap'))
  write_table(data.frame(Cell = rownames(emb), emb, subtype = obj$subtype),
              file.path(dest, paste0(nt, '_UMAP_coordinates.csv')))
  rm(obj)
}
write_session(file.path(dest, 'sessionInfo.txt'))
