suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})
source('R/common.R')
source('R/preprocessing.R')

ctx <- pipeline_context()
cfg <- ctx$cfg
out <- new_stage(ctx, '02_annotation')
obj <- readRDS(stage_path(ctx, '01_preprocessing', 'merged.rds'))
obj <- join_rna_layers(obj)
DefaultAssay(obj) <- 'RNA'
integer_counts(rna_matrix(obj))
obj <- NormalizeData(obj, normalization.method = 'LogNormalize',
                      scale.factor = cfg$global$scale_factor, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = 'vst', nfeatures = cfg$global$hvg, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj),
                  vars.to.regress = cfg$global$vars_to_regress, verbose = FALSE)
npcs <- available_pcs(obj, cfg$global$pcs)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs,
               seed.use = cfg$seed, verbose = FALSE)
choice <- select_global_pcs(Stdev(obj, reduction = 'pca'), cfg$global$pc_cap)
write_table(choice$table, file.path(out, 'PCA_contributions.csv'))
write_table(data.frame(co1 = choice$co1, co2 = choice$co2, requested_PCs = cfg$global$pcs,
                        available_PCs = npcs, selected_PCs = choice$n), file.path(out, 'PCA_selection.csv'))
p <- ElbowPlot(obj, ndims = npcs) + annotation_theme(ctx)
plot_export(p, file.path(out, 'PCA_elbow'), 8, 6, cfg$plot$dpi)
obj <- run_sample_harmony(obj, seq_len(choice$n))
obj <- FindNeighbors(obj, reduction = 'harmony', dims = seq_len(choice$n),
                      k.param = cfg$global$neighbors_k, verbose = FALSE)
for (res in cfg$global$resolutions) {
  obj <- FindClusters(obj, resolution = res, random.seed = cfg$seed, verbose = FALSE)
}
obj <- FindClusters(obj, resolution = cfg$global$resolution, random.seed = cfg$seed, verbose = FALSE)
obj$cluster <- as.character(Idents(obj))
obj <- RunUMAP(obj, reduction = 'harmony', dims = seq_len(choice$n),
                 n.neighbors = cfg$global$umap_neighbors, min.dist = cfg$global$umap_min_dist,
                 metric = cfg$global$umap_metric, seed.use = cfg$seed, verbose = FALSE)
Idents(obj) <- obj$cluster
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = cfg$markers$min_pct,
                           logfc.threshold = cfg$markers$logfc_threshold,
                           max.cells.per.ident = cfg$markers$max_cells_per_ident,
                           random.seed = cfg$seed, verbose = FALSE)
write_table(markers, file.path(out, 'cluster_markers.csv'))
clusters <- sort(unique(obj$cluster))
mapping <- tryCatch(annotation_mapping(cfg$paths$major_annotations, clusters, 'cell_type'),
                    error = function(e) {
                      warning(conditionMessage(e))
                      stats::setNames(rep('Unknown', length(clusters)), clusters)
                    })
obj$cell_type <- unname(mapping[obj$cluster])
write_table(data.frame(cluster = clusters, cell_type = unname(mapping)),
             file.path(out, 'major_annotations.csv'))
write_table(data.frame(cluster = clusters, cell_type = ''), file.path(out, 'major_annotations_template.csv'))
write_table(as.data.frame(table(cluster = obj$cluster, cell_type = obj$cell_type)),
             file.path(out, 'cluster_annotation_counts.csv'))
save_object(obj, file.path(out, 'all_nuclei.rds'))
export_embedding(obj, file.path(out, 'UMAP_coordinates.csv'))
p <- DimPlot(obj, reduction = 'umap', group.by = 'cluster', label = TRUE, repel = TRUE,
              pt.size = 0.3) + annotation_theme(ctx)
plot_export(p, file.path(out, 'UMAP_clusters'), 9, 7, cfg$plot$dpi)
export_marker_plots(ctx, obj, major_marker_panel(), 'cluster', file.path(out, 'markers_by_cluster'))
if (all(obj$cell_type != 'Unknown')) {
  save_object(obj, file.path(out, 'annotated.rds'))
  cols <- c(Excitatory = '#1E3A8A', Inhibitory = '#DC2626', Oligodendrocytes = '#7C2D12',
            OPC = '#6D28D9', Astrocytes = '#059669', Microglia = '#E76F51')
  other <- setdiff(unique(obj$cell_type), names(cols))
  if (length(other)) cols <- c(cols, stats::setNames(grDevices::hcl.colors(length(other), 'Dark 3'), other))
  p <- DimPlot(obj, reduction = 'umap', group.by = 'cell_type', label = TRUE,
                repel = TRUE, cols = cols, pt.size = 0.3) + annotation_theme(ctx)
  plot_export(p, file.path(out, 'UMAP_cell_types'), 9, 7, cfg$plot$dpi)
  export_marker_plots(ctx, obj, major_marker_panel(), 'cell_type', file.path(out, 'markers_by_cell_type'))
} else {
  message('Complete the cluster-to-cell-type table in ', cfg$paths$major_annotations,
          ' before proceeding to neuron subclustering. Unassigned clusters remain Unknown.')
}
write_session(file.path(out, 'sessionInfo.txt'))
