suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})
source('R/common.R')
source('R/preprocessing.R')

ctx <- pipeline_context()
cfg <- ctx$cfg
input <- stage_path(ctx, '02_annotation', 'annotated.rds')
if (!file.exists(input)) stop('A complete major-cell-type annotation is required: ', input)
out <- new_stage(ctx, '03_neurons')
all_nuclei <- readRDS(input)
if (!'cell_type' %in% names(all_nuclei[[]])) stop('Missing cell_type annotation')
templates <- list()
for (neuron_type in c('Excitatory', 'Inhibitory')) {
  cells <- colnames(all_nuclei)[all_nuclei$cell_type == neuron_type]
  if (length(cells) < 10L) stop('Too few nuclei for ', neuron_type)
  obj <- subset(all_nuclei, cells = cells)
  obj <- join_rna_layers(obj)
  DefaultAssay(obj) <- 'RNA'
  obj$neuron_type <- neuron_type
  obj <- NormalizeData(obj, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = 'vst', nfeatures = cfg$neurons$hvg, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  npcs <- available_pcs(obj, cfg$neurons$pcs)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs,
                 seed.use = cfg$seed, verbose = FALSE)
  harmony_dims <- seq_len(min(npcs, as.integer(cfg$neurons$harmony_pcs)))
  obj <- run_sample_harmony(obj, harmony_dims)
  dims <- seq_len(min(length(harmony_dims), as.integer(cfg$neurons$graph_dims)))
  pars <- cfg$neurons[[neuron_type]]
  obj <- FindNeighbors(obj, reduction = 'harmony', dims = dims, k.param = pars$k, verbose = FALSE)
  obj <- FindClusters(obj, resolution = cfg$neurons$resolution,
                        random.seed = cfg$seed, verbose = FALSE)
  obj$cluster <- as.character(Idents(obj))
  if (neuron_type == 'Excitatory') {
    obj <- RunUMAP(obj, reduction = 'harmony', dims = dims,
                     n.neighbors = pars$umap_neighbors, min.dist = pars$umap_min_dist,
                     spread = pars$umap_spread, repulsion.strength = pars$umap_repulsion,
                     seed.use = cfg$seed, verbose = FALSE)
  } else {
    obj <- RunUMAP(obj, reduction = 'harmony', dims = dims, seed.use = cfg$seed, verbose = FALSE)
  }
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = cfg$markers$min_pct,
                             logfc.threshold = cfg$markers$logfc_threshold,
                             max.cells.per.ident = cfg$markers$max_cells_per_ident,
                             random.seed = cfg$seed, verbose = FALSE)
  dir <- file.path(out, neuron_type)
  write_table(markers, file.path(dir, 'cluster_markers.csv'))
  clusters <- sort(unique(obj$cluster))
  mapping <- tryCatch(annotation_mapping(cfg$paths$subtype_annotations, clusters, 'subtype', neuron_type),
                      error = function(e) {
                        warning(conditionMessage(e))
                        stats::setNames(rep('Unknown', length(clusters)), clusters)
                      })
  obj$subtype <- unname(mapping[obj$cluster])
  templates[[neuron_type]] <- data.frame(neuron_type = neuron_type, cluster = clusters, subtype = '')
  write_table(data.frame(neuron_type = neuron_type, cluster = clusters, subtype = unname(mapping)),
               file.path(dir, 'subtype_annotations.csv'))
  write_table(data.frame(cell = colnames(obj), obj[[]], check.names = FALSE), file.path(dir, 'metadata.csv'))
  write_table(data.frame(neuron_type = neuron_type, available_PCs = npcs,
                          Harmony_PCs = length(harmony_dims), graph_PCs = length(dims),
                          neighbors_k = pars$k, resolution = cfg$neurons$resolution),
               file.path(dir, 'parameters.csv'))
  export_embedding(obj, file.path(dir, 'UMAP_coordinates.csv'))
  p <- DimPlot(obj, reduction = 'umap', group.by = 'cluster', label = TRUE,
                repel = TRUE, pt.size = 0.5) + annotation_theme(ctx)
  plot_export(p, file.path(dir, 'UMAP_clusters'), 8, 7, cfg$plot$dpi)
  export_marker_plots(ctx, obj, neuron_marker_panel(neuron_type), 'cluster', file.path(dir, 'markers_by_cluster'))
  if (all(obj$subtype != 'Unknown')) {
    p <- DimPlot(obj, reduction = 'umap', group.by = 'subtype', label = TRUE,
                  repel = TRUE, pt.size = 0.5) + annotation_theme(ctx)
    plot_export(p, file.path(dir, 'UMAP_subtypes'), 8, 7, cfg$plot$dpi)
    export_marker_plots(ctx, obj, neuron_marker_panel(neuron_type), 'subtype', file.path(dir, 'markers_by_subtype'))
  } else {
    message(neuron_type, ': complete ', cfg$paths$subtype_annotations,
            ' before subtype-level comparisons. Unassigned clusters remain Unknown.')
  }
  save_object(obj, file.path(out, paste0(neuron_type, '.rds')))
  rm(obj)
  invisible(gc())
}
write_table(do.call(rbind, templates), file.path(out, 'subtype_annotations_template.csv'))
write_session(file.path(out, 'sessionInfo.txt'))
