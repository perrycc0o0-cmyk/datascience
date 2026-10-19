suppressPackageStartupMessages(library(Seurat))
source('R/common.R')
source('R/preprocessing.R')
ctx <- pipeline_context()

if (ctx$level == 'major') {
  output <- stage_path(ctx, '02_annotation', 'annotated.rds')
  if (file.exists(output)) stop('Annotated object already exists: ', output)
  obj <- readRDS(stage_path(ctx, '02_annotation', 'all_nuclei.rds'))
  clusters <- sort(unique(as.character(obj$cluster)))
  mapping <- annotation_mapping(ctx$cfg$paths$major_annotations, clusters, 'cell_type')
  if (any(mapping == 'Unknown')) stop('Complete all cluster-to-cell-type labels first')
  obj$cell_type <- unname(mapping[as.character(obj$cluster)])
  save_object(obj, output)
  write_table(data.frame(cluster = clusters, cell_type = unname(mapping)),
              stage_path(ctx, '02_annotation', 'major_annotations_applied.csv'))
} else if (ctx$level == 'subtype') {
  for (nt in c('Excitatory', 'Inhibitory')) {
    output <- stage_path(ctx, '03_neurons', paste0(nt, '_annotated.rds'))
    if (file.exists(output)) stop('Annotated neuron object already exists: ', output)
  }
  objects <- list()
  for (nt in c('Excitatory', 'Inhibitory')) {
    obj <- readRDS(stage_path(ctx, '03_neurons', paste0(nt, '.rds')))
    clusters <- sort(unique(as.character(obj$cluster)))
    mapping <- annotation_mapping(ctx$cfg$paths$subtype_annotations, clusters, 'subtype', nt)
    if (any(mapping == 'Unknown')) stop('Complete all ', nt, ' subtype labels first')
    allowed <- names(neuron_marker_panel(nt))
    if (any(!mapping %in% allowed)) stop('Unexpected subtype name in ', nt, ' annotation table')
    obj$subtype <- unname(mapping[as.character(obj$cluster)])
    require_neuron_metadata(obj)
    objects[[nt]] <- obj
  }
  for (nt in names(objects)) {
    save_object(objects[[nt]], stage_path(ctx, '03_neurons', paste0(nt, '_annotated.rds')))
  }
} else stop('--level must be major or subtype')
