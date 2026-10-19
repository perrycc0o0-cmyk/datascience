suppressPackageStartupMessages(library(jsonlite))

pipeline_context <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  opt <- list(config = 'config/analysis.json', output = NULL, level = 'major')
  if (length(a) %% 2L != 0L) stop('Arguments must be --config PATH and/or --output PATH')
  if (length(a)) for (i in seq.int(1L, length(a), by = 2L)) {
    key <- sub('^--', '', a[i])
    if (!key %in% names(opt)) stop('Unknown argument: ', a[i])
    opt[key] <- list(a[i + 1L])
  }
  project <- normalizePath('.', winslash = '/', mustWork = TRUE)
  cfg <- jsonlite::read_json(opt$config, simplifyVector = TRUE)
  out <- if (is.null(opt$output)) cfg$paths$output else opt$output
  set.seed(as.integer(cfg$seed))
  list(project = project, cfg = cfg, out = out, config_path = opt$config, level = opt$level)
}

stage_path <- function(ctx, stage, ...) file.path(ctx$out, stage, ...)

neuron_object_path <- function(ctx, neuron_type) {
  candidates <- stage_path(ctx, '03_neurons', c(paste0(neuron_type, '_annotated.rds'),
                                                paste0(neuron_type, '.rds')))
  present <- candidates[file.exists(candidates)]
  if (!length(present)) stop('Missing neuron object: ', neuron_type)
  present[1L]
}

new_stage <- function(ctx, stage) {
  p <- stage_path(ctx, stage)
  if (dir.exists(p) && length(list.files(p, all.files = TRUE, no.. = TRUE))) {
    stop('Output stage is not empty: ', p, '. Choose a new --output directory.')
  }
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = 'NA')
}

save_object <- function(x, path) {
  if (file.exists(path)) stop('Refusing to overwrite object: ', path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path, compress = 'gzip')
}

write_session <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), path)
}

read_samples <- function(ctx) {
  x <- utils::read.csv(ctx$cfg$paths$samples, stringsAsFactors = FALSE)
  need <- c('sample', 'condition', 'filtered_h5', 'raw_h5')
  if (!all(need %in% names(x))) stop('Sample table requires: ', paste(need, collapse = ', '))
  if (nrow(x) != 9L || anyDuplicated(x$sample) || anyNA(x[, need])) stop('Expected nine unique libraries')
  if (!all(x$condition %in% c('WD', 'Ext', 'Ret-Ext')) ||
      !all(table(factor(x$condition, levels = c('WD', 'Ext', 'Ret-Ext'))) == 3L)) {
    stop('Expected three independent pooled samples per condition')
  }
  x
}

join_rna_layers <- function(obj) {
  if (inherits(obj[['RNA']], 'Assay5')) {
    n <- SeuratObject::Layers(obj[['RNA']], search = '^counts($|\\.)')
    if (length(n) > 1L) obj <- SeuratObject::JoinLayers(obj, assay = 'RNA')
  }
  obj
}

rna_matrix <- function(obj, layer = 'counts') {
  if (inherits(obj[['RNA']], 'Assay5')) {
    SeuratObject::LayerData(obj[['RNA']], layer = layer)
  } else {
    SeuratObject::GetAssayData(obj, assay = 'RNA', slot = layer)
  }
}

integer_counts <- function(x) {
  v <- if (inherits(x, 'sparseMatrix')) x@x else as.vector(x)
  if (any(!is.finite(v)) || any(v < 0) || any(abs(v - round(v)) > 1e-8)) {
    stop('RNA counts must be finite, nonnegative integer UMI counts')
  }
  invisible(TRUE)
}

require_neuron_metadata <- function(obj) {
  m <- obj[[]]
  need <- c('sample', 'condition', 'subtype')
  if (!all(need %in% names(m)) || anyNA(m[, need])) stop('Missing sample/condition/subtype metadata')
  if (any(!nzchar(as.character(m$subtype))) || any(as.character(m$subtype) == 'Unknown')) {
    stop('Complete manual subtype annotation before downstream analysis')
  }
  invisible(TRUE)
}

bh <- function(p) stats::p.adjust(p, method = 'BH')

plot_export <- function(p, base, width, height, dpi = 600) {
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(base, '.png'), p, width = width, height = height,
                 dpi = dpi, bg = 'white', limitsize = FALSE)
  ggplot2::ggsave(paste0(base, '.pdf'), p, width = width, height = height,
                 device = grDevices::cairo_pdf, bg = 'white', limitsize = FALSE)
}
