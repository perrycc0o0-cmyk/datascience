require_version <- function(package, expected) {
  if (!requireNamespace(package, quietly = TRUE)) stop('Required package: ', package)
  found <- as.character(utils::packageVersion(package))
  if (!identical(found, as.character(expected))) {
    stop(package, ' version ', expected, ' is required; found ', found)
  }
}

read_gene_expression_h5 <- function(path) {
  if (!file.exists(path)) stop('Missing matrix: ', path)
  x <- Seurat::Read10X_h5(path, use.names = TRUE, unique.features = TRUE)
  if (is.list(x)) {
    if (!'Gene Expression' %in% names(x)) stop('Gene Expression matrix is absent: ', path)
    x <- x[['Gene Expression']]
  }
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop('Matrix requires unique feature names and barcodes: ', path)
  }
  integer_counts(x)
  x
}

available_pcs <- function(obj, requested) {
  k <- min(as.integer(requested), ncol(obj) - 1L, length(Seurat::VariableFeatures(obj)) - 1L)
  if (!is.finite(k) || k < 2L) stop('Insufficient cells or variable genes for PCA')
  as.integer(k)
}

select_global_pcs <- function(stdev, cap) {
  if (length(stdev) < 2L || any(!is.finite(stdev)) || sum(stdev) <= 0) {
    stop('PCA standard deviations are not valid')
  }
  pct <- stdev / sum(stdev) * 100
  cumulative <- cumsum(pct)
  a <- which(cumulative > 90 & pct < 5)
  b <- which(head(pct, -1L) - tail(pct, -1L) > 0.1)
  co1 <- if (length(a)) a[1L] else NA_integer_
  co2 <- if (length(b)) max(b) + 1L else NA_integer_
  candidates <- c(co1, co2, as.integer(cap), length(stdev))
  selected <- max(2L, min(candidates[is.finite(candidates)]))
  list(n = as.integer(selected), co1 = co1, co2 = co2,
       table = data.frame(PC = seq_along(stdev), stdev = stdev,
                          contribution_percent = pct, cumulative_percent = cumulative))
}

run_sample_harmony <- function(obj, dims) {
  if (!'sample' %in% colnames(obj[[]]) || anyNA(obj$sample)) stop('Sample metadata are required')
  harmony::RunHarmony(object = obj, group.by.vars = 'sample', reduction.use = 'pca',
                      dims.use = dims, project.dim = FALSE,
                      reduction.save = 'harmony', verbose = TRUE)
}

annotation_mapping <- function(path, clusters, label, neuron_type = NULL) {
  out <- stats::setNames(rep('Unknown', length(clusters)), clusters)
  if (!file.exists(path)) return(out)
  x <- utils::read.csv(path, colClasses = 'character', stringsAsFactors = FALSE)
  required <- c('cluster', label, if (!is.null(neuron_type)) 'neuron_type')
  if (!all(required %in% names(x))) stop('Annotation columns required: ', paste(required, collapse = ', '))
  if (!is.null(neuron_type)) x <- x[x$neuron_type == neuron_type & !is.na(x$neuron_type), , drop = FALSE]
  if (!nrow(x)) return(out)
  if (anyDuplicated(x$cluster)) stop('Duplicate annotation clusters in ', path)
  if (anyNA(x$cluster) || any(!nzchar(trimws(x$cluster)))) stop('Empty annotation cluster in ', path)
  unknown_clusters <- setdiff(x$cluster, clusters)
  if (length(unknown_clusters)) stop('Annotations contain absent clusters: ', paste(unknown_clusters, collapse = ', '))
  keep <- !is.na(x[[label]]) & nzchar(trimws(x[[label]]))
  out[x$cluster[keep]] <- trimws(x[[label]][keep])
  out
}

major_marker_panel <- function() list(
  Excitatory = c('Slc17a7', 'Satb2', 'Slc30a3'),
  Inhibitory = c('Gad1', 'Gad2', 'Dlx1'),
  Oligodendrocytes = c('Mog', 'Mbp', 'St18'),
  OPC = c('Cacng4', 'Pdgfra', 'Vcan'),
  Astrocytes = c('Slco1c1', 'Gfap', 'Aqp4'),
  Microglia = c('P2ry12', 'Tmem119', 'Cx3cr1'))

neuron_marker_panel <- function(neuron_type) {
  if (neuron_type == 'Excitatory') return(list(
    GLU_L6_CLA = 'Gng2', GLU_L5_6 = c('Syt6', 'Foxp2', 'Tle4'),
    GLU_L2_L3 = 'Glis3', GLU_L5 = 'Parm1', GLU_L4_5 = c('Il1rapl2', 'Rorb'),
    GLU_ADGM = c('Rerg', 'Arhgap6')))
  list(GABA_Amb = 'Slc17a6', GABA_Etv1 = 'Etv1', GABA_Lamp5 = 'Lamp5',
       GABA_MSN = 'Gng7', GABA_Sst = 'Sst', GABA_Hs3st2 = 'Hs3st2',
       GABA_Pvalb = 'Slit2', GABA_Vip = 'Vip')
}

annotation_theme <- function(ctx) {
  ggplot2::theme_classic(base_family = ctx$cfg$plot$font, base_size = ctx$cfg$plot$base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

export_embedding <- function(obj, path, reduction = 'umap') {
  e <- as.data.frame(Seurat::Embeddings(obj, reduction))
  e$cell <- rownames(e)
  m <- obj[[]]
  keep <- intersect(c('sample', 'condition', 'cluster', 'cell_type', 'neuron_type', 'subtype'), names(m))
  e <- cbind(e[, c('cell', setdiff(names(e), 'cell')), drop = FALSE], m[rownames(e), keep, drop = FALSE])
  write_table(e, path)
}

export_marker_plots <- function(ctx, obj, panels, label, out) {
  genes <- unique(unlist(panels, use.names = FALSE))
  availability <- data.frame(gene = genes, present = genes %in% rownames(obj))
  write_table(availability, file.path(out, 'marker_availability.csv'))
  genes <- intersect(genes, rownames(obj))
  if (!length(genes)) return(invisible(NULL))
  p <- Seurat::DotPlot(obj, features = genes, group.by = label, dot.scale = 10,
                       cols = c('lightgrey', 'red'), scale = TRUE) +
    ggplot2::scale_colour_gradientn(colours = c('lightgrey', 'blue', 'red')) +
    annotation_theme(ctx) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  plot_export(p, file.path(out, 'marker_dotplot'), max(10, length(genes) * 0.45), 7,
              ctx$cfg$plot$dpi)
  for (gene in genes) {
    p <- Seurat::FeaturePlot(obj, features = gene, reduction = 'umap', order = TRUE,
                            pt.size = 0.3, min.cutoff = 'q10', max.cutoff = 'q90',
                            cols = c('lightgrey', 'blue', 'red')) + annotation_theme(ctx)
    plot_export(p, file.path(out, 'marker_UMAP', gene), 8, 7, ctx$cfg$plot$dpi)
  }
}
