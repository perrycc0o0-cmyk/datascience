read_one_to_one_orthologs <- function(path) {
  if (!file.exists(path)) {
    stop('Provide a rat-human ortholog table at ', path,
         ' with rat_symbol, human_symbol and orthology_type columns')
  }
  x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c('rat_symbol', 'human_symbol', 'orthology_type')
  if (!all(required %in% names(x))) stop('Ortholog table requires: ', paste(required, collapse = ', '))
  x <- x[x$orthology_type %in% c('one_to_one', 'ortholog_one2one'), required, drop = FALSE]
  if (!nrow(x) || anyNA(x) || any(!nzchar(x$rat_symbol)) || any(!nzchar(x$human_symbol))) {
    stop('No complete one-to-one rat-human gene-symbol mappings')
  }
  if (anyDuplicated(x$rat_symbol) || anyDuplicated(x$human_symbol)) {
    stop('One-to-one mapping requires unique rat symbols and unique human symbols')
  }
  x
}

cellchat_input <- function(obj, neuron_type, orthologs) {
  require_neuron_metadata(obj)
  obj <- join_rna_layers(obj)
  metadata <- obj[[]]
  if (!'IEGhigh' %in% names(metadata) || anyNA(metadata$IEGhigh)) {
    stop('IEGhigh metadata is required: ', neuron_type)
  }
  high <- metadata$IEGhigh
  if (is.factor(high)) high <- as.character(high)
  if (is.character(high)) {
    if (!all(high %in% c('TRUE', 'FALSE'))) stop('IEGhigh must be logical')
    high <- high == 'TRUE'
  } else if (!is.logical(high)) {
    if (!all(high %in% c(0, 1))) stop('IEGhigh must be logical')
    high <- as.logical(high)
  }
  metadata <- metadata[high, c('sample', 'condition', 'subtype'), drop = FALSE]
  metadata$neuron_type <- neuron_type
  if (!nrow(metadata)) stop('No IEG-high nuclei: ', neuron_type)
  expression <- rna_matrix(obj, 'data')
  counts <- rna_matrix(obj, 'counts')
  integer_counts(counts)
  if (is.null(rownames(expression)) || anyDuplicated(rownames(expression)) ||
      !all(rownames(metadata) %in% colnames(expression))) {
    stop('RNA expression features or cell names are not unique and aligned: ', neuron_type)
  }
  map <- orthologs[orthologs$rat_symbol %in% rownames(expression), , drop = FALSE]
  if (!nrow(map)) stop('RNA row names must match rat gene symbols in the ortholog table')
  expression <- expression[map$rat_symbol, rownames(metadata), drop = FALSE]
  rownames(expression) <- map$human_symbol
  nonzero <- if (inherits(expression, 'sparseMatrix')) expression@x else as.vector(expression)
  if (any(!is.finite(nonzero)) || any(nonzero < 0)) {
    stop('CellChat requires finite nonnegative log-normalized RNA data')
  }
  metadata$subtype <- as.character(metadata$subtype)
  metadata$condition <- as.character(metadata$condition)
  if (any(!metadata$condition %in% c('WD', 'Ext', 'Ret-Ext'))) stop('Unexpected condition label')
  list(expression = expression, metadata = metadata, mapping = map,
       input_features = nrow(counts))
}

merge_cellchat_inputs <- function(inputs) {
  features <- Reduce(intersect, lapply(inputs, function(x) rownames(x$expression)))
  if (!length(features)) stop('No shared mapped human genes between neuronal classes')
  metadata <- do.call(rbind, unname(lapply(inputs, `[[`, 'metadata')))
  if (anyDuplicated(rownames(metadata))) stop('Nucleus identifiers must be unique across neuronal classes')
  expression <- do.call(cbind, lapply(inputs, function(x) x$expression[features, , drop = FALSE]))
  expression <- expression[, rownames(metadata), drop = FALSE]
  list(expression = expression, metadata = metadata)
}

long_network_matrix <- function(x, value_name) {
  if (length(dim(x)) != 2L || !length(x)) {
    z <- data.frame(source = character(), target = character(), value = numeric(),
                    stringsAsFactors = FALSE)
  } else {
    z <- expand.grid(source = rownames(x), target = colnames(x),
                     stringsAsFactors = FALSE)
    z$value <- as.vector(x)
  }
  names(z)[names(z) == 'value'] <- value_name
  z
}

outgoing_network_summary <- function(cellchat, scope, condition) {
  weight <- cellchat@net$weight
  count <- cellchat@net$count
  groups <- rownames(weight)
  do.call(rbind, lapply(groups, function(group) {
    targets <- setdiff(colnames(weight), group)
    data.frame(scope = scope, condition = condition, subtype = group,
               n_nuclei = sum(as.character(cellchat@idents) == group),
               outgoing_weight = sum(weight[group, targets, drop = FALSE]),
               outgoing_interactions = sum(count[group, targets, drop = FALSE]),
               autocrine_weight = weight[group, group],
               autocrine_interactions = count[group, group],
               stringsAsFactors = FALSE)
  }))
}

compare_outgoing_conditions <- function(x) {
  ext <- x[x$condition == 'Ext', , drop = FALSE]
  ret <- x[x$condition == 'Ret-Ext', , drop = FALSE]
  keys <- c('scope', 'subtype')
  selected <- c(keys, 'n_nuclei', 'outgoing_weight', 'outgoing_interactions')
  z <- merge(ext[, selected, drop = FALSE], ret[, selected, drop = FALSE],
             by = keys, suffixes = c('_Ext', '_RetExt'), all = TRUE)
  for (metric in c('outgoing_weight', 'outgoing_interactions')) {
    baseline <- z[[paste0(metric, '_Ext')]]
    ret_value <- z[[paste0(metric, '_RetExt')]]
    estimable <- is.finite(baseline) & is.finite(ret_value) & baseline > 0
    z[[paste0(metric, '_difference')]] <- ret_value - baseline
    z[[paste0(metric, '_fold_change')]] <- ifelse(estimable, ret_value / baseline, NA_real_)
    z[[paste0(metric, '_percent_change')]] <- ifelse(estimable, 100 * (ret_value - baseline) / baseline,
                                                     NA_real_)
    z[[paste0(metric, '_zero_baseline')]] <- !is.na(baseline) & baseline == 0
  }
  z
}
