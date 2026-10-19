#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(Matrix)
})
source('R/common.R')
source('R/communication.R')
ctx <- pipeline_context()
out <- new_stage(ctx, '09_CellChat')
settings <- ctx$cfg$cellchat
if (!identical(settings$database, 'human') || !identical(settings$ortholog_policy, 'one_to_one')) {
  stop('CellChat requires the human database and one-to-one rat-human ortholog mapping')
}
orthologs <- read_one_to_one_orthologs(ctx$cfg$paths$orthologs)
inputs <- lapply(c('Excitatory', 'Inhibitory'), function(neuron_type) {
  cellchat_input(readRDS(stage_path(ctx, '04_ieg', paste0(neuron_type, '.rds'))),
                 neuron_type, orthologs)
})
names(inputs) <- c('Excitatory', 'Inhibitory')
samples <- read_samples(ctx)
sample_conditions <- setNames(samples$condition, samples$sample)
for (neuron_type in names(inputs)) {
  m <- inputs[[neuron_type]]$metadata
  if (any(!m$sample %in% samples$sample) ||
      any(m$condition != sample_conditions[as.character(m$sample)])) {
    stop('Metadata sample and condition disagree with the sample table: ', neuron_type)
  }
}
for (neuron_type in names(inputs)) {
  write_table(inputs[[neuron_type]]$mapping,
              file.path(out, paste0(neuron_type, '_orthologs.csv')))
}
write_table(do.call(rbind, lapply(names(inputs), function(neuron_type) {
  x <- inputs[[neuron_type]]
  data.frame(neuron_type = neuron_type, RNA_features = x$input_features,
             mapped_human_features = nrow(x$expression), IEGhigh_nuclei = ncol(x$expression))
})), file.path(out, 'mapping_summary.csv'))

scopes <- list(All = merge_cellchat_inputs(inputs),
               Excitatory = inputs$Excitatory, Inhibitory = inputs$Inhibitory)
conditions <- c('Overall', 'WD', 'Ext', 'Ret-Ext')
status <- list()
outgoing <- list()

run_network <- function(input, scope, condition, destination) {
  keep <- if (condition == 'Overall') rep(TRUE, nrow(input$metadata)) else {
    input$metadata$condition == condition
  }
  metadata <- input$metadata[keep, , drop = FALSE]
  expression <- input$expression[, rownames(metadata), drop = FALSE]
  if (!ncol(expression)) stop('No IEG-high nuclei for ', scope, ' / ', condition)
  metadata$subtype <- droplevels(factor(metadata$subtype))
  if (nlevels(metadata$subtype) < 2L) stop('Fewer than two annotated subtypes')
  write_table(data.frame(subtype = levels(metadata$subtype),
                         nuclei = as.integer(table(metadata$subtype))),
              file.path(destination, 'group_sizes.csv'))
  write_table(data.frame(cell = rownames(metadata), metadata, row.names = NULL),
              file.path(destination, 'cell_metadata.csv'))

  cellchat <- CellChat::createCellChat(object = expression, meta = metadata, group.by = 'subtype')
  database <- new.env(parent = emptyenv())
  utils::data('CellChatDB.human', package = 'CellChat', envir = database)
  cellchat@DB <- database$CellChatDB.human
  if (is.null(cellchat@DB)) stop('The human CellChat database is unavailable')
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, type = settings$probability_type,
                                         raw.use = isTRUE(settings$raw_use),
                                         population.size = isTRUE(settings$population_size),
                                         seed.use = as.integer(ctx$cfg$seed))
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = as.integer(settings$min_cells))
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  centrality_status <- 'no_inferred_pathways'
  if (length(cellchat@netP$pathways)) {
    centrality <- tryCatch(CellChat::netAnalysis_computeCentrality(cellchat, slot.name = 'netP'),
                           error = function(e) e)
    if (inherits(centrality, 'error')) {
      centrality_status <- paste('centrality_error:', conditionMessage(centrality))
    } else {
      cellchat <- centrality
      centrality_status <- 'computed'
    }
  }

  lr <- if (any(cellchat@net$prob > 0)) {
    CellChat::subsetCommunication(cellchat)
  } else data.frame(source = character(), target = character(), ligand = character(),
                    receptor = character(), prob = numeric(), pval = numeric())
  pathways <- if (length(cellchat@netP$pathways)) {
    CellChat::subsetCommunication(cellchat, slot.name = 'netP')
  } else data.frame(source = character(), target = character(), pathway_name = character(),
                    prob = numeric())
  pathway_edges <- pathways[, c('source', 'target', 'pathway_name', 'prob'), drop = FALSE]
  names(pathway_edges) <- c('source', 'target', 'pathway', 'probability')
  edges <- merge(long_network_matrix(cellchat@net$count, 'interaction_count'),
                  long_network_matrix(cellchat@net$weight, 'interaction_weight'),
                  by = c('source', 'target'), all = TRUE)
  write_table(lr, file.path(destination, 'ligand_receptor_communication.csv'))
  write_table(pathways, file.path(destination, 'pathway_communication.csv'))
  write_table(pathway_edges, file.path(destination, 'outgoing_pathways.csv'))
  write_table(edges, file.path(destination, 'edges.csv'))
  save_object(cellchat, file.path(destination, 'object.rds'))
  list(outgoing = outgoing_network_summary(cellchat, scope, condition),
       centrality = centrality_status, n_nuclei = nrow(metadata),
       n_groups = nlevels(metadata$subtype))
}

for (scope in names(scopes)) {
  for (condition in conditions) {
    destination <- file.path(out, scope, condition)
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    warnings <- character()
    answer <- tryCatch(withCallingHandlers(
      run_network(scopes[[scope]], scope, condition, destination),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart('muffleWarning')
      }
    ), error = function(e) e)
    if (inherits(answer, 'error')) {
      status[[length(status) + 1L]] <- data.frame(
        scope = scope, condition = condition, status = 'error', centrality = NA_character_,
        n_nuclei = NA_integer_, n_groups = NA_integer_,
        warnings = paste(unique(warnings), collapse = ' | '),
        error = conditionMessage(answer), stringsAsFactors = FALSE
      )
    } else {
      outgoing[[length(outgoing) + 1L]] <- answer$outgoing
      status[[length(status) + 1L]] <- data.frame(
        scope = scope, condition = condition, status = 'completed', centrality = answer$centrality,
        n_nuclei = answer$n_nuclei, n_groups = answer$n_groups,
        warnings = paste(unique(warnings), collapse = ' | '), error = '', stringsAsFactors = FALSE
      )
    }
    write_table(do.call(rbind, status), file.path(out, 'network_status.csv'))
  }
}
if (length(outgoing)) {
  outgoing <- do.call(rbind, outgoing)
  write_table(outgoing, file.path(out, 'outgoing_activity_self_excluded.csv'))
  comparison <- compare_outgoing_conditions(outgoing)
  write_table(comparison, file.path(out, 'RetExt_vs_Ext_outgoing_comparison.csv'))
  write_table(comparison[comparison$subtype %in% c('GLU_L6_CLA', 'GLU_L5_6'), , drop = FALSE],
              file.path(out, 'RetExt_vs_Ext_L6_CLA_L5_6.csv'))
}
write_session(file.path(out, 'sessionInfo.txt'))
if (any(vapply(status, function(x) x$status == 'error', logical(1)))) {
  stop('Some CellChat networks failed; see network_status.csv')
}
