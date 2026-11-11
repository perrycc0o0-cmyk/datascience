#!/usr/bin/env Rscript
# Nine-pool DESeq2 ranking and the twenty prespecified Monocle 2 signature fits.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
here <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(here, "workflow_helpers.R"))
source(file.path(here, "monocle2_compatibility.R"))
args <- commandArgs(TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript signature_trajectories.R <excitatory_payload.rds> <six_pool_exc_ret_ext_all_genes.csv> <output_dir>")
}
limit_threads()
invisible(Sys.setlocale("LC_COLLATE", "C"))
require_packages(c("DESeq2", "Matrix", "monocle", "DDRTree", "igraph", "VGAM", "Biobase", "dplyr"))
suppressPackageStartupMessages({
  library(dplyr)
  library(monocle)
  library(igraph)
  library(Matrix)
})
versions <- c(DESeq2 = "1.50.2", monocle = "2.38.0", DDRTree = "0.1.5",
              igraph = "2.2.1", dplyr = "1.1.4")
observed <- vapply(names(versions), function(p) as.character(packageVersion(p)), character(1))
if (!identical(versions, observed)) {
  stop("Use the recorded package versions for trajectory reproduction: ",
       paste(names(versions), versions, sep = "=", collapse = ", "))
}
out <- new_output_dir(args[3])
q <- read_ieg_payload(args[1], "Excitatory")
pb <- aggregate_pools(q)
baseline <- read_table(args[2])
stopifnot(identical(baseline$gene, rownames(pb)), all(baseline$neuron_type == "Excitatory"),
          all(baseline$n_pools == 6L), all(baseline$numerator == "Ret-Ext"),
          all(baseline$denominator == "Ext"), all(baseline$independentFiltering),
          all(baseline$alpha == 0.05))
keep <- as.logical(baseline$eligible)
md <- make_pool_metadata(pool_ids, c("Ext", "Ret-Ext", "WD"))
fit <- fit_native_deseq(pb, md, keep, "Ret-Ext", "Ext", contrast = c(0, 1, 0))
stopifnot(identical(DESeq2::resultsNames(fit$model),
                    c("Intercept", "condition_Ret.Ext_vs_Ext", "condition_WD_vs_Ext")))
de <- native_de_table(fit, pb, q$gene_annotation, keep, "Excitatory", "Ret-Ext", "Ext")
write_table(de, file.path(out, "nine_pool_DESeq2_ALL_genes.csv"))
write_table(data.frame(gene = rownames(pb), eligible_in_six_pool_fit = keep),
            file.path(out, "six_pool_filter_used_for_nine_pool_fit.csv"))
write_table(data.frame(md, total_counts = colSums(pb)), file.path(out, "nine_pool_metadata.csv"))

sizes <- c(10L, 20L, 30L, 40L, 50L, 100L, 200L, 300L, 400L, 500L)
valid <- de$tested & is.finite(de$log2FoldChange) & de$log2FoldChange != 0
directions <- list(F = valid & de$log2FoldChange > 0, G = valid & de$log2FoldChange < 0)
selected <- list()
for (panel in names(directions)) {
  d <- de[directions[[panel]], , drop = FALSE]
  d <- d[order(d$pvalue, -abs(d$log2FoldChange), d$gene), , drop = FALSE]
  stopifnot(nrow(d) >= max(sizes))
  d$rank_in_direction <- seq_len(nrow(d))
  for (size in sizes) {
    id <- sprintf("%s_Top%03d", panel, size)
    z <- head(d, size)
    z$case_id <- id
    z$panel <- panel
    z$N <- size
    selected[[id]] <- z
  }
}
write_table(do.call(rbind, selected), file.path(out, "all20_selected_genes.csv"))

# Signature-independent count preprocessing, computed once on all eligible nuclei.
nuclear_md <- q$metadata
nuclear_md$condition <- factor(nuclear_md$condition, levels = c("WD", "Ext", "Ret-Ext"))
rownames(nuclear_md) <- nuclear_md$Cell
features <- data.frame(gene_short_name = rownames(q$counts), row.names = rownames(q$counts))
set.seed(42)
base <- newCellDataSet(q$counts,
                       phenoData = methods::new("AnnotatedDataFrame", data = nuclear_md),
                       featureData = methods::new("AnnotatedDataFrame", data = features),
                       expressionFamily = VGAM::negbinomial.size())
base <- estimateSizeFactors(base)
base <- estimateDispersions(base, cores = 1)
apply_monocle2_compatibility()

edge_table <- function(graph) {
  x <- igraph::as_data_frame(graph, what = "edges")
  ends <- cbind(pmin(x$from, x$to), pmax(x$from, x$to))
  x$from <- ends[, 1]
  x$to <- ends[, 2]
  x <- x[order(x$from, x$to), , drop = FALSE]
  rownames(x) <- NULL
  x
}

fit_case <- function(s) {
  id <- s$case_id[1]
  panel <- s$panel[1]
  n <- s$N[1]
  dest <- file.path(out, "cases", id)
  result <- data.frame(case_id = id, panel = panel, N = n, status = "failed",
                       failure_stage = NA_character_, failure_reason = NA_character_,
                       mean_Ext = NA_real_, SEM_Ext = NA_real_,
                       mean_RetExt = NA_real_, SEM_RetExt = NA_real_,
                       delta_RetExt_minus_Ext = NA_real_, Welch_t = NA_real_,
                       Welch_df = NA_real_, Welch_P = NA_real_,
                       CI95_difference_low = NA_real_, CI95_difference_high = NA_real_,
                       secondary_exact20_pool_permutation_P = NA_real_,
                       selected_State = NA_integer_, WD_nuclei_selected_State = NA_integer_,
                       principal_root = NA_character_, projected_root = NA_character_)
  stage <- "DDRTree fit"
  tryCatch({
    cds <- setOrderingFilter(base, s$gene)
    set.seed(42)
    cds <- reduceDimension(cds, max_components = 2, reduction_method = "DDRTree")
    stage <- "Initial ordering"
    cds <- orderCells(cds)
    pre <- pData(cds)
    principal <- minSpanningTree(cds)
    graph <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree
    closest <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    stopifnot(identical(rownames(pre), colnames(cds)),
              identical(rownames(closest), colnames(cds)),
              identical(V(graph)$name, colnames(cds)), ecount(graph) == ncol(cds) - 1L)
    leaves <- names(which(degree(graph) == 1L))
    states <- sort(unique(as.integer(as.character(pre$State))))
    stopifnot(length(states) > 0L, !anyNA(pre$State), !anyNA(states))
    stage <- "Root selection"
    choices <- lapply(states, function(state) {
      root <- tryCatch(getFromNamespace("select_root_cell", "monocle")(cds, root_state = state, reverse = FALSE),
                       error = function(e) NA_character_)
      index <- which(V(principal)$name %in% root)
      valid_root <- length(root) == 1L && !is.na(root) && length(index) == 1L
      mapped <- if (valid_root) V(graph)$name[which(as.vector(closest) == index)] else character()
      terminal <- mapped[mapped %in% leaves]
      data.frame(State = state, WD_nuclei = sum(as.character(pre$State) == state & pre$condition == "WD"),
                 principal_root = if (valid_root) root else NA_character_,
                 projected_root = if (length(terminal)) terminal[1] else NA_character_,
                 eligible = valid_root && length(terminal) > 0L)
    })
    choices <- do.call(rbind, choices)
    choices$selected <- FALSE
    feasible <- which(choices$eligible)
    if (length(feasible)) {
      ix <- feasible[order(-choices$WD_nuclei[feasible], choices$State[feasible])][1]
      choices$selected[ix] <- TRUE
    }
    write_table(choices, file.path(dest, "root_state_eligibility.csv"))
    if (!length(feasible)) stop("No structurally feasible root state; no fallback root is used.")
    chosen <- choices[choices$selected, , drop = FALSE]
    reduced_s <- reducedDimS(cds)
    reduced_k <- reducedDimK(cds)
    stage <- "Ordering at selected root"
    set.seed(42)
    cds <- orderCells(cds, root_state = chosen$State)
    pd <- pData(cds)
    final_graph <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree
    root_cell <- cds@auxOrderingData[["DDRTree"]]$root_cell
    stopifnot(identical(root_cell, chosen$projected_root),
              identical(reducedDimS(cds), reduced_s), identical(reducedDimK(cds), reduced_k),
              identical(edge_table(minSpanningTree(cds)), edge_table(principal)),
              identical(edge_table(final_graph), edge_table(graph)),
              identical(pd$State, pre$State), all(is.finite(pd$Pseudotime)))
    root_distance <- as.numeric(distances(final_graph, v = root_cell, to = colnames(cds),
                                           weights = edge_attr(final_graph, "weight")))
    stopifnot(max(abs(root_distance - pd$Pseudotime)) < 1e-8)
    cells <- data.frame(Cell = colnames(cds), pool_id = pd$pool_id, condition = pd$condition,
                        State = pd$State, Pseudotime = pd$Pseudotime,
                        DDRTree_1 = reduced_s[1, ], DDRTree_2 = reduced_s[2, ])
    write_table(cells, file.path(dest, "nuclear_pseudotime.csv"))
    stage <- "Pooled-sample comparison"
    means <- aggregate(Pseudotime ~ pool_id + condition, cells, mean)
    names(means)[names(means) == "Pseudotime"] <- "sample_mean_pseudotime"
    means$n_nuclei <- as.integer(table(cells$pool_id)[as.character(means$pool_id)])
    means$case_id <- id
    means$panel <- panel
    means$N <- n
    x <- means$sample_mean_pseudotime[means$condition == "Ext"]
    y <- means$sample_mean_pseudotime[means$condition == "Ret-Ext"]
    stopifnot(length(x) == 3L, length(y) == 3L)
    tt <- t.test(y, x, var.equal = FALSE, alternative = "two.sided")
    pooled <- c(x, y)
    permuted <- apply(combn(6L, 3L), 2L, function(ix) mean(pooled[ix]) - mean(pooled[-ix]))
    result$status <- "success"
    result$mean_Ext <- mean(x); result$SEM_Ext <- sd(x) / sqrt(3)
    result$mean_RetExt <- mean(y); result$SEM_RetExt <- sd(y) / sqrt(3)
    result$delta_RetExt_minus_Ext <- mean(y) - mean(x)
    result$Welch_t <- unname(tt$statistic); result$Welch_df <- unname(tt$parameter)
    result$Welch_P <- tt$p.value
    result$CI95_difference_low <- tt$conf.int[1]; result$CI95_difference_high <- tt$conf.int[2]
    result$secondary_exact20_pool_permutation_P <-
      mean(abs(permuted) >= abs(result$delta_RetExt_minus_Ext) - 1e-12)
    result$selected_State <- chosen$State
    result$WD_nuclei_selected_State <- chosen$WD_nuclei
    result$principal_root <- chosen$principal_root
    result$projected_root <- root_cell
    write_table(means, file.path(dest, "all9pool_mean_pseudotime.csv"))
    list(test = result, means = means)
  }, error = function(e) {
    result$failure_stage <- stage
    result$failure_reason <- conditionMessage(e)
    message(id, " failed at ", stage, ": ", conditionMessage(e))
    list(test = result, means = NULL)
  })
}

fits <- lapply(selected, fit_case)
tests <- do.call(rbind, lapply(fits, `[[`, "test"))
tests$size_N <- tests$N
tests$BH20 <- p.adjust(tests$Welch_P, method = "BH", n = 20L)
tests <- tests[order(tests$N, tests$panel), , drop = FALSE]
write_table(tests, file.path(out, "tables", "all20_tests.csv"))
means <- do.call(rbind, lapply(fits, `[[`, "means"))
if (!is.null(means)) {
  means$size_N <- means$N
  write_table(means, file.path(out, "source_data", "all20_pool_means.csv"))
}
write_session(out)
