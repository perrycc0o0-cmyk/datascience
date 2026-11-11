#!/usr/bin/env Rscript
# Six original pools per contrast; integer counts and native independent filtering.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script, mustWork = TRUE)), "workflow_helpers.R"))
args <- commandArgs(TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript pseudobulk_deseq2.R <excitatory_payload.rds> <inhibitory_payload.rds> <output_dir>")
}
limit_threads()
require_packages(c("DESeq2", "edgeR", "Matrix"))
out <- new_output_dir(args[3])
types <- c("Excitatory", "Inhibitory")
pairs <- list(RetExt_vs_Ext = c("Ret-Ext", "Ext"),
              Ext_vs_WD = c("Ext", "WD"), RetExt_vs_WD = c("Ret-Ext", "WD"))
summary_rows <- list()
for (i in seq_along(types)) {
  type <- types[i]
  q <- read_ieg_payload(args[i], type)
  all_counts <- aggregate_pools(q)
  for (contrast in names(pairs)) {
    numerator <- pairs[[contrast]][1]
    denominator <- pairs[[contrast]][2]
    pools <- c(paste0(denominator, 1:3), paste0(numerator, 1:3))
    pb <- all_counts[, pools, drop = FALSE]
    md <- make_pool_metadata(pools, c(denominator, numerator))
    n <- as.integer(table(factor(q$metadata$pool_id, levels = pools)))
    keep <- expression_filter(pb, md)
    fit <- fit_native_deseq(pb, md, keep, numerator, denominator,
                           allow_convergence_retry = contrast != "RetExt_vs_Ext")
    tab <- native_de_table(fit, pb, q$gene_annotation, keep, type, numerator, denominator)
    tab$contrast <- contrast
    key <- paste(type, contrast, sep = "__")
    dest <- file.path(out, "datasets", key)
    write_table(data.frame(gene = rownames(pb), pb, check.names = FALSE),
                file.path(dest, "inputs", "counts.csv"))
    write_table(data.frame(md, nuclei = n, total_counts = colSums(pb)),
                file.path(dest, "inputs", "pool_metadata.csv"))
    write_table(data.frame(gene = rownames(pb), eligible = keep),
                file.path(dest, "source_data", "fresh_expression_filter.csv"))
    write_table(tab, file.path(dest, "tables", "ALL_genes.csv"))
    write_table(tab[tab$BH_significant, , drop = FALSE],
                file.path(dest, "tables", "BH005_genes.csv"))
    write_table(tab[tab$display_DEG, , drop = FALSE],
                file.path(dest, "tables", "BH005_FC115_genes.csv"))
    summary_rows[[key]] <- data.frame(
      neuron_type = type, contrast = contrast, n_nuclei = sum(n),
      n_pools = 6L, n_per_condition = 3L, eligible_genes = sum(keep),
      finite_P = sum(tab$tested), BH_family = sum(tab$in_BH_family),
      BH_up = sum(tab$BH_significant & tab$log2FoldChange > 0, na.rm = TRUE),
      BH_down = sum(tab$BH_significant & tab$log2FoldChange < 0, na.rm = TRUE),
      native_independent_filter_threshold = fit$filter_threshold,
      fit_type = attr(DESeq2::dispersionFunction(fit$model), "fitType"),
      same_model_convergence_retry = fit$convergence_retry)
  }
}
write_table(do.call(rbind, summary_rows), file.path(out, "DESeq2_pairwise_summary.csv"))
write_session(out)
