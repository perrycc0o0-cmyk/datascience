#!/usr/bin/env Rscript
# Preranked GSEA of the fixed rat candidate pathway panel.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(script_dir, "workflow_helpers.R"))
args <- commandArgs(TRUE)
if (length(args) != 5L) {
  stop(paste("Usage: Rscript candidate_gsea.R <deseq2_output_dir>",
             "<excitatory_symbol_entrez.csv> <inhibitory_symbol_entrez.csv>",
             "<frozen_pathway_membership.csv> <output_dir>"))
}
limit_threads()
require_packages("fgsea")
out <- new_output_dir(args[5])
candidates <- read_table(file.path(script_dir, "config", "candidate_pathways.csv"))
stopifnot(nrow(candidates) == 43L, !anyDuplicated(candidates$ID),
          identical(candidates$candidate_order, seq_len(nrow(candidates))))
membership <- read_table(args[4])
stopifnot(all(c("ID", "Entrez") %in% names(membership)))
membership$Entrez <- as.character(membership$Entrez)
membership <- membership[membership$ID %in% candidates$ID, , drop = FALSE]
stopifnot(nrow(membership) > 0L, !anyNA(membership$Entrez), all(nzchar(membership$Entrez)))
sets <- setNames(lapply(candidates$ID, function(id) unique(membership$Entrez[membership$ID == id])),
                 candidates$ID)
candidates$annotation_present <- lengths(sets) > 0L
candidates$full_set_size <- lengths(sets)
write_table(candidates, file.path(out, "source_data", "candidate_pathway_coverage.csv"))
write_table(unique(membership[, c("ID", "Entrez")]),
            file.path(out, "source_data", "union_gene_set_membership.csv"))

rank_genes <- function(z, mapping) {
  stopifnot(all(c("SYMBOL", "ENTREZID") %in% names(mapping)), !anyDuplicated(z$gene))
  eligible <- z$eligible & is.finite(z$pvalue) & is.finite(z$stat) &
    is.finite(z$log2FoldChange) & is.finite(z$lfcSE)
  by_symbol <- split(as.character(mapping$ENTREZID), as.character(mapping$SYMBOL))
  possible <- lapply(z$original_symbol, function(symbol) {
    if (is.na(symbol) || !nzchar(symbol)) return(character())
    ids <- by_symbol[[symbol]]
    sort(unique(ids[!is.na(ids) & nzchar(ids)]))
  })
  audit <- data.frame(
    gene = z$gene, symbol = z$original_symbol, eligible_native_P = eligible,
    native_BH_family = z$in_BH_family, baseMean = z$baseMean, n_entrez = lengths(possible),
    all_entrez = vapply(possible, paste, character(1), collapse = ";"),
    ENTREZID = vapply(possible, function(v) if (length(v) == 1L) v else NA_character_, character(1)))
  audit$status <- ifelse(!eligible, "excluded_no_finite_native_test",
                         ifelse(audit$n_entrez == 0L, "excluded_unmapped",
                                ifelse(audit$n_entrez > 1L, "excluded_ambiguous_symbol", "eligible_unique_symbol")))
  use <- which(eligible & audit$n_entrez == 1L)
  use <- use[order(-audit$baseMean[use], audit$gene[use])]
  repeated <- duplicated(audit$ENTREZID[use])
  audit$status[use[repeated]] <- "excluded_duplicate_entrez_lower_baseMean"
  use <- use[!repeated]
  audit$status[use] <- "ranked"
  ranks <- data.frame(audit[use, , drop = FALSE], GSEA_rank = z$stat[use],
                      log2FoldChange = z$log2FoldChange[use], lfcSE = z$lfcSE[use],
                      DE_pvalue = z$pvalue[use], DE_padj = z$padj[use])
  ranks <- ranks[order(-ranks$GSEA_rank, ranks$ENTREZID), , drop = FALSE]
  ranks$rank <- seq_len(nrow(ranks))
  stopifnot(nrow(ranks) > 0L, !anyDuplicated(ranks$ENTREZID),
            all(abs(ranks$GSEA_rank - ranks$log2FoldChange / ranks$lfcSE) < 1e-10))
  list(ranks = ranks, audit = audit)
}

summaries <- list()
types <- c("Excitatory", "Inhibitory")
for (i in seq_along(types)) {
  type <- types[i]
  key <- paste0(type, "__RetExt_vs_Ext")
  de <- read_table(file.path(args[1], "datasets", key, "tables", "ALL_genes.csv"))
  stopifnot(all(de$neuron_type == type), all(de$numerator == "Ret-Ext"), all(de$denominator == "Ext"),
            all(de$method == "DESeq2"), all(de$n_pools == 6L))
  mapped <- rank_genes(de, read_table(args[i + 1L]))
  rr <- mapped$ranks
  write_table(mapped$audit, file.path(out, "source_data", paste0(type, "_gene_mapping_audit.csv")))
  write_table(rr, file.path(out, "source_data", paste0(key, "_ranked_genes.csv")))
  stats <- setNames(rr$GSEA_rank, rr$ENTREZID)
  z <- candidates
  z$neuron_type <- type
  z$contrast <- "RetExt_vs_Ext"
  z$effective_genes <- vapply(sets, function(ids) sum(ids %in% names(stats)), integer(1))
  z$status <- ifelse(!z$annotation_present, "no_rat_annotation_in_frozen_snapshot",
                     ifelse(z$effective_genes < 10L, "below_minSize10",
                            ifelse(z$effective_genes > 1000L, "above_maxSize1000", "tested")))
  if (!any(z$status == "tested")) stop("No testable pathways for ", type)
  set.seed(42L)
  fit <- as.data.frame(fgsea::fgseaMultilevel(pathways = sets, stats = stats,
    minSize = 10L, maxSize = 1000L, scoreType = "std", gseaParam = 1,
    eps = 0, nproc = 1L, sampleSize = 101L, nPermSimple = 1000L))
  stopifnot(setequal(fit$pathway, z$ID[z$status == "tested"]))
  j <- match(z$ID, fit$pathway)
  for (column in c("ES", "NES", "size", "log2err")) z[[column]] <- fit[[column]][j]
  z$pvalue <- fit$pval[j]
  z$BH <- NA_real_
  finite <- is.finite(z$pvalue)
  z$BH[finite] <- p.adjust(z$pvalue[finite], method = "BH")
  z$P005 <- finite & z$pvalue < 0.05
  z$BH005 <- is.finite(z$BH) & z$BH < 0.05
  z$observed_direction <- ifelse(!is.finite(z$NES), "not_estimable",
                                 ifelse(z$NES > 0, "RetExt_higher", "Ext_higher"))
  z$leadingEdgeENTREZ <- z$leadingEdgeSymbols <- ""
  lead <- list()
  for (k in which(!is.na(j))) {
    ids <- as.character(fit$leadingEdge[[j[k]]])
    ii <- match(ids, rr$ENTREZID)
    stopifnot(!anyNA(ii))
    z$leadingEdgeENTREZ[k] <- paste(ids, collapse = ";")
    z$leadingEdgeSymbols[k] <- paste(rr$symbol[ii], collapse = ";")
    if (length(ii)) {
      lead[[z$ID[k]]] <- data.frame(ID = z$ID[k], Description = z$Description[k], rr[ii, , drop = FALSE],
                                     ES = z$ES[k], NES = z$NES[k], pvalue = z$pvalue[k],
                                     BH = z$BH[k], display_P005 = z$P005[k])
    }
  }
  selected <- z[z$P005, , drop = FALSE]
  selected <- selected[order(-selected$NES, selected$candidate_order), , drop = FALSE]
  write_table(z, file.path(out, "tables", paste0(key, "_ALL43_GSEA.csv")))
  write_table(selected, file.path(out, "tables", paste0(key, "_P005_displayed_GSEA.csv")))
  write_table(if (length(lead)) do.call(rbind, lead) else data.frame(),
              file.path(out, "source_data", paste0(key, "_ALL43_leading_edge_genes.csv")))
  summaries[[type]] <- data.frame(neuron_type = type, candidates = nrow(z),
    testable = sum(z$status == "tested"), ranked_genes = nrow(rr), nominal_P005 = sum(z$P005),
    candidate_BH005 = sum(z$BH005), positive_P005 = sum(z$P005 & z$NES > 0, na.rm = TRUE),
    negative_P005 = sum(z$P005 & z$NES < 0, na.rm = TRUE))
}
write_table(do.call(rbind, summaries), file.path(out, "tables", "GSEA_union_summary.csv"))
write_session(out)
