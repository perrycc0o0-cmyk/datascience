#!/usr/bin/env Rscript
# Recreate the study's accepted cluster-6 annotation change from stage-2 results.
# The input label, size and output counts are verified before writing the mapping.
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 1L) stop('Usage: Rscript 00_apply_accepted_subtype_override.R STAGE2_DIRECTORY')
root <- normalizePath(a[[1L]], mustWork = TRUE)
ann <- file.path(root, '03_annotation', 'Inh_cell_annotations_UMAP.csv.gz')
if (!file.exists(ann)) stop('Complete inhibitory subtype annotation first: ', ann)
d <- read.csv(gzfile(ann), stringsAsFactors = FALSE, check.names = FALSE)
required <- c('cell_id', 'subcluster', 'subtype_proposed', 'sample_id', 'condition_label', 'UMAP_1', 'UMAP_2')
stopifnot(all(required %in% names(d)), nrow(d) == 21693L, !anyDuplicated(d$cell_id))
before <- as.character(d$subtype_proposed)
move <- as.character(d$subcluster) == '6'
stopifnot(sum(move) == 1377L, all(before[move] == 'GABA_Lamp5'))
after <- before
after[move] <- 'GABA_Etv1'
stopifnot(sum(after == 'GABA_Etv1') == 8139L, sum(after == 'GABA_Lamp5') == 1085L,
          all(after[as.character(d$subcluster) == '8'] == 'GABA_Lamp5'), all(after[!move] == before[!move]))
mapping <- data.frame(cell_id = d$cell_id, subcluster = as.character(d$subcluster), original_label = before,
                      preview_label = after, sample_id = d$sample_id, condition_label = d$condition_label,
                      UMAP_1 = d$UMAP_1, UMAP_2 = d$UMAP_2)
out <- file.path(root, '07_cluster6_to_Etv1_preview', 'source_data', 'preview_cell_mapping.csv.gz')
if (file.exists(out)) {
  old <- read.csv(gzfile(out), stringsAsFactors = FALSE, check.names = FALSE)
  old$subcluster <- as.character(old$subcluster)
  stopifnot(isTRUE(all.equal(old, mapping, check.attributes = FALSE, tolerance = 1e-10)))
} else {
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(out, 'wt')
  tryCatch(write.csv(mapping, con, row.names = FALSE), finally = close(con))
}
message('Verified accepted inhibitory annotation mapping.')
