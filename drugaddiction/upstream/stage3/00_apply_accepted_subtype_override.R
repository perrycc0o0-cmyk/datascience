# 00 apply accepted subtype override

# 1. Packages
library(utils, include.only = character())

# 2. Paths
input_dir <- "/path/to/matrices"
work_dir <- "/path/to/work"
result_dir <- file.path(work_dir, "stage3")
figure_dir <- result_dir
fontconfig_file <- NULL
script_path <- "/path/to/00_apply_accepted_subtype_override.R"
stage1_dir <- file.path(work_dir, "stage1")
stage2_dir <- file.path(work_dir, "stage2")

# 3. Samples and conditions
sample_cols <- c("WD1", "WD2", "WD3", "Ext1", "Ext2", "Ext3", "Ret-Ext1", "Ret-Ext2", "Ret-Ext3")
metadata <- data.frame(
  pool_id = sample_cols,
  source_sample_id = c("Con1", "Con2", "Con3", "Ext1", "Ext2", "Ext3", "Ret_ext1", "Ret_ext2", "Ret_ext3"),
  model = rep("Rat IL", 9L),
  condition = factor(rep(c("WD", "Ext", "Ret-Ext"), each = 3L), levels = c("WD", "Ext", "Ret-Ext")),
  rats_per_pool = rep(5L, 9L), stringsAsFactors = FALSE)

# 4. Parameters and marker definitions
neuron_settings <- list(
  `schema_version` = "1.0",
  `run_id` = "20261101_Manuscript_matched_stage2_01",
  `stage_gate` = "neuron_subtype_review",
  `source_class_column` = "major_class_six",
  `keep_previously_forced_clusters` = TRUE,
  `do_not_refilter` = TRUE,
  `analysis_seed` = 42L,
  `runtime` = list(
    `threads` = 4L,
    `blas_threads` = 1L,
    `r_libraries` = list(),
    `fontconfig_file` = NULL,
    `font_family` = "Arial"
  ),
  `processing` = list(
    `assay` = "RNA",
    `counts_basis` = "Existing SoupX-corrected integer UMI counts; original observed counts retained separately",
    `normalization` = "LogNormalize",
    `scale_factor` = 10000L,
    `hvg_method` = "vst",
    `n_hvg` = 2000L,
    `scale_vars_to_regress` = list(),
    `pca_components_max` = 100L,
    `harmony_pca_dimensions_max` = 30L,
    `harmony_group_by` = "sample_id",
    `harmony_max_iter` = 10L,
    `harmony_project_dim` = FALSE,
    `neighbor_harmony_dimensions_max` = 19L,
    `clustering_algorithm` = 1L,
    `resolution` = 0.5,
    `umap_engine` = "uwot::umap",
    `umap_metric` = "cosine",
    `umap_repulsion_strength` = 1L,
    `umap_sgd_threads` = 1L
  ),
  `lineages` = list(
    list(
      `id` = "Exc",
      `source_class` = "Excitatory",
      `expected_nuclei` = 62991L,
      `output_directory` = "01_Exc",
      `neighbors_k` = 20L,
      `umap_n_neighbors` = 10L,
      `umap_min_dist` = 0.8,
      `umap_spread` = 6L
    ),
    list(
      `id` = "Inh",
      `source_class` = "Inhibitory",
      `expected_nuclei` = 21693L,
      `output_directory` = "02_Inh",
      `neighbors_k` = 30L,
      `umap_n_neighbors` = 30L,
      `umap_min_dist` = 0.3,
      `umap_spread` = 1L
    )
  ),
  `annotation_marker_discovery` = list(
    `test` = "wilcox",
    `min_pct` = 0.25,
    `logfc_threshold` = 0.25,
    `max_cells_per_ident` = 1000L,
    `only_pos` = FALSE,
    `return_thresh` = 1L,
    `purpose` = "Annotation markers only, not between-condition DE"
  ),
  `parent_stage1_root` = NULL,
  `input_rds` = NULL
)
neuron_settings$parent_stage1_root <- stage1_dir
neuron_settings$input_rds <- file.path(stage1_dir, "05_six_class_review", "global_celltype_six_class.rds")
neuron_settings$runtime["fontconfig_file"] <- list(fontconfig_file)
analysis_settings <- list(
  `run_id` = "20261101_Manuscript_reanalysis_08",
  `seed` = 42L,
  `threads` = 4L,
  `font` = "Arial",
  `backend` = "R",
  `formats` = list(
    "png",
    "pdf"
  ),
  `expected_counts` = list(
    `all` = 131480L,
    `Excitatory` = 62991L,
    `Inhibitory` = 21693L,
    `Inh_Etv1` = 8139L,
    `Inh_Lamp5` = 1085L
  ),
  `annotation` = "User-approved six major classes, six Exc subtypes and eight Inh working groups; retain original uncertainty and lineage-conflict flags",
  `conditions` = list(
    "WD",
    "Ext",
    "Ret-Ext"
  ),
  `libraries_per_condition` = 3L,
  `rats_per_library` = 5L,
  `IEG_genes` = list(
    "Arc",
    "Bdnf",
    "Btg2",
    "Fos",
    "Fosl2",
    "Homer1",
    "Npas4",
    "Nr4a1"
  ),
  `IEG_quantile` = 0.95,
  `IEG_quantile_type` = 7L,
  `IEG_comparison` = ">=",
  `IEG_scope` = "within each neuron class across WD,Ext,Ret-Ext, using corrected RNA/data; retain ties",
  `contrasts` = list(
    `RetExt_vs_Ext` = list(
      "Ret-Ext",
      "Ext"
    ),
    `Ext_vs_WD` = list(
      "Ext",
      "WD"
    ),
    `RetExt_vs_WD` = list(
      "Ret-Ext",
      "WD"
    )
  ),
  `base` = NULL,
  `upstream_stage2` = "stage2"
)
analysis_settings$base <- dirname(stage2_dir)
analysis_settings$upstream_stage2 <- basename(stage2_dir)

# 5. Helper functions
save_settings <- function(value, filename) {
  path <- file.path(result_dir, "provenance", filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    previous <- jsonlite::read_json(path, simplifyVector = FALSE)
    if (!isTRUE(all.equal(previous, value, tolerance = 0))) stop("Existing settings differ: ", path)
  } else {
    jsonlite::write_json(value, path, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = 17)
  }
  path
}

# 6. Analysis and output
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
config_snapshot <- save_settings(analysis_settings, "analysis_settings.json")
# Recreate the study's accepted cluster-6 annotation change from stage-2 results.
# The input label, size and output counts are verified before writing the mapping.
root <- normalizePath(stage2_dir, mustWork = TRUE)
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
