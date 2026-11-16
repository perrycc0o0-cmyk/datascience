# 03 plot neuron review

# 1. Packages
library(jsonlite, include.only = character())
library(ggplot2, include.only = character())
library(ggrepel, include.only = character())
library(patchwork, include.only = character())
library(ragg, include.only = character())
library(Matrix, include.only = character())
library(scales, include.only = character())

# 2. Paths
input_dir <- "/path/to/matrices"
work_dir <- "/path/to/work"
result_dir <- file.path(work_dir, "stage2")
figure_dir <- result_dir
fontconfig_file <- NULL
script_path <- "/path/to/03_plot_neuron_review.R"
stage1_dir <- file.path(work_dir, "stage1")
output_name <- "04_marker_review"

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
marker_settings <- list(
  `schema_version` = "1.0",
  `purpose` = "Fresh neuronal subcluster marker review; existing label vocabulary without any historical numeric cluster mapping",
  `reference_doi` = "10.1038/s41467-024-50465-z",
  `reference_notes_physical_pages` = list(
    54L,
    58L,
    59L,
    60L
  ),
  `annotation_is_reference_label_transfer` = FALSE,
  `annotation_is_probability` = FALSE,
  `display_excluded_genes` = list(
    "Reln",
    "Slco1c1",
    "Slc30a3"
  ),
  `lineage_review_markers` = list(
    "Slc17a6",
    "Gad1",
    "Gad2",
    "Slc17a7",
    "Snap25"
  ),
  `review_colors` = list(
    `Unknown` = "#6B7280"
  ),
  `lineages` = list(
    `Exc` = list(
      `name` = "Excitatory",
      `display_scheme` = "B",
      `order` = list(
        "GLU_L6_CLA",
        "GLU_L5_6",
        "GLU_L2_L3",
        "GLU_L5",
        "GLU_L4_5",
        "GLU_ADGM"
      ),
      `display_panels` = list(
        `GLU_L6_CLA` = list(
          "Gng2",
          "Cck"
        ),
        `GLU_L5_6` = list(
          "Syt6",
          "Thsd7b"
        ),
        `GLU_L2_L3` = list(
          "Glis3",
          "Plk5"
        ),
        `GLU_L5` = list(
          "Parm1",
          "Trpc7"
        ),
        `GLU_L4_5` = list(
          "Il1rapl2",
          "Rorb"
        ),
        `GLU_ADGM` = list(
          "Rerg",
          "Arhgap6"
        )
      ),
      `annotation_panels` = list(
        `GLU_L6_CLA` = list(
          `genes` = list(
            "Synpr",
            "Nr4a2",
            "Gng2"
          ),
          `use_for_identity_signature` = TRUE
        ),
        `GLU_L5_6` = list(
          `genes` = list(
            "Syt6",
            "Foxp2",
            "Tle4",
            "Thsd7b"
          ),
          `use_for_identity_signature` = TRUE
        ),
        `GLU_L2_L3` = list(
          `genes` = list(
            "Glis3",
            "Rasgrf2",
            "Ccbe1",
            "Otof"
          ),
          `use_for_identity_signature` = TRUE
        ),
        `GLU_L5` = list(
          `genes` = list(
            "Parm1",
            "Serpine2",
            "Bcl11b"
          ),
          `use_for_identity_signature` = TRUE
        ),
        `GLU_L4_5` = list(
          `genes` = list(
            "Rorb",
            "Thsd7a",
            "Il1rapl2"
          ),
          `use_for_identity_signature` = TRUE
        ),
        `GLU_ADGM` = list(
          `genes` = list(
            "Arhgap6",
            "Rerg",
            "Pkp2",
            "Rfx3"
          ),
          `use_for_identity_signature` = TRUE
        )
      ),
      `context_markers` = list(
        `GLU_L6_CLA` = list(
          "Gng2",
          "Synpr",
          "Nr4a2",
          "Lxn"
        ),
        `GLU_L5_6` = list(
          "Syt6",
          "Foxp2",
          "Tle4",
          "Thsd7b",
          "Vwc2l",
          "Tshz2"
        ),
        `GLU_L2_L3` = list(
          "Glis3",
          "Otof",
          "Rasgrf2",
          "Ccbe1",
          "Cux2"
        ),
        `GLU_L5` = list(
          "Parm1",
          "Serpine2",
          "Bcl11b"
        ),
        `GLU_L4_5` = list(
          "Il1rapl2",
          "Rorb",
          "Thsd7a"
        ),
        `GLU_ADGM` = list(
          "Rerg",
          "Arhgap6",
          "Pkp2",
          "Rfx3"
        )
      ),
      `colors` = list(
        `GLU_L6_CLA` = "#8DD3C7",
        `GLU_L5_6` = "#FFFFB3",
        `GLU_L2_L3` = "#BEBADA",
        `GLU_L5` = "#FB8072",
        `GLU_L4_5` = "#80B1D3",
        `GLU_ADGM` = "#FDB462",
        `Unknown` = "#6B7280"
      ),
      `scientific_cautions` = list(
        "Cck, Plk5 and Trpc7 are descriptive display candidates, not independent identity anchors.",
        "The L6/CLA-associated panel does not prove anatomical claustrum identity or a specific projection.",
        "Vwc2l and Tshz2 describe an adjacent broad deep-layer program; the whole GLU_L5_6 label is not automatically corticothalamic.",
        "Cux2 and Lxn are review context rather than additional author-curated primary identity seeds; missing or weak expression remains visible in source data."
      )
    ),
    `Inh` = list(
      `name` = "Inhibitory",
      `display_scheme` = "A",
      `order` = list(
        "GABA_Amb",
        "GABA_Etv1",
        "GABA_Lamp5",
        "GABA_MSN",
        "GABA_Sst",
        "GABA_Hs3st2",
        "GABA_Pvalb",
        "GABA_Vip"
      ),
      `display_panels` = list(
        `GABA_Amb` = list(
          "Slc17a6"
        ),
        `GABA_Etv1` = list(
          "Etv1"
        ),
        `GABA_Lamp5` = list(
          "Lamp5",
          "Jam2"
        ),
        `GABA_MSN` = list(
          "Rgs9",
          "Rarb"
        ),
        `GABA_Sst` = list(
          "Sst"
        ),
        `GABA_Hs3st2` = list(
          "Hs3st2"
        ),
        `GABA_Pvalb` = list(
          "Pvalb",
          "Slit2"
        ),
        `GABA_Vip` = list(
          "Vip"
        )
      ),
      `annotation_panels` = list(
        `GABA_Amb` = list(
          `genes` = list(
            "Slc17a6",
            "Glra3"
          ),
          `use_for_identity_signature` = FALSE
        ),
        `GABA_Etv1` = list(
          `genes` = list(
            "Etv1",
            "Pde3a",
            "Eya1",
            "Dpf3"
          ),
          `use_for_identity_signature` = TRUE
        ),
        `GABA_Lamp5` = list(
          `genes` = list(
            "Jam2",
            "Lamp5",
            "Sv2c",
            "Reln"
          ),
          `use_for_identity_signature` = TRUE
        ),
        `GABA_MSN` = list(
          `genes` = list(
            "Rgs9",
            "Rarb",
            "Gng7"
          ),
          `use_for_identity_signature` = TRUE
        ),
        `GABA_Sst` = list(
          `genes` = list(
            "Sst",
            "Nos1",
            "Npy",
            "Chodl"
          ),
          `use_for_identity_signature` = TRUE
        ),
        `GABA_Hs3st2` = list(
          `genes` = list(
            "Col5a2",
            "Hs3st2",
            "Vipr2"
          ),
          `use_for_identity_signature` = TRUE
        ),
        `GABA_Pvalb` = list(
          `genes` = list(
            "Pvalb",
            "Slit2",
            "Eya4",
            "Esrrg",
            "Rspo2"
          ),
          `use_for_identity_signature` = TRUE
        ),
        `GABA_Vip` = list(
          `genes` = list(
            "Vip"
          ),
          `use_for_identity_signature` = TRUE
        )
      ),
      `context_markers` = list(
        `GABA_Amb` = list(
          "Slc17a6",
          "Glra3",
          "Gad1",
          "Gad2",
          "Slc17a7"
        ),
        `GABA_Etv1` = list(
          "Etv1",
          "Pde3a",
          "Eya1",
          "Dpf3"
        ),
        `GABA_Lamp5` = list(
          "Lamp5",
          "Reln",
          "Jam2",
          "Sv2c"
        ),
        `GABA_MSN` = list(
          "Gng7",
          "Rgs9",
          "Rarb",
          "Ppp1r1b",
          "Drd2"
        ),
        `GABA_Sst` = list(
          "Sst",
          "Nos1",
          "Npy",
          "Chodl",
          "Grik3",
          "Elfn1",
          "Lypd6"
        ),
        `GABA_Hs3st2` = list(
          "Hs3st2",
          "Col5a2",
          "Vipr2"
        ),
        `GABA_Pvalb` = list(
          "Pvalb",
          "Slit2",
          "Eya4",
          "Rspo2",
          "Esrrg"
        ),
        `GABA_Vip` = list(
          "Vip"
        )
      ),
      `colors` = list(
        `GABA_Amb` = "#66C2A5",
        `GABA_Etv1` = "#FC8D62",
        `GABA_Lamp5` = "#8DA0CB",
        `GABA_MSN` = "#E78AC3",
        `GABA_Sst` = "#A6D854",
        `GABA_Hs3st2` = "#FFD92F",
        `GABA_Pvalb` = "#E5C494",
        `GABA_Vip` = "#FF6347",
        `Unknown` = "#6B7280"
      ),
      `scientific_cautions` = list(
        "GABA_Amb is an unresolved display container; Slc17a6 is not proof of inhibitory identity and this panel does not compete as a validated signature.",
        "Retain the original global cluster26 and lineage-conflict flags for explicit review after fresh subclustering.",
        "Drd2 and Ppp1r1b support only part of the reference MSN program, not every MSN sourcecluster.",
        "Grik3, Elfn1 and Lypd6 describe a neighboring Sst-associated reference program, not the full primary Sst program.",
        "Slit2 alone does not independently prove Pvalb identity; Vip has no invented second exclusive author-curated marker.",
        "Reln is omitted only from displayed panels, not removed from RNA or annotation evidence."
      )
    )
  ),
  `style` = list(
    `backend` = "R",
    `font_family` = "Arial",
    `png_dpi` = 600L,
    `font_sizes` = list(
      `base` = 8L,
      `axis` = 7.5,
      `title` = 9L,
      `legend` = 6.5,
      `tag` = 10L
    ),
    `colors` = list(
      "#EEEEEE",
      "#E5E5E5",
      "#F4ACB8",
      "#DA304F",
      "#A60F2D"
    ),
    `color_knots` = list(
      -2.5,
      0L,
      0.5,
      1.5,
      2.5
    ),
    `color_limits` = list(
      -2.5,
      2.5
    ),
    `color_space` = "Lab",
    `point_area_mapping` = "area proportional to percent expressed",
    `max_group_point_mm` = 6.5,
    `max_cluster_point_mm` = 4L,
    `hide_below_detection_percent` = 5L,
    `size_limits` = list(
      0L,
      100L
    ),
    `size_breaks` = list(
      25L,
      50L,
      75L
    ),
    `umap_label_mm` = 2.6,
    `umap_raster_pointsize` = 1.2,
    `umap_raster_pixels` = 1200L,
    `umap_vector_point_mm` = 0.09,
    `numeric_cluster_ramp` = list(
      "#1E3A8A",
      "#DC2626",
      "#059669",
      "#7C2D12",
      "#6D28D9",
      "#0F766E",
      "#BE185D",
      "#374151",
      "#D97706",
      "#0369A1",
      "#4F46E5",
      "#A16207"
    ),
    `single_panels_untagged` = TRUE,
    `combined_tags` = "uppercase"
  ),
  `scientific_cautions` = list(
    "The reference is a rat brain single-nucleus study, not a normal IL-specific truth set or executed label transfer.",
    "Copy marker vocabulary and colors only; do not transfer historical numeric cluster assignments to fresh clusters.",
    "Candidate annotation remains pending human review; retain Unknown when evidence is insufficient or conflicting.",
    "Marker selection and clustering in the same data are descriptive, not independent anatomical or projection validation."
  )
)

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
config_snapshot <- save_settings(neuron_settings, "neuron_settings.json")
marker_snapshot <- save_settings(marker_settings, "marker_settings.json")
# Render saved neuronal-subcluster annotation review only; never fit or relabel.

nr_flat <- function(x) as.character(unlist(x, use.names = FALSE))
nr_named <- function(x) unlist(x, use.names = TRUE)
nr_columns <- function(x, needed, label) {
  absent <- setdiff(needed, names(x))
  if (length(absent)) stop(label, " lacks fields: ", paste(absent, collapse = ", "))
}
nr_sort_ids <- function(x) {
  x <- unique(as.character(x))
  if (all(grepl("^[0-9]+$", x))) x[order(as.numeric(x))] else sort(x)
}
nr_warning_present <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  z <- trimws(tolower(as.character(x)))
  !is.na(z) & nzchar(z) & !z %in% c("false", "0", "none", "no", "not_flagged")
}

neuron_review_plot_main <- function(root, runtime_config = NULL, output_name = "04_marker_review") {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
  if (!is.null(runtime_config)) {
    rt <- neuron_settings
    if (!is.null(rt$runtime$r_libraries)) .libPaths(unique(c(nr_flat(rt$runtime$r_libraries), .libPaths())))
    if (!is.null(rt$runtime$fontconfig_file)) Sys.setenv(FONTCONFIG_FILE = rt$runtime$fontconfig_file)
  }
  needed <- c("ggplot2", "ggrepel", "patchwork", "ragg", "Matrix", "scales", "jsonlite")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing R plotting packages: ", paste(absent, collapse = ", "))
  if (!isTRUE(capabilities("cairo"))) stop("Cairo is required for editable PDF text")
  fc <- Sys.which("fc-match")
  if (!nzchar(fc)) stop("fc-match is required to verify Arial")
  font <- system2(fc, c("-f", shQuote("%{family}\t%{file}"), "Arial"), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(font, "status")) || !length(font) ||
      !grepl("(^|[,;])\\s*Arial\\s*([,;]|$)", strsplit(font[[1L]], "\t", fixed = TRUE)[[1L]][1L], ignore.case = TRUE))
    stop("Arial was not resolved: ", paste(font, collapse = " "))
  if (requireNamespace("showtext", quietly = TRUE)) showtext::showtext_auto(FALSE)
  config <- file.path(root, "config", "subtype_markers_and_style.json")
  cfg <- marker_settings
  style <- cfg$style
  if (!identical(style$backend, "R")) stop("This renderer is R-only")
  if (length(output_name) != 1L || !grepl("^[A-Za-z0-9_]+$", output_name)) stop("Invalid output directory")
  out <- file.path(root, output_name)
  existing <- if (dir.exists(out)) list.files(out, recursive = TRUE, full.names = TRUE, pattern = "\\.(png|pdf)$") else character()
  if (file.exists(file.path(out, "complete.json")) || length(existing))
    stop("Review exports already exist; select an explicit new output directory")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  nums <- function(x) as.numeric(unlist(x, use.names = FALSE))
  excluded <- nr_flat(cfg$display_excluded_genes)
  fs <- style$font_sizes
  palette <- nr_flat(style$colors)
  knots <- nums(style$color_knots)
  color_limits <- nums(style$color_limits)
  palette_values <- (knots - color_limits[1L]) / diff(color_limits)
  size_limits <- nums(style$size_limits)
  cutoff <- as.numeric(style$hide_below_detection_percent)
  raster <- requireNamespace("scattermore", quietly = TRUE)
  base_theme <- function(umap = FALSE) {
    p <- ggplot2::theme_classic(base_family = "Arial", base_size = as.numeric(fs$base)) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = as.numeric(fs$title), hjust = .5),
        plot.subtitle = ggplot2::element_text(size = as.numeric(fs$legend), hjust = .5),
        axis.text = ggplot2::element_text(size = as.numeric(fs$axis), color = "black"),
        axis.title = ggplot2::element_text(size = as.numeric(fs$axis)),
        axis.line = ggplot2::element_line(linewidth = .3, color = "black"),
        axis.ticks = ggplot2::element_line(linewidth = .3, color = "black"),
        legend.text = ggplot2::element_text(size = as.numeric(fs$legend)),
        legend.title = ggplot2::element_text(size = as.numeric(fs$legend), face = "bold"),
        panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.key = ggplot2::element_rect(fill = "white", color = NA),
        plot.margin = ggplot2::margin(5, 5, 5, 5), plot.tag = ggplot2::element_blank())
    if (umap) p <- p + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
    p
  }
  point_layer <- function() {
    if (raster) scattermore::geom_scattermore(pointsize = as.numeric(style$umap_raster_pointsize),
      pixels = rep(as.integer(style$umap_raster_pixels), 2L))
    else ggplot2::geom_point(size = as.numeric(style$umap_vector_point_mm), stroke = 0)
  }
  tag_theme <- ggplot2::theme(plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = as.numeric(fs$tag)))
  combine <- function(plots, ncol) {
    patchwork::wrap_plots(plots, ncol = ncol) +
      patchwork::plot_annotation(tag_levels = "A", theme = tag_theme) & tag_theme
  }
  all_status <- list()
  for (lineage in c("Exc", "Inh")) {
    spec <- cfg$lineages[[lineage]]
    input <- file.path(root, "03_annotation", paste0(lineage, "_marker_review_payload.rds"))
    if (!file.exists(input)) stop("Missing annotation review payload: ", lineage)
    protected <- c(input, config)
    hashes_before <- tools::md5sum(protected)
    cache <- readRDS(input)
    nr_columns(cache, c("expression", "coordinates", "marker_stats_by_cluster", "marker_stats_by_subtype",
      "proposals", "display_panels", "subtype_order", "colors"), paste(lineage, "payload"))
    d <- cache$coordinates
    nr_columns(d, c("cell_id", "sample_id", "condition_label", "subcluster", "subtype_proposed", "UMAP_1", "UMAP_2",
      "original_global_cluster", "initial_major_class", "forced_display_only", "vascular_review_flag", "lineage_identity_warning"),
      paste(lineage, "coordinates"))
    if (!identical(colnames(cache$expression), as.character(d$cell_id)) || anyDuplicated(d$cell_id) ||
        any(!is.finite(d$UMAP_1)) || any(!is.finite(d$UMAP_2))) stop("Invalid marker/coordinate alignment")
    if (!inherits(cache$expression, "sparseMatrix")) stop("Expected sparse normalized marker expression")
    if (any(!is.finite(cache$expression@x)) || any(cache$expression@x < 0)) stop("Invalid normalized marker values")
    original_d <- d
    xy_before <- d[, c("cell_id", "UMAP_1", "UMAP_2")]
    raw_ids <- nr_sort_ids(d$subcluster)
    prop <- cache$proposals
    nr_columns(prop, c("subcluster", "subtype_proposed", "n_cells", "annotation_status", "review_notes"), "proposals")
    if (anyDuplicated(prop$subcluster) || !setequal(as.character(prop$subcluster), raw_ids)) stop("Proposal/cluster mismatch")
    expected_labels <- as.character(prop$subtype_proposed[match(as.character(d$subcluster), as.character(prop$subcluster))])
    if (!identical(expected_labels, as.character(d$subtype_proposed))) stop("Cell/proposal labels disagree")
    desired_order <- nr_flat(spec$order)
    declared_order <- nr_flat(cache$subtype_order)
    known_present <- desired_order[desired_order %in% unique(d$subtype_proposed)]
    if (!identical(declared_order[declared_order %in% desired_order], desired_order[desired_order %in% declared_order]))
      stop("Payload changes approved subtype order")
    subtype_order <- unique(c(known_present, declared_order[declared_order %in% unique(d$subtype_proposed)]))
    subtype_order <- c(subtype_order, setdiff(unique(as.character(d$subtype_proposed)), subtype_order))
    colors <- nr_named(spec$colors)
    colors <- c(colors, nr_named(cfg$review_colors)[!names(nr_named(cfg$review_colors)) %in% names(colors)])
    payload_colors <- nr_named(cache$colors)
    shared_colors <- intersect(names(payload_colors), names(colors))
    if (length(shared_colors) && any(toupper(payload_colors[shared_colors]) != toupper(colors[shared_colors])))
      stop("Payload changes original subtype palette")
    if (length(setdiff(subtype_order, names(colors)))) stop("No declared color for subtype: ", paste(setdiff(subtype_order, names(colors)), collapse = ","))
    colors <- colors[subtype_order]
    expected_panels <- lapply(spec$display_panels, nr_flat)
    supplied_panels <- lapply(cache$display_panels, nr_flat)
    if (!identical(expected_panels, supplied_panels)) stop("Payload changes the selected B/A display panel")
    requested_display <- unique(nr_flat(spec$display_panels))
    requested_context <- unique(c(nr_flat(spec$context_markers),
      nr_flat(lapply(spec$annotation_panels, function(v) v$genes)), nr_flat(cfg$lineage_review_markers)))
    gene_names <- rownames(cache$expression)
    genes <- intersect(setdiff(requested_display, excluded), gene_names)
    context <- intersect(setdiff(requested_context, excluded), gene_names)
    if (!length(genes) || !length(context)) stop("No available display or context genes")
    dest <- file.path(out, lineage)
    for (p in c("01_UMAP", "02_DotPlot", "source_data", "provenance")) dir.create(file.path(dest, p), recursive = TRUE, showWarnings = FALSE)
    write.csv(data.frame(gene = unique(c(requested_display, requested_context)),
      available = unique(c(requested_display, requested_context)) %in% gene_names,
      excluded_from_display_only = unique(c(requested_display, requested_context)) %in% excluded),
      file.path(dest, "source_data", "marker_availability.csv"), row.names = FALSE)
    con <- gzfile(file.path(dest, "source_data", "coordinates_and_review_flags.csv.gz"), "wt")
    write.csv(original_d, con, row.names = FALSE); close(con)
    write.csv(prop, file.path(dest, "source_data", "annotation_proposals.csv"), row.names = FALSE)
    manifest <- list()
    save_plot <- function(p, stem, width, height) {
      path <- file.path(dest, stem)
      grDevices::cairo_pdf(paste0(path, ".pdf"), width = width / 25.4, height = height / 25.4, family = "Arial", bg = "white")
      tryCatch(print(p), finally = grDevices::dev.off())
      ragg::agg_png(paste0(path, ".png"), width = width, height = height, units = "mm", res = as.integer(style$png_dpi), background = "white")
      tryCatch(print(p), finally = grDevices::dev.off())
      manifest[[length(manifest) + 1L]] <<- data.frame(stem = stem, width_mm = width, height_mm = height,
        png_dpi = as.integer(style$png_dpi), pdf_editable_text = TRUE, stringsAsFactors = FALSE)
    }
    add_arrows <- function(p) {
      xr <- range(d$UMAP_1); yr <- range(d$UMAP_2); dx <- diff(xr); dy <- diff(yr)
      if (dx <= 0 || dy <= 0) stop("UMAP has a zero-range axis")
      x0 <- xr[1L] - .05 * dx; y0 <- yr[1L] - .05 * dy
      p + ggplot2::annotate("segment", x = x0, xend = x0 + .15 * dx, y = y0, yend = y0,
        linewidth = .35, color = "black", arrow = grid::arrow(length = grid::unit(1.7, "mm"), type = "closed")) +
        ggplot2::annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + .15 * dy,
          linewidth = .35, color = "black", arrow = grid::arrow(length = grid::unit(1.7, "mm"), type = "closed")) +
        ggplot2::annotate("text", x = x0 + .03 * dx, y = y0 - .03 * dy, label = "UMAP-1", hjust = 0, vjust = 1,
          family = "Arial", size = 2) +
        ggplot2::annotate("text", x = x0 - .03 * dx, y = y0 + .12 * dy, label = "UMAP-2", angle = 90, hjust = 1,
          vjust = 0, family = "Arial", size = 2) + ggplot2::coord_fixed(clip = "off")
    }
    umap <- function(group, palette, title = NULL, label = TRUE) {
      dd <- d; dd$plot_group <- factor(as.character(dd[[group]]), levels = names(palette))
      if (anyNA(dd$plot_group)) stop("Unmapped UMAP group: ", group)
      dd <- dd[order((seq_len(nrow(dd)) * 104729) %% 1000003), , drop = FALSE]
      p <- ggplot2::ggplot(dd, ggplot2::aes(UMAP_1, UMAP_2, color = plot_group)) + point_layer() +
        ggplot2::scale_color_manual(values = palette, drop = FALSE) + ggplot2::labs(title = title, color = NULL) + base_theme(TRUE)
      if (label) {
        lab <- stats::aggregate(dd[, c("UMAP_1", "UMAP_2")], list(plot_group = dd$plot_group), stats::median)
        p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = plot_group), color = "black", family = "Arial",
          size = as.numeric(style$umap_label_mm), seed = 42L, max.overlaps = Inf, min.segment.length = 0,
          segment.size = .2, show.legend = FALSE) + ggplot2::theme(legend.position = "none")
      } else p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1.7)))
      add_arrows(p)
    }
    ramp <- nr_flat(style$numeric_cluster_ramp)
    raw_colors <- stats::setNames(grDevices::colorRampPalette(ramp)(length(raw_ids)), raw_ids)
    pA <- umap("subcluster", raw_colors, paste(spec$name, "subclusters"))
    pB <- umap("subtype_proposed", colors, NULL)
    save_plot(pA, "01_UMAP/A_numeric_clusters", 110, 100)
    save_plot(pB, "01_UMAP/B_proposed_subtypes", 110, 100)
    sample_order <- c(paste0("WD", 1:3), paste0("Ext", 1:3), paste0("Ret-Ext", 1:3))
    if (!all(as.character(d$sample_id) %in% sample_order) ||
        !all(as.character(d$condition_label) %in% c("WD", "Ext", "Ret-Ext"))) stop("Canonical sample/condition names are invalid")
    sample_colors <- stats::setNames(c("#377EB8", "#66A9D3", "#A6CEE3", "#984EA3", "#C08FD0", "#DCC6E7", "#E41A1C", "#F4777F", "#FBB4AE"), sample_order)
    sample_colors <- sample_colors[sample_order %in% d$sample_id]
    condition_colors <- c(WD = "#888888", Ext = "#4BBAD3", `Ret-Ext` = "#EE4B35")
    condition_colors <- condition_colors[names(condition_colors) %in% d$condition_label]
    pE <- umap("sample_id", sample_colors, "Library", FALSE)
    pF <- umap("condition_label", condition_colors, "Condition", FALSE)
    parent_ids <- nr_sort_ids(d$original_global_cluster)
    parent_colors <- stats::setNames(grDevices::colorRampPalette(ramp)(length(parent_ids)), parent_ids)
    pG <- umap("original_global_cluster", parent_colors, "Parent global clusters")
    if (lineage == "Inh" && "26" %in% parent_ids)
      pG <- pG + ggplot2::labs(subtitle = "Global cluster 26 retained for lineage-identity review")
    save_plot(pE, "01_UMAP/E_sample", 140, 110)
    save_plot(pF, "01_UMAP/F_condition", 125, 110)
    save_plot(pG, "01_UMAP/G_parent_global_cluster", 110, 100)
    stats_c <- cache$marker_stats_by_cluster
    stats_d <- cache$marker_stats_by_subtype
    for (tt in list(stats_c, stats_d)) {
      nr_columns(tt, c("gene", "cluster", "n_cells", "pct_expressed", "avg_expression_scaled"), "marker statistics")
      if (any(!is.finite(tt$pct_expressed)) || any(tt$pct_expressed < 0 | tt$pct_expressed > 100) ||
          any(!is.finite(tt$avg_expression_scaled)) || any(tt$avg_expression_scaled < color_limits[1L] - 1e-10 | tt$avg_expression_scaled > color_limits[2L] + 1e-10))
        stop("Invalid saved marker statistic or color scale")
    }
    if (!setequal(as.character(unique(stats_c$cluster)), raw_ids) ||
        !setequal(as.character(unique(stats_d$cluster)), unique(as.character(d$subtype_proposed)))) stop("Marker statistics have different groups")
    raw_labels <- stats::setNames(paste0(prop$subcluster, " | ", prop$subtype_proposed), as.character(prop$subcluster))
    stats_c$source_group <- as.character(stats_c$cluster)
    stats_c$cluster <- unname(raw_labels[stats_c$source_group])
    stats_d$source_group <- as.character(stats_d$cluster)
    dotplot <- function(tab, gene_order, group_order, title, max_size, id) {
      dd <- tab[tab$gene %in% gene_order, , drop = FALSE]
      dd$gene <- factor(dd$gene, levels = gene_order); dd$cluster <- factor(dd$cluster, levels = group_order)
      if (anyNA(dd$gene) || anyNA(dd$cluster) || any(as.character(dd$gene) %in% excluded)) stop("DotPlot gene/group mismatch")
      dd$pct_display <- ifelse(dd$pct_expressed < cutoff, NA_real_, dd$pct_expressed)
      p <- ggplot2::ggplot(dd, ggplot2::aes(gene, cluster, size = pct_display, color = avg_expression_scaled)) +
        ggplot2::geom_point(shape = 16, stroke = 0, na.rm = TRUE) +
        ggplot2::scale_color_gradientn(colors = palette, values = palette_values, limits = color_limits,
          breaks = c(-2, -1, 0, 1, 2), name = "Average Expression") +
        ggplot2::scale_size_area(max_size = max_size, limits = size_limits, breaks = nums(style$size_breaks), name = "Percent Expressed") +
        ggplot2::labs(title = title, x = NULL, y = NULL) + base_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = .5, face = "italic"),
          axis.text.y = ggplot2::element_text(face = "bold"), legend.position = "bottom", legend.box = "horizontal") +
        ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top", direction = "horizontal",
          barwidth = grid::unit(28, "mm"), barheight = grid::unit(2.2, "mm")), size = ggplot2::guide_legend(title.position = "top", nrow = 1))
      actual <- ggplot2::ggplot_build(p)$data[[1L]]
      if (nrow(actual) != nrow(dd) ||
          !isTRUE(all.equal(as.numeric(actual$x), as.numeric(dd$gene))) ||
          !isTRUE(all.equal(as.numeric(actual$y), as.numeric(dd$cluster))) ||
          !isTRUE(all.equal(actual$size, max_size * sqrt(dd$pct_display / 100), tolerance = 1e-10)))
        stop("Actual point placement or area differs from the saved source")
      expected_color <- scales::gradient_n_pal(palette, values = palette_values)((dd$avg_expression_scaled - color_limits[1L]) / diff(color_limits))
      if (!identical(toupper(actual$colour), toupper(expected_color))) stop("Actual color mapping differs from true saved Z")
      source_points <- data.frame(plot_id = id, source_key = paste(as.character(dd$gene), dd$source_group, sep = "::"),
        gene = as.character(dd$gene), group = dd$source_group, group_label = as.character(dd$cluster),
        pct_expressed = dd$pct_expressed, avg_expression_scaled = dd$avg_expression_scaled,
        pct_display = dd$pct_display, plotted_x = as.numeric(actual$x), plotted_y = as.numeric(actual$y),
        plotted_size = actual$size, plotted_colour = actual$colour, plotted_shape = actual$shape,
        plotted_stroke = actual$stroke, stringsAsFactors = FALSE)
      write.csv(source_points, file.path(dest, "source_data", paste0(id, "_actual_points.csv")), row.names = FALSE)
      p
    }
    cluster_height <- max(110, 45 + length(raw_ids) * 5)
    pC <- dotplot(stats_c, genes, rev(unname(raw_labels[raw_ids])), "Selected markers across fresh subclusters",
      as.numeric(style$max_cluster_point_mm), "C")
    pD <- dotplot(stats_d, genes, rev(subtype_order), "Candidate subtype markers",
      as.numeric(style$max_group_point_mm), "D")
    pH <- dotplot(stats_c, context, rev(unname(raw_labels[raw_ids])), "Subtype context and lineage review",
      as.numeric(style$max_cluster_point_mm), "H")
    save_plot(pC, "02_DotPlot/C_by_cluster", 183, cluster_height)
    save_plot(pD, "02_DotPlot/D_by_subtype", 183, if (lineage == "Exc") 100 else 115)
    save_plot(pH, "02_DotPlot/H_context", max(240, 40 + length(context) * 4.5), cluster_height)
    save_plot(combine(list(pA, pB), 2), "01_UMAP/Combined_numeric_and_subtypes", 220, 100)
    save_plot(combine(list(pE, pF), 2), "01_UMAP/Combined_sample_condition", 240, 110)
    write.csv(cache$marker_stats_by_cluster, file.path(dest, "source_data", "saved_marker_stats_by_cluster.csv"), row.names = FALSE)
    write.csv(cache$marker_stats_by_subtype, file.path(dest, "source_data", "saved_marker_stats_by_subtype.csv"), row.names = FALSE)
    write.csv(do.call(rbind, manifest), file.path(dest, "provenance", "figure_export_manifest.csv"), row.names = FALSE)
    if (!identical(original_d, cache$coordinates) || !identical(xy_before, d[, c("cell_id", "UMAP_1", "UMAP_2")]) ||
        !identical(hashes_before, tools::md5sum(protected))) stop("Plotting changed source values, metadata or coordinates")
    status <- list(status = "complete", completed_UTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      lineage = lineage, n_cells = nrow(d), n_subclusters = length(raw_ids), candidate_subtypes = subtype_order,
      annotation_pending = TRUE, annotation_approved = FALSE, plotting_only = TRUE,
      source_values_labels_coordinates_unchanged = TRUE, source_md5 = as.list(hashes_before),
      font = font, backend = "R", png_dpi = as.integer(style$png_dpi), single_panels_untagged = TRUE, combined_labels = "uppercase",
      displayed_markers = genes, displayed_context = context, display_excluded_genes = excluded,
      missing_display_genes = setdiff(setdiff(requested_display, excluded), gene_names),
      missing_context_genes = setdiff(setdiff(requested_context, excluded), gene_names),
      point_size = "Point area proportional to unchanged detection percentage; below 5 percent not displayed",
      color = "Saved gene-wise scaled expression, true limits -2.5 to 2.5; common continuous palette",
      subclustering_recomputed_by_renderer = FALSE, UMAP_recomputed_by_renderer = FALSE,
      parent_global_cluster26_retained = lineage == "Inh" && "26" %in% parent_ids,
      retained_identity_warning_nuclei = sum(nr_warning_present(d$lineage_identity_warning)),
      downstream_started = FALSE, stop_before_downstream = TRUE,
      waiting_for = "User review of fresh Exc/Inh UMAP and marker DotPlots before downstream analysis")
    jsonlite::write_json(status, file.path(dest, "complete.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
    all_status[[lineage]] <- status
    cat(lineage, "marker review plots completed; stop before downstream.\n")
  }
  final <- list(status = "awaiting_user_subtype_marker_review", completed_UTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    lineages = all_status, annotation_approved = FALSE, downstream_started = FALSE, stop_before_downstream = TRUE,
    session_info = capture.output(utils::sessionInfo()))
  jsonlite::write_json(final, file.path(out, "complete.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
  invisible(final)
}

neuron_review_plot_main(result_dir, runtime_config = config_snapshot, output_name = output_name)
