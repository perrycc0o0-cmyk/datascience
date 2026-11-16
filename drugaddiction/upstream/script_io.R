# Command-line paths and configuration for individually runnable upstream scripts.
# This file defines utilities only; it never starts an analysis stage.

.upstream_main <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
.upstream_directory <- dirname(dirname(normalizePath(.upstream_main, mustWork = TRUE)))

upstream_script <- function(stage, name) file.path(.upstream_directory, stage, name)

upstream_option <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, a)
  if (is.na(i)) return(default)
  if (i == length(a) || startsWith(a[[i + 1L]], "--")) stop("Missing value for ", flag)
  a[[i + 1L]]
}

upstream_work_directory <- function(path) {
  if (is.null(path) || !nzchar(path)) stop("An output directory is required; see --help.")
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  code <- normalizePath(dirname(.upstream_directory), winslash = "/", mustWork = TRUE)
  if (tolower(path) == tolower(code) || startsWith(tolower(path), paste0(tolower(code), "/")))
    stop("Keep analysis outputs outside the public code directory.")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

upstream_read_configuration <- function(root, name, supplied = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required to read configuration.")
  saved <- file.path(root, "config", name)
  path <- if (!is.null(supplied)) supplied else if (file.exists(saved)) saved else
    file.path(.upstream_directory, "config", name)
  jsonlite::read_json(path, simplifyVector = FALSE)
}

upstream_save_configuration <- function(root, name, cfg) {
  path <- file.path(root, "config", name)
  if (file.exists(path)) {
    previous <- jsonlite::read_json(path, simplifyVector = FALSE)
    if (!isTRUE(all.equal(previous, cfg, tolerance = 0))) stop("Existing run configuration differs: ", path,
      "; use the original paths and settings or a new output directory.")
  } else {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(cfg, path, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = 17)
  }
  path
}

upstream_auxiliary_configuration <- function(root, name) {
  path <- file.path(root, "config", name)
  if (!file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(.upstream_directory, "config", name), path, overwrite = FALSE))
      stop("Could not initialize configuration: ", path)
  }
  path
}

upstream_runtime_paths <- function(cfg) {
  font <- upstream_option("--fontconfig-file")
  if (!is.null(font)) cfg$runtime$fontconfig_file <- normalizePath(font, mustWork = TRUE)
  cfg
}

upstream_processing_configuration <- function(root, supplied = NULL) {
  cfg <- upstream_read_configuration(root, "processing_config.json", supplied)
  input <- upstream_option("--input-dir")
  if (!is.null(input)) cfg$input_root <- normalizePath(input, winslash = "/", mustWork = TRUE)
  if (is.null(cfg$input_root) || !nzchar(cfg$input_root))
    stop("Supply --input-dir with the original nine-library matrix directory.")
  cfg <- upstream_runtime_paths(cfg)
  path <- upstream_save_configuration(root, "processing_config.json", cfg)
  upstream_auxiliary_configuration(root, "sample_naming.json")
  path
}

upstream_neuron_configuration <- function(root, supplied = NULL) {
  cfg <- upstream_read_configuration(root, "stage2_config.json", supplied)
  parent <- upstream_option("--stage1", cfg$parent_stage1_root)
  if (is.null(parent)) parent <- file.path(dirname(root), "stage1")
  cfg$parent_stage1_root <- normalizePath(parent, winslash = "/", mustWork = TRUE)
  if (!is.null(upstream_option("--stage1")) || is.null(cfg$input_rds))
    cfg$input_rds <- file.path(cfg$parent_stage1_root, "05_six_class_review", "global_celltype_six_class.rds")
  cfg <- upstream_runtime_paths(cfg)
  path <- upstream_save_configuration(root, "stage2_config.json", cfg)
  upstream_auxiliary_configuration(root, "subtype_markers_and_style.json")
  path
}

upstream_analysis_configuration <- function(root, supplied = NULL) {
  cfg <- upstream_read_configuration(root, "analysis_contract.json", supplied)
  parent <- upstream_option("--stage2")
  if (!is.null(parent)) {
    parent <- normalizePath(parent, winslash = "/", mustWork = TRUE)
    cfg$base <- dirname(parent)
    cfg$upstream_stage2 <- basename(parent)
  } else if (is.null(cfg$base)) {
    cfg$base <- dirname(root)
  }
  upstream_save_configuration(root, "analysis_contract.json", cfg)
}
