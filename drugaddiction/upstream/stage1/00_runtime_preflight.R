#!/usr/bin/env Rscript
if (sys.nframe() == 0L && any(commandArgs(TRUE) %in% c("--help", "-h"))) {
  cat("Usage: Rscript upstream/stage1/00_runtime_preflight.R --root STAGE1_DIR --input-dir MATRIX_DIR [--config JSON] [--fontconfig-file FILE]\n")
  cat("Run this script directly; supply its saved inputs. No workflow launcher is required.\n")
  quit(save = "no", status = 0L)
}
.upstream_entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(dirname(normalizePath(.upstream_entry, mustWork = TRUE))), "script_io.R"))

args <- commandArgs(trailingOnly = TRUE)
root_arg <- upstream_option('--root', if (length(args) && !startsWith(args[[1L]], '--')) args[[1L]] else NULL)
root <- upstream_work_directory(root_arg)
cfgpath <- upstream_processing_configuration(root, upstream_option('--config'))
if (!requireNamespace('jsonlite', quietly = TRUE)) stop('jsonlite unavailable')
cfg <- jsonlite::read_json(cfgpath, simplifyVector = FALSE)
if (!identical(cfg$stage_gate, 'global_celltype_review')) stop('Stage gate is missing')
libs <- as.character(unlist(cfg$runtime$r_libraries))
.libPaths(unique(c(libs, .libPaths())))
Sys.setenv(OMP_NUM_THREADS = '4', OMP_THREAD_LIMIT = '4', OPENBLAS_NUM_THREADS = '1', MKL_NUM_THREADS = '1')
packages <- c('Seurat','SeuratObject','SoupX','scDblFinder','SingleCellExperiment','Matrix','hdf5r','harmony','jsonlite','ggplot2','ggrepel','patchwork','ragg','systemfonts','uwot','presto','xgboost','future')
ok <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
versions <- setNames(lapply(packages[ok], function(p) as.character(utils::packageVersion(p))), packages[ok])
samples <- vapply(cfg$samples, function(s) s$sample_id, character(1))
inputs <- expand.grid(sample = samples, kind = c('filtered_feature_bc_matrix.h5','raw_feature_bc_matrix.h5'), stringsAsFactors = FALSE)
inputs$path <- file.path(cfg$input_root, inputs$sample, inputs$kind)
inputs$exists <- file.exists(inputs$path)
inputs$bytes <- file.info(inputs$path)$size
parsecheck <- lapply(c('01_process_libraries.R','02_global_celltype.R','03_plot_global_review.R'), function(p) {
  f <- upstream_script('stage1', p)
  if (!file.exists(f)) stop('Script missing: ', p)
  parse(f)
  list(file = p, parsed = TRUE, md5 = unname(tools::md5sum(f)))
})
fonts_conf <- cfg$runtime$fontconfig_file
if (!is.null(fonts_conf)) {
  stopifnot(file.exists(fonts_conf))
  Sys.setenv(FONTCONFIG_FILE = fonts_conf)
}
fontmatch <- system2('fc-match', c('-f', shQuote('%{family}'), shQuote('Arial')), stdout = TRUE, stderr = TRUE)
font_ok <- any(grepl('Arial', fontmatch, fixed = TRUE))
report <- list(status = if (all(ok) && all(inputs$exists) && font_ok) 'pass' else 'failed',
  checked = format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z'), R = R.version.string,
  package_versions = versions, missing_packages = packages[!ok], input_files = inputs,
  parsed_scripts = parsecheck, library_paths = .libPaths(), Arial_match = fontmatch,
  font_ok = font_ok, stage_gate = cfg$stage_gate, downstream_authorized = FALSE)
dir.create(file.path(root,'qa'), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(report, file.path(root,'qa','runtime_preflight.json'), auto_unbox = TRUE, pretty = TRUE, na = 'null')
writeLines(capture.output(sessionInfo()), file.path(root,'qa','preflight_sessionInfo.txt'))
cat(jsonlite::toJSON(report, auto_unbox = TRUE, pretty = TRUE, na = 'null'), '\n')
if (!identical(report$status, 'pass')) quit(status = 2L)
