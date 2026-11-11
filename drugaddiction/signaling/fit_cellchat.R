# Fit abundance-weighted neuronal CellChat networks.
# Entry: Rscript --vanilla fit_cellchat.R --root <new run root>
# source() defines functions only. All result writes remain inside 06_CellChat.

.cc_stop <- function(...) stop(..., call. = FALSE)
.cc_atomic_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".incomplete_", Sys.getpid())
  jsonlite::write_json(x, tmp, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = 17)
  if (!file.rename(tmp, path)) .cc_stop("Atomic JSON rename failed: ", path)
}
.cc_atomic_rds <- function(x, path) {
  if (file.exists(path)) .cc_stop("Refusing to overwrite completed RDS: ", path)
  tmp <- paste0(path, ".incomplete_", Sys.getpid())
  saveRDS(x, tmp, compress = "gzip")
  if (!file.rename(tmp, path)) .cc_stop("Atomic RDS rename failed: ", path)
}
.cc_csv <- function(x, path, row.names = FALSE) {
  utils::write.csv(x, path, row.names = row.names, na = "NA")
}
.cc_hash <- function(x, file = FALSE) digest::digest(x, algo = "sha256", serialize = !file, file = file)
.cc_utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
.cc_defaults <- function(nms) {
  setNames(lapply(nms, function(nm) {
    f <- formals(get(nm, envir = asNamespace("CellChat")))
    setNames(lapply(seq_along(f), function(i) paste(deparse(f[[i]]), collapse = " ")), names(f))
  }), nms)
}
.cc_runtime <- function() {
  Sys.setenv(OMP_NUM_THREADS = "4", OMP_THREAD_LIMIT = "4", RCPP_PARALLEL_NUM_THREADS = "4",
             OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
             NUMEXPR_NUM_THREADS = "1")
  options(mc.cores = 1L)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(4L)
  }
  future::plan(future::sequential)
  set.seed(42L)
}

# Map only frozen, verified one-to-one orthologs. Any same-gene duplicate
# features are aggregated from counts using the full original RNA library size.
.cc_prepare_payload <- function(x, lineage, dest, orthology) {
  need <- c("normalized", "counts", "metadata", "gene_annotation", "IEG_threshold")
  if (!all(need %in% names(x))) .cc_stop("Missing payload fields for ", lineage)
  z <- x$normalized; counts <- x$counts; md <- x$metadata; fm <- x$gene_annotation
  if (!inherits(z, "sparseMatrix") || !inherits(counts, "sparseMatrix")) .cc_stop("Sparse matrices required")
  if (!identical(dim(z), dim(counts)) || !identical(dimnames(z), dimnames(counts))) .cc_stop("Count/data dimnames differ")
  if (any(!is.finite(z@x)) || any(z@x < 0) || any(!is.finite(counts@x)) || any(counts@x < 0)) .cc_stop("Invalid expression values")
  if (any(abs(counts@x - round(counts@x)) > 1e-8)) .cc_stop("Corrected counts must be integer-valued")
  if (!all(c("seurat_feature", "original_symbol", "gene_id") %in% names(fm))) .cc_stop("Gene mapping fields missing")
  if (anyDuplicated(fm$seurat_feature)) .cc_stop("Duplicate seurat_feature in annotation")
  map_ix <- match(rownames(z), fm$seurat_feature)
  if (anyNA(map_ix)) .cc_stop("Feature mapping incomplete")
  fm <- fm[map_ix, , drop = FALSE]
  if (!all(c("Cell", "neuron_type", "pool_id", "condition", "subtype", "IEGhigh") %in% names(md))) .cc_stop("Metadata fields missing")
  if (anyDuplicated(md$Cell) || anyNA(match(colnames(z), md$Cell))) .cc_stop("Cell join incomplete")
  md <- md[match(colnames(z), md$Cell), , drop = FALSE]
  rownames(md) <- md$Cell
  if (!identical(rownames(md), colnames(z))) .cc_stop("Cell order mismatch")
  if (anyNA(md$IEGhigh) || !all(md$IEGhigh)) .cc_stop("Payload must contain exclusively IEG-high nuclei")
  if (!all(md$neuron_type == lineage)) .cc_stop("Lineage mismatch")
  if (anyNA(md$condition) || !all(md$condition %in% c("WD", "Ext", "Ret-Ext"))) .cc_stop("Unexpected condition names")
  expected <- if (lineage == "Excitatory") c("GLU_L2_L3", "GLU_L4_5", "GLU_L5", "GLU_L5_6", "GLU_L6_CLA", "GLU_ADGM") else
    c("GABA_Amb", "GABA_Etv1", "GABA_Lamp5", "GABA_MSN", "GABA_Sst", "GABA_Hs3st2", "GABA_Pvalb", "GABA_Vip")
  if (anyNA(md$subtype) || !all(md$subtype %in% expected)) .cc_stop("Unexpected frozen working labels")
  lib <- Matrix::colSums(counts)
  if (any(!is.finite(lib) | lib <= 0)) .cc_stop("Full corrected-count library sizes must be positive")
  source_normalized_hash <- .cc_hash(z)
  source_counts_hash <- .cc_hash(counts)
  lookup <- orthology$accepted
  raw <- orthology$raw
  if (anyDuplicated(lookup$rat_gene_id) || anyDuplicated(lookup$human_gene_id) || anyDuplicated(lookup$human_symbol) ||
      !all(lookup$orthology_type == "ortholog_one2one")) .cc_stop("Invalid frozen one-to-one map")
  fm$rat_gene_id_stripped <- sub("\\.[0-9]+$", "", as.character(fm$gene_id))
  ix <- match(fm$rat_gene_id_stripped, lookup$rat_gene_id)
  known_id <- fm$rat_gene_id_stripped %in% raw$rat_gene_id
  fm$mapping_method <- ifelse(!is.na(ix), "Ensembl_rat_gene_id", "unmapped")
  unique_raw <- unique(raw[, c("rat_gene_id", "rat_symbol")])
  symbol_n_ids <- vapply(split(unique_raw$rat_gene_id, unique_raw$rat_symbol), function(v) length(unique(v)), integer(1))
  symbol_n_rows <- table(lookup$rat_symbol)
  unique_symbol <- !is.na(lookup$rat_symbol) & nzchar(lookup$rat_symbol) &
    symbol_n_ids[lookup$rat_symbol] == 1L & symbol_n_rows[lookup$rat_symbol] == 1L
  symbol_lookup <- lookup[which(unique_symbol), , drop = FALSE]
  fallback <- which(is.na(ix) & !known_id & !is.na(fm$original_symbol) & nzchar(as.character(fm$original_symbol)))
  fi <- match(as.character(fm$original_symbol[fallback]), symbol_lookup$rat_symbol)
  resolved <- fallback[!is.na(fi)]
  ix[resolved] <- match(symbol_lookup$rat_gene_id[fi[!is.na(fi)]], lookup$rat_gene_id)
  fm$mapping_method[resolved] <- "unique_original_rat_symbol_fallback"
  fm$mapping_status <- ifelse(!is.na(ix), "verified_unique_one2one", ifelse(known_id,
    "known_rat_id_without_usable_unique_one2one", "unmapped_rat_id_and_no_unique_symbol_fallback"))
  fm$orthology_rat_gene_id <- lookup$rat_gene_id[ix]
  fm$human_gene_id <- lookup$human_gene_id[ix]
  fm$human_symbol <- lookup$human_symbol[ix]
  fm$orthology_type <- lookup$orthology_type[ix]
  fm$orthology_confidence <- lookup$orthology_confidence[ix]
  valid <- !is.na(ix)
  if (!any(valid)) .cc_stop("No verified human orthologs in this payload")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  .cc_csv(fm, file.path(dest, "feature_orthology_map_all_features.csv"))
  .cc_csv(as.data.frame(table(fm$mapping_status, fm$mapping_method)), file.path(dest, "mapping_status_counts.csv"))
  n_input_features <- nrow(z)
  z <- z[valid, , drop = FALSE]
  mapped_counts <- counts[valid, , drop = FALSE]
  fm <- fm[valid, , drop = FALSE]
  symbol <- as.character(fm$human_symbol)
  symbols <- unique(symbol)
  first <- match(symbols, symbol)
  n_members <- tabulate(match(symbol, symbols), nbins = length(symbols))
  mat <- z[first, , drop = FALSE]
  rownames(mat) <- symbols
  dup <- which(n_members > 1L)
  detail <- vector("list", length(dup))
  for (k in seq_along(dup)) {
    i <- dup[k]; rows <- which(symbol == symbols[i])
    if (length(unique(fm$orthology_rat_gene_id[rows])) != 1L || length(unique(fm$human_gene_id[rows])) != 1L) .cc_stop("Cannot collapse unrelated genes sharing a symbol")
    combined <- Matrix::colSums(mapped_counts[rows, , drop = FALSE])
    new_value <- log1p(1e4 * combined / lib)
    old_first <- as.numeric(z[rows[1L], ])
    mat[i, ] <- new_value
    detail[[k]] <- data.frame(lineage = lineage, human_symbol = symbols[i], n_input_features = length(rows),
                             seurat_features = paste(fm$seurat_feature[rows], collapse = ";"),
                             original_symbols = paste(fm$original_symbol[rows], collapse = ";"),
                             gene_ids = paste(fm$gene_id[rows], collapse = ";"),
                             n_cells_nonzero_combined_counts = sum(combined > 0),
                             n_cells_changed_from_first_row = sum(abs(new_value - old_first) > 1e-12),
                             max_abs_change_from_first_row = max(abs(new_value - old_first)), stringsAsFactors = FALSE)
  }
  mat <- methods::as(mat, "dgCMatrix")
  fm$human_symbol <- symbol
  fm$n_features_same_symbol <- n_members[match(symbol, symbols)]
  fm$normalization_policy <- ifelse(fm$n_features_same_symbol > 1L,
    "corrected_counts_sum_then_log1p_1e4_over_full_corrected_library_size", "original_normalized_row_unchanged")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  .cc_csv(fm, file.path(dest, "mapped_feature_orthology_map.csv"))
  duplicate_audit <- if (length(detail)) do.call(rbind, detail) else data.frame(
    lineage = character(), human_symbol = character(), n_input_features = integer(), seurat_features = character(),
    original_symbols = character(), gene_ids = character(), n_cells_nonzero_combined_counts = integer(),
    n_cells_changed_from_first_row = integer(), max_abs_change_from_first_row = numeric())
  .cc_csv(duplicate_audit, file.path(dest, "duplicate_symbol_aggregation.csv"))
  .cc_csv(data.frame(cell_id = md$Cell, lineage = lineage, full_corrected_library_size = lib), file.path(dest, "library_sizes.csv"))
  .cc_csv(md, file.path(dest, "IEGhigh_metadata.csv"))
  info <- list(lineage = lineage, n_cells = ncol(z), n_input_features = n_input_features, n_mapped_features = sum(valid), n_unmapped_features = sum(!valid), n_symbols = nrow(mat),
               n_duplicate_symbol_groups = length(dup), n_features_in_duplicate_groups = sum(n_members[dup]),
               IEG_threshold = x$IEG_threshold, normalized_unique_rows_preserved = TRUE,
               duplicate_policy = "sum corrected counts for multiple source features of the same verified rat-human gene, log1p(1e4*sum/FULL ORIGINAL RNA corrected library size)",
               orthology_policy = "Ensembl116 unique one-to-one; exact ratID preferred; original-symbol fallback only for unknown/absent IDs and globally unique rat symbols; known non-one2one IDs never bypassed",
               source_normalized_matrix_hash = source_normalized_hash, source_counts_matrix_hash = source_counts_hash,
               full_RNA_library_denominator_preserved = TRUE, metadata_columns_all_retained = TRUE,
               input_metadata_hash = .cc_hash(x$metadata), output_metadata_hash = .cc_hash(md))
  .cc_atomic_json(info, file.path(dest, "preparation.json"))
  list(mat = mat, meta = md, info = info)
}

.cc_group_counts <- function(md) {
  out <- aggregate(rep(1L, nrow(md)), by = list(group = as.character(md$celltype),
                     lineage = as.character(md$neuron_type), subtype = as.character(md$subtype)), FUN = sum)
  names(out)[4L] <- "n_cells"
  out$eligible_min_cells_10 <- out$n_cells > 10L
  out[order(out$group), , drop = FALSE]
}
.cc_array_csv <- function(prob, pval, path, index_name) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  if (is.null(prob) || length(prob) == 0L || length(dim(prob)) != 3L || any(dim(prob) == 0L)) {
    empty <- data.frame(source = character(), target = character(), entry = character(), probability = numeric())
    names(empty)[3L] <- index_name
    if (!is.null(pval)) empty$p_value <- numeric()
    write.table(empty, con, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
    return(invisible(0L))
  }
  dn <- dimnames(prob)
  if (any(vapply(dn, is.null, logical(1)))) .cc_stop("Probability array needs named dimensions")
  grid <- expand.grid(source = dn[[1L]], target = dn[[2L]], stringsAsFactors = FALSE)
  n_written <- 0L
  for (k in seq_len(dim(prob)[3L])) {
    d <- grid
    d[[index_name]] <- dn[[3L]][k]
    d$probability <- as.vector(prob[, , k, drop = FALSE])
    if (!is.null(pval)) d$p_value <- as.vector(pval[, , k, drop = FALSE])
    write.table(d, con, sep = ",", row.names = FALSE, col.names = k == 1L, quote = TRUE)
    n_written <- n_written + nrow(d)
  }
  invisible(n_written)
}
.cc_export_outgoing <- function(cellchat, run_dir, network, condition) {
  w <- cellchat@net$weight
  if (!is.matrix(w) || any(!is.finite(w))) .cc_stop("Invalid net weight")
  self <- numeric(nrow(w)); names(self) <- rownames(w)
  j <- match(rownames(w), colnames(w)); ok <- !is.na(j)
  self[ok] <- w[cbind(which(ok), j[ok])]
  d <- data.frame(network = network, condition = condition, source = rownames(w),
                  outgoing_including_self = rowSums(w), self_interaction = self,
                  outgoing_excluding_self = rowSums(w) - self, stringsAsFactors = FALSE)
  .cc_csv(d, file.path(run_dir, "total_outgoing.csv"))
  pp <- cellchat@netP$prob
  pd <- list()
  if (!is.null(pp) && length(pp) && length(dim(pp)) == 3L) for (k in seq_len(dim(pp)[3L])) {
    m <- pp[, , k, drop = FALSE]; dim(m) <- dim(pp)[1:2]; dimnames(m) <- dimnames(pp)[1:2]
    ss <- numeric(nrow(m)); jj <- match(rownames(m), colnames(m)); yes <- !is.na(jj)
    ss[yes] <- m[cbind(which(yes), jj[yes])]
    pd[[k]] <- data.frame(network = network, condition = condition, source = rownames(m),
                          pathway = dimnames(pp)[[3L]][k], outgoing_including_self = rowSums(m),
                          self_interaction = ss, outgoing_excluding_self = rowSums(m) - ss)
  }
  pathway <- if (length(pd)) do.call(rbind, pd) else data.frame(network = character(), condition = character(),
    source = character(), pathway = character(), outgoing_including_self = numeric(), self_interaction = numeric(), outgoing_excluding_self = numeric())
  .cc_csv(pathway, file.path(run_dir, "pathway_outgoing.csv"))
  invisible(list(total = d, pathway = pathway))
}

.cc_run_one <- function(input, network, condition, output_root, db, contract, contract_hash) {
  parent <- file.path(output_root, network, condition)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  complete_path <- file.path(parent, "complete.json")
  if (file.exists(complete_path)) {
    old <- jsonlite::read_json(complete_path, simplifyVector = TRUE)
    if (identical(old$status, "complete") && identical(old$contract_hash, contract_hash) &&
        file.exists(file.path(old$attempt_dir, "cellchat_object.rds"))) {
      if (is.null(old$artifact_sha256) || !length(old$artifact_sha256)) .cc_stop("Complete run lacks artifact hash audit")
      for (nm in names(old$artifact_sha256)) {
        path <- file.path(old$attempt_dir, nm)
        if (!file.exists(path) || !identical(.cc_hash(path, file = TRUE), old$artifact_sha256[[nm]]))
          .cc_stop("Complete artifact missing/hash mismatch: ", path)
      }
      cat("SKIP matching complete run: ", network, " / ", condition, "\n", sep = "")
      return(old)
    }
    .cc_stop("Existing complete run has a different contract or missing object: ", parent)
  }
  run_dir <- file.path(parent, paste0("attempt_", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "_", Sys.getpid()))
  if (dir.exists(run_dir)) .cc_stop("Attempt directory collision")
  dir.create(run_dir)
  started <- .cc_utc()
  tryCatch({
    keep <- if (condition == "Overall") rep(TRUE, nrow(input$meta)) else input$meta$condition == condition
    md <- input$meta[keep, , drop = FALSE]
    mat <- input$mat[, keep, drop = FALSE]
    if (!ncol(mat)) .cc_stop("No IEG-high cells for run")
    md$celltype <- if (network == "CrossType") paste0(md$neuron_type, "_", md$subtype) else as.character(md$subtype)
    md$celltype <- factor(md$celltype, levels = sort(unique(md$celltype)))
    rownames(md) <- md$Cell
    group_counts <- .cc_group_counts(md)
    .cc_csv(group_counts, file.path(run_dir, "actual_group_counts.csv"))
    exc_labels <- c("GLU_L2_L3", "GLU_L4_5", "GLU_L5", "GLU_L5_6", "GLU_L6_CLA", "GLU_ADGM")
    inh_labels <- c("GABA_Amb", "GABA_Etv1", "GABA_Lamp5", "GABA_MSN", "GABA_Sst", "GABA_Hs3st2", "GABA_Pvalb", "GABA_Vip")
    declared <- switch(network, CrossType = c(paste0("Excitatory_", exc_labels), paste0("Inhibitory_", inh_labels)), Excitatory = exc_labels, Inhibitory = inh_labels)
    eligibility <- data.frame(group = declared, n_cells = as.integer(table(factor(as.character(md$celltype), levels = declared))), stringsAsFactors = FALSE)
    eligibility$status <- ifelse(eligibility$n_cells == 0L, "absent_not_modeled", ifelse(eligibility$n_cells <= 10L, "rare_native_min_cells_mask", "eligible"))
    .cc_csv(eligibility, file.path(run_dir, "declared_group_eligibility.csv"))
    .cc_csv(as.data.frame(table(md$pool_id, md$condition, md$celltype)), file.path(run_dir, "pool_condition_group_counts.csv"))
    .cc_csv(md, file.path(run_dir, "input_cell_metadata.csv"))
    if (sum(group_counts$eligible_min_cells_10) < 2L) .cc_stop("Fewer than two actual min.cells-eligible groups")
    set.seed(42L)
    object <- CellChat::createCellChat(object = mat, meta = md, group.by = "celltype", datatype = "RNA", do.sparse = TRUE)
    object@DB <- db
    object <- CellChat::subsetData(object, features = NULL)
    native_genes <- rownames(object@data.signaling)
    signal_coverage <- data.frame(human_symbol = contract$db_expanded_genes,
                                  present_in_input = contract$db_expanded_genes %in% rownames(mat),
                                  retained_native_subsetData = contract$db_expanded_genes %in% native_genes)
    detected <- Matrix::rowSums(object@data.signaling != 0)
    signal_coverage$n_cells_detected <- unname(detected[match(signal_coverage$human_symbol, names(detected))])
    signal_coverage$n_cells_detected[is.na(signal_coverage$n_cells_detected)] <- 0L
    .cc_csv(signal_coverage, file.path(run_dir, "native_signaling_gene_coverage.csv"))
    if (any(signal_coverage$present_in_input != signal_coverage$retained_native_subsetData)) .cc_stop("Native DB expansion mismatch")
    object <- CellChat::identifyOverExpressedGenes(object, data.use = NULL, group.by = NULL, idents.use = NULL,
      invert = FALSE, group.dataset = NULL, pos.dataset = NULL, features.name = "features", only.pos = TRUE,
      features = NULL, return.object = TRUE, thresh.pc = 0, thresh.fc = 0, thresh.p = 0.05)
    object <- CellChat::identifyOverExpressedInteractions(object, features.name = "features", features = NULL, return.object = TRUE)
    .cc_csv(object@LR$LRsig, file.path(run_dir, "LR_interactions_used.csv"))
    if (!nrow(object@LR$LRsig)) .cc_stop("No overexpressed LR interactions; cannot infer a network")
    object <- CellChat::computeCommunProb(object, type = "triMean", trim = 0.1, LR.use = NULL,
      raw.use = TRUE, population.size = TRUE, distance.use = TRUE, interaction.length = 200,
      scale.distance = 0.01, k.min = 10, nboot = 100, seed.use = 42L, Kh = 0.5, n = 1)
    .cc_array_csv(object@net$prob, object@net$pval, file.path(run_dir, "all_LR_probabilities_before_min_cells.csv.gz"), "interaction_name")
    object <- CellChat::filterCommunication(object, min.cells = 10L)
    .cc_array_csv(object@net$prob, object@net$pval, file.path(run_dir, "all_LR_probabilities.csv.gz"), "interaction_name")
    object <- CellChat::computeCommunProbPathway(object, net = NULL, pairLR.use = NULL, thresh = 0.05)
    object <- CellChat::aggregateNet(object, sources.use = NULL, targets.use = NULL, signaling = NULL,
      pairLR.use = NULL, remove.isolate = TRUE, thresh = 0.05, return.object = TRUE)
    object <- CellChat::netAnalysis_computeCentrality(object, slot.name = "netP", net = NULL, net.name = NULL, thresh = 0.05)
    eligibility$present_in_aggregate_source <- eligibility$group %in% rownames(object@net$weight)
    eligibility$present_in_aggregate_target <- eligibility$group %in% colnames(object@net$weight)
    .cc_csv(eligibility, file.path(run_dir, "native_target_universe_audit.csv"))
    .cc_array_csv(object@netP$prob, NULL, file.path(run_dir, "all_pathway_probabilities.csv.gz"), "pathway")
    .cc_csv(as.data.frame(object@net$count), file.path(run_dir, "net_count.csv"), row.names = TRUE)
    .cc_csv(as.data.frame(object@net$weight), file.path(run_dir, "net_weight.csv"), row.names = TRUE)
    .cc_csv(CellChat::subsetCommunication(object, slot.name = "net", thresh = 0.05), file.path(run_dir, "significant_LR_pairs_native_P0p05.csv"))
    .cc_export_outgoing(object, run_dir, network, condition)
    .cc_atomic_rds(object, file.path(run_dir, "cellchat_object.rds"))
    artifacts <- list.files(run_dir, pattern = "\\.(csv|csv\\.gz|rds)$", full.names = TRUE)
    artifact_hashes <- setNames(lapply(artifacts, .cc_hash, file = TRUE), basename(artifacts))
    result <- list(status = "complete", network = network, condition = condition, attempt_dir = run_dir,
                   started_UTC = started, finished_UTC = .cc_utc(), contract_hash = contract_hash,
                   n_input_cells = ncol(mat), n_actual_groups = nrow(group_counts),
                   n_min_cells_eligible_groups = sum(group_counts$eligible_min_cells_10),
                   rare_groups = as.character(group_counts$group[!group_counts$eligible_min_cells_10]),
                   n_input_symbols = nrow(mat), n_native_signaling_genes = length(native_genes),
                   n_LR_used = nrow(object@LR$LRsig), n_pathways = length(object@netP$pathways),
                   object_sha256 = .cc_hash(file.path(run_dir, "cellchat_object.rds"), file = TRUE),
                   artifact_sha256 = artifact_hashes,
                   eligible_declared_targets = as.character(eligibility$group[eligibility$status == "eligible"]),
                   absent_or_rare_declared_groups = eligibility[eligibility$status != "eligible", , drop = FALSE],
                   actual_aggregate_target_groups = colnames(object@net$weight),
                   target_policy = "Native condition-specific groups and min.cells10 mask; outgoing sums all native included targets excluding self; no hidden common-target intersection or replacement model",
                   set_seed = 42L, computeCommunProb_seed_use = 42L,
                   inferential_scope = "descriptive predicted communication; native ligand-receptor bootstrap P is not biological-pool condition inference")
    .cc_atomic_json(result, file.path(run_dir, "run_complete.json"))
    .cc_atomic_json(result, complete_path)
    rm(object); gc()
    result
  }, error = function(e) {
    .cc_atomic_json(list(status = "failed", network = network, condition = condition, attempt_dir = run_dir,
                        started_UTC = started, failed_UTC = .cc_utc(), contract_hash = contract_hash,
                        error = conditionMessage(e)), file.path(run_dir, "failure.json"))
    stop(e)
  })
}

.cc_comparisons <- function(results, dest) {
  for (network in c("CrossType", "Excitatory", "Inhibitory")) {
    keys <- paste(network, c("Ext", "Ret-Ext"), sep = "__")
    if (!all(keys %in% names(results))) .cc_stop("Missing contrast runs")
    for (kind in c("total", "pathway")) {
      files <- vapply(results[keys], function(x) file.path(x$attempt_dir, paste0(kind, "_outgoing.csv")), character(1))
      a <- read.csv(files[1L], stringsAsFactors = FALSE, check.names = FALSE)
      b <- read.csv(files[2L], stringsAsFactors = FALSE, check.names = FALSE)
      join <- c("network", "source", if (kind == "pathway") "pathway")
      a$condition <- NULL; b$condition <- NULL
      d <- merge(a, b, by = join, all = TRUE, suffixes = c("_Ext", "_RetExt"))
      count_a <- read.csv(file.path(results[[keys[1L]]]$attempt_dir, "actual_group_counts.csv"), stringsAsFactors = FALSE)
      count_b <- read.csv(file.path(results[[keys[2L]]]$attempt_dir, "actual_group_counts.csv"), stringsAsFactors = FALSE)
      d$estimable_source_Ext <- d$source %in% count_a$group[count_a$eligible_min_cells_10]
      d$estimable_source_RetExt <- d$source %in% count_b$group[count_b$eligible_min_cells_10]
      # Only a pathway absent for an actual eligible source is a modeled zero.
      # Missing or rare (<10 cells) sources remain unavailable, not zero.
      for (suffix in c("Ext", "RetExt")) {
        estimable <- d[[paste0("estimable_source_", suffix)]]
        for (metric in c("outgoing_including_self", "self_interaction", "outgoing_excluding_self")) {
          nm <- paste0(metric, "_", suffix)
          if (kind == "pathway") d[[nm]][is.na(d[[nm]]) & estimable] <- 0
          d[[nm]][!estimable] <- NA_real_
        }
      }
      d$delta_RetExt_minus_Ext_without_self <- d$outgoing_excluding_self_RetExt - d$outgoing_excluding_self_Ext
      d$delta_Ext_minus_RetExt_without_self <- -d$delta_RetExt_minus_Ext_without_self
      d$percent_change_without_self <- ifelse(d$outgoing_excluding_self_Ext > 0,
        100 * d$delta_RetExt_minus_Ext_without_self / d$outgoing_excluding_self_Ext, NA_real_)
      d$percent_status <- ifelse(!d$estimable_source_Ext | !d$estimable_source_RetExt |
        is.na(d$outgoing_excluding_self_Ext) | is.na(d$outgoing_excluding_self_RetExt),
        "unavailable_or_rare_group", ifelse(d$outgoing_excluding_self_Ext == 0, "undefined_zero_Ext_baseline", "defined"))
      d$delta_RetExt_minus_Ext_with_self <- d$outgoing_including_self_RetExt - d$outgoing_including_self_Ext
      d$percent_change_with_self <- ifelse(d$outgoing_including_self_Ext > 0,
        100 * d$delta_RetExt_minus_Ext_with_self / d$outgoing_including_self_Ext, NA_real_)
      .cc_csv(d, file.path(dest, paste0(network, "_", kind, "_RetExt_vs_Ext.csv")))
    }
  }
}

run_clean_cellchat <- function(root) {
  for (pkg in c("CellChat", "Matrix", "digest", "jsonlite", "future")) if (!requireNamespace(pkg, quietly = TRUE)) .cc_stop("Missing package: ", pkg)
  if (as.character(utils::packageVersion("CellChat")) != "1.6.1") .cc_stop("CellChat1.6.1 required; no silent version switch")
  .cc_runtime()
  root <- normalizePath(root, mustWork = TRUE)
  output <- file.path(root, "06_CellChat")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  orthology_dir <- file.path(output, "orthology")
  orthology_complete <- jsonlite::read_json(file.path(orthology_dir, "complete.json"), simplifyVector = TRUE)
  if (!identical(orthology_complete$status, "complete") || as.integer(orthology_complete$release) != 116L) .cc_stop("Frozen Ensembl116 orthology is required")
  orthology_path <- file.path(orthology_dir, "orthology_lookup.rds")
  orthology <- readRDS(orthology_path)
  api_names <- c("createCellChat", "subsetData", "identifyOverExpressedGenes", "identifyOverExpressedInteractions",
    "computeCommunProb", "filterCommunication", "computeCommunProbPathway", "aggregateNet", "netAnalysis_computeCentrality", "subsetCommunication")
  db_env <- new.env(parent = emptyenv())
  utils::data(list = "CellChatDB.human", package = "CellChat", envir = db_env)
  db <- db_env$CellChatDB.human
  if (is.null(db) || nrow(db$interaction) != 1939L) .cc_stop("Unexpected human DB")
  expanded <- get("extractGene", envir = asNamespace("CellChat"))(db)
  db_path <- file.path(output, paste0("CellChatDB_human_", substr(.cc_hash(db), 1L, 12L), ".rds"))
  if (file.exists(db_path)) {
    if (!identical(.cc_hash(readRDS(db_path)), .cc_hash(db))) .cc_stop("Existing frozen DB hash mismatch")
  } else .cc_atomic_rds(db, db_path)
  input_paths <- setNames(file.path(root, "01_inputs", c("Excitatory", "Inhibitory"), "IEGhigh_payload.rds"), c("Excitatory", "Inhibitory"))
  if (!all(file.exists(input_paths))) .cc_stop("Both frozen IEG-high payloads are required")
  input_hashes <- lapply(input_paths, .cc_hash, file = TRUE)
  contract <- list(schema_version = "1.0", R_version = as.character(getRversion()), CellChat_version = "1.6.1",
    script_sha256 = .cc_hash(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]), file = TRUE),
    DB = "CellChatDB.human", DB_sha256_serialized = .cc_hash(db), db_expanded_genes = expanded,
    input_paths = as.list(input_paths), input_sha256 = input_hashes, api_defaults = .cc_defaults(api_names),
    actual_parameters = list(OE_only_pos = TRUE, OE_thresh_pc = 0, OE_thresh_fc = 0, OE_thresh_p = 0.05,
      probability_type = "triMean", trim = 0.1, raw_use = TRUE, population_size = TRUE, distance_use = TRUE,
      interaction_length = 200, scale_distance = 0.01, k_min = 10L, nboot = 100L, seed_use = 42L, Kh = 0.5, hill_n = 1,
      min_cells = 10L, pathway_threshold = 0.05, aggregate_threshold = 0.05, remove_isolate = TRUE, centrality_threshold = 0.05),
    seed_note = "global seed42 and computeCommunProb seed.use42; the latter explicitly differs from API default seed.use1",
    symbol_policy = "Verified Ensembl116 rat-human one-to-one orthology, ratID preferred and unique exact original-symbol fallback audited; official human symbols, no case conversion",
    orthology = list(release = 116L, lookup_path = orthology_path, lookup_sha256 = .cc_hash(orthology_path, file = TRUE), public_source = orthology_complete$actual_mart_host,
      raw_table_sha256 = .cc_hash(file.path(orthology_dir, "rat_human_Ensembl116_raw.rds"), file = TRUE), n_one2one = nrow(orthology$accepted)),
    duplicate_policy = "unique mapped input rows unchanged; same-verified-gene duplicate corrected counts summed then log1p(1e4*sum/FULL ORIGINAL RNA corrected library size)",
    native_complex_cofactor = "native subsetData/extractGene retained; no extra direct ligand/receptor name filter",
    biological_scope = "frozen user-approved working display groups; uncertain review flags preserved, not newly validated anatomical types",
    method_differences = c("new clean cells/annotations", "IEG-high scope rather than old full-human script lacking high filter",
      "frozen Ensembl116 true unique one-to-one orthology", "native complex/cofactor expansion retained", "explicit seed.use42 instead of API default1",
      "12 Overall/condition x lineage/joint networks", "self-excluded outgoing totals plus including-self audit"),
    resources = list(future = "sequential", BLAS_threads = 1L, OMP_max = 4L, GPU = FALSE))
  contract_hash <- .cc_hash(contract)
  .cc_atomic_json(contract, file.path(output, "analysis_contract.json"))
  .cc_csv(db$interaction, file.path(output, "CellChatDB_human_interactions.csv"), row.names = TRUE)
  .cc_csv(db$complex, file.path(output, "CellChatDB_human_complexes.csv"), row.names = TRUE)
  .cc_csv(db$cofactor, file.path(output, "CellChatDB_human_cofactors.csv"), row.names = TRUE)
  prep <- lapply(names(input_paths), function(lineage) {
    x <- readRDS(input_paths[[lineage]])
    result <- .cc_prepare_payload(x, lineage, file.path(output, "input_audit", lineage), orthology = orthology)
    rm(x); gc(); result
  })
  names(prep) <- names(input_paths)
  if (!identical(rownames(prep$Excitatory$mat), rownames(prep$Inhibitory$mat))) .cc_stop("Lineage symbol universes differ; no silent intersect")
  joint_meta <- rbind(prep$Excitatory$meta, prep$Inhibitory$meta)
  joint_mat <- cbind(prep$Excitatory$mat, prep$Inhibitory$mat)
  if (anyDuplicated(colnames(joint_mat)) || !identical(colnames(joint_mat), joint_meta$Cell)) .cc_stop("CrossType cell join failed")
  inputs <- list(CrossType = list(mat = joint_mat, meta = joint_meta), Excitatory = prep$Excitatory, Inhibitory = prep$Inhibitory)
  .cc_csv(data.frame(human_symbol = rownames(joint_mat), matched_expanded_DB = rownames(joint_mat) %in% expanded),
          file.path(output, "input_audit", "human_DB_symbol_coverage.csv"))
  started <- .cc_utc(); results <- list()
  for (network in names(inputs)) for (condition in c("Overall", "WD", "Ext", "Ret-Ext")) {
    key <- paste(network, condition, sep = "__")
    cat("\nBEGIN ", key, " ", .cc_utc(), "\n", sep = "")
    results[[key]] <- .cc_run_one(inputs[[network]], network, condition, output, db, contract, contract_hash)
    .cc_atomic_json(list(status = "running", started_UTC = started, updated_UTC = .cc_utc(), contract_hash = contract_hash,
                        n_complete = length(results), runs = results), file.path(output, "run_manifest.json"))
    gc()
  }
  .cc_comparisons(results, output)
  capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
  manifest <- list(status = "complete", started_UTC = started, finished_UTC = .cc_utc(), contract_hash = contract_hash,
    n_complete = length(results), runs = results, plots_created = FALSE,
    inference = "descriptive predicted communication; no biological-pool differential communication hypothesis test",
    input_preparation = lapply(prep, `[[`, "info"))
  .cc_atomic_json(manifest, file.path(output, "run_manifest.json"))
  invisible(manifest)
}

.cc_main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  at <- match("--root", args)
  if (is.na(at) || at == length(args)) .cc_stop("Usage: fit_cellchat.R --root <new run root>")
  root <- args[at + 1L]
  tryCatch(run_clean_cellchat(root), error = function(e) {
    if (requireNamespace("jsonlite", quietly = TRUE)) .cc_atomic_json(
      list(status = "failed", failed_UTC = .cc_utc(), error = conditionMessage(e)),
      file.path(root, "06_CellChat", "execution_failure.json"))
    stop(e)
  })
}
if (sys.nframe() == 0L) .cc_main()
