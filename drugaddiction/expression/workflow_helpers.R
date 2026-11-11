# Shared input checks and native DESeq2 calculations for the expression workflows.
require_packages <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install the required R packages: ", paste(absent, collapse = ", "))
}

limit_threads <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
  }
}

new_output_dir <- function(path) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE))) {
    stop("Choose an empty output directory: ", path)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

read_table <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           fileEncoding = "UTF-8-BOM")
}

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

write_session <- function(out) {
  writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
}

pool_ids <- c(paste0("WD", 1:3), paste0("Ext", 1:3), paste0("Ret-Ext", 1:3))
pool_conditions <- setNames(rep(c("WD", "Ext", "Ret-Ext"), each = 3L), pool_ids)

read_ieg_payload <- function(path, neuronal_class) {
  require_packages("Matrix")
  q <- readRDS(path)
  stopifnot(all(c("counts", "metadata", "gene_annotation") %in% names(q)))
  x <- q$counts
  md <- q$metadata
  ann <- q$gene_annotation
  stopifnot(inherits(x, "sparseMatrix"),
            all(c("Cell", "pool_id", "condition", "IEGhigh") %in% names(md)),
            all(c("seurat_feature", "gene_id", "original_symbol") %in% names(ann)),
            identical(colnames(x), as.character(md$Cell)),
            identical(rownames(x), as.character(ann$seurat_feature)),
            !anyDuplicated(md$Cell), !anyDuplicated(rownames(x)),
            all(md$IEGhigh), setequal(md$pool_id, pool_ids),
            identical(as.character(md$condition), unname(pool_conditions[as.character(md$pool_id)])),
            all(is.finite(x@x)), all(x@x >= 0), all(x@x == round(x@x)))
  if ("neuron_type" %in% names(md)) stopifnot(all(md$neuron_type == neuronal_class))
  q
}

aggregate_pools <- function(q) {
  group <- Matrix::sparseMatrix(i = seq_len(ncol(q$counts)),
                               j = match(q$metadata$pool_id, pool_ids),
                               x = 1, dims = c(ncol(q$counts), length(pool_ids)))
  counts <- as.matrix(q$counts %*% group)
  stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts <= .Machine$integer.max))
  dimnames(counts) <- list(rownames(q$counts), pool_ids)
  storage.mode(counts) <- "integer"
  counts
}

make_pool_metadata <- function(pools, levels) {
  data.frame(pool_id = pools, condition = factor(unname(pool_conditions[pools]), levels),
             row.names = pools)
}

expression_filter <- function(counts, metadata) {
  unname(edgeR::filterByExpr(edgeR::DGEList(counts, group = metadata$condition),
                            design = model.matrix(~ condition, metadata)))
}

fit_native_deseq <- function(counts, metadata, eligible, numerator, denominator,
                             allow_convergence_retry = FALSE,
                             contrast = c("condition", numerator, denominator)) {
  stopifnot(identical(colnames(counts), rownames(metadata)), any(eligible),
            length(eligible) == nrow(counts), !anyNA(eligible))
  dds <- DESeq2::DESeqDataSetFromMatrix(counts[eligible, , drop = FALSE],
                                       metadata, design = ~ condition)
  dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric", sfType = "ratio",
                       betaPrior = FALSE, minReplicatesForReplace = 7L,
                       quiet = TRUE, parallel = FALSE)
  initial <- S4Vectors::mcols(dds)$betaConv
  retried <- any(is.na(initial) | !initial)
  if (retried && allow_convergence_retry) {
    dds <- DESeq2::nbinomWaldTest(dds, betaPrior = FALSE, maxit = 1000, quiet = TRUE)
  }
  stopifnot(all(S4Vectors::mcols(dds)$betaConv),
            all(DESeq2::counts(dds) == counts[eligible, , drop = FALSE]),
            !("replaceCounts" %in% SummarizedExperiment::assayNames(dds)))
  res <- DESeq2::results(dds, contrast = contrast,
                         alpha = 0.05, independentFiltering = TRUE,
                         pAdjustMethod = "BH", parallel = FALSE)
  threshold <- as.numeric(S4Vectors::metadata(res)$filterThreshold)
  family <- is.finite(res$pvalue) & res$baseMean >= threshold
  stopifnot(length(threshold) == 1L, is.finite(threshold),
            identical(unname(is.finite(res$padj)), unname(family)))
  if (any(family)) {
    stopifnot(all(abs(res$padj[family] - p.adjust(res$pvalue[family], "BH")) < 1e-12))
  }
  list(model = dds, result = res, filter_threshold = threshold,
       convergence_retry = retried && allow_convergence_retry)
}

native_de_table <- function(fit, counts, annotation, eligible, neuronal_class,
                            numerator, denominator) {
  r <- as.data.frame(fit$result)
  idx <- match(rownames(counts), rownames(r))
  mc <- S4Vectors::mcols(fit$model)
  cooks <- SummarizedExperiment::assays(fit$model)[["cooks"]]
  max_cooks <- apply(cooks, 1, function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_)
  z <- data.frame(
    gene = rownames(counts), gene_id = annotation$gene_id,
    original_symbol = annotation$original_symbol, eligible = eligible,
    tested = eligible & is.finite(r$pvalue[idx]),
    in_BH_family = eligible & is.finite(r$padj[idx]),
    baseMean = r$baseMean[idx], log2FoldChange = r$log2FoldChange[idx],
    lfcSE = r$lfcSE[idx], stat = r$stat[idx], pvalue = r$pvalue[idx],
    padj = r$padj[idx], maxCooks = max_cooks[idx], betaConv = mc$betaConv[idx],
    neuron_type = neuronal_class, numerator = numerator, denominator = denominator,
    method = "DESeq2", n_pools = ncol(counts), n_per_condition = 3L,
    independentFiltering = TRUE, alpha = 0.05, stringsAsFactors = FALSE)
  z$status <- ifelse(!eligible, "filtered_low_expression_before_fit",
                     ifelse(!z$tested, "native_P_withheld",
                            ifelse(!z$in_BH_family, "independent_filter_low_mean", "tested_with_native_BH")))
  z$BH_significant <- z$in_BH_family & z$padj < 0.05
  z$display_DEG <- z$BH_significant & abs(z$log2FoldChange) >= log2(1.15)
  z$DE_direction <- ifelse(z$display_DEG, ifelse(z$log2FoldChange > 0, "Up", "Down"),
                           ifelse(z$in_BH_family, "Not significant", "Not in BH family"))
  z
}
