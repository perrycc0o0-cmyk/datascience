# Summarize the Ext and Ret-Ext conditions from the same joint neuronal model.
# Usage: Rscript summarize_joint_cellchat.R RUN_ROOT
suppressPackageStartupMessages(library(CellChat))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript summarize_joint_cellchat.R RUN_ROOT")
root <- normalizePath(args[1], mustWork = TRUE)
model_dir <- file.path(root, "06_CellChat")
out <- file.path(model_dir, "manuscript_tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(out, pattern = "\\.csv$"))) {
  stop("Joint-model manuscript tables already exist. Archive the existing tables before rebuilding.")
}
write_table <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE,
                                         fileEncoding = "UTF-8")
manifest <- jsonlite::read_json(file.path(model_dir, "run_manifest.json"))
if (!identical(manifest$status, "complete")) stop("Complete the CellChat fit before summarizing the joint model.")

conditions <- c("Ext", "Ret-Ext")
models <- setNames(lapply(conditions, function(condition) {
  entry <- manifest$runs[[paste0("CrossType__", condition)]]
  if (!identical(entry$status, "complete")) stop("Missing completed joint model: ", condition)
  path <- file.path(entry$attempt_dir, "cellchat_object.rds")
  if (!file.exists(path)) {
    path <- file.path(model_dir, "CrossType", condition, basename(entry$attempt_dir), "cellchat_object.rds")
  }
  if (!identical(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE), entry$object_sha256)) {
    stop("CellChat model checksum does not match its manifest: ", path)
  }
  readRDS(path)
}), conditions)
subtypes <- c("GLU_L2_L3", "GLU_L4_5", "GLU_L5", "GLU_L5_6", "GLU_L6_CLA", "GLU_ADGM")
pretty <- setNames(c("GLU-L2/3", "GLU-L4/5", "GLU-L5", "GLU-L5/6", "GLU-L6 CLA", "GLU-ADGM"), subtypes)
senders <- paste0("Excitatory_", subtypes)
scopes <- c("Excitatory", "All")
pathways <- c("NRXN", "PTN", "CADM", "NRG", "CDH", "NCAM", "LAMININ",
              "CNTN", "EPHA", "NEGR", "SEMA5", "NGL", "EPHB")
p_source <- c("GLU_L6_CLA", "GLU_L6_CLA", "GLU_L5", "GLU_L5_6", "GLU_L5_6")
p_scope <- c("All", "Excitatory", "All", "Excitatory", "All")
p_receiver <- c("All neuronal targets", "Exc targets", "All other subtypes",
                "Exc targets", "All neuronal targets")
scope_label <- function(scope) ifelse(scope == "All", "All", "Exc")
mask_limit <- 10L

summarize_one <- function(x, condition) {
  probability <- x@net$prob
  groups <- dimnames(probability)[[1]]
  stopifnot(identical(groups, dimnames(probability)[[2]]))
  counts <- setNames(as.integer(table(factor(as.character(x@idents), levels = groups))), groups)
  keep <- groups[counts > mask_limit]
  total <- probability
  total[x@net$pval >= 0.05] <- 0
  weight <- apply(total, c(1L, 2L), sum)
  w <- x@net$weight
  stopifnot(all(abs(weight[rownames(w), colnames(w), drop = FALSE] - w) < 1e-12))
  probability[x@net$pval > 0.05] <- 0
  lr <- x@LR$LRsig
  tested <- dimnames(probability)[[3]]
  idx <- match(tested, lr$interaction_name)
  if (anyNA(idx)) idx <- match(tested, rownames(lr))
  stopifnot(!anyNA(idx))
  lr_pathway <- as.character(lr$pathway_name[idx])
  pathway_weights <- setNames(lapply(pathways, function(p) {
    selected <- which(lr_pathway == p)
    if (!length(selected)) return(NULL)
    apply(probability[, , selected, drop = FALSE], c(1L, 2L), sum)
  }), pathways)

  grid <- expand.grid(source = subtypes, target_scope = scopes, stringsAsFactors = FALSE)
  out <- lapply(seq_len(nrow(grid)), function(i) {
    source <- grid$source[i]
    scope <- grid$target_scope[i]
    sender <- paste0("Excitatory_", source)
    n <- if (sender %in% groups) counts[[sender]] else 0L
    candidate <- if (scope == "All") keep else keep[startsWith(keep, "Excitatory_")]
    receivers <- setdiff(candidate, sender)
    state <- if (n == 0L) "absent_not_modeled" else if (n <= mask_limit) {
      "native_mask_n_le_10"
    } else if (!length(candidate)) "no_eligible_targets" else "estimable"
    value <- function(mat) {
      if (state != "estimable" || is.null(mat)) return(NA_real_)
      sum(mat[sender, receivers, drop = FALSE])
    }
    totals <- data.frame(source, target_scope = scope, condition, source_n = n,
                         n_receivers = length(receivers), outgoing = value(weight), status = state)
    vals <- vapply(pathway_weights, value, numeric(1))
    status <- ifelse(vapply(pathway_weights, is.null, logical(1)), "pathway_not_tested",
                     ifelse(is.finite(vals) & vals == 0, "tested_filtered_zero", state))
    paths <- data.frame(source, target_scope = scope, condition, source_n = n,
                        pathway = pathways, outgoing = unname(vals), status = unname(status))
    list(totals = totals, paths = paths)
  })
  list(totals = do.call(rbind, lapply(out, `[[`, "totals")),
       paths = do.call(rbind, lapply(out, `[[`, "paths")), weight = weight, counts = counts)
}

fits <- Map(summarize_one, models, conditions)
compare <- function(which, keys) {
  a <- fits$Ext[[which]]
  b <- fits[["Ret-Ext"]][[which]]
  a$condition <- b$condition <- NULL
  d <- merge(a, b, by = keys, suffixes = c("_Ext", "_RetExt"), all = TRUE)
  d$delta_RetExt_minus_Ext <- d$outgoing_RetExt - d$outgoing_Ext
  d$delta_Ext_minus_RetExt <- -d$delta_RetExt_minus_Ext
  d$percent_change <- ifelse(is.finite(d$outgoing_Ext) & d$outgoing_Ext > 0,
                             100 * d$delta_RetExt_minus_Ext / d$outgoing_Ext, NA_real_)
  d
}
totals <- compare("totals", c("source", "target_scope"))
path_totals <- compare("paths", c("source", "target_scope", "pathway"))
write_table(totals, "joint_outgoing_all_excitatory_senders")
write_table(path_totals, "joint_outgoing_candidate_pathways")

main <- totals[match(paste(p_source, p_scope), with(totals, paste(source, target_scope))), ]
stopifnot(nrow(main) == 5L, !anyNA(main$source))
main$display_order <- seq_len(nrow(main))
main$display_label <- paste(unname(pretty[main$source]), "\u2192", p_receiver)
write_table(main, "outgoing_selected_comparisons")

# ADGM remains in the audit table; it is not a displayed zero.
display <- totals[totals$source != "GLU_ADGM", ]
display$display_order <- match(display$source, subtypes)
display$display_label <- unname(pretty[display$source])
write_table(display, "outgoing_by_receiver_scope")

nodes <- setdiff(subtypes, "GLU_ADGM")
edges <- subset(expand.grid(source = nodes, target = nodes, stringsAsFactors = FALSE), source != target)
weights <- function(model) {
  src <- paste0("Excitatory_", edges$source)
  dst <- paste0("Excitatory_", edges$target)
  s <- match(src, rownames(model$weight))
  t <- match(dst, colnames(model$weight))
  valid <- !is.na(s) & !is.na(t) & model$counts[src] > mask_limit & model$counts[dst] > mask_limit
  valid[is.na(valid)] <- FALSE
  result <- rep(NA_real_, nrow(edges))
  result[valid] <- model$weight[cbind(s[valid], t[valid])]
  result
}
edges$weight_Ext <- weights(fits$Ext)
edges$weight_RetExt <- weights(fits[["Ret-Ext"]])
edges$delta_RetExt_minus_Ext <- edges$weight_RetExt - edges$weight_Ext
edges$selected_for_display <- edges$source %in% c("GLU_L6_CLA", "GLU_L5_6", "GLU_L5") &
                              is.finite(edges$delta_RetExt_minus_Ext) & edges$delta_RetExt_minus_Ext != 0
edges$source_label <- unname(pretty[edges$source])
edges$target_label <- unname(pretty[edges$target])
write_table(edges, "joint_differential_edges")

# C and D rank directional effect sizes across both scopes, without significance testing.
directional <- function(data, senders, sign, n) {
  value <- sign * data$delta_RetExt_minus_Ext
  eligible <- is.finite(value) & value > 0 & data$source %in% senders
  d <- data[eligible, ]
  d$plot_value <- value[eligible]
  d <- head(d[order(-d$plot_value, d$source, d$target_scope, d$pathway), ], n)
  d$display_order <- seq_len(nrow(d))
  d$display_label <- paste(d$pathway, scope_label(d$target_scope), sep = " | ")
  if (sign < 0) d$display_label <- paste(sub("^GLU-", "", pretty[d$source]), d$display_label, sep = " | ")
  d
}
write_table(directional(path_totals, "GLU_L6_CLA", 1, 9L), "outgoing_pathways_increased")
write_table(directional(path_totals, c("GLU_L5", "GLU_L5_6"), -1, 12L), "outgoing_pathways_decreased")
heat <- path_totals[paste(path_totals$source, path_totals$target_scope) %in% paste(p_source, p_scope), ]
heat$row_order <- match(paste(heat$source, heat$target_scope), paste(p_source, p_scope))
heat$column_order <- match(heat$pathway, pathways)
heat$row_label <- paste(pretty[heat$source], ifelse(heat$source == "GLU_L5", "All other subtypes",
                       scope_label(heat$target_scope)), sep = " | ")
heat$display_value <- 1000 * heat$delta_RetExt_minus_Ext
heat$color_value <- asinh(heat$display_value / 5)
write_table(heat[order(heat$row_order, heat$column_order), ], "outgoing_pathways_heatmap")
