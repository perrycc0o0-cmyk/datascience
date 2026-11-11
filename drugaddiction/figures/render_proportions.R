# Draw observed pooled-sample percentages using an external, finalized statistics table.
# Usage: Rscript render_proportions.R NEURON_METADATA.csv EXC_STATS.csv INH_STATS.csv OUTPUT_DIR
suppressPackageStartupMessages({library(ggplot2); library(grid)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript render_proportions.R NEURON_METADATA.csv EXC_STATS.csv INH_STATS.csv OUTPUT_DIR")
}
metadata <- read.csv(args[1], fileEncoding = "UTF-8-BOM")
stopifnot(all(c("Cell", "pool_id", "condition", "neuron_type", "subtype", "IEGhigh") %in% names(metadata)),
          !anyDuplicated(metadata$Cell), is.logical(metadata$IEGhigh), !anyNA(metadata$IEGhigh))
out <- args[4]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stems <- c("excitatory_recruitment", "inhibitory_recruitment")
if (any(file.exists(outer(file.path(out, stems),
                          c(".pdf", ".png", ".svg", "_sample_points.csv", "_observed_ranges.csv", "_brackets.csv"), paste0)))) {
  stop("A recruitment output already exists. Choose a new output directory.")
}
subtypes <- list(
  Excitatory = c("GLU_L2_L3", "GLU_L4_5", "GLU_L5", "GLU_L5_6", "GLU_L6_CLA", "GLU_ADGM"),
  Inhibitory = c("GABA_Amb", "GABA_Etv1", "GABA_Lamp5", "GABA_MSN", "GABA_Sst", "GABA_Hs3st2", "GABA_Pvalb", "GABA_Vip")
)
pools <- c(paste0("Ext", 1:3), paste0("Ret-Ext", 1:3))
expected <- c(rep("Ext", 3), rep("Ret-Ext", 3))
observed <- unique(metadata[metadata$pool_id %in% pools, c("pool_id", "condition")])
stopifnot(!anyDuplicated(observed$pool_id), identical(observed$condition[match(pools, observed$pool_id)], expected))

read_stats <- function(path, order) {
  d <- read.csv(path, na.strings = c("", "NA"), fileEncoding = "UTF-8-BOM")
  required <- c("subtype", "P", "BH_adjusted_P", "Label")
  if (!all(required %in% names(d)) || anyDuplicated(d$subtype) || !setequal(d$subtype, order)) {
    stop("Statistics table must contain exactly one row per subtype and columns: ", paste(required, collapse = ", "))
  }
  d <- d[match(order, d$subtype), required]
  if (!is.numeric(d$P) || !is.numeric(d$BH_adjusted_P) ||
      !identical(is.na(d$P), is.na(d$BH_adjusted_P)) ||
      any(!is.finite(d$P[!is.na(d$P)]) | d$P[!is.na(d$P)] < 0 | d$P[!is.na(d$P)] > 1)) {
    stop("P and BH-adjusted P must be valid numeric values or paired NA.")
  }
  q <- p.adjust(d$P, method = "BH", n = length(order))
  if (any(abs(q - d$BH_adjusted_P) > 5.1e-6, na.rm = TRUE)) {
    stop("Supplied adjusted P values do not agree with the full subtype BH family.")
  }
  labels <- ifelse(is.na(d$BH_adjusted_P), NA_character_,
    ifelse(d$BH_adjusted_P < .0001, "****", ifelse(d$BH_adjusted_P < .001, "***",
      ifelse(d$BH_adjusted_P < .01, "**", ifelse(d$BH_adjusted_P < .05, "*", "ns")))))
  if (!identical(d$Label, labels)) stop("Supplied significance labels do not match the adjusted P values.")
  d
}

observed_points <- function(md, lineage) {
  ord <- subtypes[[lineage]]
  d <- subset(md, neuron_type == lineage & pool_id %in% pools)
  if (anyNA(d$subtype) || !all(d$subtype %in% ord)) stop("Unexpected subtype labels: ", lineage)
  result <- do.call(rbind, lapply(ord, function(s) {
    do.call(rbind, lapply(pools, function(p) {
      m <- d[d$pool_id == p, ]
      n <- sum(m$IEGhigh & m$subtype == s)
      den <- if (lineage == "Excitatory") sum(m$IEGhigh) else sum(m$subtype == s)
      data.frame(subtype = s, pool_id = p, condition = expected[match(p, pools)],
                 numerator = n, denominator = den, value = if (den > 0) 100 * n / den else NA_real_)
    }))
  }))
  result <- result[is.finite(result$value), ]
  result$condition <- factor(result$condition, levels = c("Ret-Ext", "Ext"))
  result$x <- match(result$subtype, ord) + ifelse(result$condition == "Ext", -.13, .13)
  set.seed(42)
  result$xj <- result$x + runif(nrow(result), -.035, .035)
  result
}

save_plot <- function(plot, name, width, height) {
  stem <- file.path(out, name)
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width, height, family = "Arial", bg = "white")
  tryCatch(print(plot), finally = dev.off())
  ragg::agg_png(paste0(stem, ".png"), width, height, units = "in", res = 600, background = "white")
  tryCatch(print(plot), finally = dev.off())
  svglite::svglite(paste0(stem, ".svg"), width, height, bg = "white")
  tryCatch(print(plot), finally = dev.off())
}

for (lineage in names(subtypes)) {
  inhibitory <- lineage == "Inhibitory"
  ord <- subtypes[[lineage]]
  z <- observed_points(metadata, lineage)
  st <- read_stats(args[if (inhibitory) 3 else 2], ord)
  sig <- subset(st, is.finite(BH_adjusted_P) & BH_adjusted_P < .05)
  if (!nrow(z) || !all(sig$subtype %in% z$subtype)) stop("A displayed result has no observed pool points.")
  ranges <- do.call(rbind, lapply(split(z, list(z$subtype, z$condition), drop = TRUE), function(t) {
    data.frame(subtype = t$subtype[1], condition = t$condition[1], x = t$x[1], lo = min(t$value), hi = max(t$value))
  }))
  upper <- max(20, ceiling(max(z$value) * 1.2 / 20) * 20)
  size <- if (inhibitory) 12 else 8
  p <- ggplot() +
    geom_linerange(data = ranges, aes(x, ymin = lo, ymax = hi, color = condition),
                  linewidth = if (inhibitory) .65 else .4, show.legend = FALSE) +
    geom_point(data = z, aes(xj, value, color = condition), shape = 21, fill = "white",
               size = if (inhibitory) 2.8 else 2.2, stroke = if (inhibitory) .8 else .55) +
    scale_color_manual(values = c(Ext = "#4D4DFF", "Ret-Ext" = "#FF7A7A"), name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = seq_along(ord), labels = gsub("_", " ", ord, fixed = TRUE),
                       limits = c(.5, length(ord) + .5), expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = seq(0, upper, by = if (upper <= 20) 5 else 20),
                       limits = c(-upper * .08, upper), expand = expansion(mult = 0)) +
    labs(x = NULL, y = if (inhibitory) "IEG-high nuclei within subtype (%)" else
         "Percentage of\nIEG-high excitatory neurons (%)") +
    theme_classic(base_family = "Arial", base_size = size) +
    theme(text = element_text(family = "Arial", face = "bold"), axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_text(size = size + 1),
          axis.line.x = element_line(linewidth = if (inhibitory) .96 else .64), axis.line.y = element_blank(),
          axis.ticks = element_line(linewidth = if (inhibitory) .88 else .55),
          legend.position = "inside", legend.position.inside = if (inhibitory) c(.06, .98) else c(.43, .97),
          legend.justification = if (inhibitory) c(0, 1) else c(.5, 1),
          legend.direction = if (inhibitory) "vertical" else "horizontal",
          legend.key.height = unit(if (inhibitory) 4 else 3, "mm"),
          legend.key.width = unit(3, "mm"), legend.spacing.x = unit(1, "mm"),
          plot.margin = if (inhibitory) margin(21, 13, 25, 38) else margin(21, 13, 21, 28)) +
    guides(color = guide_legend(ncol = if (inhibitory) 1 else NULL,
             override.aes = list(shape = 22, size = if (inhibitory) 3 else 2.5,
                                 fill = "white", stroke = if (inhibitory) .85 else .7))) +
    annotate("segment", x = .5, xend = .5, y = 0, yend = upper, linewidth = if (inhibitory) .96 else .64)
  sig$x <- match(sig$subtype, ord)
  sig$y <- vapply(sig$subtype, function(s) max(z$value[z$subtype == s]), numeric(1)) + upper * .045
  sig$h <- rep(upper * .035, nrow(sig))
  for (i in seq_len(nrow(sig))) {
    x <- sig$x[i]; y <- sig$y[i]; h <- sig$h[i]
    p <- p + annotate("segment", x = c(x - .23, x - .23, x + .23),
      xend = c(x - .23, x + .23, x + .23), y = c(y, y + h, y + h), yend = c(y + h, y + h, y),
      linewidth = if (inhibitory) .6 else .35) +
      annotate("text", x = x, y = y + h, label = sig$Label[i],
               size = (if (inhibitory) 9 else 8) / 2.845, vjust = -.15, family = "Arial", fontface = "bold")
  }
  stem <- if (inhibitory) "inhibitory_recruitment" else "excitatory_recruitment"
  write.csv(z, file.path(out, paste0(stem, "_sample_points.csv")), row.names = FALSE)
  write.csv(ranges, file.path(out, paste0(stem, "_observed_ranges.csv")), row.names = FALSE)
  write.csv(sig, file.path(out, paste0(stem, "_brackets.csv")), row.names = FALSE)
  save_plot(p, stem, if (inhibitory) 6.46 else 4.47, if (inhibitory) 4.04 else 2.84)
}
