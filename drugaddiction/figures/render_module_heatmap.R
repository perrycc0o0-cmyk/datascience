#!/usr/bin/env Rscript
# Fixed-list expression display. The panel is not a differential-expression test.
args <- commandArgs(TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript render_module_heatmap.R EXCITATORY_IEGHIGH_PAYLOAD.rds OUTPUT_DIR")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
config <- file.path(dirname(normalizePath(script, mustWork = TRUE)), "config")
suppressPackageStartupMessages({library(Matrix); library(grid); library(ragg)})
if (!capabilities("cairo")) stop("Cairo support is required for PDF export.")
genes <- read.csv(file.path(config, "fixed_module_genes.csv"), stringsAsFactors = FALSE)
palette <- read.csv(file.path(config, "heatmap_palette.csv"), stringsAsFactors = FALSE)
stopifnot(identical(genes$gene_order, seq_len(44L)), !anyDuplicated(genes$gene),
  identical(as.integer(table(factor(genes$module_id, levels = paste0("M", 1:6)))),
            c(4L, 4L, 11L, 7L, 9L, 9L)),
  nrow(palette) == 61L, all(diff(palette$position) > 0))
p <- readRDS(args[1])
stopifnot(all(c("normalized", "metadata") %in% names(p)),
  identical(colnames(p$normalized), as.character(p$metadata$Cell)),
  all(p$metadata$IEGhigh), all(p$metadata$neuron_type == "Excitatory"),
  setequal(unique(p$metadata$condition), c("WD", "Ext", "Ret-Ext")))
if (!all(genes$gene %in% rownames(p$normalized))) {
  stop("Missing fixed-module genes: ", paste(setdiff(genes$gene, rownames(p$normalized)), collapse = ", "))
}
x <- as.matrix(p$normalized[genes$gene, , drop = FALSE])
stopifnot(all(is.finite(x)), all(x >= 0))
conditions <- c("Ret-Ext", "Ext", "WD")
means <- vapply(conditions, function(k)
  rowMeans(x[, p$metadata$condition == k, drop = FALSE]), numeric(nrow(genes)))
rownames(means) <- genes$gene
gene_sd <- apply(means, 1, sd)
z <- sweep(sweep(means, 1, rowMeans(means), "-"), 1, gene_sd, "/")
z[!is.finite(gene_sd) | gene_sd == 0, ] <- NA_real_
mat <- z[, c("Ret-Ext", "Ext"), drop = FALSE]

module_names <- c("Secretory / dense-core vesicle pathway",
  "Activity-dependent transcriptional response", "Synaptic adhesion / circuit remodeling",
  "Endocannabinoid-linked neuromodulation", "Drug-evoked Glu/GABA synaptic plasticity",
  "Addiction-associated cAMP signaling")
module_colors <- c(M1 = "#F96495", M2 = "#FF8000", M3 = "#C0C000",
                   M4 = "#55A0FB", M5 = "#05BE78", M6 = "#F94040")
w <- 15
h <- w * 266/687
lims <- c(-1.5, 1.5)
ramp <- colorRamp(palette$hex, space = "Lab")
map_color <- function(v) {
  if (!is.finite(v)) return("#F2F2F2")
  rgb(ramp((min(lims[2], max(lims[1], v)) - lims[1])/diff(lims)), maxColorValue = 255)
}
draw <- function() {
  grid.newpage()
  pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontface = "plain", col = "black")))
  on.exit(popViewport(), add = TRUE)
  x0 <- 67/687; x1 <- 665/687; cw <- (x1-x0)/44
  yc <- c("Ret-Ext" = .786, Ext = .696); ch <- .09
  for (i in seq_len(44L)) for (k in names(yc)) {
    grid.rect(x0+(i-.5)*cw, yc[k], cw, ch, gp = gpar(fill = map_color(mat[i,k]), col = NA))
  }
  grid.rect((x0+x1)/2, mean(yc), x1-x0, 2*ch, gp = gpar(fill = NA, col = "#333333", lwd = 2))
  for (k in names(yc)) grid.text(k, x0-.012, yc[k], just = "right", gp = gpar(fontsize = 18))
  for (i in seq_len(44L)) {
    grid.text(genes$gene[i], x0+(i-.5)*cw, .637, rot = 45, just = c("right","top"),
              gp = gpar(fontsize = 12))
    grid.rect(x0+(i-.5)*cw, .888, cw, .04,
              gp = gpar(fill = module_colors[genes$module_id[i]], col = NA))
  }
  grid.roundrect(294/687, .276, 522/687, .372, r = unit(9*w/687, "in"),
                 gp = gpar(fill = "#E0E0E0", col = "#B0B0B0", lwd = 2))
  grid.text("Module", 47/687, .276, rot = 90, gp = gpar(fontsize = 18))
  rows <- 1-c(162.5,191.5,223.5)/266
  for (j in seq_len(6L)) {
    x <- if (j <= 3L) 80/687 else 327/687
    tx <- if (j <= 3L) 95/687 else 343/687
    y <- rows[(j-1L)%%3L+1L]
    grid.rect(x, y, unit(17*w/687,"in"), unit(18*w/687,"in"),
              gp = gpar(fill = module_colors[j], col = "#555555", lwd = 2))
    grid.text(paste0("M",j,"  ",module_names[j]), tx, y, just = "left", gp = gpar(fontsize = 14))
  }
  by <- 1-224/266; bh <- (224-163)/266
  for (j in seq_len(150L)) {
    v <- lims[1]+(j-.5)/150*diff(lims)
    grid.rect(611/687, by+(j-.5)/150*bh, 20/687, bh/150+.0003,
              gp = gpar(fill = map_color(v), col = NA))
  }
  grid.rect(611/687, by+bh/2, 20/687, bh, gp = gpar(fill = NA, col = "#333333", lwd = 2))
  for (v in c(-1.5, 0, 1.5)) {
    grid.text(format(v), 628/687, by+(v-lims[1])/diff(lims)*bh, just = "left", gp = gpar(fontsize = 14))
  }
  grid.text("Relative Expression", 578/687, by+bh/2, rot = 90, gp = gpar(fontsize = 16))
}
out <- args[2]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop("Use an empty output directory.")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(genes, means, check.names = FALSE),
          file.path(out, "condition_mean_normalized_expression.csv"), row.names = FALSE)
write.csv(data.frame(genes, z, check.names = FALSE),
          file.path(out, "gene_standardized_expression.csv"), row.names = FALSE)
stem <- file.path(out, "module_expression_heatmap")
cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = "Arial", bg = "white")
draw(); dev.off()
agg_png(paste0(stem, ".png"), width = w, height = h, units = "in", res = 600, background = "white")
draw(); dev.off()
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h, bg = "white")
  draw(); dev.off()
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
