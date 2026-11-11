# Render the manuscript global UMAP, major-class marker dots and marker feature panels.
# Usage: Rscript render_major_markers.R SIX_CLASS_OBJECT.rds OUTPUT_DIR
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(grid); library(scattermore)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript render_major_markers.R SIX_CLASS_OBJECT.rds OUTPUT_DIR")
obj <- readRDS(args[1])
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stems <- c("global_celltype_umap", "major_marker_dotplot", "major_marker_feature_umap")
if (any(file.exists(outer(file.path(out, stems), c(".pdf", ".png", ".svg"), paste0))) ||
    file.exists(file.path(out, "major_marker_dotplot_source.csv"))) {
  stop("A marker output already exists. Choose a new output directory.")
}
classes <- c("Excitatory", "Inhibitory", "Oligodendrocytes", "OPC", "Astrocytes", "Microglia")
display_names <- setNames(c("Excitatory Neuron", "Inhibitory Neuron", "Oligodendrocyte",
                    "OPC", "Astrocyte", "Microglia"), classes)
palette <- setNames(c("#163F7A", "#A51B28", "#684018", "#542487", "#007044", "#B84A10"), classes)
dot_genes <- c("Slc17a7", "Satb2", "Neurod6", "Gad1", "Gad2", "Dlx1", "Mog", "Mbp", "St18",
               "Cacng4", "Pdgfra", "Vcan", "Gfap", "Aqp4", "P2ry12", "Tmem119", "Cx3cr1")
feature_genes <- c("Slc17a7", "Satb2", "Gad1", "Gad2", "Mog", "Mbp", "Cacng4", "Pdgfra",
                   "Aqp4", "Gfap", "P2ry12", "Tmem119")
xy <- SeuratObject::Embeddings(obj[["umap.global"]])
md <- obj[[]]
cells <- rownames(xy)
stopifnot(ncol(xy) >= 2L, all(cells %in% rownames(md)), !anyDuplicated(cells))
d <- data.frame(cell_id = cells, UMAP_1 = xy[, 1], UMAP_2 = xy[, 2],
                group = factor(md[cells, "major_class_six"], levels = classes))
if (anyNA(d$group) || any(table(d$group) == 0L)) stop("Missing accepted broad-class assignments.")
genes <- unique(c(dot_genes, feature_genes))
if (!all(genes %in% rownames(obj[["RNA"]]))) {
  stop("Missing manuscript marker features: ", paste(setdiff(genes, rownames(obj[["RNA"]])), collapse = ", "))
}
expression <- SeuratObject::LayerData(obj[["RNA"]], layer = "data", features = genes)
expression <- expression[genes, cells, drop = FALSE]
stopifnot(identical(colnames(expression), d$cell_id))
rm(obj); invisible(gc())

save_panel <- function(draw, name, width, height) {
  stem <- file.path(out, name)
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width, height, family = "Arial", bg = "white")
  tryCatch(draw(), finally = dev.off())
  ragg::agg_png(paste0(stem, ".png"), width, height, units = "in", res = 600, background = "white")
  tryCatch(draw(), finally = dev.off())
  svglite::svglite(paste0(stem, ".svg"), width, height, bg = "white")
  tryCatch(draw(), finally = dev.off())
}

# Deterministic point order and anchors; never refit the embedding.
lab <- aggregate(d[, c("UMAP_1", "UMAP_2")], list(group = d$group), median)
lab$displayed_name <- unname(display_names[as.character(lab$group)])
draw_order <- order((seq_len(nrow(d)) * 104729) %% 1000003)
xr <- range(d$UMAP_1); yr <- range(d$UMAP_2)
dx <- diff(xr); dy <- diff(yr)
x0 <- xr[1] - .04 * dx; y0 <- yr[1] - .04 * dy
extra_mm <- 1 / (72.27 / 25.4)
umap <- ggplot(d[draw_order, ], aes(UMAP_1, UMAP_2, color = group)) +
  geom_scattermore(pointsize = 2.4, pixels = c(1600, 1600), alpha = 1) +
  scale_color_manual(values = palette, drop = FALSE) +
  ggrepel::geom_text_repel(data = lab, aes(label = displayed_name), family = "Arial", color = "black",
    size = 1.95 + extra_mm, seed = 42, max.overlaps = Inf, segment.color = NA,
    box.padding = .08, point.padding = .05, min.segment.length = Inf, show.legend = FALSE,
    bg.colour = "white", bg.r = .08) +
  annotate("segment", x = x0, xend = x0 + .15 * dx, y = y0, yend = y0, linewidth = .4,
           arrow = arrow(length = unit(1, "mm"), type = "closed")) +
  annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + .15 * dy, linewidth = .4,
           arrow = arrow(length = unit(1, "mm"), type = "closed")) +
  annotate("text", x = x0 + .07 * dx, y = y0 - .03 * dy, label = "UMAP-1",
           size = 1.7 + extra_mm, vjust = 1, family = "Arial") +
  annotate("text", x = x0 - .03 * dx, y = y0 + .07 * dy, label = "UMAP-2",
           angle = 90, size = 1.7 + extra_mm, vjust = 0, family = "Arial") +
  coord_fixed(clip = "off") +
  theme_void(base_family = "Arial") +
  theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2),
        plot.background = element_rect(fill = "white", color = NA))
save_panel(function() print(umap), "global_celltype_umap", 2.92, 2.67)

# Group-average expression is calculated on the linear normalized scale, then log-transformed and standardized by gene.
dots <- do.call(rbind, lapply(dot_genes, function(g) {
  raw <- as.numeric(expression[g, ])
  z <- do.call(rbind, lapply(classes, function(cl) {
    v <- raw[d$group == cl]
    data.frame(gene = g, group = cl, n = length(v), n_detected = sum(v > 0),
               pct_expressed = 100 * mean(v > 0), mean_log_normalized = mean(v),
               mean_linear_normalized = mean(expm1(v)))
  }))
  z$log1p_mean_linear <- log1p(z$mean_linear_normalized)
  z$avg_expression_scaled <- pmax(-2.5, pmin(2.5, as.numeric(scale(z$log1p_mean_linear))))
  z
}))
dots$gene <- factor(dots$gene, levels = dot_genes)
dots$group <- factor(dots$group, levels = classes)
dots$pct_display <- ifelse(dots$pct_expressed < 5, NA_real_, dots$pct_expressed)
write.csv(dots, file.path(out, "major_marker_dotplot_source.csv"), row.names = FALSE)
size <- 7.6
dot <- ggplot(dots, aes(gene, group, size = pct_display, color = avg_expression_scaled)) +
  geom_point(shape = 16, stroke = 0, na.rm = TRUE) +
  scale_y_discrete(labels = function(x) unname(display_names[x])) +
  scale_color_gradientn(colors = c("#D3D3D3", "#EF0038"), limits = c(0, 2), oob = scales::squish,
                        breaks = 0:2, name = "Average Expression") +
  scale_size_area(max_size = 3.9, limits = c(0, 100), breaks = c(25, 50, 75), name = "Percent Expressed") +
  labs(x = NULL, y = NULL) + theme_classic(base_family = "Arial", base_size = size) +
  theme(text = element_text(family = "Arial", face = "plain"),
        axis.text = element_text(size = size, color = "black", face = "plain"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.line = element_line(linewidth = .65, color = "black"),
        axis.ticks = element_line(linewidth = .585, color = "black"),
        legend.position = "bottom", legend.box = "horizontal", legend.spacing.x = unit(1, "mm"),
        legend.key.width = unit(3, "mm"), legend.key.height = unit(3, "mm"),
        legend.text = element_text(size = size - 1), legend.title = element_text(size = size - 1),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(30, 3, 16, 16)) +
  guides(color = guide_colorbar(title.position = "top", direction = "horizontal",
                                barwidth = unit(22, "mm"), barheight = unit(2, "mm"), order = 1),
         size = guide_legend(nrow = 1, title.position = "left", order = 2, override.aes = list(color = "black")))
save_panel(function() print(dot), "major_marker_dotplot", 5.35, 2.87)

feature <- function(gene) {
  raw <- as.numeric(expression[gene, ])
  pos <- which(raw > 0)
  if (!length(pos)) stop("No positive nuclei for displayed marker: ", gene)
  cap <- quantile(raw[pos], .95, type = 7, names = FALSE)
  order <- pos[order(raw[pos], pos)]
  foreground <- data.frame(d[order, c("UMAP_1", "UMAP_2")], expression = raw[order])
  ggplot() +
    geom_scattermore(data = d, aes(UMAP_1, UMAP_2), color = "#9E9E9E", pointsize = 1.1, pixels = c(800, 800)) +
    geom_scattermore(data = foreground, aes(UMAP_1, UMAP_2, color = expression),
                     pointsize = 1.1, pixels = c(800, 800)) +
    scale_color_gradientn(colors = c("#9E9E9E", "#9155B6", "#640097"), limits = c(0, cap),
                          oob = scales::squish, breaks = c(0, cap / 2, cap),
                          labels = c("0", format(signif(cap / 2, 2), trim = TRUE),
                                     paste0("\u2265", format(signif(cap, 2), trim = TRUE))), name = NULL) +
    coord_fixed() + theme_void(base_family = "Arial") + labs(title = gene) +
    theme(text = element_text(face = "bold"), plot.title = element_text(size = 10, face = "bold", hjust = .5),
          legend.position = "right", legend.text = element_text(size = 5, face = "bold"),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
          legend.spacing.x = unit(.2, "mm"), plot.margin = margin(1, 1, 1, 1)) +
    guides(color = guide_colorbar(direction = "vertical", barheight = unit(19, "mm"),
                                  barwidth = unit(2, "mm"), ticks = FALSE))
}
feature_grid <- patchwork::wrap_plots(lapply(feature_genes, feature), ncol = 4)
draw_feature_grid <- function() {
  grid.newpage()
  grid.draw(patchwork::patchworkGrob(feature_grid))
  gp <- gpar(col = "black", fill = "black", lwd = .56 * 72.27 / 25.4, lineend = "butt")
  grid.segments(unit(.11, "in"), unit(.12, "in"), unit(.73, "in"), unit(.12, "in"),
                arrow = arrow(length = unit(1.3, "mm"), type = "closed"), gp = gp)
  grid.segments(unit(.11, "in"), unit(.12, "in"), unit(.11, "in"), unit(.77, "in"),
                arrow = arrow(length = unit(1.3, "mm"), type = "closed"), gp = gp)
  tp <- gpar(fontfamily = "Arial", fontface = "bold", fontsize = 6.3, col = "black")
  grid.text("UMAP-1", x = unit(.42, "in"), y = unit(.047, "in"), gp = tp)
  grid.text("UMAP-2", x = unit(.045, "in"), y = unit(.445, "in"), rot = 90, gp = tp)
}
save_panel(draw_feature_grid, "major_marker_feature_umap", 7, 4.9)
