# Render condition-pooled descriptive composition of major classes and IEG-high subtypes.
# Usage: Rscript render_composition.R GLOBAL_ANNOTATION.csv.gz NEURON_METADATA.csv OUTPUT_DIR
suppressPackageStartupMessages({library(ggplot2); library(grid)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript render_composition.R GLOBAL_ANNOTATION.csv.gz NEURON_METADATA.csv OUTPUT_DIR")
global <- read.csv(args[1])
neuron <- read.csv(args[2])
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stems <- c("major_class_composition", "excitatory_ieghigh_composition", "inhibitory_ieghigh_composition")
if (any(file.exists(outer(file.path(out, stems),
                          c(".pdf", ".png", ".svg", "_all_conditions.csv", "_displayed.csv"), paste0)))) {
  stop("A composition output already exists. Choose a new output directory.")
}
stopifnot(all(c("cell_id", "condition_label", "proposed_major_class") %in% names(global)),
          all(c("Cell", "condition", "neuron_type", "subtype", "IEGhigh") %in% names(neuron)),
          !anyDuplicated(global$cell_id), !anyDuplicated(neuron$Cell),
          is.logical(neuron$IEGhigh), !anyNA(neuron$IEGhigh))
major <- c(Excitatory = "#3B6FBF", Inhibitory = "#E85D5D", Oligodendrocytes = "#B5895A",
           OPC = "#9A6FD1", Astrocytes = "#4DB894", Microglia = "#F0986D")
major_names <- c("Excitatory Neuron", "Inhibitory Neuron", "Oligodendrocyte", "OPC", "Astrocyte", "Microglia")
exc <- c(GLU_L6_CLA = "#05BE78", GLU_L5_6 = "#FFFFA0", GLU_L2_L3 = "#F7ABE8", GLU_L5 = "#FF8080",
         GLU_L4_5 = "#90BFF9", GLU_ADGM = "#FF8000")
inh <- c(GABA_Amb = "#66C2A5", GABA_Etv1 = "#FC8D62", GABA_Lamp5 = "#8DA0CB", GABA_MSN = "#E78AC3",
         GABA_Sst = "#A6D854", GABA_Hs3st2 = "#FFD92F", GABA_Pvalb = "#E5C494", GABA_Vip = "#FF6347")
composition <- function(condition, labels, palette) {
  conditions <- c("WD", "Ext", "Ret-Ext")
  stopifnot(all(condition %in% conditions), all(labels %in% names(palette)))
  d <- as.data.frame(table(condition = factor(condition, levels = conditions),
                           label = factor(labels, levels = names(palette))), stringsAsFactors = FALSE)
  names(d)[3] <- "n"
  total <- ave(d$n, d$condition, FUN = sum)
  if (any(total == 0L)) stop("A condition has no eligible nuclei.")
  d$percent <- 100 * d$n / total
  d
}
thin_key <- function(data, params, size) {
  grid::rectGrob(width = unit(.85, "npc"), height = unit(.18, "npc"), gp = gpar(fill = data$fill, col = NA))
}
stack <- function(d, palette, labels, ylabel, reverse, size, plain, margin, legend_height, border) {
  d <- d[d$condition %in% c("Ext", "Ret-Ext"), ]
  d$condition <- factor(d$condition, levels = c("Ext", "Ret-Ext"))
  d$label <- factor(d$label, levels = names(palette))
  stopifnot(all(abs(tapply(d$percent, d$condition, sum) - 100) < 1e-7))
  ggplot(d, aes(condition, percent, fill = label)) +
    geom_col(width = .65, position = position_stack(reverse = reverse), key_glyph = thin_key) +
    scale_fill_manual(values = palette, labels = labels, drop = FALSE, name = "Subtype") +
    scale_y_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = 0)) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(x = NULL, y = ylabel) +
    theme_classic(base_family = "Arial", base_size = size) +
    theme(text = element_text(family = "Arial", face = if (plain) "plain" else "bold"),
          axis.text = element_text(size = size, color = "black"), axis.title = element_text(size = size + 1),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(linewidth = .65, color = "black"),
          axis.ticks = element_line(linewidth = .585, color = "black"),
          panel.border = element_rect(fill = NA, linewidth = border, color = "black"),
          legend.key.height = unit(legend_height, "mm"), legend.key.width = unit(4, "mm"),
          legend.text = element_text(size = size - 1), legend.margin = margin(0, 0, 0, 0),
          legend.justification = "top", plot.margin = margin)
}
save_panel <- function(p, data, name, width, height) {
  stem <- file.path(out, name)
  write.csv(data, paste0(stem, "_all_conditions.csv"), row.names = FALSE)
  write.csv(p$data, paste0(stem, "_displayed.csv"), row.names = FALSE)
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width, height, family = "Arial", bg = "white")
  tryCatch(print(p), finally = dev.off())
  ragg::agg_png(paste0(stem, ".png"), width, height, units = "in", res = 600, background = "white")
  tryCatch(print(p), finally = dev.off())
  svglite::svglite(paste0(stem, ".svg"), width, height, bg = "white")
  tryCatch(print(p), finally = dev.off())
}
d <- composition(global$condition_label, global$proposed_major_class, major)
p <- stack(d, major, major_names, "Percentage (%)", TRUE, 8, TRUE, margin(21, 13, 23, 21), 4.8, 1.10)
save_panel(p, d, "major_class_composition", 3.18, 2.82)
for (lineage in c("Excitatory", "Inhibitory")) {
  is_exc <- lineage == "Excitatory"
  palette <- if (is_exc) exc else inh
  md <- neuron[neuron$neuron_type == lineage & neuron$IEGhigh, ]
  d <- composition(md$condition, md$subtype, palette)
  p <- stack(d, palette, gsub("_", " ", names(palette), fixed = TRUE),
             "Percentage of IEG-high Cells (%)", is_exc,
             if (is_exc) 8.5 else 10, FALSE,
             if (is_exc) margin(26, 20, 16.5, 34) else margin(46, 42, 26, 58), 3, .72)
  save_panel(p, d, paste0(tolower(lineage), "_ieghigh_composition"),
             if (is_exc) 3.29 else 4.71, if (is_exc) 2.84 else 4.04)
}
