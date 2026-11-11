# Shared Arial figure export. Fonts are supplied by the operating system.
load_figure_packages <- function(extra = character()) {
  require_packages(unique(c("grid", "ragg", "svglite", "systemfonts", extra)))
  stopifnot(capabilities("cairo"))
  regular <- systemfonts::match_fonts("Arial")
  bold <- systemfonts::match_fonts("Arial", weight = "bold")
  paths <- tolower(basename(c(regular$path, bold$path)))
  if (length(paths) != 2L || !all(grepl("arial", paths))) {
    stop("Install Arial Regular and Arial Bold before rendering.")
  }
}

export_plot <- function(stem, width, height, draw) {
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width, height, family = "Arial", bg = "white")
  tryCatch(draw(), finally = grDevices::dev.off())
  ragg::agg_png(paste0(stem, ".png"), width, height, units = "in", res = 600, background = "white")
  tryCatch(draw(), finally = grDevices::dev.off())
  svglite::svglite(paste0(stem, ".svg"), width, height, bg = "white")
  tryCatch(draw(), finally = grDevices::dev.off())
}
