# ============================================================
# R/99_vih_figure_export_helpers.R
# ============================================================
# Helpers for Value in Health / ISPOR submission-ready figures.
#
# Implemented rules:
# - each figure is exported as a separate image file;
# - no title, subtitle, or caption is burned into the image;
# - PNG by default, TIFF available via --tiff;
# - 600 DPI by default via config$output$fig_dpi;
# - axis labels use explicit units;
# - figure titles and explanations are written to a separate Figure Legends file;
# - a manifest/checklist is written for submission QA.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# IMPORTANT: this function intentionally does not add scale_x_continuous()
# or scale_y_continuous(). Axes may be discrete in several plots; scaling
# is applied only in figure-specific code where the axis type is known.
vih_clean_plot <- function(
  p,
  base_size = 10,
  legend_position = "bottom",
  remove_titles = TRUE
) {
  p2 <- p

  if (isTRUE(remove_titles)) {
    p2 <- p2 + labs(title = NULL, subtitle = NULL, caption = NULL)
  }

  p2 +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = legend_position,
      legend.title = element_blank(),
      legend.box = "horizontal",
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = max(base_size - 1, 7)),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = margin(6, 6, 6, 6)
    )
}

vih_safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

vih_billion_labels <- function(x) {
  # European/biomedical-friendly: numeric values in billions, no scientific notation.
  scales::label_number(accuracy = 0.1, big.mark = ",", trim = TRUE)(x / 1e9)
}

vih_million_labels <- function(x) {
  scales::label_number(accuracy = 0.1, big.mark = ",", trim = TRUE)(x / 1e6)
}

vih_number_labels <- function(x) {
  scales::label_number(big.mark = ",", trim = TRUE)(x)
}

vih_apply_points_billion_axis <- function(p, y_label = "Effect on total settlement points (billion)", include_zero = TRUE) {
  p +
    scale_y_continuous(
      labels = vih_billion_labels,
      limits = if (isTRUE(include_zero)) function(x) range(c(0, x), na.rm = TRUE) else waiver()
    ) +
    labs(y = y_label)
}

vih_apply_points_million_axis <- function(p, y_label = "Root mean squared error (million settlement points)", include_zero = TRUE) {
  p +
    scale_y_continuous(
      labels = vih_million_labels,
      limits = if (isTRUE(include_zero)) function(x) range(c(0, x), na.rm = TRUE) else waiver()
    ) +
    labs(y = y_label)
}

vih_save_figure <- function(
  plot,
  figure_no,
  stem,
  width = 6.5,
  height = 4.2,
  dpi = 600,
  path = "submission/figures",
  device = c("png", "tiff", "jpg"),
  base_size = 10,
  legend_position = "bottom",
  prefix = "Figure"
) {
  device <- match.arg(device)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  filename <- sprintf(
    "%s_%02d_%s.%s",
    prefix,
    as.integer(figure_no),
    vih_safe_filename(stem),
    device
  )
  out_file <- file.path(path, filename)

  p_clean <- vih_clean_plot(
    plot,
    base_size = base_size,
    legend_position = legend_position
  )

  ggsave(
    filename = out_file,
    plot = p_clean,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    device = device,
    limitsize = FALSE
  )

  message("Saved submission figure: ", out_file)
  invisible(out_file)
}

vih_save_supplement_figure <- function(
  plot,
  supp_no,
  stem,
  width = 6.5,
  height = 4.2,
  dpi = 600,
  path = "submission/figures/supplement",
  device = c("png", "tiff", "jpg"),
  base_size = 10,
  legend_position = "bottom"
) {
  vih_save_figure(
    plot = plot,
    figure_no = supp_no,
    stem = stem,
    width = width,
    height = height,
    dpi = dpi,
    path = path,
    device = device,
    base_size = base_size,
    legend_position = legend_position,
    prefix = "Supplementary_Figure"
  )
}

vih_write_manifest <- function(manifest, path = "submission/figures/FIGURE_MANIFEST.csv") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(manifest, path)
  message("Saved figure manifest: ", path)
  invisible(path)
}

vih_write_checklist <- function(
  path = "submission/figures/VIH_FIGURE_QA_CHECKLIST.csv",
  dpi = 600,
  device = "png",
  main_figures_n = NA_integer_,
  supplementary_figures_n = NA_integer_
) {
  checklist <- data.table::data.table(
    item = c(
      "Each main figure is a separate image file",
      "No figure is embedded in the main manuscript by this script",
      "No title, subtitle, or caption is burned into exported images",
      "Figure legends are written separately",
      "Main-figure numbering uses Arabic numerals",
      "Resolution is at least 300 DPI and set to configured DPI",
      "Default output format is journal-acceptable PNG; TIFF available with --tiff",
      "Main figures limited to the manuscript display-item plan",
      "Supplementary figures placed in a separate supplement folder"
    ),
    value = c(
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      dpi >= 300,
      device %in% c("png", "tiff", "jpg"),
      main_figures_n <= 3L,
      TRUE
    ),
    detail = c(
      "submission/figures/Figure_*.png or .tiff",
      "Main manuscript should contain text callouts only: Figure 1, Figure 2, Figure 3",
      "Titles are removed by vih_clean_plot(); titles belong in Figure_Legends_VIH.*",
      "Use scripts/run_98b_make_figure_legends.R",
      "Figure_01, Figure_02, Figure_03",
      paste0("DPI = ", dpi),
      paste0("Device = ", device),
      paste0("Main figures exported = ", main_figures_n),
      paste0("Supplementary figures exported = ", supplementary_figures_n)
    )
  )
  data.table::fwrite(checklist, path)
  message("Saved figure QA checklist: ", path)
  invisible(path)
}
