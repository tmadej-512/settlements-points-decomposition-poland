#!/usr/bin/env Rscript
# ============================================================
# scripts/run_97_make_diagnostics_validation_plots.R
# ------------------------------------------------------------
# One-stop script for the plots that are easy to miss:
#   1. Main-model residual diagnostics, including DHARMa if installed
#   2. Bootstrap interval plot
#   3. Sequential-order sensitivity plots
#
# Run after the main compute step:
#   source("scripts/run_97_make_diagnostics_validation_plots.R")
# ============================================================

source("R/99_load_all.R")

cfg <- load_config("config.R")
tdir <- cfg$output$tables_dir
fdir <- cfg$output$figures_dir
ensure_dir(tdir)
ensure_dir(fdir)
ensure_dir(file.path(fdir, "model_diagnostics"))
ensure_dir(file.path(fdir, "validation"))
ensure_dir(file.path("submission", "figures", "supplement", "diagnostics_validation"))

copy_if_exists <- function(path, dest_dir = file.path("submission", "figures", "supplement", "diagnostics_validation")) {
  if (file.exists(path)) {
    file.copy(path, file.path(dest_dir, basename(path)), overwrite = TRUE)
    TRUE
  } else {
    FALSE
  }
}

manifest <- tibble::tibble(
  artifact = character(),
  path = character(),
  status = character(),
  note = character()
)

add_manifest <- function(artifact, path, status, note = "") {
  manifest <<- dplyr::bind_rows(
    manifest,
    tibble::tibble(artifact = artifact, path = path, status = status, note = note)
  )
}

message("== 1. Main-model diagnostics ==")
fit_path <- file.path("output", "rds", "model_fit.rds")
dat_path <- file.path("output", "rds", "aggregated_data.rds")
if (file.exists(fit_path) && file.exists(dat_path)) {
  fit <- readRDS(fit_path)
  dat <- readRDS(dat_path)

  resid_df <- extract_main_model_residuals(fit, dat)
  diag_summary <- summarise_main_model_diagnostics(fit, resid_df)
  readr::write_csv(resid_df, file.path(tdir, "model_residuals_long.csv"))
  readr::write_csv(diag_summary, file.path(tdir, "model_residual_diagnostics_summary.csv"))

  p1 <- plot_main_model_residuals_vs_fitted(resid_df)
  p2 <- plot_main_model_residual_histogram(resid_df)
  p3 <- plot_main_model_residual_qq(resid_df)
  p4 <- plot_main_model_residuals_by_year(resid_df)

  diag_dir <- file.path(fdir, "model_diagnostics")
  diag_files <- c(
    residuals_vs_fitted = file.path(diag_dir, "model_residuals_vs_fitted.png"),
    residual_histogram  = file.path(diag_dir, "model_residual_histogram.png"),
    residual_qq         = file.path(diag_dir, "model_residual_qq.png"),
    residuals_by_year   = file.path(diag_dir, "model_residuals_by_year.png")
  )

  safe_ggsave(p1, diag_files[["residuals_vs_fitted"]], cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
  safe_ggsave(p2, diag_files[["residual_histogram"]],  cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
  safe_ggsave(p3, diag_files[["residual_qq"]],         cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
  safe_ggsave(p4, diag_files[["residuals_by_year"]],   cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)

  for (nm in names(diag_files)) {
    copy_if_exists(diag_files[[nm]])
    add_manifest(nm, diag_files[[nm]], "created")
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    suppressPackageStartupMessages(library(patchwork))
    p_combined <- (p1 | p3) / (p2 | p4) + patchwork::plot_annotation(tag_levels = "A")
    combined_path <- file.path(diag_dir, "model_residual_diagnostics_combined.png")
    safe_ggsave(p_combined, combined_path, cfg$output$fig_width + 4, cfg$output$fig_height + 2, cfg$output$fig_dpi)
    copy_if_exists(combined_path)
    add_manifest("model_residual_diagnostics_combined", combined_path, "created")
  }

  dharma_dir <- file.path(diag_dir, "DHARMa")
  dharma_summary <- make_dharma_diagnostics(
    fit = fit,
    dat = dat,
    out_dir = dharma_dir,
    nsim = cfg$validation$dharma_nsim %||% 250,
    width = 1800,
    height = 1400,
    res = 200
  )
  if (!is.null(dharma_summary)) {
    readr::write_csv(dharma_summary, file.path(tdir, "model_residual_dharma_summary.csv"))
    dharma_files <- list.files(dharma_dir, pattern = "\\.png$", full.names = TRUE)
    for (p in dharma_files) {
      copy_if_exists(p)
      add_manifest(paste0("DHARMa_", tools::file_path_sans_ext(basename(p))), p, "created")
    }
  } else {
    add_manifest("DHARMa", dharma_dir, "skipped", "Package DHARMa not installed or diagnostics returned NULL.")
  }
} else {
  add_manifest("main_model_diagnostics", file.path(fdir, "model_diagnostics"), "skipped",
               "Missing output/rds/model_fit.rds or output/rds/aggregated_data.rds. Run scripts/run_01_fit_model.R first.")
}

message("== 2. Bootstrap interval plot ==")
boot_path <- file.path(tdir, "bootstrap_summary.csv")
method_path <- file.path(tdir, "method_comparison.csv")
if (file.exists(boot_path) && file.exists(method_path)) {
  boot_summary <- readr::read_csv(boot_path, show_col_types = FALSE)
  method_comparison <- readr::read_csv(method_path, show_col_types = FALSE)
  p_boot <- plot_bootstrap_intervals(boot_summary, method_comparison)
  boot_fig <- file.path(fdir, "bootstrap_intervals.png")
  safe_ggsave(p_boot, boot_fig, cfg$output$fig_width, cfg$output$fig_height, cfg$output$fig_dpi)
  copy_if_exists(boot_fig)
  add_manifest("bootstrap_intervals", boot_fig, "created")
} else {
  add_manifest("bootstrap_intervals", file.path(fdir, "bootstrap_intervals.png"), "skipped",
               paste("Need", boot_path, "and", method_path, ". Run scripts/run_06_bootstrap.R and scripts/run_02_compare_methods.R."))
}

message("== 3. Sequential-order sensitivity plots ==")
order_summary_path <- file.path(tdir, "order_sensitivity_summary.csv")
all_orders_path <- file.path(tdir, "sequential_all_24_orders.csv")
if (file.exists(order_summary_path)) {
  order_summary <- readr::read_csv(order_summary_path, show_col_types = FALSE)
  p_order <- plot_order_sensitivity(order_summary)
  order_fig <- file.path(fdir, "order_sensitivity.png")
  safe_ggsave(p_order, order_fig, cfg$output$fig_width, cfg$output$fig_height, cfg$output$fig_dpi)
  copy_if_exists(order_fig)
  add_manifest("order_sensitivity_range", order_fig, "created")
} else {
  add_manifest("order_sensitivity_range", file.path(fdir, "order_sensitivity.png"), "skipped",
               paste("Need", order_summary_path, ". Run scripts/run_03_order_sensitivity.R."))
}

if (file.exists(all_orders_path)) {
  all_orders <- readr::read_csv(all_orders_path, show_col_types = FALSE)
  p_all_orders <- plot_sequential_all_orders(all_orders)
  all_orders_fig <- file.path(fdir, "sequential_all_24_orders.png")
  safe_ggsave(p_all_orders, all_orders_fig, cfg$output$fig_width + 3, cfg$output$fig_height + 1, cfg$output$fig_dpi)
  copy_if_exists(all_orders_fig)
  add_manifest("sequential_all_24_orders", all_orders_fig, "created")
} else {
  add_manifest("sequential_all_24_orders", file.path(fdir, "sequential_all_24_orders.png"), "skipped",
               paste("Need", all_orders_path, ". Run scripts/run_03_order_sensitivity.R."))
}

manifest_path <- file.path(tdir, "diagnostics_validation_plot_manifest.csv")
readr::write_csv(manifest, manifest_path)
message("Saved manifest: ", manifest_path)
message("Copied created figures to: submission/figures/supplement/diagnostics_validation")
print(manifest)
