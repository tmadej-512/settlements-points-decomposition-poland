#!/usr/bin/env Rscript
# ============================================================
# scripts/run_01b_residual_diagnostics.R
# ------------------------------------------------------------
# Residual diagnostics for the main fitted models.
# Run this after scripts/run_01_fit_model.R.
# ============================================================

source("R/99_load_all.R")

cfg <- load_config("config.R")
ensure_dir(cfg$output$tables_dir)
ensure_dir(cfg$output$figures_dir)
ensure_dir(file.path(cfg$output$figures_dir, "model_diagnostics"))

fit_path <- "output/rds/model_fit.rds"
dat_path <- "output/rds/aggregated_data.rds"
if (!file.exists(fit_path) || !file.exists(dat_path)) {
  stop("Run scripts/run_01_fit_model.R first (need output/rds/model_fit.rds and output/rds/aggregated_data.rds).", call. = FALSE)
}

fit <- readRDS(fit_path)
dat <- readRDS(dat_path)

message("== Residual diagnostics for main models ==")
resid_df <- extract_main_model_residuals(fit, dat)
diag_summary <- summarise_main_model_diagnostics(fit, resid_df)
print(diag_summary)

readr::write_csv(resid_df, file.path(cfg$output$tables_dir, "model_residuals_long.csv"))
readr::write_csv(diag_summary, file.path(cfg$output$tables_dir, "model_residual_diagnostics_summary.csv"))

p1 <- plot_main_model_residuals_vs_fitted(resid_df)
p2 <- plot_main_model_residual_histogram(resid_df)
p3 <- plot_main_model_residual_qq(resid_df)
p4 <- plot_main_model_residuals_by_year(resid_df)

out_dir <- file.path(cfg$output$figures_dir, "model_diagnostics")
safe_ggsave(p1, file.path(out_dir, "model_residuals_vs_fitted.png"), cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
safe_ggsave(p2, file.path(out_dir, "model_residual_histogram.png"), cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
safe_ggsave(p3, file.path(out_dir, "model_residual_qq.png"), cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
safe_ggsave(p4, file.path(out_dir, "model_residuals_by_year.png"), cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)

if (requireNamespace("patchwork", quietly = TRUE)) {
  suppressPackageStartupMessages(library(patchwork))
  p_combined <- (p1 | p3) / (p2 | p4) + patchwork::plot_annotation(tag_levels = "A")
  safe_ggsave(p_combined, file.path(out_dir, "model_residual_diagnostics_combined.png"), cfg$output$fig_width + 4, cfg$output$fig_height + 2, cfg$output$fig_dpi)
}

message("Saved residual diagnostics tables:")
message("  ", file.path(cfg$output$tables_dir, "model_residuals_long.csv"))
message("  ", file.path(cfg$output$tables_dir, "model_residual_diagnostics_summary.csv"))

message("Saved residual diagnostics figures in: ", out_dir)

dharma_out_dir <- file.path(out_dir, "DHARMa")
dharma_summary <- make_dharma_diagnostics(
  fit = fit,
  dat = dat,
  out_dir = dharma_out_dir,
  nsim = 250,
  width = 1800,
  height = 1400,
  res = 200
)
if (!is.null(dharma_summary)) {
  readr::write_csv(dharma_summary, file.path(cfg$output$tables_dir, "model_residual_dharma_summary.csv"))
  message("Saved DHARMa diagnostics table:")
  message("  ", file.path(cfg$output$tables_dir, "model_residual_dharma_summary.csv"))
  message("Saved DHARMa figures in: ", dharma_out_dir)
} else {
  message("DHARMa diagnostics skipped (package not installed).")
}
