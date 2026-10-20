#!/usr/bin/env Rscript
# ============================================================
# scripts/run_03_order_sensitivity.R
# ------------------------------------------------------------
# Runs the sequential decomposition under all 24 possible
# component orderings and compares the resulting range to the
# order-independent Shapley value. This is the empirical evidence
# motivating Shapley over a single arbitrary sequential order.
#
# Run from the project root, after run_01_fit_model.R:
#   Rscript scripts/run_03_order_sensitivity.R
# ============================================================

source("R/99_load_all.R")

cfg <- load_config("config.R")
make_plots <- !(exists(".NFZ_SKIP_PLOTS") && isTRUE(.NFZ_SKIP_PLOTS)) && !("--no-plots" %in% commandArgs(trailingOnly = TRUE))
dat <- readRDS("output/rds/aggregated_data.rds")
fit <- readRDS("output/rds/model_fit.rds")

baseline_year   <- cfg$data$baseline_year
comparison_year <- max(dat$year)

df_pair <- dat %>% dplyr::filter(year %in% c(baseline_year, comparison_year))

message("Running all 24 sequential orderings + Shapley for ",
        baseline_year, " vs ", comparison_year, " ...")
result <- compare_order_sensitivity(df_pair, fit, baseline_year, comparison_year)

message("\nRange of the sequential effect across 24 orderings vs. Shapley value:")
print(result$range_summary)

ensure_dir(cfg$output$tables_dir)
readr::write_csv(result$all_orders, file.path(cfg$output$tables_dir, "sequential_all_24_orders.csv"))
readr::write_csv(result$range_summary, file.path(cfg$output$tables_dir, "order_sensitivity_summary.csv"))

if (make_plots) {
  p <- plot_order_sensitivity(result$range_summary)
  ensure_dir(cfg$output$figures_dir)
  safe_ggsave(p, file.path(cfg$output$figures_dir, "order_sensitivity.png"),
              cfg$output$fig_width, cfg$output$fig_height, cfg$output$fig_dpi)
} else {
  message("Skipping order-sensitivity plot because .NFZ_SKIP_PLOTS=TRUE or --no-plots was used.")
}

message("\nSaved: ", file.path(cfg$output$tables_dir, "sequential_all_24_orders.csv"))
message("Saved: ", file.path(cfg$output$tables_dir, "order_sensitivity_summary.csv"))
if (make_plots) message("Saved: ", file.path(cfg$output$figures_dir, "order_sensitivity.png"))
message("Done.")
