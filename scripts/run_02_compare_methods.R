#!/usr/bin/env Rscript
# ============================================================
# scripts/run_02_compare_methods.R
# ------------------------------------------------------------
# Loads the fitted model chain (from run_01) and runs all four
# decomposition methods -- naive, OFAT, sequential, Shapley -- for
# the configured baseline-vs-comparison year pair. Saves a combined
# results table and a comparison plot.
#
# Run from the project root, after run_01_fit_model.R:
#   Rscript scripts/run_02_compare_methods.R
# ============================================================

source("R/99_load_all.R")

cfg <- load_config("config.R")
make_plots <- !(exists(".NFZ_SKIP_PLOTS") && isTRUE(.NFZ_SKIP_PLOTS)) && !("--no-plots" %in% commandArgs(trailingOnly = TRUE))
dat <- readRDS("output/rds/aggregated_data.rds")
fit <- readRDS("output/rds/model_fit.rds")

baseline_year   <- cfg$data$baseline_year
comparison_year <- max(dat$year)

disease_label <- unique(as.character(dat$icd10))[1]
message("Comparing ", baseline_year, " vs ", comparison_year,
        " for disease group '", disease_label, "'")

df_pair <- dat %>% dplyr::filter(year %in% c(baseline_year, comparison_year))

message("\n== Naive (log-share) ==")
res_naive <- decompose_naive_2y(df_pair, NULL, baseline_year, comparison_year)
print(res_naive)

message("\n== OFAT ==")
res_ofat <- decompose_ofat_2y(df_pair, fit, baseline_year, comparison_year)
print(res_ofat)

message("\n== Sequential (single order: ", paste(cfg$decomposition$sequential_order, collapse = " -> "), ") ==")
res_sequential <- decompose_sequential_2y(
  df_pair, fit, baseline_year, comparison_year,
  order = cfg$decomposition$sequential_order
)
print(res_sequential)

message("\n== Shapley (average of all 24 orderings) ==")
res_shapley <- decompose_shapley_2y(df_pair, fit, baseline_year, comparison_year)
print(res_shapley)

stacked <- stack_methods(res_naive, res_ofat, res_sequential, res_shapley)

# Keep observed raw totals next to the decomposition totals. For naive the
# totals are observed aggregate totals; for OFAT/sequential/Shapley they are
# model-implied V(empty)/V(full) totals. This avoids the common manuscript
# ambiguity between "observed total change" and "model-implied total change".
raw_identity_pair <- make_raw_identity_table(df_pair, years = c(baseline_year, comparison_year))
model_vs_observed <- make_model_vs_observed_table(stacked, raw_identity_pair)
stacked <- stacked %>%
  dplyr::left_join(
    model_vs_observed %>%
      dplyr::select(method, decomposition_total_type, observed_baseline_total,
                    observed_comparison_total, observed_total_change,
                    change_difference_vs_observed, change_difference_pct_observed),
    by = "method"
  )

ensure_dir(cfg$output$tables_dir)
readr::write_csv(stacked, file.path(cfg$output$tables_dir, "method_comparison.csv"))
readr::write_csv(raw_identity_pair, file.path(cfg$output$tables_dir, "SUPP_raw_identity_baseline_comparison.csv"))
readr::write_csv(model_vs_observed, file.path(cfg$output$tables_dir, "TABLE_model_vs_observed_change.csv"))

if (make_plots) {
  p1 <- plot_method_comparison_bars(stacked)
  p2 <- plot_residual_comparison(stacked)

  ensure_dir(cfg$output$figures_dir)
  safe_ggsave(p1, file.path(cfg$output$figures_dir, "method_comparison_bars.png"),
              cfg$output$fig_width, cfg$output$fig_height, cfg$output$fig_dpi)
  safe_ggsave(p2, file.path(cfg$output$figures_dir, "residual_comparison.png"),
              cfg$output$fig_width, cfg$output$fig_height, cfg$output$fig_dpi)
} else {
  message("Skipping method-comparison plots because .NFZ_SKIP_PLOTS=TRUE or --no-plots was used.")
}

message("\nSaved: ", file.path(cfg$output$tables_dir, "method_comparison.csv"))
if (make_plots) {
  message("Saved: ", file.path(cfg$output$figures_dir, "method_comparison_bars.png"))
  message("Saved: ", file.path(cfg$output$figures_dir, "residual_comparison.png"))
}
message("Done.")
