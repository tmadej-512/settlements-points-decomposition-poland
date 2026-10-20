#!/usr/bin/env Rscript
# ============================================================
# scripts/run_10_value_in_health_outputs.R
# ------------------------------------------------------------
# Extra manuscript-ready outputs for a Value-in-Health-style methods paper.
#
# Run after:
#   scripts/run_01_fit_model.R
#   scripts/run_02_compare_methods.R
#   scripts/run_05b_simulation_two_scenarios.R
#
# It writes:
#   TABLE_methods_feature_comparison.csv
#   TABLE_method_decision_guide.csv
#   TABLE_simulation_performance_overall.csv
#   TABLE_shapley_vs_lmdi_simulation.csv
#   TABLE_empirical_shapley_vs_lmdi.csv
#   TABLE_empirical_multiyear_lmdi_shapley.csv
# and figures:
#   FIGURE_simulation_overall_rmse.png
#   FIGURE_shapley_vs_lmdi_difference.png
#   FIGURE_empirical_multiyear_lmdi_shapley.png
# ============================================================

source("R/99_load_all.R")

cfg <- load_config("config.R")
make_plots <- !(exists(".NFZ_SKIP_PLOTS") && isTRUE(.NFZ_SKIP_PLOTS)) && !("--no-plots" %in% commandArgs(trailingOnly = TRUE))
tdir <- cfg$output$tables_dir
fdir <- cfg$output$figures_dir
ensure_dir(tdir)
ensure_dir(fdir)

# ---- 1. Conceptual/methodological tables ----
methods_feature <- make_methods_feature_table()
decision_guide <- make_method_decision_table()

readr::write_csv(methods_feature, file.path(tdir, "TABLE_methods_feature_comparison.csv"))
readr::write_csv(decision_guide, file.path(tdir, "TABLE_method_decision_guide.csv"))

# ---- 2. Simulation performance tables and figures ----
sim_path <- file.path(tdir, "simulation_three_scenarios_raw.csv")
if (!file.exists(sim_path)) {
  stop("Missing ", sim_path, ". Run scripts/run_05b_simulation_two_scenarios.R first.")
}
sim_raw <- readr::read_csv(sim_path, show_col_types = FALSE)

# If older outputs lack scenario_label, create a readable fallback.
if (!"scenario_label" %in% names(sim_raw)) {
  sim_raw <- sim_raw %>%
    dplyr::mutate(scenario_label = .data$scenario)
}

sim_perf <- make_simulation_performance_table(sim_raw)
shapley_vs_lmdi_sim <- make_shapley_vs_lmdi_simulation_table(sim_raw)

readr::write_csv(sim_perf, file.path(tdir, "TABLE_simulation_performance_overall.csv"))
readr::write_csv(shapley_vs_lmdi_sim, file.path(tdir, "TABLE_shapley_vs_lmdi_simulation.csv"))

if (make_plots) {
  p_rmse <- plot_simulation_overall_rmse(sim_perf)
  p_diff <- plot_shapley_vs_lmdi_difference(sim_raw)

  safe_ggsave(p_rmse, file.path(fdir, "FIGURE_simulation_overall_rmse.png"),
              cfg$output$fig_width + 3, cfg$output$fig_height + 1, cfg$output$fig_dpi)
  safe_ggsave(p_diff, file.path(fdir, "FIGURE_shapley_vs_lmdi_difference.png"),
              cfg$output$fig_width + 3, cfg$output$fig_height + 1, cfg$output$fig_dpi)
} else {
  message("Skipping Value-in-Health simulation figures because .NFZ_SKIP_PLOTS=TRUE or --no-plots was used.")
}

# ---- 3. Empirical Shapley vs LMDI-style comparison ----
method_comparison_path <- file.path(tdir, "method_comparison.csv")
if (file.exists(method_comparison_path)) {
  method_comparison <- readr::read_csv(method_comparison_path, show_col_types = FALSE)
  emp_comp <- make_empirical_shapley_vs_lmdi_table(method_comparison)
  readr::write_csv(emp_comp, file.path(tdir, "TABLE_empirical_shapley_vs_lmdi.csv"))
} else {
  warning("Missing method_comparison.csv; skipping empirical Shapley-vs-LMDI table.")
}

# ---- 4. Multi-year empirical illustration: baseline vs every later year ----
dat_path <- "output/rds/aggregated_data.rds"
fit_path <- "output/rds/model_fit.rds"
if (file.exists(dat_path) && file.exists(fit_path)) {
  dat <- readRDS(dat_path)
  fit <- readRDS(fit_path)
  baseline_year <- cfg$data$baseline_year

  empirical_multiyear <- run_empirical_multiyear_decomposition(
    dat = dat,
    fit = fit,
    baseline_year = baseline_year,
    methods = c("naive", "shapley")
  )

  readr::write_csv(empirical_multiyear, file.path(tdir, "TABLE_empirical_multiyear_lmdi_shapley.csv"))

  if (make_plots) {
    p_multi <- plot_empirical_multiyear_components(empirical_multiyear)
    safe_ggsave(p_multi, file.path(fdir, "FIGURE_empirical_multiyear_lmdi_shapley.png"),
                cfg$output$fig_width + 3, cfg$output$fig_height + 1, cfg$output$fig_dpi)
  } else {
    message("Skipping empirical multi-year figure because .NFZ_SKIP_PLOTS=TRUE or --no-plots was used.")
  }
} else {
  warning("Missing output/rds/aggregated_data.rds or output/rds/model_fit.rds; skipping multi-year empirical outputs.")
}

message("Value-in-Health-style outputs written:")
message("  ", file.path(tdir, "TABLE_methods_feature_comparison.csv"))
message("  ", file.path(tdir, "TABLE_method_decision_guide.csv"))
message("  ", file.path(tdir, "TABLE_simulation_performance_overall.csv"))
message("  ", file.path(tdir, "TABLE_shapley_vs_lmdi_simulation.csv"))
if (make_plots) {
  message("  ", file.path(fdir, "FIGURE_simulation_overall_rmse.png"))
  message("  ", file.path(fdir, "FIGURE_shapley_vs_lmdi_difference.png"))
}
message("Done.")
