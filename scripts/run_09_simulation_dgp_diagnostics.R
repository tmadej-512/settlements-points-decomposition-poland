#!/usr/bin/env Rscript
# ============================================================
# scripts/run_09_simulation_dgp_diagnostics.R
# ------------------------------------------------------------
# Regenerates synthetic data-generating-process diagnostic figures
# without re-running the full simulation study.
# ============================================================

source("R/99_load_all.R")
cfg <- load_config("config.R")

sim_years <- cfg$validation$simulation_years %||% (cfg$data$baseline_year + 0:5)
sim_n_units <- cfg$validation$simulation_n_units %||% 30
sim_sex_levels <- cfg$validation$simulation_sex_levels %||% c("F", "M")
sim_age_levels <- cfg$validation$simulation_age_levels %||% c("0-17", "18-44", "45-64", "65+")
sim_heterogeneous_slopes <- cfg$validation$simulation_heterogeneous_slopes %||% TRUE
sim_regional_random_slopes <- cfg$validation$simulation_regional_random_slopes %||% TRUE
sim_regional_slope_sd <- cfg$validation$simulation_regional_slope_sd %||% list(N = 0.003, PI = 0.008, LAMBDA = 0.008, KAPPA = 0.006)
sim_true_growth <- cfg$validation$simulation_true_growth %||% list(N = 0.01, PI = 0.03, LAMBDA = -0.01, KAPPA = 0.02)

plot_simulation_dgp_diagnostics(
  seed = cfg$validation$simulation_dgp_seed %||% cfg$validation$simulation_seed,
  out_dir = file.path(cfg$output$figures_dir, "simulation_dgp"),
  years = sim_years,
  n_units = sim_n_units,
  sex_levels = sim_sex_levels,
  age_levels = sim_age_levels,
  true_growth = sim_true_growth,
  heterogeneous_slopes = sim_heterogeneous_slopes,
  regional_random_slopes = sim_regional_random_slopes,
  regional_slope_sd = sim_regional_slope_sd,
  width = cfg$output$fig_width + 1,
  height = cfg$output$fig_height,
  dpi = cfg$output$fig_dpi
)

message("Saved DGP diagnostics to: ", file.path(cfg$output$figures_dir, "simulation_dgp"))
