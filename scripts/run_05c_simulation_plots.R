#!/usr/bin/env Rscript
# ============================================================
# scripts/run_05c_simulation_plots.R
# ------------------------------------------------------------
# Plot-only script for simulation outputs.
#
# It reads existing simulation CSVs from output/tables/ and regenerates
# simulation figures in output/figures/. It does NOT refit models and
# does NOT rerun simulation replicates.
#
# Run after:
#   source("scripts/run_05b_simulation_generate.R")
# or after any full run that has created:
#   output/tables/simulation_three_scenarios_raw.csv
# ============================================================

source("R/99_load_all.R")
cfg <- load_config("config.R")
args <- commandArgs(trailingOnly = TRUE)
make_dgp <- isTRUE(cfg$validation$simulation_make_dgp_plots %||% TRUE) && !("--no-dgp" %in% args)

sim_path <- file.path(cfg$output$tables_dir, "simulation_three_scenarios_raw.csv")
if (!file.exists(sim_path)) {
  stop("Missing ", sim_path, ". Run scripts/run_05b_simulation_generate.R first.", call. = FALSE)
}

ensure_dir(cfg$output$figures_dir)
sim_raw <- readr::read_csv(sim_path, show_col_types = FALSE)
scenarios <- get_simulation_scenarios(cfg)

# Use scenarios present in the raw data, but preserve configured order.
available <- unique(sim_raw$scenario)
scenarios <- scenarios[names(scenarios) %in% available]
if (length(scenarios) == 0) {
  stop("No configured scenarios found in ", sim_path, call. = FALSE)
}

message("Rendering simulation plots from existing CSV: ", sim_path)

sim_plot_limits <- compute_simulation_plot_limits(sim_raw)
scenario_plot_objects <- list()
for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  message("  scenario: ", sc$scenario)

  p_truth <- plot_simulation_recovery_vs_truth_by_scenario(sim_raw, sc$scenario, y_limits = sim_plot_limits$truth)
  p_error <- plot_simulation_error_by_scenario(sim_raw, sc$scenario, y_limits = sim_plot_limits$error)

  scenario_plot_objects[[paste0(sc$scenario, "_truth")]] <- p_truth
  scenario_plot_objects[[paste0(sc$scenario, "_error")]] <- p_error

  safe_ggsave(
    p_truth,
    file.path(cfg$output$figures_dir, paste0("simulation_recovery_vs_truth_", sc$scenario, ".png")),
    cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi
  )
  safe_ggsave(
    p_error,
    file.path(cfg$output$figures_dir, paste0("simulation_error_", sc$scenario, ".png")),
    cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi
  )
}

# Backward-compatible main-text figure names for clean_recovery.
if ("clean_recovery" %in% available) {
  p_clean_truth <- scenario_plot_objects[["clean_recovery_truth"]]
  p_clean_error <- scenario_plot_objects[["clean_recovery_error"]]
  safe_ggsave(p_clean_error, file.path(cfg$output$figures_dir, "simulation_recovery.png"),
              cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
  safe_ggsave(p_clean_truth, file.path(cfg$output$figures_dir, "simulation_recovery_vs_truth.png"),
              cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
}

if (make_dgp) {
  message("Rendering DGP diagnostic plots. Use --no-dgp to skip them.")
  for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]
    plot_simulation_dgp_diagnostics(
      seed = cfg$validation$simulation_dgp_seed %||% sc$seed,
      out_dir = file.path(cfg$output$figures_dir, "simulation_dgp", sc$scenario),
      years = sc$years,
      n_units = sc$n_units,
      sex_levels = sc$sex_levels,
      age_levels = sc$age_levels,
      true_growth = sc$true_growth,
      heterogeneous_slopes = sc$heterogeneous_slopes,
      regional_random_slopes = sc$regional_random_slopes,
      regional_slope_sd = sc$regional_slope_sd,
      width = cfg$output$fig_width + 1,
      height = cfg$output$fig_height,
      dpi = cfg$output$fig_dpi
    )
  }
} else {
  message("Skipping DGP diagnostic plots.")
}

message("Simulation plots regenerated in: ", cfg$output$figures_dir)
