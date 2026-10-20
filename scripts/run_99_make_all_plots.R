#!/usr/bin/env Rscript
# ============================================================
# scripts/run_99_make_all_plots.R
# ------------------------------------------------------------
# Regenerates all figures from existing tables/RDS outputs.
#
# This script does NOT rerun expensive simulations or bootstrap.
# It is meant for iterative figure polishing.
#
# Typical workflow:
#   source("scripts/run_00_compute_all.R")       # once, expensive
#   source("scripts/run_99_make_all_plots.R")    # many times, cheap
#
# Optional:
#   Rscript scripts/run_99_make_all_plots.R --no-dgp
# ============================================================

source("R/99_load_all.R")
cfg <- load_config("config.R")
args <- commandArgs(trailingOnly = TRUE)
make_dgp <- !("--no-dgp" %in% args)

tdir <- cfg$output$tables_dir
fdir <- cfg$output$figures_dir
ensure_dir(fdir)

need_file <- function(path, hint) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path, "\n", hint, call. = FALSE)
  }
  invisible(path)
}

message("Regenerating all figures from existing outputs.")

# ------------------------------------------------------------
# 1. Main method-comparison figures
# ------------------------------------------------------------
method_comparison_path <- file.path(tdir, "method_comparison.csv")
need_file(method_comparison_path, "Run scripts/run_02_compare_methods.R first, or scripts/run_00_compute_all.R.")
method_comparison <- readr::read_csv(method_comparison_path, show_col_types = FALSE)

p_method <- plot_method_comparison_bars(method_comparison)
p_resid <- plot_residual_comparison(method_comparison)
safe_ggsave(p_method, file.path(fdir, "method_comparison_bars.png"),
            cfg$output$fig_width, cfg$output$fig_height, cfg$output$fig_dpi)
safe_ggsave(p_resid, file.path(fdir, "residual_comparison.png"),
            cfg$output$fig_width, cfg$output$fig_height, cfg$output$fig_dpi)

# ------------------------------------------------------------
# 2. Order sensitivity
# ------------------------------------------------------------
order_path <- file.path(tdir, "order_sensitivity_summary.csv")
if (file.exists(order_path)) {
  order_sensitivity <- readr::read_csv(order_path, show_col_types = FALSE)
  p_order <- plot_order_sensitivity(order_sensitivity)
  safe_ggsave(p_order, file.path(fdir, "order_sensitivity.png"),
              cfg$output$fig_width, cfg$output$fig_height, cfg$output$fig_dpi)
} else {
  warning("Skipping order_sensitivity.png; missing ", order_path)
}

# ------------------------------------------------------------
# 3. Bootstrap intervals
# ------------------------------------------------------------
boot_path <- file.path(tdir, "bootstrap_summary.csv")
if (file.exists(boot_path)) {
  boot_summary <- readr::read_csv(boot_path, show_col_types = FALSE)
  p_boot <- plot_bootstrap_intervals(boot_summary, method_comparison)
  safe_ggsave(p_boot, file.path(fdir, "bootstrap_intervals.png"),
              cfg$output$fig_width, cfg$output$fig_height, cfg$output$fig_dpi)
} else {
  warning("Skipping bootstrap_intervals.png; missing ", boot_path)
}

# ------------------------------------------------------------
# 4. Simulation plots from existing simulation CSV
# ------------------------------------------------------------
sim_path <- file.path(tdir, "simulation_three_scenarios_raw.csv")
if (file.exists(sim_path)) {
  sim_raw <- readr::read_csv(sim_path, show_col_types = FALSE)
  scenarios <- get_simulation_scenarios(cfg)
  available <- unique(sim_raw$scenario)
  scenarios <- scenarios[names(scenarios) %in% available]

  sim_plot_limits <- compute_simulation_plot_limits(sim_raw)
  scenario_plot_objects <- list()
  for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]
    p_truth <- plot_simulation_recovery_vs_truth_by_scenario(sim_raw, sc$scenario, y_limits = sim_plot_limits$truth)
    p_error <- plot_simulation_error_by_scenario(sim_raw, sc$scenario, y_limits = sim_plot_limits$error)

    scenario_plot_objects[[paste0(sc$scenario, "_truth")]] <- p_truth
    scenario_plot_objects[[paste0(sc$scenario, "_error")]] <- p_error

    safe_ggsave(p_truth, file.path(fdir, paste0("simulation_recovery_vs_truth_", sc$scenario, ".png")),
                cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
    safe_ggsave(p_error, file.path(fdir, paste0("simulation_error_", sc$scenario, ".png")),
                cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
  }

  if ("clean_recovery" %in% available) {
    safe_ggsave(scenario_plot_objects[["clean_recovery_error"]], file.path(fdir, "simulation_recovery.png"),
                cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
    safe_ggsave(scenario_plot_objects[["clean_recovery_truth"]], file.path(fdir, "simulation_recovery_vs_truth.png"),
                cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
  }

  if (make_dgp && isTRUE(cfg$validation$simulation_make_dgp_plots %||% TRUE)) {
    for (sc_name in names(scenarios)) {
      sc <- scenarios[[sc_name]]
      plot_simulation_dgp_diagnostics(
        seed = cfg$validation$simulation_dgp_seed %||% sc$seed,
        out_dir = file.path(fdir, "simulation_dgp", sc$scenario),
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
  }
} else {
  warning("Skipping simulation plots; missing ", sim_path)
}

# ------------------------------------------------------------
# 5. Combined validation figure
# ------------------------------------------------------------
if (file.exists(order_path) && file.exists(boot_path) &&
    file.exists(file.path(tdir, "simulation_raw.csv"))) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    suppressPackageStartupMessages(library(patchwork))
    order_sensitivity <- readr::read_csv(order_path, show_col_types = FALSE)
    boot_summary <- readr::read_csv(boot_path, show_col_types = FALSE)
    sim_raw_clean <- readr::read_csv(file.path(tdir, "simulation_raw.csv"), show_col_types = FALSE)

    combined <- (plot_method_comparison_bars(method_comparison) + ggplot2::labs(title = "A. Method comparison") +
                 plot_order_sensitivity(order_sensitivity) + ggplot2::labs(title = "B. Order sensitivity")) /
                (plot_simulation_recovery_vs_truth(sim_raw_clean) + ggplot2::labs(title = "C. Simulation recovery") +
                 plot_bootstrap_intervals(boot_summary, method_comparison) + ggplot2::labs(title = "D. Bootstrap intervals"))

    ggplot2::ggsave(file.path(fdir, "FIGURE_validation_combined.png"),
                    combined,
                    width = cfg$output$fig_width + 6,
                    height = (cfg$output$fig_height + 1) * 2,
                    dpi = cfg$output$fig_dpi)
  } else {
    warning("Skipping FIGURE_validation_combined.png; package {patchwork} is not installed.")
  }
} else {
  warning("Skipping combined validation figure; missing one or more required tables.")
}

# ------------------------------------------------------------
# 6. Value-in-Health figures from existing tables
# ------------------------------------------------------------
if (file.exists(sim_path)) {
  sim_raw <- readr::read_csv(sim_path, show_col_types = FALSE)
  if (!"scenario_label" %in% names(sim_raw)) {
    sim_raw <- sim_raw %>% dplyr::mutate(scenario_label = .data$scenario)
  }
  sim_perf <- make_simulation_performance_table(sim_raw)
  p_rmse <- plot_simulation_overall_rmse(sim_perf)
  p_diff <- plot_shapley_vs_lmdi_difference(sim_raw)
  safe_ggsave(p_rmse, file.path(fdir, "FIGURE_simulation_overall_rmse.png"),
              cfg$output$fig_width + 3, cfg$output$fig_height + 1, cfg$output$fig_dpi)
  safe_ggsave(p_diff, file.path(fdir, "FIGURE_shapley_vs_lmdi_difference.png"),
              cfg$output$fig_width + 3, cfg$output$fig_height + 1, cfg$output$fig_dpi)
}

multi_path <- file.path(tdir, "TABLE_empirical_multiyear_lmdi_shapley.csv")
if (file.exists(multi_path)) {
  empirical_multiyear <- readr::read_csv(multi_path, show_col_types = FALSE)
  p_multi <- plot_empirical_multiyear_components(empirical_multiyear)
  safe_ggsave(p_multi, file.path(fdir, "FIGURE_empirical_multiyear_lmdi_shapley.png"),
              cfg$output$fig_width + 3, cfg$output$fig_height + 1, cfg$output$fig_dpi)
} else {
  warning("Skipping FIGURE_empirical_multiyear_lmdi_shapley.png; missing ", multi_path)
}

message("All available figures regenerated in: ", fdir)
