#!/usr/bin/env Rscript
# ============================================================
# scripts/run_05_simulation_study.R
# ------------------------------------------------------------
# Simulation study with RDS cache.
#
# Expensive model-fitting results are saved to
# cfg$validation$simulation_cache_file. Re-running this script will
# reuse the cache when simulation_use_cache = TRUE and
# simulation_force_rerun = FALSE, provided the cache matches the
# requested number of reps/seed and contains all four methods.
#
# This lets you change plots/tables without recalculating 200-500
# simulation replicates.
# ============================================================

source("R/99_load_all.R")

cfg <- load_config("config.R")

# Reps can be overridden two ways:
#   .NFZ_REPS_OVERRIDE <- 20; source("scripts/run_05_simulation_study.R")
#   Rscript scripts/run_05_simulation_study.R --reps 20
reps <- cfg$validation$simulation_reps
args <- commandArgs(trailingOnly = TRUE)
if (exists(".NFZ_REPS_OVERRIDE")) {
  reps <- .NFZ_REPS_OVERRIDE
} else if ("--reps" %in% args) {
  reps <- as.integer(args[which(args == "--reps") + 1])
}

sim_cache_file <- cfg$validation$simulation_cache_file %||% "output/cache/simulation_results.rds"
sim_use_cache <- cfg$validation$simulation_use_cache %||% TRUE
sim_force_rerun <- cfg$validation$simulation_force_rerun %||% FALSE
if ("--force" %in% args || "--force-simulation" %in% args ||
    (exists(".NFZ_FORCE_RERUN") && isTRUE(.NFZ_FORCE_RERUN)) ||
    (exists(".NFZ_SIMULATION_FORCE_RERUN") && isTRUE(.NFZ_SIMULATION_FORCE_RERUN))) {
  sim_force_rerun <- TRUE
}

required_methods <- c("naive", "ofat", "sequential", "shapley")
required_components <- c("N", "PI", "LAMBDA", "KAPPA")

sim_years <- cfg$validation$simulation_years %||% (cfg$data$baseline_year + 0:5)
sim_comparison_year <- cfg$validation$simulation_comparison_year %||% max(sim_years)
sim_n_units <- cfg$validation$simulation_n_units %||% 30
sim_sex_levels <- cfg$validation$simulation_sex_levels %||% c("F", "M")
sim_age_levels <- cfg$validation$simulation_age_levels %||% c("0-17", "18-44", "45-64", "65+")
sim_heterogeneous_slopes <- cfg$validation$simulation_heterogeneous_slopes %||% TRUE
sim_regional_random_slopes <- cfg$validation$simulation_regional_random_slopes %||% TRUE
sim_regional_slope_sd <- cfg$validation$simulation_regional_slope_sd %||% list(N = 0.003, PI = 0.008, LAMBDA = 0.008, KAPPA = 0.006)
sim_true_growth <- cfg$validation$simulation_true_growth %||% list(N = 0.01, PI = 0.03, LAMBDA = -0.01, KAPPA = 0.02)
sim_cache_version <- "simulation-cache-v3-regional-random-slopes"

ensure_dir(cfg$output$tables_dir)
ensure_dir(cfg$output$figures_dir)
ensure_dir_for_file(sim_cache_file)

cache_ok <- FALSE
if (should_use_cache(sim_use_cache, sim_force_rerun, sim_cache_file)) {
  message("Checking cached simulation results: ", sim_cache_file)
  cached <- tryCatch(cache_read(sim_cache_file), error = function(e) NULL)
  cache_ok <- !is.null(cached) &&
    cache_metadata_matches(cached, expected_reps = reps, expected_seed = cfg$validation$simulation_seed) &&
    identical(cached$metadata$cache_version %||% NA_character_, sim_cache_version) &&
    identical(cached$metadata$n_units %||% NA_integer_, as.integer(sim_n_units)) &&
    identical(cached$metadata$heterogeneous_slopes %||% NA, sim_heterogeneous_slopes) &&
    identical(cached$metadata$regional_random_slopes %||% NA, sim_regional_random_slopes) &&
    identical(cached$metadata$regional_slope_sd %||% list(), sim_regional_slope_sd) &&
    identical(cached$metadata$years %||% NA_integer_, as.integer(sim_years)) &&
    simulation_has_required_methods(cached$raw, required_methods, required_components)

  if (cache_ok) {
    message("Using cached simulation results: ", sim_cache_file)
    sim_result <- cached
  } else {
    message("Simulation cache exists but is stale/incomplete or uses different reps/seed; recomputing.")
  }
}

if (!cache_ok) {
  message("Running simulation study: ", reps, " replicates ",
          "(each refits the 3-stage GAM/GLM chain)...")
  message("This validates the FULL pipeline (model fit + counterfactual ",
          "prediction + attribution), not just the attribution rule.\n")

  sim_result <- run_simulation_study(
    model_cfg       = cfg$model,
    baseline_year   = cfg$data$baseline_year,
    comparison_year = sim_comparison_year,
    years           = sim_years,
    n_units         = sim_n_units,
    sex_levels      = sim_sex_levels,
    age_levels      = sim_age_levels,
    heterogeneous_slopes = sim_heterogeneous_slopes,
    regional_random_slopes = sim_regional_random_slopes,
    regional_slope_sd = sim_regional_slope_sd,
    true_growth     = sim_true_growth,
    reps            = reps,
    seed            = cfg$validation$simulation_seed,
    methods         = required_methods,
    stop_on_fit_failure = TRUE,
    stop_on_method_failure = TRUE
  )

  sim_result$audit <- sim_result$raw %>%
    dplyr::count(method, component, name = "n_rows") %>%
    tidyr::complete(
      method = required_methods,
      component = required_components,
      fill = list(n_rows = 0)
    )

  sim_result$metadata <- list(
    step = "simulation",
    reps = as.integer(reps),
    seed = as.integer(cfg$validation$simulation_seed),
    methods = required_methods,
    components = required_components,
    created_at = as.character(Sys.time()),
    cache_version = sim_cache_version,
    n_units = as.integer(sim_n_units),
    years = as.integer(sim_years),
    comparison_year = as.integer(sim_comparison_year),
    age_levels = sim_age_levels,
    sex_levels = sim_sex_levels,
    heterogeneous_slopes = sim_heterogeneous_slopes,
    regional_random_slopes = sim_regional_random_slopes,
    regional_slope_sd = sim_regional_slope_sd,
    true_growth = sim_true_growth
  )

  cache_write(sim_result, sim_cache_file)
  message("Saved simulation cache: ", sim_cache_file)
}

# Defensive audit: do not silently create a naive-only plot.
simulation_method_audit <- sim_result$raw %>%
  dplyr::count(method, component, name = "n_rows") %>%
  tidyr::complete(
    method = required_methods,
    component = required_components,
    fill = list(n_rows = 0)
  )
missing_cells <- simulation_method_audit %>% dplyr::filter(n_rows == 0)
if (nrow(missing_cells) > 0) {
  stop(
    "Simulation results are incomplete. Missing method/component cells:\n",
    paste(paste(missing_cells$method, missing_cells$component, sep = " / "), collapse = "\n"),
    "\nSet cfg$validation$simulation_force_rerun <- TRUE after fixing the simulation code.",
    call. = FALSE
  )
}

message("\n== Simulation summary: bias and RMSE vs. known ground truth ==")
print(sim_result$summary)

readr::write_csv(sim_result$summary, file.path(cfg$output$tables_dir, "simulation_summary.csv"))
readr::write_csv(sim_result$raw, file.path(cfg$output$tables_dir, "simulation_raw.csv"))
readr::write_csv(simulation_method_audit, file.path(cfg$output$tables_dir, "simulation_method_audit.csv"))

p <- plot_simulation_recovery(sim_result$raw)
safe_ggsave(p, file.path(cfg$output$figures_dir, "simulation_recovery.png"),
            cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)

p_truth <- plot_simulation_recovery_vs_truth(sim_result$raw)
safe_ggsave(p_truth, file.path(cfg$output$figures_dir, "simulation_recovery_vs_truth.png"),
            cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)

if (isTRUE(cfg$validation$simulation_make_dgp_plots %||% TRUE)) {
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
}

message("\nSaved: ", file.path(cfg$output$tables_dir, "simulation_summary.csv"))
message("Saved: ", file.path(cfg$output$tables_dir, "simulation_raw.csv"))
message("Saved: ", file.path(cfg$output$tables_dir, "simulation_method_audit.csv"))
message("Saved: ", file.path(cfg$output$figures_dir, "simulation_recovery.png"))
message("Saved: ", file.path(cfg$output$figures_dir, "simulation_recovery_vs_truth.png"))
if (isTRUE(cfg$validation$simulation_make_dgp_plots %||% TRUE)) {
  message("Saved DGP diagnostics to: ", file.path(cfg$output$figures_dir, "simulation_dgp"))
}
message("Done.")
