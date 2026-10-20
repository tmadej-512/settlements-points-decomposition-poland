#!/usr/bin/env Rscript
# ============================================================
# scripts/run_05b_simulation_two_scenarios.R
# ------------------------------------------------------------
# Runs methodological validation scenarios:
#   0. oracle_decomposition    -> known value function, no model fit
#   1. clean_recovery          -> main-text model-based recovery evidence
#   2. heterogeneous_stress    -> supplementary stress test
#   3. cell_composition_shift  -> explicit cell-mix shift stress test
#   4. model_misspecification  -> deliberately unmodelled nonlinear stress test
#
# Backward-compatible outputs are also written:
#   simulation_raw.csv, simulation_summary.csv, simulation_recovery*.png
# use the clean_recovery scenario, so run_07_make_paper_outputs.R still works.
# ============================================================

source("R/99_load_all.R")
cfg <- load_config("config.R")

args <- commandArgs(trailingOnly = TRUE)
make_plots <- !(exists(".NFZ_SKIP_PLOTS") && isTRUE(.NFZ_SKIP_PLOTS)) && !("--no-plots" %in% args)
reps <- cfg$validation$simulation_reps
if (exists(".NFZ_REPS_OVERRIDE")) {
  reps <- .NFZ_REPS_OVERRIDE
} else if ("--reps" %in% args) {
  reps <- as.integer(args[which(args == "--reps") + 1])
}

sim_cache_file <- cfg$validation$simulation_two_scenarios_cache_file %||%
  "output/cache/simulation_two_scenarios_results.rds"
sim_use_cache <- cfg$validation$simulation_use_cache %||% TRUE
sim_force_rerun <- cfg$validation$simulation_force_rerun %||% FALSE
if ("--force" %in% args || "--force-simulation" %in% args ||
    (exists(".NFZ_FORCE_RERUN") && isTRUE(.NFZ_FORCE_RERUN)) ||
    (exists(".NFZ_SIMULATION_FORCE_RERUN") && isTRUE(.NFZ_SIMULATION_FORCE_RERUN))) {
  sim_force_rerun <- TRUE
}

required_methods <- c("naive", "ofat", "sequential", "shapley")
required_components <- c("N", "PI", "LAMBDA", "KAPPA")
cache_version <- "simulation-cache-v6-oracle-clean-stress-composition-misspec"

scenarios <- get_simulation_scenarios(cfg)
required_scenarios <- names(scenarios)

ensure_dir(cfg$output$tables_dir)
ensure_dir(cfg$output$figures_dir)
ensure_dir_for_file(sim_cache_file)

# Apply reps override after construction.
for (nm in names(scenarios)) scenarios[[nm]]$reps <- reps

cache_ok <- FALSE
if (should_use_cache(sim_use_cache, sim_force_rerun, sim_cache_file)) {
  message("Checking cached simulation results: ", sim_cache_file)
  cached <- tryCatch(cache_read(sim_cache_file), error = function(e) NULL)
  cache_ok <- !is.null(cached) &&
    cache_metadata_matches(cached, expected_reps = reps, expected_seed = cfg$validation$simulation_seed) &&
    identical(cached$metadata$cache_version %||% NA_character_, cache_version) &&
    simulation_has_required_methods(cached$raw, required_methods, required_components) &&
    all(required_scenarios %in% unique(cached$raw$scenario))

  if (cache_ok) {
    message("Using cached simulation results: ", sim_cache_file)
    sim_result <- cached
  } else {
    message("Simulation cache is stale/incomplete; recomputing.")
  }
}

if (!cache_ok) {
  message("Running multi-scenario simulation validation: ", reps, " replicates per scenario.")
  message("Scenario 0 isolates attribution with an oracle value function.")
  message("Scenario 1 checks full model-based recovery under controlled conditions.")
  message("Additional scenarios test cell-composition shift and model misspecification.\n")

  scenario_results <- list()

  for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]
    message("\n", strrep("=", 60))
    message("Simulation scenario: ", sc$label)
    message(strrep("=", 60))
    message(sc$description)

    if (identical(sc$type %||% "model_based", "oracle")) {
      one <- run_oracle_decomposition_study(
        baseline_year   = sc$baseline_year,
        comparison_year = sc$comparison_year,
        years           = sc$years,
        n_units         = sc$n_units,
        sex_levels      = sc$sex_levels,
        age_levels      = sc$age_levels,
        heterogeneous_slopes = sc$heterogeneous_slopes,
        regional_random_slopes = sc$regional_random_slopes,
        regional_slope_sd = sc$regional_slope_sd,
        true_growth     = sc$true_growth,
        reps            = sc$reps,
        seed            = sc$seed,
        methods         = required_methods,
        pandemic_shock  = sc$pandemic_shock,
        noise_scale     = sc$noise_scale
      )
    } else {
      one <- run_simulation_study(
        model_cfg       = cfg$model,
        baseline_year   = sc$baseline_year,
        comparison_year = sc$comparison_year,
        years           = sc$years,
        n_units         = sc$n_units,
        sex_levels      = sc$sex_levels,
        age_levels      = sc$age_levels,
        heterogeneous_slopes = sc$heterogeneous_slopes,
        regional_random_slopes = sc$regional_random_slopes,
        regional_slope_sd = sc$regional_slope_sd,
        true_growth     = sc$true_growth,
        reps            = sc$reps,
        seed            = sc$seed,
        methods         = required_methods,
        stop_on_fit_failure = TRUE,
        stop_on_method_failure = TRUE,
        pandemic_shock  = sc$pandemic_shock,
        noise_scale     = sc$noise_scale,
        composition_shift = sc$composition_shift %||% FALSE,
        composition_shift_strength = sc$composition_shift_strength %||% 0.25,
        misspecification_stress = sc$misspecification_stress %||% FALSE,
        misspecification_strength = sc$misspecification_strength %||% 0.25
      )
    }

    one$raw <- one$raw %>%
      dplyr::mutate(
        scenario = sc$scenario,
        scenario_label = sc$label,
        validation_level = sc$type %||% "model_based",
        .before = 1
      )

    one$summary <- one$summary %>%
      dplyr::mutate(
        scenario = sc$scenario,
        scenario_label = sc$label,
        validation_level = sc$type %||% "model_based",
        .before = 1
      )

    scenario_results[[sc_name]] <- one
  }

  sim_raw <- dplyr::bind_rows(lapply(scenario_results, `[[`, "raw"))
  metrics_by_component <- compute_simulation_recovery_metrics(sim_raw)
  metrics_overall <- compute_simulation_recovery_metrics_overall(metrics_by_component)

  sim_result <- list(
    raw = sim_raw,
    summary = metrics_by_component,
    overall = metrics_overall,
    audit = sim_raw %>%
      dplyr::count(scenario, method, component, name = "n_rows") %>%
      tidyr::complete(
        scenario = names(scenarios),
        method = required_methods,
        component = required_components,
        fill = list(n_rows = 0)
      ),
    metadata = list(
      step = "simulation_three_scenarios",
      reps = as.integer(reps),
      seed = as.integer(cfg$validation$simulation_seed),
      methods = required_methods,
      components = required_components,
      scenarios = names(scenarios),
      created_at = as.character(Sys.time()),
      cache_version = cache_version
    )
  )

  cache_write(sim_result, sim_cache_file)
  message("Saved simulation cache: ", sim_cache_file)
}

# Defensive audit.
missing_cells <- sim_result$audit %>% dplyr::filter(n_rows == 0)
if (nrow(missing_cells) > 0) {
  stop(
    "Simulation results are incomplete. Missing cells:\n",
    paste(paste(missing_cells$scenario, missing_cells$method, missing_cells$component, sep = " / "), collapse = "\n"),
    call. = FALSE
  )
}

# Full scenario outputs.
readr::write_csv(sim_result$raw, file.path(cfg$output$tables_dir, "simulation_three_scenarios_raw.csv"))
readr::write_csv(sim_result$summary, file.path(cfg$output$tables_dir, "simulation_three_scenarios_metrics_by_component.csv"))
readr::write_csv(sim_result$overall, file.path(cfg$output$tables_dir, "simulation_three_scenarios_metrics_overall.csv"))
readr::write_csv(sim_result$audit, file.path(cfg$output$tables_dir, "simulation_three_scenarios_method_audit.csv"))

# Backward-compatible names retained for older scripts/manuscript assembly.
readr::write_csv(sim_result$raw, file.path(cfg$output$tables_dir, "simulation_two_scenarios_raw.csv"))
readr::write_csv(sim_result$summary, file.path(cfg$output$tables_dir, "simulation_two_scenarios_metrics_by_component.csv"))
readr::write_csv(sim_result$overall, file.path(cfg$output$tables_dir, "simulation_two_scenarios_metrics_overall.csv"))
readr::write_csv(sim_result$audit, file.path(cfg$output$tables_dir, "simulation_two_scenarios_method_audit.csv"))

# Backward-compatible main-text model-based outputs: clean scenario only.
clean_raw <- sim_result$raw %>% dplyr::filter(.data$scenario == "clean_recovery")
clean_summary <- sim_result$summary %>% dplyr::filter(.data$scenario == "clean_recovery") %>%
  dplyr::select(-dplyr::any_of(c("scenario", "scenario_label", "validation_level")))

# Separate oracle outputs for the main methodological validation.
oracle_raw <- sim_result$raw %>% dplyr::filter(.data$scenario == "oracle_decomposition")
oracle_summary <- sim_result$summary %>% dplyr::filter(.data$scenario == "oracle_decomposition")

readr::write_csv(clean_summary, file.path(cfg$output$tables_dir, "simulation_summary.csv"))
readr::write_csv(clean_raw %>% dplyr::select(-dplyr::any_of(c("scenario", "scenario_label", "validation_level"))),
                 file.path(cfg$output$tables_dir, "simulation_raw.csv"))
readr::write_csv(oracle_summary, file.path(cfg$output$tables_dir, "simulation_oracle_summary.csv"))
readr::write_csv(oracle_raw, file.path(cfg$output$tables_dir, "simulation_oracle_raw.csv"))

message("\n== Oracle metrics, by method/component ==")
print(oracle_summary)
message("\n== Clean model-based recovery metrics, by method/component ==")
print(clean_summary)
message("\n== Overall metrics, all validation scenarios ==")
print(sim_result$overall)

if (make_plots) {
  # Plots: one recovery-vs-truth and one error plot for every configured scenario.
  sim_plot_limits <- compute_simulation_plot_limits(sim_result$raw)
  scenario_plot_objects <- list()
  for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]
    p_truth <- plot_simulation_recovery_vs_truth_by_scenario(sim_result$raw, sc$scenario, y_limits = sim_plot_limits$truth)
    p_error <- plot_simulation_error_by_scenario(sim_result$raw, sc$scenario, y_limits = sim_plot_limits$error)

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

  # Backward-compatible figure names for downstream assembly: clean scenario.
  p_clean_truth <- scenario_plot_objects[["clean_recovery_truth"]]
  p_clean_error <- scenario_plot_objects[["clean_recovery_error"]]
  safe_ggsave(p_clean_error, file.path(cfg$output$figures_dir, "simulation_recovery.png"),
              cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)
  safe_ggsave(p_clean_truth, file.path(cfg$output$figures_dir, "simulation_recovery_vs_truth.png"),
              cfg$output$fig_width + 2, cfg$output$fig_height, cfg$output$fig_dpi)

  # DGP diagnostics: write separate folders for each scenario.
  if (isTRUE(cfg$validation$simulation_make_dgp_plots %||% TRUE)) {
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
  }
} else {
  message("Skipping all simulation plots because .NFZ_SKIP_PLOTS=TRUE or --no-plots was used.")
}

message("\nMulti-scenario simulation validation complete.")
message("Main methodological validation outputs:")
message("  ", file.path(cfg$output$tables_dir, "simulation_oracle_summary.csv"))
message("  ", file.path(cfg$output$figures_dir, "simulation_recovery_vs_truth_oracle_decomposition.png"))
message("Main model-based recovery outputs:")
message("  ", file.path(cfg$output$tables_dir, "simulation_summary.csv"))
message("  ", file.path(cfg$output$figures_dir, "simulation_recovery_vs_truth.png"))
message("Full multi-scenario outputs:")
message("  ", file.path(cfg$output$tables_dir, "simulation_three_scenarios_metrics_by_component.csv"))
message("  ", file.path(cfg$output$tables_dir, "simulation_three_scenarios_metrics_overall.csv"))
