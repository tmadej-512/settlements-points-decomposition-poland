# ============================================================
# R/10b_simulation_scenarios.R
# ------------------------------------------------------------
# Five validation scenarios:
#   0) oracle_decomposition: known value function, no model fit.
#   1) clean_recovery: controlled model-based recovery check.
#   2) heterogeneous_stress: harder DGP for sensitivity/supplement.
#   3) cell_composition_shift: explicit cell-mix shift correlated with high-value cells.
#   4) model_misspecification: deliberately unmodelled nonlinear/interaction stress.
# ============================================================

get_simulation_scenarios <- function(cfg) {
  base_reps <- cfg$validation$simulation_reps %||% 200
  if (exists(".NFZ_REPS_OVERRIDE", envir = .GlobalEnv)) {
    base_reps <- get(".NFZ_REPS_OVERRIDE", envir = .GlobalEnv)
  }
  base_seed <- cfg$validation$simulation_seed %||% 1
  base_years <- cfg$validation$simulation_years %||% (cfg$data$baseline_year + 0:6)
  base_n_units <- cfg$validation$simulation_n_units %||% 80
  base_sex <- cfg$validation$simulation_sex_levels %||% c("F", "M")
  base_age <- cfg$validation$simulation_age_levels %||% c("0-17", "18-44", "45-64", "65+")
  base_growth <- cfg$validation$simulation_true_growth %||% list(N = 0.005, PI = 0.025, LAMBDA = 0.035, KAPPA = -0.005)
  comparison_year <- cfg$validation$simulation_comparison_year %||% max(base_years)

  list(
    oracle_decomposition = list(
      scenario = "oracle_decomposition",
      label = "Oracle decomposition",
      description = paste(
        "Oracle validation: uses the known deterministic DGP value function directly;",
        "no noisy observations and no GAM/GLM refit. This isolates the attribution rule."
      ),
      type = "oracle",
      reps = base_reps,
      seed = base_seed + 500000,
      years = base_years,
      baseline_year = cfg$data$baseline_year,
      comparison_year = comparison_year,
      n_units = base_n_units,
      sex_levels = base_sex,
      age_levels = base_age,
      true_growth = base_growth,
      heterogeneous_slopes = TRUE,
      regional_random_slopes = FALSE,
      regional_slope_sd = list(N = 0, PI = 0, LAMBDA = 0, KAPPA = 0),
      pandemic_shock = FALSE,
      noise_scale = 1.00
    ),

    clean_recovery = list(
      scenario = "clean_recovery",
      label = "Clean model-based recovery",
      description = paste(
        "Controlled model-based recovery scenario: moderate noise, age-sex-specific slopes,",
        "no regional random slopes and no pandemic shock."
      ),
      type = "model_based",
      reps = base_reps,
      seed = base_seed,
      years = base_years,
      baseline_year = cfg$data$baseline_year,
      comparison_year = comparison_year,
      n_units = base_n_units,
      sex_levels = base_sex,
      age_levels = base_age,
      true_growth = base_growth,
      heterogeneous_slopes = TRUE,
      regional_random_slopes = FALSE,
      regional_slope_sd = list(N = 0, PI = 0, LAMBDA = 0, KAPPA = 0),
      pandemic_shock = FALSE,
      noise_scale = 0.65
    ),

    heterogeneous_stress = list(
      scenario = "heterogeneous_stress",
      label = "Heterogeneous stress test",
      description = paste(
        "Stress-test scenario: regional random slopes, stronger noise,",
        "and a temporary pandemic-like service disruption."
      ),
      type = "model_based",
      reps = base_reps,
      seed = base_seed + 100000,
      years = base_years,
      baseline_year = cfg$data$baseline_year,
      comparison_year = comparison_year,
      n_units = base_n_units,
      sex_levels = base_sex,
      age_levels = base_age,
      true_growth = base_growth,
      heterogeneous_slopes = TRUE,
      regional_random_slopes = TRUE,
      regional_slope_sd = cfg$validation$simulation_regional_slope_sd %||% list(N = 0.003, PI = 0.008, LAMBDA = 0.008, KAPPA = 0.006),
      pandemic_shock = TRUE,
      noise_scale = 1.35,
      composition_shift = FALSE,
      misspecification_stress = FALSE
    ),

    cell_composition_shift = list(
      scenario = "cell_composition_shift",
      label = "Cell-composition shift",
      description = paste(
        "Explicit cell-mix stress test: population weight shifts toward older/high-use cells,",
        "so aggregate LMDI-style allocation can remain close on totals but lose cell-level attribution."
      ),
      type = "model_based",
      reps = base_reps,
      seed = base_seed + 200000,
      years = base_years,
      baseline_year = cfg$data$baseline_year,
      comparison_year = comparison_year,
      n_units = base_n_units,
      sex_levels = base_sex,
      age_levels = base_age,
      true_growth = base_growth,
      heterogeneous_slopes = TRUE,
      regional_random_slopes = TRUE,
      regional_slope_sd = cfg$validation$simulation_regional_slope_sd %||% list(N = 0.004, PI = 0.010, LAMBDA = 0.010, KAPPA = 0.008),
      pandemic_shock = FALSE,
      noise_scale = 1.10,
      composition_shift = TRUE,
      composition_shift_strength = cfg$validation$simulation_composition_shift_strength %||% 0.32,
      misspecification_stress = FALSE
    ),

    model_misspecification = list(
      scenario = "model_misspecification",
      label = "Model-misspecification stress",
      description = paste(
        "Stress test with a nonlinear age-by-sex-by-time shock deliberately not aligned",
        "with the simulation-safe fitting formula."
      ),
      type = "model_based",
      reps = base_reps,
      seed = base_seed + 300000,
      years = base_years,
      baseline_year = cfg$data$baseline_year,
      comparison_year = comparison_year,
      n_units = base_n_units,
      sex_levels = base_sex,
      age_levels = base_age,
      true_growth = base_growth,
      heterogeneous_slopes = TRUE,
      regional_random_slopes = TRUE,
      regional_slope_sd = cfg$validation$simulation_regional_slope_sd %||% list(N = 0.003, PI = 0.010, LAMBDA = 0.010, KAPPA = 0.010),
      pandemic_shock = TRUE,
      noise_scale = 1.45,
      composition_shift = TRUE,
      composition_shift_strength = cfg$validation$simulation_composition_shift_strength %||% 0.22,
      misspecification_stress = TRUE,
      misspecification_strength = cfg$validation$simulation_misspecification_strength %||% 0.30
    )
  )
}
