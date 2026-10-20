# ============================================================
# R/99_load_all.R
# ------------------------------------------------------------
# Sources every file in R/ in the correct dependency order.
# Each scripts/run_*.R file starts by calling source("R/99_load_all.R")
# rather than sourcing individual files, so the load order is defined
# in exactly one place.
# ============================================================

.r_files_in_order <- c(
  "R/00_utils.R",
  "R/00_safe_fwrite.R",
  "R/01_data_schema.R",
  "R/02_models.R",
  "R/03_predict_helpers.R",
  "R/04_decomposition_naive.R",
  "R/05_decomposition_ofat.R",
  "R/06_decomposition_sequential.R",
  "R/07_decomposition_shapley.R",
  "R/08_validation_axioms.R",
  "R/09_validation_bootstrap.R",
  "R/10_simulation_ground_truth.R",
  "R/10e_oracle_decomposition.R",
  "R/10d_simulation_recovery_plots.R",
  "R/10c_simulation_recovery_metrics.R",
  "R/10b_simulation_scenarios.R",
  "R/11_order_sensitivity.R",
  "R/12_plots.R",
  "R/13_distribution_diagnostics.R",
  "R/14_value_in_health_outputs.R",
  "R/15_model_residual_diagnostics.R"
)

for (f in .r_files_in_order) {
  if (!file.exists(f)) {
    stop("Expected source file not found: ", f,
         "\nRun this from the project root (the directory containing config.R).")
  }
  source(f)
}

rm(.r_files_in_order, f)
