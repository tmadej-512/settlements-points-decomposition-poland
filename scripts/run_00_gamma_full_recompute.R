#!/usr/bin/env Rscript
# ============================================================
# scripts/run_00_gamma_full_recompute.R
# ------------------------------------------------------------
# Full recomputation for the diagnostic-updated main specification:
# value-stage Gamma log-link model with offset log(n_provision).
#
# This script does NOT require editing config.R. It sets a global override
# for this R session and forces recomputation of cached bootstrap/simulation.
# ============================================================

if (!file.exists("config.R")) {
  stop("config.R not found. Run from the project root.", call. = FALSE)
}

.NFZ_VALUE_FAMILY_OVERRIDE <- "gamma"
.NFZ_FORCE_RERUN <- TRUE
if (exists(".NFZ_QUICK", envir = .GlobalEnv, inherits = FALSE)) {
  rm(.NFZ_QUICK, envir = .GlobalEnv)
}

message("Running full recomputation with .NFZ_VALUE_FAMILY_OVERRIDE = 'gamma'.")
message("Forcing simulation/bootstrap cache refresh via .NFZ_FORCE_RERUN = TRUE.")

source("scripts/run_00_compute_all.R")
source("scripts/run_01b_residual_diagnostics.R")
source("scripts/run_01c_quick_value_family_test.R")
source("scripts/run_99_make_all_plots.R")
source("scripts/run_97_make_diagnostics_validation_plots.R")
source("scripts/run_98_export_vih_figures.R")
source("scripts/run_98b_make_figure_legends.R")

message("Gamma main-specification recomputation complete.")
message("Review output/tables/model_quality.csv and output/tables/quick_value_family_test.csv.")
message("Review submission/figures and submission/Figure_Legends_VIH.*.")
