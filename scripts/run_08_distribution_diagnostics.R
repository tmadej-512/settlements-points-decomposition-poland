#!/usr/bin/env Rscript
# ============================================================
# scripts/run_08_distribution_diagnostics.R
# ------------------------------------------------------------
# Empirically checks whether the families assumed in R/02_models.R
# (negative binomial for patients/provision, Tweedie for value) are
# well-supported by the actual disease-group data, BEFORE trusting
# the fitted model chain's counterfactual predictions.
#
# This is diagnostic only: it does not change the model fit. If a
# stage's recommendation differs from the current default, update
# build_formula()'s family argument in R/02_models.R accordingly,
# then re-run scripts/run_01_fit_model.R onward.
#
# Run from the project root, ideally BEFORE run_01_fit_model.R (no
# dependency on its output -- this script reads directly from the
# configured data path):
#   Rscript scripts/run_08_distribution_diagnostics.R
# ============================================================

source("R/99_load_all.R")

cfg <- load_config("config.R")
ensure_dir(cfg$output$tables_dir)

message("Loading and aggregating data for distribution diagnostics...")
df_raw <- read_input_data(cfg$data$path)
dat0   <- normalize_schema(df_raw, value_var = cfg$data$value_var)
dat    <- filter_and_aggregate(dat0, disease_filter = cfg$data$disease_filter)
disease_label <- attr(dat, "disease_filter") %||% cfg$data$disease_filter %||% unique(dat$icd10)[1]

message("Disease group: '", disease_label, "', ", nrow(dat), " cells, years ",
        paste(range(dat$year), collapse = "-"), "\n")

diag <- run_distribution_diagnostics(dat)

# ---- Save all sub-tables for inclusion in the paper's supplementary material ----
readr::write_csv(diag$patients$moment_check,    file.path(cfg$output$tables_dir, "diag_patients_moments.csv"))
readr::write_csv(diag$patients$aic_comparison,   file.path(cfg$output$tables_dir, "diag_patients_aic.csv"))
readr::write_csv(diag$patients$zero_check,       file.path(cfg$output$tables_dir, "diag_patients_zeros.csv"))

readr::write_csv(diag$provision$moment_check,    file.path(cfg$output$tables_dir, "diag_provision_moments.csv"))
readr::write_csv(diag$provision$aic_comparison,  file.path(cfg$output$tables_dir, "diag_provision_aic.csv"))
readr::write_csv(diag$provision$zero_check,      file.path(cfg$output$tables_dir, "diag_provision_zeros.csv"))

readr::write_csv(diag$value$zero_mass_check,             file.path(cfg$output$tables_dir, "diag_value_zero_mass.csv"))
readr::write_csv(diag$value$aic_comparison,                file.path(cfg$output$tables_dir, "diag_value_aic.csv"))
readr::write_csv(diag$value$mean_variance_relationship,   file.path(cfg$output$tables_dir, "diag_value_power_law.csv"))

message("\nSaved 9 diagnostic tables to ", cfg$output$tables_dir, "/diag_*.csv")
message("Done.")
