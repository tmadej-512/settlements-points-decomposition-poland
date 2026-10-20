#!/usr/bin/env Rscript
# ============================================================
# scripts/run_01_fit_model.R
# ------------------------------------------------------------
# Loads the raw NFZ extract, filters to a single disease group
# (config$data$disease_filter, e.g. "Cardiovascular diseases"),
# aggregates to analysis cells, and fits the three-stage GAM/GLM
# model chain. Saves the fitted chain and the aggregated data to
# output/ for use by all later scripts.
#
# Run from the project root:
#   Rscript scripts/run_01_fit_model.R
# ============================================================

source("R/99_load_all.R")

cfg <- load_config("config.R")
ensure_dir(cfg$output$tables_dir)
ensure_dir(cfg$output$figures_dir)
ensure_dir("output/rds")

message("== Step 1: load and normalise data ==")
df_raw <- read_input_data(cfg$data$path)
dat0   <- normalize_schema(df_raw, value_var = cfg$data$value_var)

message("Disease groups available in this extract:")
print(sort(unique(dat0$icd10)))

message("\n== Step 2: filter/aggregate disease group and aggregate to cells ==")
dat <- filter_and_aggregate(dat0, disease_filter = cfg$data$disease_filter)
disease_label <- attr(dat, "disease_filter") %||% cfg$data$disease_filter %||% unique(dat$icd10)[1]
message("Selected disease group: '", disease_label, "'")

message("Aggregated data: ", nrow(dat), " rows, years ",
        paste(range(dat$year), collapse = "-"), ", ",
        length(unique(dat$code)), " administrative units.")

# Raw descriptive identity table for the supplement. This is intentionally
# computed before model fitting: it shows the observed aggregate population,
# patients, services, value, and the three ratios pi/lambda/kappa.
raw_identity_all <- make_raw_identity_table(dat)
raw_identity_main <- make_raw_identity_table(dat, years = c(cfg$data$baseline_year, max(dat$year)))
readr::write_csv(raw_identity_all, file.path(cfg$output$tables_dir, "SUPP_raw_identity_all_years.csv"))
readr::write_csv(raw_identity_main, file.path(cfg$output$tables_dir, "SUPP_raw_identity_baseline_comparison.csv"))
message("Saved raw identity checks: SUPP_raw_identity_all_years.csv and SUPP_raw_identity_baseline_comparison.csv")

message("\n== Step 3: fit the 3-stage model chain (patients -> provision -> value) ==")
fit <- fit_model_chain(dat, cfg$model)

q <- quality_summary(fit, disease_label)
message("\nModel quality (deviance explained):")
print(q)

saveRDS(dat, "output/rds/aggregated_data.rds")
saveRDS(fit, "output/rds/model_fit.rds")
readr::write_csv(q, file.path(cfg$output$tables_dir, "model_quality.csv"))

message("\nSaved: output/rds/aggregated_data.rds, output/rds/model_fit.rds")
message("Saved: ", file.path(cfg$output$tables_dir, "model_quality.csv"))
message("Done.")
