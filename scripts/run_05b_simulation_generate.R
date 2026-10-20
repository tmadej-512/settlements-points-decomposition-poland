#!/usr/bin/env Rscript
# ============================================================
# scripts/run_05b_simulation_generate.R
# ------------------------------------------------------------
# Data-only wrapper for the multi-scenario simulation.
#
# Use this when the expensive simulation should be generated once,
# without rendering any figures. It writes the simulation CSVs and RDS
# cache. Later regenerate plots only with:
#
#   source("scripts/run_05c_simulation_plots.R")
#
# Optional:
#   .NFZ_QUICK <- TRUE
#   source("scripts/run_05b_simulation_generate.R")
# ============================================================

.NFZ_SKIP_PLOTS <- TRUE
source("scripts/run_05b_simulation_two_scenarios.R")
rm(.NFZ_SKIP_PLOTS, envir = .GlobalEnv)
