#!/usr/bin/env Rscript
# ============================================================
# scripts/run_00_compute_all.R
# ------------------------------------------------------------
# Runs the full computational pipeline but suppresses all figure
# rendering. Use this when you want to generate/cache data, tables,
# simulations and bootstrap once.
#
# Then regenerate figures without recomputing with:
#   source("scripts/run_99_make_all_plots.R")
#
# Quick smoke test:
#   .NFZ_QUICK <- TRUE
#   source("scripts/run_00_compute_all.R")
# ============================================================

.NFZ_SKIP_PLOTS <- TRUE
source("scripts/run_00_all.R")
rm(.NFZ_SKIP_PLOTS, envir = .GlobalEnv)
