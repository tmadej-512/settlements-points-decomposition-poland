# ============================================================
# scripts/run_00_all.R
# ------------------------------------------------------------
# Runs the entire pipeline end to end, in order, using source()
# for every step. This works IDENTICALLY whether you call it from
# the command line (Rscript) or from inside R/RStudio (source()) --
# there is no dependency on spawning a separate Rscript process.
#
# FROM R / RSTUDIO:
#   setwd("path/to/nfz-shapley-decomposition")   # must contain config.R
#   source("scripts/run_00_all.R")
#
#   # Quick smoke-test (small bootstrap/simulation reps, to confirm
#   # everything runs before committing to the full study):
#   .NFZ_QUICK <- TRUE
#   source("scripts/run_00_all.R")
#
# FROM THE COMMAND LINE:
#   Rscript scripts/run_00_all.R
#   Rscript scripts/run_00_all.R --quick
#
# Either way, must be run from the project root (the directory
# containing config.R and the R/ folder).
# ============================================================

if (!file.exists("config.R")) {
  stop(
    "config.R not found in the current working directory.\n",
    "Set your working directory to the project root first, e.g.:\n",
    '  setwd("path/to/nfz-shapley-decomposition")\n',
    "then re-run source(\"scripts/run_00_all.R\")."
  )
}

# Quick mode can be requested two ways, both supported:
#   1. From R/RStudio: set .NFZ_QUICK <- TRUE before sourcing this file.
#   2. From the command line: pass --quick (handled via commandArgs()).
quick <- FALSE
force_expensive <- FALSE
args <- commandArgs(trailingOnly = TRUE)
if (exists(".NFZ_QUICK") && isTRUE(.NFZ_QUICK)) {
  quick <- TRUE
} else {
  if ("--quick" %in% args) quick <- TRUE
}
if (exists(".NFZ_FORCE_RERUN") && isTRUE(.NFZ_FORCE_RERUN)) {
  force_expensive <- TRUE
} else {
  if ("--force-expensive" %in% args || "--force-cache" %in% args || "--force" %in% args) {
    force_expensive <- TRUE
  }
}
if (force_expensive) {
  .NFZ_FORCE_RERUN <- TRUE
  message("Forcing recomputation of expensive cached steps: simulation and bootstrap.\n")
}

if (quick) {
  message("Running in QUICK mode: simulation study and bootstrap will use ",
          "20 replicates instead of the config.R default, for a fast ",
          "end-to-end smoke test.\n")
  .NFZ_REPS_OVERRIDE <- 20
}

# Report cache settings before expensive steps are reached.
source("R/99_load_all.R")
cfg_for_cache_report <- load_config("config.R")
message("Cache settings:")
message("  simulation: use_cache=", cfg_for_cache_report$validation$simulation_use_cache,
        ", force_rerun=", cfg_for_cache_report$validation$simulation_force_rerun,
        ", file=", cfg_for_cache_report$validation$simulation_cache_file)
message("  bootstrap:  use_cache=", cfg_for_cache_report$validation$bootstrap_use_cache,
        ", force_rerun=", cfg_for_cache_report$validation$bootstrap_force_rerun,
        ", file=", cfg_for_cache_report$validation$bootstrap_cache_file)
rm(cfg_for_cache_report)

run_step <- function(label, path) {
  message("\n", strrep("=", 60))
  message(label)
  message(strrep("=", 60))
  source(path)
}

run_step("STEP 1/7: fit model chain",            "scripts/run_01_fit_model.R")
run_step("STEP 2/7: compare 4 methods",          "scripts/run_02_compare_methods.R")
run_step("STEP 3/7: order-sensitivity analysis", "scripts/run_03_order_sensitivity.R")
run_step("STEP 4/7: axiom validation",           "scripts/run_04_validate_axioms.R")
run_step("STEP 5/7: simulation study",           "scripts/run_05_simulation_study.R")
run_step("STEP 6/7: bootstrap",                  "scripts/run_06_bootstrap.R")
run_step("STEP 7/7: assemble paper outputs",     "scripts/run_07_make_paper_outputs.R")

if (exists(".NFZ_REPS_OVERRIDE")) rm(.NFZ_REPS_OVERRIDE, envir = .GlobalEnv)

message("\n", strrep("=", 60))
message("PIPELINE COMPLETE. See output/tables/ and output/figures/.")
message(strrep("=", 60))
