#!/usr/bin/env Rscript
# ============================================================
# scripts/run_04_validate_axioms.R
# ------------------------------------------------------------
# Runs the formal Shapley axiom checks (efficiency, null player,
# symmetry) on toy games with known analytic answers, plus the
# empirical efficiency check on the real fitted model's Shapley
# decomposition. This is the "does the algorithm do what it claims"
# evidence, distinct from the simulation study's "does the whole
# estimation pipeline recover the truth" evidence (run_05).
#
# Run from the project root, after run_01_fit_model.R and
# run_02_compare_methods.R:
#   Rscript scripts/run_04_validate_axioms.R
# ============================================================

source("R/99_load_all.R")

cfg <- load_config("config.R")

message("== Toy-game axiom checks (known analytic ground truth) ==\n")
axiom_results <- run_axiom_validation_suite()

ensure_dir(cfg$output$tables_dir)
readr::write_csv(axiom_results, file.path(cfg$output$tables_dir, "axiom_checks_toy_games.csv"))

message("\n== Empirical efficiency check on the real fitted model ==")
method_comparison_path <- file.path(cfg$output$tables_dir, "method_comparison.csv")
if (!file.exists(method_comparison_path)) {
  stop("Run scripts/run_02_compare_methods.R first (need: ", method_comparison_path, ")")
}
mc <- readr::read_csv(method_comparison_path, show_col_types = FALSE)
shap_row <- mc %>% dplyr::filter(method == "shapley")

eff_check <- check_efficiency_empirical(shap_row)
print(eff_check)

readr::write_csv(eff_check, file.path(cfg$output$tables_dir, "axiom_efficiency_empirical.csv"))

message("\nSaved: ", file.path(cfg$output$tables_dir, "axiom_checks_toy_games.csv"))
message("Saved: ", file.path(cfg$output$tables_dir, "axiom_efficiency_empirical.csv"))
message("Done.")
