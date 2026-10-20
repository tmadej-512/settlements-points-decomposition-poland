#!/usr/bin/env Rscript
# ============================================================
# scripts/run_07_make_paper_outputs.R
# ------------------------------------------------------------
# Assembles the tables and figures referenced in the paper's
# Results and Validation sections from the outputs of scripts
# 01-06. Produces a single manuscript-ready table (Table 1: method
# comparison with bootstrap CIs) and a combined validation figure.
#
# Run from the project root, after running scripts 01-06 (or after
# scripts/run_00_all.R):
#   Rscript scripts/run_07_make_paper_outputs.R
# ============================================================

source("R/99_load_all.R")

cfg <- load_config("config.R")
make_plots <- !(exists(".NFZ_SKIP_PLOTS") && isTRUE(.NFZ_SKIP_PLOTS)) && !("--no-plots" %in% commandArgs(trailingOnly = TRUE))
tdir <- cfg$output$tables_dir
fdir <- cfg$output$figures_dir

need <- c(
  "method_comparison.csv", "bootstrap_summary.csv",
  "order_sensitivity_summary.csv", "simulation_summary.csv",
  "axiom_checks_toy_games.csv", "axiom_efficiency_empirical.csv",
  "model_quality.csv", "SUPP_raw_identity_baseline_comparison.csv",
  "TABLE_model_vs_observed_change.csv"
)
missing <- need[!file.exists(file.path(tdir, need))]
if (length(missing) > 0) {
  stop("Missing required input table(s): ", paste(missing, collapse = ", "),
       "\nRun scripts/run_00_all.R (or the individual run_01..06 scripts) first.")
}

method_comparison <- readr::read_csv(file.path(tdir, "method_comparison.csv"), show_col_types = FALSE)
boot_summary       <- readr::read_csv(file.path(tdir, "bootstrap_summary.csv"), show_col_types = FALSE)
order_sensitivity   <- readr::read_csv(file.path(tdir, "order_sensitivity_summary.csv"), show_col_types = FALSE)
sim_summary         <- readr::read_csv(file.path(tdir, "simulation_summary.csv"), show_col_types = FALSE)
axiom_toy           <- readr::read_csv(file.path(tdir, "axiom_checks_toy_games.csv"), show_col_types = FALSE)
axiom_empirical      <- readr::read_csv(file.path(tdir, "axiom_efficiency_empirical.csv"), show_col_types = FALSE)
model_quality        <- readr::read_csv(file.path(tdir, "model_quality.csv"), show_col_types = FALSE)
raw_identity_main    <- readr::read_csv(file.path(tdir, "SUPP_raw_identity_baseline_comparison.csv"), show_col_types = FALSE)
model_vs_observed    <- readr::read_csv(file.path(tdir, "TABLE_model_vs_observed_change.csv"), show_col_types = FALSE)

# ---- Table 0 / supplement: raw identity and model-vs-observed alignment ----
readr::write_csv(model_quality, file.path(tdir, "TABLE0_model_quality.csv"))
readr::write_csv(raw_identity_main, file.path(tdir, "SUPP_TABLE_raw_identity_check.csv"))
readr::write_csv(model_vs_observed, file.path(tdir, "SUPP_TABLE_model_vs_observed_change.csv"))

# ---- Table 1: point estimates + bootstrap CIs, all 4 methods, all 4 components ----
table1 <- method_comparison %>%
  tidyr::pivot_longer(
    cols = c(N_effect, PI_effect, LAMBDA_effect, KAPPA_effect),
    names_to = "component", values_to = "point_estimate"
  ) %>%
  mutate(component = sub("_effect$", "", component)) %>%
  dplyr::select(method, decomposition_total_type, component, point_estimate, residual,
                total_change, observed_total_change, change_difference_vs_observed,
                change_difference_pct_observed) %>%
  left_join(boot_summary, by = c("method", "component")) %>%
  mutate(
    component = factor(component, levels = c("N", "PI", "LAMBDA", "KAPPA")),
    method = factor(method, levels = c("naive", "ofat", "sequential", "shapley"))
  ) %>%
  arrange(component, method)

readr::write_csv(table1, file.path(tdir, "TABLE1_method_comparison_with_CI.csv"))

# ---- Table 2: order-sensitivity vs Shapley ----
readr::write_csv(order_sensitivity, file.path(tdir, "TABLE2_order_sensitivity.csv"))

# ---- Table 3: simulation recovery (bias / RMSE vs ground truth) ----
readr::write_csv(sim_summary, file.path(tdir, "TABLE3_simulation_recovery.csv"))

# ---- Table 4: axiom validation summary ----
table4 <- bind_rows(
  axiom_toy %>% mutate(check_type = "toy game (analytic ground truth)"),
  axiom_empirical %>%
    transmute(
      axiom = "Efficiency (empirical)",
      passed = passed,
      detail = sprintf("residual = %.4g on real fitted model, %s vs %s",
                        residual, baseline_year, comparison_year),
      check_type = "real fitted model"
    )
)
readr::write_csv(table4, file.path(tdir, "TABLE4_axiom_validation.csv"))

message("Paper-ready tables written to ", tdir, ":")
message("  TABLE0_model_quality.csv")
message("  SUPP_TABLE_raw_identity_check.csv")
message("  SUPP_TABLE_model_vs_observed_change.csv")
message("  TABLE1_method_comparison_with_CI.csv")
message("  TABLE2_order_sensitivity.csv")
message("  TABLE3_simulation_recovery.csv")
message("  TABLE4_axiom_validation.csv")

# ---- Combined validation figure (4-panel) ----
if (make_plots) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package {patchwork} is not installed; skipping FIGURE_validation_combined.png")
  } else {
    suppressPackageStartupMessages(library(patchwork))

    p1 <- plot_method_comparison_bars(method_comparison) + ggplot2::labs(title = "A. Method comparison")
    p2 <- plot_order_sensitivity(order_sensitivity) + ggplot2::labs(title = "B. Order sensitivity")
    p3 <- plot_simulation_recovery_vs_truth(
      readr::read_csv(file.path(tdir, "simulation_raw.csv"), show_col_types = FALSE)
    ) + ggplot2::labs(title = "C. Simulation recovery")

    p4 <- plot_bootstrap_intervals(
      boot_summary,
      method_comparison
    ) + ggplot2::labs(title = "D. Bootstrap intervals")

    combined <- (p1 + p2) / (p3 + p4)

    ensure_dir(fdir)
    ggplot2::ggsave(
      file.path(fdir, "FIGURE_validation_combined.png"),
      combined, width = cfg$output$fig_width + 6, height = (cfg$output$fig_height + 1) * 2,
      dpi = cfg$output$fig_dpi
    )

    message("\nSaved: ", file.path(fdir, "FIGURE_validation_combined.png"))
  }
} else {
  message("Skipping combined validation figure because .NFZ_SKIP_PLOTS=TRUE or --no-plots was used.")
}
message("Done.")
