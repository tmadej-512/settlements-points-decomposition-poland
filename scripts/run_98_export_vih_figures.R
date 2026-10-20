#!/usr/bin/env Rscript
# ============================================================
# scripts/run_98_export_vih_figures.R
# ============================================================
# Exports the curated Value in Health submission figure set.
#
# Default behaviour in this version:
#   - exports all main figures;
#   - exports all supplementary figures we have in the manuscript/supplement;
#   - exports each figure as a separate image file;
#   - removes title/subtitle/caption from the image itself;
#   - writes a manifest and QA checklist.
#
# Main figures:
#   Figure 1. Empirical method comparison by component, 2018-2024
#   Figure 2. Multi-year empirical LMDI-style and Shapley decompositions
#   Figure 3. Overall simulation RMSE by method and scenario
#
# Supplementary figures (curated manuscript set):
#   S1. Order sensitivity (former S1)
#   S2. Bootstrap intervals (former S2)
#   S3. Shapley vs LMDI differences (former S5)
# Former S3 additivity residual, S4 validation panel and S6 residual
# diagnostics are not submission figures. They remain available in analysis outputs.
# Moran's I is a separate statistical diagnostic, never a numbered figure.
# The exporter does NOT fit models or rerun simulations.
#
# Optional flags:
#   --tiff                 export TIFF instead of PNG
#   --jpg                  export JPG instead of PNG
#   --main-only            export only main manuscript figures
#   --no-scenario-specific do not export scenario-specific recovery plots
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

source("R/99_load_all.R")
source("R/99_vih_figure_export_helpers.R")

cfg <- load_config("config.R")
args <- commandArgs(trailingOnly = TRUE)

fig_device <- if ("--tiff" %in% args) "tiff" else if ("--jpg" %in% args) "jpg" else "png"
export_supp <- !("--main-only" %in% args)
export_scenario_specific <- "--scenario-specific" %in% args && !("--no-scenario-specific" %in% args)

fig_dpi <- cfg$output$fig_dpi %||% 600
fig_w <- 6.5
fig_h <- 4.2

tdir <- cfg$output$tables_dir
out_main <- file.path("submission", "figures")
out_supp <- file.path("submission", "figures", "supplement")
out_scen <- file.path(out_supp, "scenario_specific")

dir.create(out_main, recursive = TRUE, showWarnings = FALSE)
if (export_supp) dir.create(out_supp, recursive = TRUE, showWarnings = FALSE)
if (export_supp && export_scenario_specific) dir.create(out_scen, recursive = TRUE, showWarnings = FALSE)

read_csv_dt <- function(path) {
  if (!file.exists(path)) return(NULL)
  data.table::fread(path)
}

manifest <- data.table::data.table(
  figure = character(),
  file = character(),
  title_for_legend = character(),
  manuscript_status = character(),
  notes = character()
)

add_manifest <- function(figure, file, title, manuscript_status = "main", notes = "") {
  manifest <<- data.table::rbindlist(list(
    manifest,
    data.table::data.table(
      figure = figure,
      file = file,
      title_for_legend = title,
      manuscript_status = manuscript_status,
      notes = notes
    )
  ), fill = TRUE)
}

message("Exporting curated Value in Health submission figures (main 1-3, supplementary S1-S3)")
message("Device: ", fig_device, "; DPI: ", fig_dpi)
if (!export_supp) message("--main-only supplied: supplementary figures will not be exported.")

# ============================================================
# Main Figure 1. Empirical method comparison
# ============================================================

method_comparison <- read_csv_dt(file.path(tdir, "method_comparison.csv"))

if (!is.null(method_comparison)) {
  p_fig1 <- plot_method_comparison_bars(method_comparison) +
    labs(x = NULL) +
    geom_hline(yintercept = 0, linewidth = 0.35)
  p_fig1 <- vih_apply_points_billion_axis(
    p_fig1,
    y_label = "Effect on total settlement points (billion)",
    include_zero = TRUE
  )

  f1 <- vih_save_figure(
    p_fig1,
    figure_no = 1,
    stem = "empirical_method_comparison",
    width = fig_w,
    height = fig_h,
    dpi = fig_dpi,
    path = out_main,
    device = fig_device,
    base_size = 10
  )
  add_manifest(
    "Figure 1", f1, "Empirical method comparison by component, 2018-2024",
    "main", "Axis values exported in billions of settlement points."
  )
} else {
  warning("Skipping Figure 1; missing method_comparison.csv in ", tdir)
}

# ============================================================
# Main Figure 2. Empirical multi-year decomposition
# ============================================================

empirical_multiyear <- read_csv_dt(file.path(tdir, "TABLE_empirical_multiyear_lmdi_shapley.csv"))

if (!is.null(empirical_multiyear)) {
  p_fig2 <- plot_empirical_multiyear_components(empirical_multiyear) +
    labs(x = "Comparison year") +
    geom_hline(yintercept = 0, linewidth = 0.35)
  p_fig2 <- vih_apply_points_billion_axis(
    p_fig2,
    y_label = "Component contribution (billion settlement points)",
    include_zero = TRUE
  )

  f2 <- vih_save_figure(
    p_fig2,
    figure_no = 2,
    stem = "empirical_multiyear_decomposition",
    width = fig_w + 1.0,
    height = fig_h + 0.4,
    dpi = fig_dpi,
    path = out_main,
    device = fig_device,
    base_size = 10
  )
  add_manifest(
    "Figure 2", f2, "Multi-year empirical decomposition from 2018 baseline",
    "main", "Axis values exported in billions of settlement points."
  )
} else {
  warning("Skipping Figure 2; missing TABLE_empirical_multiyear_lmdi_shapley.csv in ", tdir)
}

# ============================================================
# Main Figure 3. Simulation RMSE
# ============================================================

sim_raw <- read_csv_dt(file.path(tdir, "simulation_three_scenarios_raw.csv"))
if (is.null(sim_raw)) sim_raw <- read_csv_dt(file.path(tdir, "simulation_raw.csv"))

if (!is.null(sim_raw)) {
  if (!"scenario_label" %in% names(sim_raw)) sim_raw[, scenario_label := scenario]
  sim_perf <- make_simulation_performance_table(sim_raw)

  p_fig3 <- plot_simulation_overall_rmse(sim_perf) +
    labs(x = NULL)
  p_fig3 <- vih_apply_points_million_axis(
    p_fig3,
    y_label = "Root mean squared error (million settlement points)",
    include_zero = TRUE
  )

  f3 <- vih_save_figure(
    p_fig3,
    figure_no = 3,
    stem = "simulation_rmse",
    width = fig_w + 1.5,
    height = fig_h + 0.5,
    dpi = fig_dpi,
    path = out_main,
    device = fig_device,
    base_size = 9
  )
  add_manifest(
    "Figure 3", f3, "Overall simulation RMSE by method and scenario",
    "main", "Axis values exported in millions of settlement points."
  )
} else {
  warning("Skipping Figure 3; missing simulation raw CSV in ", tdir)
}

# ============================================================
# Supplementary figures: exported by default in this all-figures package
# ============================================================

supp_count <- 0L

if (export_supp) {
  supp_no <- 1L

  # Supplementary Figure S1. Order sensitivity
  order_sensitivity <- read_csv_dt(file.path(tdir, "order_sensitivity_summary.csv"))
  if (!is.null(order_sensitivity)) {
    p_s1 <- plot_order_sensitivity(order_sensitivity) +
      labs(x = NULL) +
      geom_hline(yintercept = 0, linewidth = 0.35)
    p_s1 <- vih_apply_points_billion_axis(
      p_s1,
      y_label = "Effect on total settlement points (billion)",
      include_zero = TRUE
    )
    fs1 <- vih_save_supplement_figure(p_s1, supp_no, "order_sensitivity", fig_w, fig_h, fig_dpi, out_supp, fig_device)
    add_manifest("Supplementary Figure S1", fs1, "Order sensitivity of sequential substitution", "supplement", "Axis values exported in billions.")
    supp_no <- supp_no + 1L
  } else {
    warning("Skipping Supplementary Figure S1; missing order_sensitivity_summary.csv")
  }

  # Supplementary Figure S2. Bootstrap intervals
  boot_summary <- read_csv_dt(file.path(tdir, "bootstrap_summary.csv"))
  if (!is.null(boot_summary) && !is.null(method_comparison)) {
    p_s2 <- plot_bootstrap_intervals(boot_summary, method_comparison) +
      labs(x = NULL) +
      geom_hline(yintercept = 0, linewidth = 0.35)
    p_s2 <- vih_apply_points_billion_axis(
      p_s2,
      y_label = "Effect on total settlement points (billion)",
      include_zero = TRUE
    )
    fs2 <- vih_save_supplement_figure(p_s2, supp_no, "bootstrap_intervals", fig_w, fig_h, fig_dpi, out_supp, fig_device)
    add_manifest("Supplementary Figure S2", fs2, "Bootstrap intervals for empirical component contributions", "supplement", "Axis values exported in billions.")
    supp_no <- supp_no + 1L
  } else {
    warning("Skipping Supplementary Figure S2; missing bootstrap_summary.csv or method_comparison.csv")
  }

  # Supplementary Figure S3 (previously S5). Shapley vs LMDI-style differences across scenarios
  if (!is.null(sim_raw)) {
    p_sdiff <- plot_shapley_vs_lmdi_difference(sim_raw) +
      labs(x = NULL) +
      geom_hline(yintercept = 0, linewidth = 0.35)
    p_sdiff <- vih_apply_points_million_axis(
      p_sdiff,
      y_label = "Difference in estimated contribution (million settlement points)",
      include_zero = TRUE
    )
    fsd <- vih_save_supplement_figure(p_sdiff, supp_no, "shapley_vs_lmdi_difference", fig_w + 1.5, fig_h + 0.5, fig_dpi, out_supp, fig_device, base_size = 9)
    add_manifest("Supplementary Figure S3", fsd, "Difference between model-based Shapley and aggregate LMDI-style estimates", "supplement", "Axis values exported in millions.")
    supp_no <- supp_no + 1L
  }

  # Optional development-only scenario-specific plots; never numbered as manuscript supplement.
  if (export_scenario_specific && !is.null(sim_raw)) {
    scenarios <- sort(unique(as.character(sim_raw$scenario)))
    sim_plot_limits <- compute_simulation_plot_limits(sim_raw)
    sc_no <- 1L
    for (sc in scenarios) {
      p_truth <- plot_simulation_recovery_vs_truth_by_scenario(sim_raw, sc, y_limits = sim_plot_limits$truth) +
        labs(x = NULL, y = "Estimated contribution")
      fs <- vih_save_figure(
        p_truth, sc_no, paste0("simulation_recovery_vs_truth_", sc),
        width = fig_w + 1.5, height = fig_h, dpi = fig_dpi,
        path = out_scen, device = fig_device, base_size = 9,
        prefix = "Scenario_Figure"
      )
      add_manifest(paste0("Scenario-specific Figure ", sc_no), fs, paste0("Simulation recovery versus truth: ", sc), "additional_supplement", "Exported as additional file; not counted in main display-item limit.")
      sc_no <- sc_no + 1L
    }
  }

  supp_count <- supp_no - 1L
  if (supp_count != 3L) warning("Expected exactly three supplementary figures S1-S3, produced ", supp_count)

}

vih_write_manifest(manifest, file.path(out_main, "FIGURE_MANIFEST.csv"))
main_count <- sum(manifest$manuscript_status == "main")
all_supp_count <- sum(manifest$manuscript_status %in% c("supplement", "additional_supplement"))
vih_write_checklist(
  path = file.path(out_main, "VIH_FIGURE_QA_CHECKLIST.csv"),
  dpi = fig_dpi,
  device = fig_device,
  main_figures_n = main_count,
  supplementary_figures_n = all_supp_count
)

cat("\n============================================================\n")
cat("VIH ALL-FIGURE EXPORT FINISHED\n")
cat("Main figures: ", out_main, "\n", sep = "")
if (export_supp) cat("Supplementary figures: ", out_supp, "\n", sep = "")
if (export_supp && export_scenario_specific) cat("Scenario-specific figures: ", out_scen, "\n", sep = "")
cat("Manifest: ", file.path(out_main, "FIGURE_MANIFEST.csv"), "\n", sep = "")
cat("QA checklist: ", file.path(out_main, "VIH_FIGURE_QA_CHECKLIST.csv"), "\n", sep = "")
cat("============================================================\n")
