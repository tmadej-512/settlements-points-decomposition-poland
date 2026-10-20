# ============================================================
# R/14_value_in_health_outputs.R
# ------------------------------------------------------------
# Extra tables/figures for a Value-in-Health-style methods paper.
#
# These outputs are designed to make the methodological story explicit:
#   1) aggregate LMDI-style benchmark is a strong comparator,
#   2) model-based Shapley adds value under cell heterogeneity, cell-mix shift,
#      order dependence, and modelled counterfactual scenarios,
#   3) the main estimand is not "a new Shapley value", but an applied
#      cell-level counterfactual attribution framework.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

method_display_name <- function(x) {
  dplyr::recode(
    x,
    naive = "Aggregate LMDI-style",
    ofat = "OFAT",
    sequential = "Sequential",
    shapley = "Model-based Shapley",
    .default = x
  )
}

method_axis_name <- function(x) {
  dplyr::recode(
    x,
    "Aggregate LMDI-style" = "LMDI-style",
    "OFAT" = "OFAT",
    "Sequential" = "Sequential",
    "Model-based Shapley" = "Shapley",
    naive = "LMDI-style",
    ofat = "OFAT",
    sequential = "Sequential",
    shapley = "Shapley",
    .default = x
  )
}

component_symbol <- function(x) {
  dplyr::recode(
    x,
    N = "N",
    PI = "π",
    LAMBDA = "λ",
    KAPPA = "κ",
    .default = x
  )
}

#' Manuscript table: methods compared.
make_methods_feature_table <- function() {
  tibble::tribble(
    ~method, ~scale, ~uses_cell_structure, ~order_independent, ~requires_model, ~exact_for_aggregate_identity, ~main_estimand, ~main_limitation,
    "Raw additive", "Levels", "No", "Yes", "No", "No", "Negative-control comparison: treats a product as a sum", "Not valid for multiplicative identities",
    "Sequential substitution", "Counterfactual levels", "Yes", "No", "Yes", "Yes for chosen order", "Marginal effects under one replacement sequence", "Depends on arbitrary component order",
    "Aggregate LMDI-style", "Aggregate log-ratio / relative change", "No", "Yes", "No", "Yes", "Aggregate decomposition of observed multiplicative identity", "Cannot adjust cell-level composition or modelled counterfactuals",
    "Model-based Shapley", "Cell-level counterfactual levels", "Yes", "Yes", "Yes", "Yes for model value function", "Shapley-averaged contribution to modelled temporal change", "Depends on model specification and cell definition"
  )
}

#' Manuscript table: when each method is sufficient/preferred.
make_method_decision_table <- function() {
  tibble::tribble(
    ~analytic_goal, ~preferred_approach, ~reason,
    "Fast aggregate monitoring of a clean multiplicative identity", "Aggregate LMDI-style", "Exact, transparent, no model estimation error",
    "Attribution under a single explicitly chosen policy sequence", "Sequential substitution", "Represents that sequence, but should not be presented as order-independent",
    "Order-independent attribution of modelled counterfactual change", "Model-based Shapley", "Averages marginal contributions over all component replacement orders",
    "Presence of strong cell heterogeneity or cell-composition shift", "Model-based Shapley, benchmarked against LMDI-style", "Aggregate ratios may hide changes in cell mix",
    "Simple descriptive reporting only", "Aggregate LMDI-style plus raw identity table", "Lower variance and easier communication",
    "Policy-facing decomposition of reimbursed value by component", "Model-based Shapley plus LMDI-style sensitivity", "Separates cell-adjusted counterfactual attribution from aggregate identity decomposition"
  )
}

#' Shapley vs aggregate LMDI-style by scenario/component.
make_shapley_vs_lmdi_simulation_table <- function(sim_raw) {
  required <- c("scenario", "scenario_label", "rep_id", "method", "component", "estimate", "truth")
  missing <- setdiff(required, names(sim_raw))
  if (length(missing) > 0) {
    stop("make_shapley_vs_lmdi_simulation_table(): missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  sim_raw %>%
    dplyr::filter(.data$method %in% c("naive", "shapley")) %>%
    dplyr::select(scenario, scenario_label, rep_id, component, method, estimate, truth) %>%
    tidyr::pivot_wider(names_from = method, values_from = estimate) %>%
    dplyr::mutate(
      lmdi_style = .data$naive,
      difference_shapley_minus_lmdi = .data$shapley - .data$lmdi_style,
      abs_difference = abs(.data$difference_shapley_minus_lmdi),
      lmdi_error = .data$lmdi_style - .data$truth,
      shapley_error = .data$shapley - .data$truth,
      abs_lmdi_error = abs(.data$lmdi_error),
      abs_shapley_error = abs(.data$shapley_error)
    ) %>%
    dplyr::group_by(.data$scenario, .data$scenario_label, .data$component) %>%
    dplyr::summarise(
      n_reps = dplyr::n(),
      truth_mean = mean(.data$truth, na.rm = TRUE),
      lmdi_mean = mean(.data$lmdi_style, na.rm = TRUE),
      shapley_mean = mean(.data$shapley, na.rm = TRUE),
      mean_difference_shapley_minus_lmdi = mean(.data$difference_shapley_minus_lmdi, na.rm = TRUE),
      mean_abs_difference = mean(.data$abs_difference, na.rm = TRUE),
      lmdi_mae = mean(.data$abs_lmdi_error, na.rm = TRUE),
      shapley_mae = mean(.data$abs_shapley_error, na.rm = TRUE),
      lmdi_rmse = sqrt(mean(.data$lmdi_error^2, na.rm = TRUE)),
      shapley_rmse = sqrt(mean(.data$shapley_error^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      rmse_ratio_lmdi_over_shapley = .data$lmdi_rmse / pmax(.data$shapley_rmse, 1e-12)
    ) %>%
    dplyr::arrange(.data$scenario, .data$component)
}

#' Overall simulation performance table: scenario x method.
make_simulation_performance_table <- function(sim_raw) {
  sim_raw %>%
    dplyr::mutate(
      error = .data$estimate - .data$truth,
      abs_error = abs(.data$error),
      sq_error = .data$error^2,
      method_label = method_display_name(as.character(.data$method))
    ) %>%
    dplyr::group_by(.data$scenario, .data$scenario_label, .data$method, .data$method_label) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      mean_bias = mean(.data$error, na.rm = TRUE),
      mean_abs_error = mean(.data$abs_error, na.rm = TRUE),
      rmse = sqrt(mean(.data$sq_error, na.rm = TRUE)),
      median_abs_error = median(.data$abs_error, na.rm = TRUE),
      rank_correlation_estimate_truth = suppressWarnings(stats::cor(.data$estimate, .data$truth, method = "spearman", use = "complete.obs")),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$scenario) %>%
    dplyr::mutate(rmse_rank = dplyr::min_rank(.data$rmse)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$scenario, .data$rmse_rank, .data$method)
}

plot_simulation_overall_rmse <- function(performance_table) {
  performance_table %>%
    dplyr::mutate(
      method_label = factor(.data$method_label,
                            levels = c("Aggregate LMDI-style", "OFAT", "Sequential", "Model-based Shapley")),
      method_axis = factor(method_axis_name(as.character(.data$method_label)),
                           levels = c("LMDI-style", "OFAT", "Sequential", "Shapley")),
      method_key = factor(dplyr::recode(as.character(.data$method_label),
                                        "Aggregate LMDI-style" = "naive",
                                        "OFAT" = "ofat",
                                        "Sequential" = "sequential",
                                        "Model-based Shapley" = "shapley"),
                          levels = c("naive", "ofat", "sequential", "shapley"))
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$method_axis, y = .data$rmse, fill = .data$method_key)) +
    ggplot2::geom_col(width = 0.65, show.legend = FALSE) +
    ggplot2::facet_wrap(~ scenario_label, scales = "fixed") +
    ggplot2::scale_fill_manual(values = METHOD_COLORS[c("naive", "ofat", "sequential", "shapley")]) +
    ggplot2::scale_y_continuous(labels = label_number_pl()) +
    ggplot2::labs(
      title = "Simulation performance by validation scenario",
      subtitle = "Lower RMSE means closer recovery of known component-level truth",
      x = NULL,
      y = "RMSE vs known truth"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1, size = 7.8),
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_shapley_vs_lmdi_difference <- function(sim_raw) {
  diff_df <- sim_raw %>%
    dplyr::filter(.data$method %in% c("naive", "shapley")) %>%
    dplyr::select(scenario, scenario_label, rep_id, component, method, estimate) %>%
    tidyr::pivot_wider(names_from = method, values_from = estimate) %>%
    dplyr::mutate(
      difference = .data$shapley - .data$naive,
      component = factor(component_symbol(.data$component), levels = c("N", "π", "λ", "κ"))
    )

  ggplot2::ggplot(diff_df, ggplot2::aes(x = .data$component, y = .data$difference)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_boxplot(outlier.alpha = 0.25, fill = METHOD_COLORS[["shapley"]], alpha = 0.55) +
    ggplot2::facet_wrap(~ scenario_label, scales = "fixed") +
    ggplot2::scale_y_continuous(labels = label_number_pl()) +
    ggplot2::labs(
      title = "Where model-based Shapley differs from aggregate LMDI-style attribution",
      subtitle = "Difference = Shapley estimate − aggregate LMDI-style estimate",
      x = NULL,
      y = "Difference in component contribution"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Empirical Shapley vs aggregate LMDI-style table from method_comparison.csv.
make_empirical_shapley_vs_lmdi_table <- function(method_comparison) {
  long <- method_comparison %>%
    tidyr::pivot_longer(
      cols = c("N_effect", "PI_effect", "LAMBDA_effect", "KAPPA_effect"),
      names_to = "component",
      values_to = "estimate"
    ) %>%
    dplyr::mutate(component = sub("_effect$", "", .data$component)) %>%
    dplyr::filter(.data$method %in% c("naive", "shapley")) %>%
    dplyr::select(baseline_year, comparison_year, component, method, estimate, dplyr::any_of(c("total_change", "observed_total_change")))

  # Different methods decompose DIFFERENT totals (observed vs model-implied).
  # Do not use total_change as a pivot key: that silently creates 8 partial rows.
  estimates <- long %>%
    dplyr::select(baseline_year, comparison_year, component, method, estimate) %>%
    tidyr::pivot_wider(
      id_cols = c(baseline_year, comparison_year, component),
      names_from = method, values_from = estimate
    ) %>%
    dplyr::rename(lmdi_style = naive) %>%
    dplyr::mutate(
      difference_shapley_minus_lmdi = .data$shapley - .data$lmdi_style,
      abs_difference = abs(.data$difference_shapley_minus_lmdi)
    ) %>%
    dplyr::arrange(.data$component)
  if (nrow(estimates) != 4L || anyNA(estimates$lmdi_style) || anyNA(estimates$shapley))
    stop("Expected four complete paired LMDI/Shapley rows")
  estimates

}

#' Multi-year empirical trends: baseline vs each later year for LMDI-style and Shapley.
run_empirical_multiyear_decomposition <- function(dat, fit, baseline_year, methods = c("naive", "shapley")) {
  years <- sort(unique(dat$year))
  years <- years[years > baseline_year]
  out <- list()

  for (yr in years) {
    df_pair <- dat %>% dplyr::filter(.data$year %in% c(baseline_year, yr))
    rows <- list()
    if ("naive" %in% methods) {
      rows$naive <- decompose_naive_2y(df_pair, NULL, baseline_year, yr) %>% dplyr::mutate(method = "naive")
    }
    if ("shapley" %in% methods) {
      rows$shapley <- decompose_shapley_2y(df_pair, fit, baseline_year, yr) %>% dplyr::mutate(method = "shapley")
    }
    out[[as.character(yr)]] <- dplyr::bind_rows(rows)
  }

  dplyr::bind_rows(out) %>%
    tidyr::pivot_longer(
      cols = c("N_effect", "PI_effect", "LAMBDA_effect", "KAPPA_effect"),
      names_to = "component",
      values_to = "effect"
    ) %>%
    dplyr::mutate(
      component = sub("_effect$", "", .data$component),
      method_label = method_display_name(as.character(.data$method))
    )
}

plot_empirical_multiyear_components <- function(multiyear_long) {
  multiyear_long %>%
    dplyr::mutate(
      component = factor(component_symbol(.data$component), levels = c("N", "π", "λ", "κ")),
      method_key = factor(dplyr::recode(as.character(.data$method), naive = "naive", shapley = "shapley"),
                          levels = c("naive", "shapley")),
      method_label = factor(.data$method_label,
                            levels = c("Aggregate LMDI-style", "Model-based Shapley"))
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$comparison_year, y = .data$effect, color = .data$method_key, linetype = .data$method_key)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.0) +
    ggplot2::facet_wrap(~ component, scales = "fixed") +
    ggplot2::scale_color_manual(values = METHOD_COLORS[c("naive", "shapley")],
                                labels = c(naive = "Aggregate LMDI-style", shapley = "Model-based Shapley")) +
    ggplot2::scale_linetype_manual(values = c(naive = "solid", shapley = "22"),
                                   labels = c(naive = "Aggregate LMDI-style", shapley = "Model-based Shapley")) +
    ggplot2::scale_y_continuous(labels = label_number_pl()) +
    ggplot2::labs(
      title = "Empirical multi-year decomposition from a common baseline",
      subtitle = "Aggregate LMDI-style benchmark vs model-based Shapley",
      x = "Comparison year",
      y = "Component contribution",
      color = NULL,
      linetype = NULL
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(order = 1), linetype = ggplot2::guide_legend(order = 1)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
}
