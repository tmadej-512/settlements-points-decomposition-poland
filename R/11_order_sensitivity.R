# ============================================================
# R/11_order_sensitivity.R
# ------------------------------------------------------------
# Quantifies how much the sequential decomposition's attribution
# would change under different substitution orderings, and compares
# it directly to the order-independent Shapley result.
#
# This is the empirical motivation for using Shapley rather than a
# single arbitrary sequential order: if the 24 orderings gave nearly
# identical attributions, order-independence would be a theoretical
# nicety with no practical consequence. If they diverge substantially
# (as is typical whenever components interact), this is direct
# evidence that the choice of method matters for substantive
# conclusions, and that Shapley's averaging is doing real work
# rather than coinciding with any one arbitrary order.
# ============================================================

suppressPackageStartupMessages(library(dplyr))

#' Compare all 24 sequential orderings against the Shapley average
#' for one (baseline, comparison) year pair.
#'
#' @param df_2y Row-level data containing exactly two distinct years.
#' @param fit Output of fit_model_chain().
#' @param baseline_year,comparison_year The two years to compare.
#' @return list(
#'   all_orders = <24-row tibble, one row per ordering, from
#'                 decompose_sequential_all_orders_2y()>,
#'   shapley    = <1-row tibble, from decompose_shapley_2y()>,
#'   range_summary = <one row per component: min/max/range/SD of the
#'                    effect across the 24 orderings, plus the Shapley
#'                    value and where it falls within that range>
#' )
compare_order_sensitivity <- function(df_2y, fit,
                                       baseline_year = min(df_2y$year),
                                       comparison_year = max(df_2y$year)) {
  all_orders <- decompose_sequential_all_orders_2y(df_2y, fit, baseline_year, comparison_year)
  shap <- decompose_shapley_2y(df_2y, fit, baseline_year, comparison_year)

  comp_cols <- c("N_effect", "PI_effect", "LAMBDA_effect", "KAPPA_effect")

  range_summary <- purrr::map_dfr(comp_cols, function(cc) {
    vals <- all_orders[[cc]]
    tibble::tibble(
      component   = sub("_effect$", "", cc),
      seq_min      = min(vals),
      seq_max      = max(vals),
      seq_range    = max(vals) - min(vals),
      seq_sd       = sd(vals),
      shapley_value = shap[[cc]],
      shapley_within_range = shap[[cc]] >= min(vals) - 1e-9 & shap[[cc]] <= max(vals) + 1e-9
    )
  })

  list(all_orders = all_orders, shapley = shap, range_summary = range_summary)
}
