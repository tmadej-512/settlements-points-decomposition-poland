# ============================================================
# R/06_decomposition_sequential.R
# ------------------------------------------------------------
# Sequential ("chained substitution") decomposition.
#
# Components are switched on one at a time, in a SINGLE fixed
# order (default N -> PI -> LAMBDA -> KAPPA, configurable via
# config$decomposition$sequential_order). Each component's effect
# is the marginal change in predicted total caused by switching it
# on, GIVEN that all earlier components in the order are already
# switched on. This is additive by construction:
#
#   V0 -> V1 -> V2 -> V3 -> V4,   effect_k = V_k - V_{k-1}
#   sum(effects) = V4 - V0 = model-implied total_change   exactly.
#
# THE ORDER-DEPENDENCE PROBLEM:
# Because components generally interact (the model chain is
# non-linear), a different ordering would, in general, attribute
# a different amount to each component. This function reports
# exactly one (arbitrary) ordering. The companion file
# 12_order_sensitivity.R quantifies how much the attribution
# would change under the other 23 possible orderings -- this is
# the empirical motivation for the order-independent Shapley
# decomposition in 07_decomposition_shapley.R.
# ============================================================

suppressPackageStartupMessages(library(dplyr))

#' Sequential decomposition for one (baseline, comparison) year pair,
#' using a single fixed order of components.
#'
#' @param df_2y Row-level data containing exactly two distinct years.
#' @param fit Output of fit_model_chain().
#' @param baseline_year,comparison_year The two years to compare.
#' @param order Character vector giving the substitution order, a
#'   permutation of c("N","PI","LAMBDA","KAPPA"). Defaults to the
#'   natural causal-chain order N -> PI -> LAMBDA -> KAPPA.
#' @return One-row tibble: N_effect, PI_effect, LAMBDA_effect, KAPPA_effect,
#'   residual (numerically ~0 by construction; reported for transparency).
decompose_sequential_2y <- function(df_2y, fit,
                                     baseline_year = min(df_2y$year),
                                     comparison_year = max(df_2y$year),
                                     order = c("N", "PI", "LAMBDA", "KAPPA")) {
  stopifnot(length(unique(df_2y$year)) == 2)
  stopifnot(setequal(order, c("N", "PI", "LAMBDA", "KAPPA")))

  df_2y <- df_2y %>% filter(year %in% c(baseline_year, comparison_year))
  p <- prep_pair(df_2y, baseline_year, comparison_year)
  tgt <- p$target
  bas <- p$base

  active <- character(0)
  V_prev <- scenario_from_set(tgt, bas, fit, baseline_year, comparison_year, active)
  V0 <- V_prev

  effects <- setNames(numeric(4), order)
  for (comp in order) {
    active <- c(active, comp)
    V_curr <- scenario_from_set(tgt, bas, fit, baseline_year, comparison_year, active)
    effects[comp] <- V_curr - V_prev
    V_prev <- V_curr
  }
  V4 <- V_prev
  total_change <- V4 - V0

  tibble::tibble(
    baseline_year    = baseline_year,
    comparison_year  = comparison_year,
    baseline_total   = V0,
    comparison_total = V4,
    total_change     = total_change,
    N_effect         = effects["N"],
    PI_effect        = effects["PI"],
    LAMBDA_effect    = effects["LAMBDA"],
    KAPPA_effect     = effects["KAPPA"],
    residual         = total_change - sum(effects),
    sequential_order = paste(order, collapse = " -> ")
  )
}

#' Run the sequential decomposition under ALL 24 possible orderings
#'
#' This is the empirical evidence for order-dependence referenced in
#' the paper: it shows how much each component's attributed effect
#' varies depending on which order it is substituted in, directly
#' motivating the order-independent Shapley average in 07.
#'
#' @return A tibble with one row per ordering (24 rows), plus the
#'   range and standard deviation of each component's effect across
#'   orderings, attached as attributes for convenience.
decompose_sequential_all_orders_2y <- function(df_2y, fit,
                                                baseline_year = min(df_2y$year),
                                                comparison_year = max(df_2y$year)) {
  comps <- c("N", "PI", "LAMBDA", "KAPPA")
  orders <- combinat_permutations(comps)

  out <- purrr::map_dfr(orders, function(ord) {
    decompose_sequential_2y(df_2y, fit, baseline_year, comparison_year, order = ord)
  })
  out
}

#' All 24 permutations of a length-4 vector (no external dependency on
#' the `combinat` package; implemented directly for portability)
combinat_permutations <- function(v) {
  stopifnot(length(v) == 4)
  res <- list()
  for (a in v) {
    for (b in setdiff(v, a)) {
      for (cc in setdiff(v, c(a, b))) {
        for (d in setdiff(v, c(a, b, cc))) {
          res <- c(res, list(c(a, b, cc, d)))
        }
      }
    }
  }
  res
}
