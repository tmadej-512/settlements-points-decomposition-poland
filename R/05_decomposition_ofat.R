# ============================================================
# R/05_decomposition_ofat.R
# ------------------------------------------------------------
# One-Factor-At-a-Time (OFAT) decomposition.
#
# For each component, build the counterfactual scenario in which
# ONLY that component is switched from baseline to comparison-year
# levels, with all THREE other components held at baseline. This is
# the simplest model-based attribution: each component's effect is
# read off in isolation, against the same baseline starting point.
#
# IMPORTANT: this method does NOT, in general, satisfy additivity --
# the four single-component effects need not sum to the model-implied
# total_change, because it ignores interactions between components.
# The `residual` column reports this gap explicitly; a paper
# reporting this method should show the residual, not hide it.
# (Contrast with sequential and Shapley, both of which are
# additive by construction -- see 06 and 07.)
# ============================================================

suppressPackageStartupMessages(library(dplyr))

#' OFAT decomposition for one (baseline, comparison) year pair
#'
#' @param df_2y Row-level data containing exactly two distinct years.
#' @param fit Output of fit_model_chain().
#' @param baseline_year,comparison_year The two years to compare.
#' @return One-row tibble: N_effect, PI_effect, LAMBDA_effect, KAPPA_effect,
#'   residual (= total_change - sum of the four effects; not exactly 0).
decompose_ofat_2y <- function(df_2y, fit,
                               baseline_year = min(df_2y$year),
                               comparison_year = max(df_2y$year)) {
  stopifnot(length(unique(df_2y$year)) == 2)
  df_2y <- df_2y %>% filter(year %in% c(baseline_year, comparison_year))
  p <- prep_pair(df_2y, baseline_year, comparison_year)
  tgt <- p$target
  bas <- p$base

  V0 <- scenario_from_set(tgt, bas, fit, baseline_year, comparison_year, character(0))
  V4 <- scenario_from_set(tgt, bas, fit, baseline_year, comparison_year, c("N", "PI", "LAMBDA", "KAPPA"))

  V_N  <- scenario_from_set(tgt, bas, fit, baseline_year, comparison_year, "N")
  V_PI <- scenario_from_set(tgt, bas, fit, baseline_year, comparison_year, "PI")
  V_LA <- scenario_from_set(tgt, bas, fit, baseline_year, comparison_year, "LAMBDA")
  V_KA <- scenario_from_set(tgt, bas, fit, baseline_year, comparison_year, "KAPPA")

  N_effect      <- V_N  - V0
  PI_effect     <- V_PI - V0
  LAMBDA_effect <- V_LA - V0
  KAPPA_effect  <- V_KA - V0
  total_change  <- V4 - V0

  tibble::tibble(
    baseline_year    = baseline_year,
    comparison_year  = comparison_year,
    baseline_total   = V0,
    comparison_total = V4,
    total_change     = total_change,
    N_effect         = N_effect,
    PI_effect        = PI_effect,
    LAMBDA_effect    = LAMBDA_effect,
    KAPPA_effect     = KAPPA_effect,
    residual         = total_change - (N_effect + PI_effect + LAMBDA_effect + KAPPA_effect)
  )
}
