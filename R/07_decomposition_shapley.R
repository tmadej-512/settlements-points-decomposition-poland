# ============================================================
# R/07_decomposition_shapley.R
# ------------------------------------------------------------
# Shapley-value decomposition.
#
# For each of the 4! = 24 possible orderings of {N, PI, LAMBDA, KAPPA},
# compute the marginal contribution of each component when it enters
# at its position in that ordering (exactly as in the sequential
# method, R/06). The Shapley value of each component is the AVERAGE
# of its marginal contribution across all 24 orderings.
#
# This is, by construction:
#   (i)   additive:   sum(shapley effects) = model-implied total_change, exactly;
#   (ii)  order-independent: unlike the sequential method (R/06),
#         it does not privilege one substitution order over another;
#   (iii) the unique allocation rule satisfying the Shapley axioms
#         (efficiency, symmetry, dummy/null player, additivity)
#         applied to the cooperative game whose "players" are the
#         four components and whose coalition value function is
#         scenario_from_set().
#
# This file is the primary method used in the paper. R/04 (naive)
# and R/05 (OFAT) are reported as comparators; R/06 (sequential,
# single order) is reported as the auxiliary analysis motivating
# why an order-independent method is needed.
# ============================================================

suppressPackageStartupMessages(library(dplyr))

#' Shapley decomposition for one (baseline, comparison) year pair
#'
#' @param df_2y Row-level data containing exactly two distinct years.
#' @param fit Output of fit_model_chain().
#' @param baseline_year,comparison_year The two years to compare.
#' @return One-row tibble: N_effect, PI_effect, LAMBDA_effect, KAPPA_effect
#'   (each a Shapley value, i.e. an average marginal contribution over
#'   24 orderings), residual (numerically ~0 by construction, reported
#'   for transparency / as a numerical-precision check), and
#'   contrib_matrix (24 x 4 matrix of all marginal contributions,
#'   attached as an attribute, used by the order-sensitivity diagnostics
#'   and by the axiom-check tests in R/08).
decompose_shapley_2y <- function(df_2y, fit,
                                  baseline_year = min(df_2y$year),
                                  comparison_year = max(df_2y$year)) {
  stopifnot(length(unique(df_2y$year)) == 2)
  df_2y <- df_2y %>% filter(year %in% c(baseline_year, comparison_year))
  p <- prep_pair(df_2y, baseline_year, comparison_year)
  tgt <- p$target
  bas <- p$base

  comps <- c("N", "PI", "LAMBDA", "KAPPA")
  perms <- combinat_permutations(comps)
  contrib_mat <- matrix(0, nrow = length(perms), ncol = 4, dimnames = list(NULL, comps))

  for (pi_idx in seq_along(perms)) {
    perm   <- perms[[pi_idx]]
    active <- character(0)
    V_prev <- scenario_from_set(tgt, bas, fit, baseline_year, comparison_year, active)
    for (comp in perm) {
      active <- c(active, comp)
      V_curr <- scenario_from_set(tgt, bas, fit, baseline_year, comparison_year, active)
      contrib_mat[pi_idx, comp] <- V_curr - V_prev
      V_prev <- V_curr
    }
  }

  sv <- colMeans(contrib_mat)
  V0 <- scenario_from_set(tgt, bas, fit, baseline_year, comparison_year, character(0))
  V4 <- scenario_from_set(tgt, bas, fit, baseline_year, comparison_year, comps)
  total_change <- V4 - V0

  out <- tibble::tibble(
    baseline_year    = baseline_year,
    comparison_year  = comparison_year,
    baseline_total   = V0,
    comparison_total = V4,
    total_change     = total_change,
    N_effect         = unname(sv["N"]),
    PI_effect        = unname(sv["PI"]),
    LAMBDA_effect    = unname(sv["LAMBDA"]),
    KAPPA_effect     = unname(sv["KAPPA"]),
    residual         = total_change - sum(sv)
  )
  attr(out, "contrib_matrix") <- contrib_mat
  out
}
