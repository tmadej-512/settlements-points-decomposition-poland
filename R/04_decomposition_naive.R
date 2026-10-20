# ============================================================
# R/04_decomposition_naive.R
# ------------------------------------------------------------
# Aggregate LMDI-style/log-ratio decomposition: no model, no counterfactual
# prediction. It just multiplies out the observed rate identity
#
#   Value = N * (n_patients / N) * (n_provision / n_patients) * (value / n_provision)
#         = N * pi * lambda * kappa
#
# at the aggregate (disease-total) level, and attributes the change
# between two periods using log-change shares (a standard "additive
# decomposition of a multiplicative identity" approach, in the
# spirit of classical demographic standardisation).
#
# PURPOSE IN THIS PROJECT:
# This is a sanity check, not a candidate method for the paper's
# main results. It answers: "if we ignore all covariate structure,
# all heterogeneity across cells, and all of the GAM/GLM machinery,
# and just look at the four aggregate ratios, what story do the raw
# numbers tell?" It should be in the same ballpark as OFAT/sequential/
# Shapley if the model chain is doing its job, and large, unexplained
# divergence from it is itself diagnostic.
#
# Method: log-change (Sato-Vartia-style) shares.
#   Let r_N = log(N1/N0), r_PI = log(pi1/pi0), r_LA = log(lambda1/lambda0),
#       r_KA = log(kappa1/kappa0). Since log(Value1/Value0) = r_N + r_PI
#   + r_LA + r_KA exactly (the identity is multiplicative), each
#   component's share of the *level* change is allocated proportionally
#   to its share of the *log* change:
#       effect_j = (r_j / (r_N + r_PI + r_LA + r_KA)) * (Value1 - Value0)
#   This guarantees exact additivity by construction, just like the
#   other three methods, but with no model and no row-level heterogeneity.
# ============================================================

suppressPackageStartupMessages(library(dplyr))

#' Aggregate LMDI-style log-ratio decomposition for one (baseline, comparison) pair
#'
#' @param df_2y Row-level data containing exactly baseline_year and comparison_year.
#' @param fit Unused; accepted only so this function has the same signature
#'   as the model-based methods and can be swapped in interchangeably by
#'   the baseline/rolling/trend runners.
#' @param baseline_year,comparison_year The two years to compare.
#' @return A one-row tibble with the same column layout as the other
#'   three decomposition methods (N_effect, PI_effect, LAMBDA_effect,
#'   KAPPA_effect, residual, ...), so it can be stacked with their output.
decompose_naive_2y <- function(df_2y, fit = NULL,
                                baseline_year = min(df_2y$year),
                                comparison_year = max(df_2y$year)) {
  stopifnot(length(unique(df_2y$year)) == 2)

  agg <- df_2y %>%
    filter(year %in% c(baseline_year, comparison_year)) %>%
    group_by(year) %>%
    summarise(
      N    = sum(pop, na.rm = TRUE),
      PAT  = sum(n_patients, na.rm = TRUE),
      PROV = sum(n_provision, na.rm = TRUE),
      VAL  = sum(value, na.rm = TRUE),
      .groups = "drop"
    )

  a0 <- agg %>% filter(year == baseline_year)
  a1 <- agg %>% filter(year == comparison_year)

  # Guard against degenerate cells (zero patients/provision in either period):
  # the multiplicative identity is undefined there, so we floor at a tiny
  # epsilon purely to keep the log-ratios finite; this never affects the
  # model-based methods, only this diagnostic one.
  eps <- 1e-9
  N0 <- a0$N; N1 <- a1$N
  pi0  <- a0$PAT  / max(a0$N,    eps); pi1  <- a1$PAT  / max(a1$N,    eps)
  lam0 <- a0$PROV / max(a0$PAT,  eps); lam1 <- a1$PROV / max(a1$PAT,  eps)
  kap0 <- a0$VAL  / max(a0$PROV, eps); kap1 <- a1$VAL  / max(a1$PROV, eps)

  r_N  <- log(max(N1, eps)   / max(N0, eps))
  r_PI <- log(max(pi1, eps)  / max(pi0, eps))
  r_LA <- log(max(lam1, eps) / max(lam0, eps))
  r_KA <- log(max(kap1, eps) / max(kap0, eps))
  r_tot <- r_N + r_PI + r_LA + r_KA

  total_change <- a1$VAL - a0$VAL

  # If r_tot is ~0 (e.g. no change at all), fall back to equal shares to
  # avoid a 0/0 division; the resulting effects are then themselves ~0.
  shares <- if (abs(r_tot) < 1e-12) {
    c(N = 0.25, PI = 0.25, LAMBDA = 0.25, KAPPA = 0.25)
  } else {
    c(N = r_N, PI = r_PI, LAMBDA = r_LA, KAPPA = r_KA) / r_tot
  }

  tibble::tibble(
    baseline_year    = baseline_year,
    comparison_year  = comparison_year,
    baseline_total    = a0$VAL,
    comparison_total  = a1$VAL,
    total_change      = total_change,
    N_effect          = shares["N"]      * total_change,
    PI_effect         = shares["PI"]     * total_change,
    LAMBDA_effect     = shares["LAMBDA"] * total_change,
    KAPPA_effect      = shares["KAPPA"]  * total_change,
    residual          = total_change - (
      shares["N"] * total_change + shares["PI"] * total_change +
      shares["LAMBDA"] * total_change + shares["KAPPA"] * total_change
    )
  )
}
