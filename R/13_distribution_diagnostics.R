# ============================================================
# R/13_distribution_diagnostics.R
# ------------------------------------------------------------
# Empirically checks whether the distributional choices made in
# R/02_models.R (negative binomial for patients/provision, Tweedie
# for value) are actually well-supported by the data, rather than
# assumed by convention.
#
# Three independent questions are addressed, one per chain stage:
#
#   1. PATIENTS (count, offset = pop):
#      Is there overdispersion relative to Poisson? Is there excess
#      zero mass beyond what NB/Poisson would predict (suggesting a
#      zero-inflated or hurdle alternative)?
#
#   2. PROVISION (count, offset = n_patients):
#      Same two questions, but provision-per-patient is typically
#      more skewed/long-tailed than patients-per-population, so NB's
#      adequacy is checked more carefully here.
#
#   3. VALUE (continuous, offset = n_provision):
#      Given that value > 0 whenever n_provision > 0 (confirmed for
#      this dataset), is Tweedie's point mass at zero actually doing
#      any work, or would a strictly-positive continuous family
#      (Gamma, inverse Gaussian, log-normal) fit at least as well
#      with one fewer parameter to estimate? Compared by AIC/BIC on
#      matched GLMs (not the full GAM, for speed and comparability
#      across families that mgcv supports differently).
#
# This module produces DIAGNOSTIC OUTPUT ONLY -- it does not change
# which family R/02_models.R actually fits. If the diagnostics
# recommend a different family, update build_formula()/fit_model_chain()
# in R/02_models.R accordingly and re-run the full pipeline.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(MASS)     # glm.nb
  library(mgcv)      # tw(), Gamma family fitting via gam/bam
})

# ------------------------------------------------------------
# 1. COUNT-STAGE DIAGNOSTICS (patients, provision)
# ------------------------------------------------------------

#' Diagnose a count outcome's distributional fit: overdispersion + zero mass
#'
#' @param y Integer vector of observed counts (e.g. n_patients or n_provision).
#' @param offset_var Numeric vector, the offset denominator (e.g. pop or
#'   n_patients), strictly positive where y > 0 is meaningfully defined.
#' @param label A short string for this stage, used in printed output.
#' @return A list with:
#'   - moment_check: the Pearson dispersion statistic from a Poisson fit
#'     (sum of squared Pearson residuals / residual df; >1 indicates
#'     overdispersion relative to Poisson, evaluated correctly at each
#'     observation's own fitted mean rather than via a naive variance of
#'     y/offset, which would conflate offset heterogeneity with genuine
#'     overdispersion)
#'   - aic_comparison: AIC for Poisson vs. negative binomial GLMs with the
#'     same offset and an intercept-only mean structure (a fast proxy for
#'     "does the full GAM need NB", without refitting the full spline model)
#'   - zero_check: observed zero count/proportion vs. the zero count/proportion
#'     each fitted model (Poisson, NB) predicts under its own distribution,
#'     to flag whether a zero-inflated or hurdle extension is warranted
diagnose_count_stage <- function(y, offset_var, label = "count stage") {
  stopifnot(length(y) == length(offset_var))
  keep <- offset_var > 0 & !is.na(y) & !is.na(offset_var)
  y <- y[keep]; offset_var <- offset_var[keep]

  df_fit <- data.frame(y = y, off = offset_var)
  m_pois <- glm(y ~ 1, offset = log(off), family = poisson(), data = df_fit)
  m_nb   <- tryCatch(
    MASS::glm.nb(y ~ 1 + offset(log(off)), data = df_fit),
    error = function(e) NULL
  )

  # CORRECT overdispersion check: the Pearson dispersion statistic,
  # sum((y - mu_hat)^2 / mu_hat) / df_residual, evaluated at the Poisson
  # fit's own fitted means. This is NOT the same as var(y/offset)/mean(y/offset):
  # that naive "rate" calculation conflates genuine overdispersion with
  # heterogeneity in the offset itself (here, population ranges from
  # ~1,800 to ~10,700 across cells, which alone inflates var(rate) in a
  # way unrelated to whether the count process is overdispersed given its
  # own mean). A dispersion statistic > 1 indicates overdispersion relative
  # to Poisson; this is the standard diagnostic (e.g. Cameron and Trivedi
  # 2013) and is what the AIC comparison below is actually testing.
  mu_hat_pois <- fitted(m_pois)
  pearson_resid_sq <- (y - mu_hat_pois)^2 / mu_hat_pois
  dispersion_stat <- sum(pearson_resid_sq) / m_pois$df.residual

  moment_check <- tibble::tibble(
    stage = label,
    mean_count = mean(y),
    pearson_dispersion_statistic = dispersion_stat,
    interpretation = ifelse(dispersion_stat > 1.5,
                            "overdispersed (Pearson statistic > 1.5): NB strongly preferred over Poisson",
                            ifelse(dispersion_stat > 1.05,
                                   "mildly overdispersed (Pearson statistic > 1.05): NB preferred",
                                   ifelse(dispersion_stat < 0.95,
                                          "underdispersed (Pearson statistic < 0.95): NB will not help; consider quasi-Poisson with dispersion<1 or check for upper-bounded/near-deterministic counts",
                                          "near-Poisson (Pearson statistic \u2248 1): NB unlikely to help much")))
  )

  aic_comparison <- tibble::tibble(
    stage = label,
    aic_poisson = AIC(m_pois),
    aic_nb = if (!is.null(m_nb)) AIC(m_nb) else NA_real_,
    delta_aic_nb_minus_poisson = if (!is.null(m_nb)) AIC(m_nb) - AIC(m_pois) else NA_real_,
    nb_theta = if (!is.null(m_nb)) m_nb$theta else NA_real_,
    interpretation = if (!is.null(m_nb) && AIC(m_nb) < AIC(m_pois) - 2) {
      "NB fits decisively better than Poisson (\u0394AIC < -2): keep NB"
    } else if (!is.null(m_nb)) {
      "NB does not improve materially on Poisson by AIC: Poisson may suffice"
    } else {
      "glm.nb() failed to converge; inspect data for separation/sparsity"
    }
  )

  mu_hat_pois <- fitted(m_pois)
  mu_hat_nb   <- if (!is.null(m_nb)) fitted(m_nb) else NA
  theta_hat   <- if (!is.null(m_nb)) m_nb$theta else NA

  p_zero_pois <- mean(dpois(0, lambda = mu_hat_pois))
  p_zero_nb   <- if (!is.null(m_nb)) mean(dnbinom(0, mu = mu_hat_nb, size = theta_hat)) else NA

  zero_check <- tibble::tibble(
    stage = label,
    n_obs = length(y),
    observed_zero_prop = mean(y == 0),
    poisson_predicted_zero_prop = p_zero_pois,
    nb_predicted_zero_prop = p_zero_nb,
    excess_zero_vs_nb = mean(y == 0) - p_zero_nb,
    interpretation = if (!is.na(p_zero_nb) && (mean(y == 0) - p_zero_nb) > 0.05) {
      "observed zeros exceed NB's own prediction by >5pp: consider zero-inflated/hurdle NB"
    } else {
      "NB's implied zero mass is consistent with observed zeros: plain NB adequate"
    }
  )

  list(moment_check = moment_check, aic_comparison = aic_comparison, zero_check = zero_check,
       fitted_poisson = m_pois, fitted_nb = m_nb)
}

# ------------------------------------------------------------
# 2. CONTINUOUS-STAGE DIAGNOSTICS (value / points)
# ------------------------------------------------------------

#' Diagnose the value/points stage: is Tweedie's zero mass doing any work,
#' or does a strictly-positive family fit at least as well?
#'
#' Fits four candidate GLMs with the same offset and intercept-only mean
#' structure (fast proxy for the full GAM, for comparability across
#' families mgcv/stats support differently), and compares them by AIC.
#' If `value` is confirmed strictly positive whenever the offset is
#' positive (the expected case for settlement points/value), the
#' zero-mass rows of this comparison are diagnostic context, not a
#' live concern -- the relevant comparison is purely among the
#' strictly-positive continuous families.
#'
#' @param value Numeric vector of observed value/points (e.g. sum_points).
#' @param offset_var Numeric vector, the offset denominator (n_provision).
#' @param label A short string for this stage, used in printed output.
#' @return A list with:
#'   - zero_mass_check: whether value > 0 whenever offset_var > 0 (confirms
#'     or refutes the premise that a zero-mass family like Tweedie is needed)
#'   - aic_comparison: AIC/BIC for Gamma, inverse Gaussian, Tweedie (p~1.5),
#'     and log-normal (via OLS on log(value)) GLMs with the same offset
#'   - mean_variance_relationship: estimated power-law exponent relating
#'     variance to the mean across offset-binned groups (Tweedie's p
#'     parameter is exactly this exponent; p=1 is Poisson-like, p=2 is
#'     Gamma-like, p=3 is inverse-Gaussian-like) -- read directly off the
#'     data rather than assumed, via Taylor's power law on binned means/variances
diagnose_value_stage <- function(value, offset_var, label = "value stage") {
  stopifnot(length(value) == length(offset_var))
  keep <- offset_var > 0 & !is.na(value) & !is.na(offset_var)
  value <- value[keep]; offset_var <- offset_var[keep]

  zero_mass_check <- tibble::tibble(
    stage = label,
    n_obs = length(value),
    n_offset_positive = sum(offset_var > 0),
    n_value_zero_given_offset_positive = sum(value == 0 & offset_var > 0),
    prop_zero_given_offset_positive = mean(value == 0 & offset_var > 0) ,
    interpretation = if (sum(value == 0 & offset_var > 0) == 0) {
      "value is strictly positive whenever offset > 0: a zero-mass family (Tweedie p<2, or ZI) is NOT needed on distributional grounds; prefer Gamma/inverse-Gaussian/log-normal"
    } else {
      "value has genuine zero mass even when offset > 0: Tweedie (or a two-part hurdle model) is appropriate"
    }
  )

  df_fit <- data.frame(value = value, off = offset_var)
  df_fit_pos <- df_fit[df_fit$value > 0, ]

  m_gamma <- tryCatch(
    glm(value ~ 1 + offset(log(off)), family = Gamma(link = "log"), data = df_fit_pos),
    error = function(e) NULL
  )
  m_ig <- tryCatch(
    glm(value ~ 1 + offset(log(off)), family = inverse.gaussian(link = "log"), data = df_fit_pos),
    error = function(e) NULL
  )
  m_tw <- tryCatch(
    mgcv::gam(value ~ 1 + offset(log(off)), family = tw(link = "log"), data = df_fit_pos, method = "REML"),
    error = function(e) NULL
  )
  # Log-normal proxy: OLS on log(value), offset absorbed as log(off) on the RHS;
  # AIC from an OLS log-response model is on a different scale (log(value) vs
  # value) and is NOT directly comparable to the GLM AICs above without a
  # Jacobian correction; we apply that correction so all four AICs are on the
  # value (not log-value) scale.
  m_lognorm <- tryCatch(lm(log(value) ~ 1 + offset(log(off)), data = df_fit_pos), error = function(e) NULL)
  aic_lognorm_corrected <- if (!is.null(m_lognorm)) {
    AIC(m_lognorm) + 2 * sum(log(df_fit_pos$value))  # Jacobian correction: d(log v)/dv = 1/v
  } else NA_real_

  aic_comparison <- tibble::tibble(
    stage = label,
    family = c("Gamma", "Inverse Gaussian", "Tweedie", "Log-normal (Jacobian-corrected)"),
    aic = c(
      if (!is.null(m_gamma)) AIC(m_gamma) else NA_real_,
      if (!is.null(m_ig)) AIC(m_ig) else NA_real_,
      if (!is.null(m_tw)) AIC(m_tw) else NA_real_,
      aic_lognorm_corrected
    )
  ) %>%
    dplyr::mutate(delta_aic = aic - min(aic, na.rm = TRUE)) %>%
    dplyr::arrange(delta_aic)

  # Taylor's power law: log(group variance) ~ a + p * log(group mean), fit
  # on offset-quantile bins, to read off the implied Tweedie power p directly
  # from the data rather than assuming p = 1.5 as in R/02_models.R's default.
  df_fit_pos$rate <- df_fit_pos$value / df_fit_pos$off
  df_fit_pos$bin <- cut(df_fit_pos$off, breaks = quantile(df_fit_pos$off, probs = seq(0, 1, 0.1)),
                        include.lowest = TRUE)
  bin_stats <- df_fit_pos %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(mean_rate = mean(rate), var_rate = var(rate), n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n >= 5, mean_rate > 0, var_rate > 0)

  power_law_fit <- if (nrow(bin_stats) >= 3) {
    lm(log(var_rate) ~ log(mean_rate), data = bin_stats)
  } else NULL

  # The power-law slope is only informative if the binned means actually
  # vary enough to estimate a slope against (Taylor's power law requires
  # spread in log(mean_rate) across bins). If kappa = value/provision does
  # not depend systematically on the offset (n_provision) -- as is the
  # case, by construction, whenever the value stage's mean is modeled as
  # depending only on sex/age/year/unit and not on provision volume itself
  # -- the binned means will be nearly flat, and the fitted slope is noise,
  # not signal. We flag this explicitly rather than reporting a precise-
  # looking but meaningless coefficient.
  mean_rate_range_ratio <- if (nrow(bin_stats) >= 3) {
    max(bin_stats$mean_rate) / min(bin_stats$mean_rate)
  } else NA_real_
  power_law_reliable <- !is.na(mean_rate_range_ratio) && mean_rate_range_ratio > 1.5

  mean_variance_relationship <- tibble::tibble(
    stage = label,
    implied_tweedie_power_p = if (!is.null(power_law_fit)) coef(power_law_fit)[2] else NA_real_,
    n_bins_used = nrow(bin_stats),
    mean_rate_range_ratio = mean_rate_range_ratio,
    reliable = power_law_reliable,
    interpretation = if (is.null(power_law_fit)) {
      "insufficient bins with positive variance to estimate power law"
    } else if (!power_law_reliable) {
      sprintf(
        paste0("UNRELIABLE: binned means vary by only %.2fx across offset deciles ",
               "(need >1.5x for a stable power-law estimate) -- the rate does not ",
               "depend systematically on the offset in this data, so the fitted ",
               "exponent (%.2f) reflects estimation noise, not a real mean-variance ",
               "relationship; rely on the AIC comparison above instead"),
        mean_rate_range_ratio, coef(power_law_fit)[2]
      )
    } else {
      p_hat <- coef(power_law_fit)[2]
      dplyr::case_when(
        p_hat < 1.2 ~ "p close to 1: Poisson-like variance, consider quasi-Poisson-style GLM",
        p_hat >= 1.2 & p_hat < 1.8 ~ "p close to 1.5: compound-Poisson-Gamma (Tweedie) variance structure supported, IF zero mass is also present",
        p_hat >= 1.8 & p_hat < 2.3 ~ "p close to 2: Gamma variance structure (Var \u221d mean^2) supported",
        p_hat >= 2.3 ~ "p > 2.3: heavier-tailed than Gamma; inverse Gaussian or log-normal may fit better"
      )
    }
  )

  list(
    zero_mass_check = zero_mass_check,
    aic_comparison = aic_comparison,
    mean_variance_relationship = mean_variance_relationship,
    fitted_gamma = m_gamma, fitted_ig = m_ig, fitted_tweedie = m_tw, fitted_lognorm = m_lognorm
  )
}

# ------------------------------------------------------------
# 3. RUN ALL THREE STAGES AND PRINT A READABLE REPORT
# ------------------------------------------------------------

#' Run the full distribution-diagnostic suite on one disease group's
#' aggregated cell-level data (the same data passed to fit_model_chain()).
#'
#' @param dat Output of filter_and_aggregate() (one disease group, all years).
#' @return A list with the three diagnose_*() outputs (patients, provision,
#'   value), each containing its own sub-tables; also printed to the console
#'   as a readable report.
run_distribution_diagnostics <- function(dat) {
  message("==================================================")
  message("DISTRIBUTION DIAGNOSTICS")
  message("==================================================\n")

  message("---- Stage 1: patients ~ offset(pop) ----")
  diag_patients <- diagnose_count_stage(dat$n_patients, dat$pop, label = "patients")
  print(diag_patients$moment_check)
  print(diag_patients$aic_comparison %>% dplyr::select(-stage))
  print(diag_patients$zero_check %>% dplyr::select(-stage))

  message("\n---- Stage 2: provision ~ offset(n_patients) ----")
  diag_provision <- diagnose_count_stage(dat$n_provision, dat$n_patients, label = "provision")
  print(diag_provision$moment_check)
  print(diag_provision$aic_comparison %>% dplyr::select(-stage))
  print(diag_provision$zero_check %>% dplyr::select(-stage))

  message("\n---- Stage 3: value ~ offset(n_provision) ----")
  diag_value <- diagnose_value_stage(dat$value, dat$n_provision, label = "value")
  print(diag_value$zero_mass_check %>% dplyr::select(-stage))
  print(diag_value$aic_comparison %>% dplyr::select(-stage))
  print(diag_value$mean_variance_relationship %>% dplyr::select(-stage))

  message("\n==================================================")
  message("END OF REPORT. See output/tables/distribution_diagnostics_*.csv for")
  message("the full tables, and update R/02_models.R's build_formula() families")
  message("if a stage's recommendation differs from the current default")
  message("(NB for patients/provision, Tweedie for value).")
  message("==================================================")

  list(patients = diag_patients, provision = diag_provision, value = diag_value)
}
