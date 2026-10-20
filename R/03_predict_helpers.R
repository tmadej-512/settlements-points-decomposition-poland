# ============================================================
# R/03_predict_helpers.R
# ------------------------------------------------------------
# Shared machinery for building counterfactual `newdata` frames
# and chaining predictions through the three fitted models
# (patients -> provision -> value).
#
# Every decomposition method (naive, OFAT, sequential, Shapley)
# is built ONLY on top of `predict_scenario()` below. This keeps
# the four methods directly comparable: they differ only in which
# (year_pat, year_prov, year_val, pop_vec) combinations they query,
# never in how a single counterfactual prediction is computed.
# ============================================================

suppressPackageStartupMessages(library(dplyr))

#' Build newdata for the patient-count model with a chosen year and pop level
mk_newdata_pat <- function(df, year_scalar, pop_vec) {
  yr <- if (length(year_scalar) == 1) rep(year_scalar, nrow(df)) else year_scalar
  df %>% mutate(
    year     = yr,
    year_num = as.numeric(yr),
    sex_age  = droplevels(interaction(sex, age_group, drop = TRUE)),
    pop      = pop_vec
  )
}

#' Build newdata for the provision-count model with a chosen year and patient level
mk_newdata_prov <- function(df, year_scalar, pat_vec) {
  yr <- if (length(year_scalar) == 1) rep(year_scalar, nrow(df)) else year_scalar
  df %>% mutate(
    year       = yr,
    year_num   = as.numeric(yr),
    sex_age    = droplevels(interaction(sex, age_group, drop = TRUE)),
    n_patients = pat_vec
  )
}

#' Build newdata for the value model with a chosen year and provision level
mk_newdata_val <- function(df, year_scalar, prov_vec) {
  yr <- if (length(year_scalar) == 1) rep(year_scalar, nrow(df)) else year_scalar
  df %>% mutate(
    year        = yr,
    year_num    = as.numeric(yr),
    sex_age     = droplevels(interaction(sex, age_group, drop = TRUE)),
    n_provision = prov_vec
  )
}

#' Predict a single counterfactual scenario's total value
#'
#' Runs the full chain patients -> provision -> value, where each
#' stage can independently be evaluated "as if" it were a different
#' year (year_pat / year_prov / year_val), and the population vector
#' feeding the first stage can independently be swapped (pop_vec).
#' This four-way independent control is what lets every decomposition
#' method below be expressed as a query into the same function.
#'
#' @param df_target The row-level data frame the scenario is evaluated on
#'   (i.e. which administrative cells / sex / age groups are summed).
#' @param fit Output of fit_model_chain().
#' @param year_pat Year to plug into the patient model (controls pi).
#' @param year_prov Year to plug into the provision model (controls lambda).
#' @param year_val Year to plug into the value model (controls kappa).
#' @param pop_vec Population vector to plug in (controls N).
#' @return A single numeric: sum of predicted value across df_target's rows.
predict_scenario <- function(df_target, fit, year_pat, year_prov, year_val, pop_vec) {
  c_hat <- as.numeric(predict(
    fit$m_pat, newdata = mk_newdata_pat(df_target, year_pat, pop_vec), type = "response"
  ))
  q_hat <- as.numeric(predict(
    fit$m_prov, newdata = mk_newdata_prov(df_target, year_prov, c_hat), type = "response"
  ))
  v_hat <- as.numeric(predict(
    fit$m_val, newdata = mk_newdata_val(df_target, year_val, q_hat), type = "response"
  ))
  v_hat <- ifelse(q_hat <= 0, 0, v_hat)
  sum(v_hat, na.rm = TRUE)
}

#' Predict the full chain at each row's own observed year/pop (the "actual" fit)
predict_chain_actual <- function(df_target, fit) {
  c_hat <- as.numeric(predict(
    fit$m_pat, newdata = mk_newdata_pat(df_target, df_target$year, df_target$pop), type = "response"
  ))
  q_hat <- as.numeric(predict(
    fit$m_prov, newdata = mk_newdata_prov(df_target, df_target$year, c_hat), type = "response"
  ))
  v_hat <- as.numeric(predict(
    fit$m_val, newdata = mk_newdata_val(df_target, df_target$year, q_hat), type = "response"
  ))
  v_hat <- ifelse(q_hat <= 0, 0, v_hat)
  tibble::tibble(c_hat = c_hat, q_hat = q_hat, v_hat = v_hat)
}

#' Prepare the (target, base) pair of row-level frames for a 2-year comparison
#'
#' `target` gets a `pop_base` column attached (the baseline-year population
#' for the same code/sex/age_group cell), needed by every counterfactual
#' scenario that "switches off" the N component.
#'
#' @param df_2y A data frame containing exactly two distinct years.
#' @param baseline_year,comparison_year The two years to compare.
#' @return list(target = ..., base = ...)
prep_pair <- function(df_2y, baseline_year, comparison_year) {
  keys <- c("code", "sex", "age_group")
  base_map <- df_2y %>%
    filter(year == baseline_year) %>%
    dplyr::select(all_of(keys), pop_base = pop)

  target_df <- df_2y %>%
    filter(year == comparison_year) %>%
    # Simulation DGP frames may already carry a helper column named
    # `pop_base`. Drop it before joining the baseline-year population map,
    # otherwise dplyr creates pop_base.x / pop_base.y and the bare
    # `pop_base` column is no longer available in mutate().
    dplyr::select(-dplyr::any_of("pop_base")) %>%
    left_join(base_map, by = keys) %>%
    mutate(pop_base = dplyr::coalesce(.data$pop_base, 0))

  base_df <- df_2y %>% filter(year == baseline_year)
  list(target = target_df, base = base_df)
}

#' Evaluate predict_scenario() for a named subset of "active" (switched-on) components
#'
#' This is the building block shared by the sequential and Shapley methods:
#' both work by toggling subsets of {N, PI, LAMBDA, KAPPA} on or off and
#' reading off predict_scenario() for that subset.
#'
#' @param tgt,bas Output of prep_pair().
#' @param fit Model chain.
#' @param baseline_year,comparison_year The two years being compared.
#' @param active_set Character vector, subset of c("N","PI","LAMBDA","KAPPA"),
#'   naming which components are evaluated at comparison_year levels.
#'   Components NOT in active_set are held at baseline_year levels.
scenario_from_set <- function(tgt, bas, fit, baseline_year, comparison_year, active_set) {
  use_N  <- "N"      %in% active_set
  use_PI <- "PI"     %in% active_set
  use_LA <- "LAMBDA" %in% active_set
  use_KA <- "KAPPA"  %in% active_set

  predict_scenario(
    df_target = if (use_N) tgt else bas,
    fit       = fit,
    year_pat  = if (use_PI) comparison_year else baseline_year,
    year_prov = if (use_LA) comparison_year else baseline_year,
    year_val  = if (use_KA) comparison_year else baseline_year,
    pop_vec   = if (use_N) tgt$pop else bas$pop
  )
}
