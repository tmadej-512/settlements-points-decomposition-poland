# ============================================================
# R/10e_oracle_decomposition.R
# ------------------------------------------------------------
# Oracle validation for the simulation study.
#
# This scenario uses the known deterministic DGP value function directly.
# No GAM/GLM model is fitted and no noisy observed counts are used for the
# model-based decompositions. It answers the clean attribution question:
#
#   If the counterfactual value function is known exactly, does the
#   decomposition recover the known counterfactual allocation?
#
# In this setting, the Shapley estimate is mathematically identical to the
# known Shapley ground truth up to numerical precision. The OFAT, sequential,
# and naive rows are retained as comparator rules, not as estimators of the
# Shapley target.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

oracle_state_table <- function(sim, baseline_year, comparison_year) {
  required_cols <- c(
    "code", "sex", "age_group", "year",
    "pop_mean", "pi_mean", "lambda_mean", "kappa_mean"
  )
  missing_cols <- setdiff(required_cols, names(sim$data))
  if (length(missing_cols) > 0) {
    stop(
      "oracle_state_table() requires deterministic DGP columns missing from sim$data: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  cell_keys <- c("code", "sex", "age_group")

  base <- sim$data %>%
    dplyr::filter(.data$year == baseline_year) %>%
    dplyr::select(
      dplyr::all_of(cell_keys),
      N_base = .data$pop_mean,
      PI_base = .data$pi_mean,
      LAMBDA_base = .data$lambda_mean,
      KAPPA_base = .data$kappa_mean
    )

  comp <- sim$data %>%
    dplyr::filter(.data$year == comparison_year) %>%
    dplyr::select(
      dplyr::all_of(cell_keys),
      N_comp = .data$pop_mean,
      PI_comp = .data$pi_mean,
      LAMBDA_comp = .data$lambda_mean,
      KAPPA_comp = .data$kappa_mean
    )

  state_tbl <- dplyr::inner_join(base, comp, by = cell_keys)
  if (nrow(state_tbl) == 0) {
    stop("No matched baseline/comparison cells found for oracle decomposition.", call. = FALSE)
  }

  state_tbl
}

oracle_expected_total <- function(state_tbl, active = character(0)) {
  use <- function(nm) nm %in% active

  N      <- if (use("N"))      state_tbl$N_comp      else state_tbl$N_base
  PI     <- if (use("PI"))     state_tbl$PI_comp     else state_tbl$PI_base
  LAMBDA <- if (use("LAMBDA")) state_tbl$LAMBDA_comp else state_tbl$LAMBDA_base
  KAPPA  <- if (use("KAPPA"))  state_tbl$KAPPA_comp  else state_tbl$KAPPA_base

  sum(N * PI * LAMBDA * KAPPA, na.rm = TRUE)
}

oracle_latent_panel <- function(sim, baseline_year, comparison_year) {
  sim$data %>%
    dplyr::filter(.data$year %in% c(baseline_year, comparison_year)) %>%
    dplyr::transmute(
      icd10 = "SIMULATED",
      code = .data$code,
      powiat = .data$code,
      x = NA_real_,
      y = NA_real_,
      sex = .data$sex,
      age_group = .data$age_group,
      year = as.integer(.data$year),
      pop = .data$pop_mean,
      n_patients = .data$pop_mean * .data$pi_mean,
      n_provision = .data$pop_mean * .data$pi_mean * .data$lambda_mean,
      value = .data$pop_mean * .data$pi_mean * .data$lambda_mean * .data$kappa_mean
    )
}

oracle_decompose_ofat <- function(state_tbl, baseline_year, comparison_year) {
  comps <- c("N", "PI", "LAMBDA", "KAPPA")
  V0 <- oracle_expected_total(state_tbl, character(0))
  V4 <- oracle_expected_total(state_tbl, comps)
  effects <- stats::setNames(vapply(comps, function(z) {
    oracle_expected_total(state_tbl, z) - V0
  }, numeric(1)), comps)

  tibble::tibble(
    baseline_year = baseline_year,
    comparison_year = comparison_year,
    baseline_total = V0,
    comparison_total = V4,
    total_change = V4 - V0,
    N_effect = unname(effects["N"]),
    PI_effect = unname(effects["PI"]),
    LAMBDA_effect = unname(effects["LAMBDA"]),
    KAPPA_effect = unname(effects["KAPPA"]),
    residual = (V4 - V0) - sum(effects)
  )
}

oracle_decompose_sequential <- function(state_tbl,
                                        baseline_year,
                                        comparison_year,
                                        order = c("N", "PI", "LAMBDA", "KAPPA")) {
  stopifnot(setequal(order, c("N", "PI", "LAMBDA", "KAPPA")))

  active <- character(0)
  V0 <- oracle_expected_total(state_tbl, active)
  V_prev <- V0
  effects <- stats::setNames(numeric(4), c("N", "PI", "LAMBDA", "KAPPA"))

  for (z in order) {
    active <- c(active, z)
    V_curr <- oracle_expected_total(state_tbl, active)
    effects[z] <- V_curr - V_prev
    V_prev <- V_curr
  }

  V4 <- oracle_expected_total(state_tbl, c("N", "PI", "LAMBDA", "KAPPA"))

  tibble::tibble(
    baseline_year = baseline_year,
    comparison_year = comparison_year,
    baseline_total = V0,
    comparison_total = V4,
    total_change = V4 - V0,
    N_effect = unname(effects["N"]),
    PI_effect = unname(effects["PI"]),
    LAMBDA_effect = unname(effects["LAMBDA"]),
    KAPPA_effect = unname(effects["KAPPA"]),
    residual = (V4 - V0) - sum(effects),
    sequential_order = paste(order, collapse = " -> ")
  )
}

oracle_decompose_shapley <- function(state_tbl, baseline_year, comparison_year) {
  comps <- c("N", "PI", "LAMBDA", "KAPPA")
  perms <- combinat_permutations(comps)
  contrib_mat <- matrix(0, nrow = length(perms), ncol = length(comps), dimnames = list(NULL, comps))

  for (i in seq_along(perms)) {
    active <- character(0)
    V_prev <- oracle_expected_total(state_tbl, active)
    for (z in perms[[i]]) {
      active <- c(active, z)
      V_curr <- oracle_expected_total(state_tbl, active)
      contrib_mat[i, z] <- V_curr - V_prev
      V_prev <- V_curr
    }
  }

  sv <- colMeans(contrib_mat)
  V0 <- oracle_expected_total(state_tbl, character(0))
  V4 <- oracle_expected_total(state_tbl, comps)

  out <- tibble::tibble(
    baseline_year = baseline_year,
    comparison_year = comparison_year,
    baseline_total = V0,
    comparison_total = V4,
    total_change = V4 - V0,
    N_effect = unname(sv["N"]),
    PI_effect = unname(sv["PI"]),
    LAMBDA_effect = unname(sv["LAMBDA"]),
    KAPPA_effect = unname(sv["KAPPA"]),
    residual = (V4 - V0) - sum(sv)
  )
  attr(out, "contrib_matrix") <- contrib_mat
  out
}

oracle_decompose_one_panel <- function(sim,
                                       baseline_year,
                                       comparison_year,
                                       methods = c("naive", "ofat", "sequential", "shapley")) {
  state_tbl <- oracle_state_table(sim, baseline_year, comparison_year)
  latent_df <- oracle_latent_panel(sim, baseline_year, comparison_year)

  rows <- list()

  if ("naive" %in% methods) {
    rows$naive <- decompose_naive_2y(latent_df, NULL, baseline_year, comparison_year) %>%
      dplyr::mutate(method = "naive")
  }

  if ("ofat" %in% methods) {
    rows$ofat <- oracle_decompose_ofat(state_tbl, baseline_year, comparison_year) %>%
      dplyr::mutate(method = "ofat")
  }

  if ("sequential" %in% methods) {
    rows$sequential <- oracle_decompose_sequential(state_tbl, baseline_year, comparison_year) %>%
      dplyr::mutate(method = "sequential")
  }

  if ("shapley" %in% methods) {
    rows$shapley <- oracle_decompose_shapley(state_tbl, baseline_year, comparison_year) %>%
      dplyr::mutate(method = "shapley")
  }

  dplyr::bind_rows(rows)
}

run_oracle_decomposition_study <- function(true_growth = list(N = 0.01, PI = 0.03, LAMBDA = -0.01, KAPPA = 0.02),
                                           baseline_year = 2018,
                                           comparison_year = 2024,
                                           years = 2018:2024,
                                           n_units = 80,
                                           sex_levels = c("F", "M"),
                                           age_levels = c("0-17", "18-44", "45-64", "65+"),
                                           heterogeneous_slopes = TRUE,
                                           regional_random_slopes = FALSE,
                                           regional_slope_sd = list(N = 0, PI = 0, LAMBDA = 0, KAPPA = 0),
                                           pandemic_shock = FALSE,
                                           noise_scale = 1,
                                           reps = 200,
                                           seed = 1,
                                           methods = c("naive", "ofat", "sequential", "shapley")) {
  pb <- txtProgressBar(min = 0, max = reps, style = 3)
  out <- vector("list", reps)

  for (r in seq_len(reps)) {
    sim <- simulate_disease_panel(
      n_units = n_units,
      years = years,
      sex_levels = sex_levels,
      age_levels = age_levels,
      true_growth = true_growth,
      heterogeneous_slopes = heterogeneous_slopes,
      regional_random_slopes = regional_random_slopes,
      regional_slope_sd = regional_slope_sd,
      pandemic_shock = pandemic_shock,
      noise_scale = noise_scale,
      seed = seed + r
    )

    truth <- truth_from_panel(sim, baseline_year, comparison_year)

    rep_tbl <- oracle_decompose_one_panel(
      sim = sim,
      baseline_year = baseline_year,
      comparison_year = comparison_year,
      methods = methods
    )

    rep_tbl$rep_id <- r
    rep_tbl$truth_N <- truth$N_effect
    rep_tbl$truth_PI <- truth$PI_effect
    rep_tbl$truth_LAMBDA <- truth$LAMBDA_effect
    rep_tbl$truth_KAPPA <- truth$KAPPA_effect

    out[[r]] <- rep_tbl
    setTxtProgressBar(pb, r)
  }

  close(pb)

  raw <- dplyr::bind_rows(out) %>%
    tidyr::pivot_longer(
      cols = c(N_effect, PI_effect, LAMBDA_effect, KAPPA_effect),
      names_to = "component",
      values_to = "estimate"
    ) %>%
    dplyr::mutate(
      component = sub("_effect$", "", .data$component),
      truth = dplyr::case_when(
        .data$component == "N" ~ .data$truth_N,
        .data$component == "PI" ~ .data$truth_PI,
        .data$component == "LAMBDA" ~ .data$truth_LAMBDA,
        .data$component == "KAPPA" ~ .data$truth_KAPPA,
        TRUE ~ NA_real_
      ),
      error = .data$estimate - .data$truth
    ) %>%
    dplyr::select(rep_id, method, component, estimate, truth, error)

  summary_tbl <- raw %>%
    dplyr::group_by(.data$method, .data$component) %>%
    dplyr::summarise(
      truth = mean(.data$truth, na.rm = TRUE),
      bias = mean(.data$error, na.rm = TRUE),
      rmse = sqrt(mean(.data$error^2, na.rm = TRUE)),
      sd_error = stats::sd(.data$error, na.rm = TRUE),
      n_reps_ok = sum(!is.na(.data$error)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      bias_pct_of_truth = dplyr::if_else(abs(.data$truth) > 0, .data$bias / abs(.data$truth), NA_real_),
      rmse_pct_of_truth = dplyr::if_else(abs(.data$truth) > 0, .data$rmse / abs(.data$truth), NA_real_)
    ) %>%
    dplyr::group_by(.data$component) %>%
    dplyr::mutate(rmse_rank_within_component = dplyr::min_rank(.data$rmse)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$component, .data$method)

  list(raw = raw, summary = summary_tbl)
}
