# ============================================================
# R/10_simulation_ground_truth.R
# ------------------------------------------------------------
# Simulation study: the strongest validation evidence in the paper.
#
# We generate synthetic data from a known data-generating process
# (DGP) where the TRUE contribution of each of N, pi, lambda, kappa
# to the total change is known exactly by construction (we set it
# ourselves when writing the DGP). We then run the full pipeline --
# fit the GAM/GLM model chain on the synthetic data exactly as we
# would on real NFZ data, then apply naive / OFAT / sequential /
# Shapley decomposition -- and check whether each method recovers
# the known truth.
#
# This answers a question the axiom checks (R/08) and the real-data
# efficiency check cannot: not "is the allocation rule mathematically
# self-consistent" but "does the whole estimation pipeline (model
# fitting + counterfactual prediction + attribution) recover the
# right answer when we already know what the right answer is".
#
# DESIGN OF THE DGP:
# For a single synthetic "disease", over `n_years` years and a panel
# of `n_units` administrative units x sex x age_group cells, we
# construct:
#   pop_it                 : grows smoothly over time (the N driver)
#   pi_it = patients/pop          : grows/shrinks at a chosen rate (pi driver)
#   lambda_it = provision/patients: grows/shrinks at a chosen rate (lambda driver)
#   kappa_it = value/provision    : grows/shrinks at a chosen rate (kappa driver)
# Each ratio's time trend is a smooth (sinusoidal-perturbed) function
# with sex/age_group-specific levels, age-sex-specific slopes, and optional
# unit-level random slopes, so the fitted GAM has genuine heterogeneous
# temporal structure to recover -- this is NOT a linear toy case.
# Counts are then drawn from the same families used for real fitting
# (negative binomial for patients/provision, Tweedie-like for value),
# so the simulation also stresses the same observation-level noise the
# real pipeline has to contend with.
#
# GROUND TRUTH CONTRIBUTION:
# Because the DGP is multiplicative in expectation (E[value] =
# pop * pi * lambda * kappa, summed over sex/age/code cells), the
# TRUE Shapley decomposition of E[value] has a closed form computed
# directly from the deterministic DGP means retained in the simulated panel
# (see truth_from_panel()), with NO sampling noise. This is the target the four methods are judged
# against after refitting on the noisy simulated counts -- isolating
# "estimation + attribution error" from "have we defined the right
# ground truth".
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# Helpers for the synthetic data-generating process (DGP)
# -----------------------------------------------------------------------------

simulation_cell_baseline <- function(sex_levels, age_levels) {
  expand.grid(sex = sex_levels, age_group = age_levels, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    mutate(
      age_old = grepl("65|75", age_group),
      age_very_old = grepl("75", age_group),
      age_child = grepl("0-4|5-12|13-17|0-17", age_group),
      age_mid = grepl("45-64|30-44", age_group),
      sex_m = sex == "M",
      pop_base    = 4200 + 1300 * age_old + 900 * age_mid - 900 * age_child - 600 * age_very_old,
      pi_base     = 0.035 + 0.035 * age_mid + 0.070 * age_old + 0.015 * sex_m,
      lambda_base = 2.4 + 1.2 * age_mid + 2.2 * age_old + 0.3 * sex_m,
      kappa_base  = 85 + 18 * age_old + 12 * sex_m,
      pop_base    = pmax(pop_base, 800),
      pi_base     = pmax(pi_base, 0.005),
      lambda_base = pmax(lambda_base, 0.5),
      kappa_base  = pmax(kappa_base, 5)
    ) %>%
    dplyr::select(-age_old, -age_very_old, -age_child, -age_mid, -sex_m)
}

add_simulation_component_trajectories <- function(dat,
                                                  true_growth = list(N = 0.01, PI = 0.03, LAMBDA = -0.01, KAPPA = 0.02),
                                                  heterogeneous_slopes = TRUE) {
  # Backward-compatible defaults: if older test fixtures or callers do not
  # contain unit-specific slope columns, they are treated as zero.
  if (!"N_unit_slope" %in% names(dat)) dat$N_unit_slope <- 0
  if (!"PI_unit_slope" %in% names(dat)) dat$PI_unit_slope <- 0
  if (!"LAMBDA_unit_slope" %in% names(dat)) dat$LAMBDA_unit_slope <- 0
  if (!"KAPPA_unit_slope" %in% names(dat)) dat$KAPPA_unit_slope <- 0

  dat %>%
    mutate(
      age_old = grepl("65|75", age_group),
      age_very_old = grepl("75", age_group),
      age_child = grepl("0-4|5-12|13-17|0-17", age_group),
      age_mid = grepl("45-64|30-44", age_group),
      sex_m = sex == "M",

      # Baseline version: all sex-age cells share a common component slope.
      # Heterogeneous version: slopes differ by age and sex. Unit random
      # slopes then add region-specific trends around those age-sex trends.
      N_slope = true_growth$N + ifelse(isTRUE(heterogeneous_slopes),
        0.004 * age_old - 0.003 * age_child - 0.002 * age_very_old + 0.001 * sex_m, 0) + N_unit_slope,
      PI_slope = true_growth$PI + ifelse(isTRUE(heterogeneous_slopes),
        0.012 * age_old + 0.004 * sex_m - 0.006 * age_child, 0) + PI_unit_slope,
      LAMBDA_slope = true_growth$LAMBDA + ifelse(isTRUE(heterogeneous_slopes),
        0.011 * age_old + 0.004 * age_mid - 0.003 * sex_m - 0.004 * age_child, 0) + LAMBDA_unit_slope,
      KAPPA_slope = true_growth$KAPPA + ifelse(isTRUE(heterogeneous_slopes),
        -0.006 * age_old + 0.004 * sex_m + 0.002 * age_child, 0) + KAPPA_unit_slope,

      N_wiggle = 0.10 * sin(t / 2) + ifelse(isTRUE(heterogeneous_slopes),
        0.030 * age_old * sin(t / 2.5) - 0.020 * age_child * cos(t / 3), 0),
      PI_wiggle = 0.05 * sin(t / 3) + ifelse(isTRUE(heterogeneous_slopes),
        0.020 * age_old * cos(t / 2) + 0.010 * sex_m * sin(t / 2.7), 0),
      LAMBDA_wiggle = 0.05 * cos(t / 3) + ifelse(isTRUE(heterogeneous_slopes),
        0.020 * sex_m * sin(t / 2) - 0.015 * age_child * cos(t / 2.5), 0),
      KAPPA_wiggle = 0.05 * sin(t / 4) + ifelse(isTRUE(heterogeneous_slopes),
        -0.015 * age_old * sin(t / 3) + 0.010 * sex_m * cos(t / 3), 0),

      N_mult      = exp(N_slope      * t + N_wiggle + unit_re),
      PI_mult     = exp(PI_slope     * t + PI_wiggle),
      LAMBDA_mult = exp(LAMBDA_slope * t + LAMBDA_wiggle),
      KAPPA_mult  = exp(KAPPA_slope  * t + KAPPA_wiggle),
      pop_mean    = pop_base * N_mult,
      pi_mean     = pi_base * PI_mult,
      lambda_mean = lambda_base * LAMBDA_mult,
      kappa_mean  = kappa_base * KAPPA_mult
    )
}

#' Generate one synthetic disease panel with known component trajectories
#'
#' @param n_units Number of administrative units (e.g. powiats).
#' @param years Integer vector of years, e.g. 2018:2023.
#' @param sex_levels,age_levels Factor levels to cross with units.
#' @param true_growth A named list with the *baseline annual log-growth rate*
#'   of each driver: list(N = ..., PI = ..., LAMBDA = ..., KAPPA = ...).
#' @param heterogeneous_slopes If TRUE, the DGP adds age-sex-specific slope
#'   modifiers and nonlinear wiggles for N, PI, LAMBDA, and KAPPA.
#' @param regional_random_slopes If TRUE, administrative units receive random
#'   slope deviations for each component, creating region-specific trends.
#' @param regional_slope_sd Named list of SDs for the unit-level random slopes.
#' @param unit_re_sd Standard deviation of the unit-level random effect.
#' @param overdispersion_patients,overdispersion_provision Negative-binomial
#'   size parameters (smaller = more overdispersed) for the two count stages.
#' @param tweedie_power,tweedie_phi Tweedie shape/dispersion for the value stage.
#' @param seed RNG seed.
#' @return A list with `data`, `true_growth`, and `dgp_params`.
simulate_disease_panel <- function(n_units = 30,
                                    years = 2018:2023,
                                    sex_levels = c("F", "M"),
                                    age_levels = c("0-17", "18-44", "45-64", "65+"),
                                    true_growth = list(N = 0.01, PI = 0.03, LAMBDA = -0.01, KAPPA = 0.02),
                                    heterogeneous_slopes = TRUE,
                                    regional_random_slopes = TRUE,
                                    regional_slope_sd = list(N = 0.003, PI = 0.008, LAMBDA = 0.008, KAPPA = 0.006),
                                    pandemic_shock = FALSE,
                                    noise_scale = 1,
                                    composition_shift = FALSE,
                                    composition_shift_strength = 0.25,
                                    misspecification_stress = FALSE,
                                    misspecification_strength = 0.25,
                                    unit_re_sd = 0.3,
                                    overdispersion_patients = 8,
                                    overdispersion_provision = 6,
                                    tweedie_power = 1.5,
                                    tweedie_phi = 2,
                                    seed = 1) {
  set.seed(seed)
  t0 <- min(years)

  grid <- expand.grid(
    code = paste0("U", seq_len(n_units)),
    sex = sex_levels,
    age_group = age_levels,
    year = years,
    stringsAsFactors = FALSE
  ) %>% as_tibble()

  slope_sd <- modifyList(list(N = 0.003, PI = 0.008, LAMBDA = 0.008, KAPPA = 0.006), regional_slope_sd)
  unit_effects <- tibble(
    code = paste0("U", seq_len(n_units)),
    unit_re = rnorm(n_units, mean = 0, sd = unit_re_sd),
    N_unit_slope = if (isTRUE(regional_random_slopes)) rnorm(n_units, 0, slope_sd$N) else 0,
    PI_unit_slope = if (isTRUE(regional_random_slopes)) rnorm(n_units, 0, slope_sd$PI) else 0,
    LAMBDA_unit_slope = if (isTRUE(regional_random_slopes)) rnorm(n_units, 0, slope_sd$LAMBDA) else 0,
    KAPPA_unit_slope = if (isTRUE(regional_random_slopes)) rnorm(n_units, 0, slope_sd$KAPPA) else 0
  )

  cell_base <- simulation_cell_baseline(sex_levels, age_levels)

  dat <- grid %>%
    left_join(unit_effects, by = "code") %>%
    left_join(cell_base, by = c("sex", "age_group")) %>%
    mutate(t = year - t0) %>%
    add_simulation_component_trajectories(
      true_growth = true_growth,
      heterogeneous_slopes = heterogeneous_slopes
    )

  # Optional temporary service disruption used only in stress-test simulations.
  # It affects intermediate-year fitting difficulty, not the definition of the
  # endpoint accounting identity. If the comparison year is post-shock, the
  # fitted model still has to learn through the disrupted years.
  if (isTRUE(pandemic_shock)) {
    dat <- dat %>%
      dplyr::mutate(
        pandemic_lambda_mult = dplyr::case_when(
          .data$year == 2020L ~ 0.72,
          .data$year == 2021L ~ 0.88,
          TRUE ~ 1.00
        ),
        lambda_mean = .data$lambda_mean * .data$pandemic_lambda_mult
      )
  } else {
    dat$pandemic_lambda_mult <- 1
  }

  # Explicit cell-composition shift: population mass moves over time toward
  # older/high-use/high-value cells. Within each administrative unit and year,
  # total population is preserved, but its age-sex mix changes. This is the
  # key stress test for aggregate LMDI-style decomposition, because aggregate
  # ratios can remain plausible while cell-level attribution changes.
  if (isTRUE(composition_shift)) {
    max_t <- max(dat$t, na.rm = TRUE)
    dat <- dat %>%
      dplyr::mutate(
        t_scaled = ifelse(max_t > 0, .data$t / max_t, 0),
        composition_score = as.numeric(scale(
          0.85 * grepl("65|75", .data$age_group) +
            0.35 * (.data$sex == "M") +
            0.25 * .data$unit_re
        )),
        pop_mean_pre_shift = .data$pop_mean,
        composition_mult = exp(composition_shift_strength * .data$t_scaled * .data$composition_score),
        pop_mean_shifted_raw = .data$pop_mean_pre_shift * .data$composition_mult
      ) %>%
      dplyr::group_by(.data$code, .data$year) %>%
      dplyr::mutate(
        pop_mean = .data$pop_mean_shifted_raw *
          sum(.data$pop_mean_pre_shift, na.rm = TRUE) /
          pmax(sum(.data$pop_mean_shifted_raw, na.rm = TRUE), 1e-9)
      ) %>%
      dplyr::ungroup()
  } else {
    dat$composition_score <- 0
    dat$composition_mult <- 1
  }

  # Deliberately unmodelled nonlinear/interaction shock. The model chain can
  # fit smooth time and broad age-sex structure, but this term creates a
  # localised age-by-sex-by-time deviation that is harder for the simulation-safe
  # formula to recover. This is an honest stress test for model dependence.
  if (isTRUE(misspecification_stress)) {
    max_t <- max(dat$t, na.rm = TRUE)
    mid_t <- max_t / 2
    dat <- dat %>%
      dplyr::mutate(
        t_scaled = ifelse(max_t > 0, .data$t / max_t, 0),
        shock_shape = exp(-0.5 * ((.data$t - mid_t) / pmax(max_t / 5, 1e-6))^2),
        high_risk_cell = as.numeric(grepl("65|75", .data$age_group) & .data$sex == "M"),
        misspecification_mult_pi = exp(misspecification_strength * .data$shock_shape * .data$high_risk_cell),
        misspecification_mult_lambda = exp(-0.55 * misspecification_strength * .data$shock_shape * .data$high_risk_cell),
        misspecification_mult_kappa = exp(0.75 * misspecification_strength * .data$t_scaled * .data$high_risk_cell),
        pi_mean = pmin(0.95, pmax(1e-6, .data$pi_mean * .data$misspecification_mult_pi)),
        lambda_mean = pmax(1e-6, .data$lambda_mean * .data$misspecification_mult_lambda),
        kappa_mean = pmax(1e-6, .data$kappa_mean * .data$misspecification_mult_kappa)
      )
  } else {
    dat$misspecification_mult_pi <- 1
    dat$misspecification_mult_lambda <- 1
    dat$misspecification_mult_kappa <- 1
  }

  # noise_scale > 1 deliberately makes the stress-test DGP noisier.
  # For negative-binomial draws, smaller size means more overdispersion.
  noise_scale <- max(as.numeric(noise_scale %||% 1), 1e-6)
  patient_size <- max(overdispersion_patients / (noise_scale^2), 0.25)
  provision_size <- max(overdispersion_provision / (noise_scale^2), 0.25)
  value_phi <- tweedie_phi * (noise_scale^2)

  # Draw observed counts from the same broad families used at fitting time.
  dat <- dat %>%
    mutate(
      pop = round(pop_mean),
      patients_mean  = pmax(pi_mean * pop, 1e-6),
      n_patients     = rnbinom(n(), mu = patients_mean, size = patient_size),
      provision_mean = pmax(lambda_mean * pmax(n_patients, 1), 1e-6),
      n_provision    = rnbinom(n(), mu = provision_mean, size = provision_size),
      value_mean     = pmax(kappa_mean * pmax(n_provision, 1), 1e-6)
    )

  dat$value <- simulate_tweedie_like(dat$value_mean, dat$n_provision,
                                      power = tweedie_power, phi = value_phi)

  list(
    data = dat %>%
      transmute(
        icd10 = "SIMULATED",
        code = factor(code), powiat = code, x = NA_real_, y = NA_real_,
        sex = factor(sex, levels = sex_levels),
        age_group = factor(age_group, levels = age_levels),
        year = as.integer(year),
        pop, n_patients, n_provision, value,
        # deterministic means and rates are retained for DGP diagnostics;
        # downstream model/decomposition functions ignore extra columns.
        pop_mean, pi_mean, lambda_mean, kappa_mean,
        patients_mean, provision_mean, value_mean,
        pop_base, pi_base, lambda_base, kappa_base,
        unit_re, N_unit_slope, PI_unit_slope, LAMBDA_unit_slope, KAPPA_unit_slope,
        N_slope, PI_slope, LAMBDA_slope, KAPPA_slope,
        pandemic_lambda_mult, composition_score, composition_mult,
        misspecification_mult_pi, misspecification_mult_lambda, misspecification_mult_kappa
      ),
    true_growth = true_growth,
    dgp_params = list(
      n_units = n_units,
      years = years,
      sex_levels = sex_levels,
      age_levels = age_levels,
      heterogeneous_slopes = heterogeneous_slopes,
      regional_random_slopes = regional_random_slopes,
      regional_slope_sd = slope_sd,
      pandemic_shock = pandemic_shock,
      noise_scale = noise_scale,
      composition_shift = composition_shift,
      composition_shift_strength = composition_shift_strength,
      misspecification_stress = misspecification_stress,
      misspecification_strength = misspecification_strength,
      unit_re_sd = unit_re_sd,
      overdispersion_patients = overdispersion_patients,
      overdispersion_provision = overdispersion_provision,
      tweedie_power = tweedie_power,
      tweedie_phi = tweedie_phi
    )
  )
}

#' Draw Tweedie-like (compound Poisson-Gamma) values without a hard
#' dependency on the {tweedie} package: zero when n_provision == 0,
#' otherwise Gamma noise around the mean with a Tweedie-consistent
#' variance function Var = phi * mean^power.
simulate_tweedie_like <- function(mean_vec, provision_vec, power = 1.5, phi = 2) {
  out <- numeric(length(mean_vec))
  pos <- provision_vec > 0
  mu <- pmax(mean_vec[pos], 1e-6)
  var_target <- phi * mu^power
  shape <- mu^2 / var_target
  scale <- var_target / mu
  out[pos] <- rgamma(sum(pos), shape = shape, rate = 1 / scale)
  out[!pos] <- 0
  out
}

#' Exact ground-truth Shapley contribution for the simulated DGP
#'
#' The DGP is sequential in how data are generated (population -> patients
#' -> services -> value), but the validation target is the order-independent
#' allocation of the known counterfactual value function. The value function
#' is computed from the deterministic DGP means retained in the simulated
#' panel with observation-level NB/Gamma noise switched off. This preserves
#' the exact age-sex slopes, unit random intercepts, and unit random slopes
#' drawn for that replicate.
#'
#' @param sim Output of simulate_disease_panel().
#' @param baseline_year,comparison_year Years to compare.
#' @return A one-row tibble in the same layout as decompose_*_2y().
truth_from_panel <- function(sim, baseline_year, comparison_year, ...) {
  # The truth is computed from deterministic DGP means retained in sim$data,
  # not from noisy observed counts. This is robust to age-sex-specific slopes,
  # regional random intercepts, and regional random slopes because all of those
  # are already embedded in the baseline/comparison deterministic trajectories.
  comps <- c("N", "PI", "LAMBDA", "KAPPA")
  perms <- combinat_permutations(comps)

  required_cols <- c("code", "sex", "age_group", "year", "pop_mean", "pi_mean", "lambda_mean", "kappa_mean")
  missing_cols <- setdiff(required_cols, names(sim$data))
  if (length(missing_cols) > 0) {
    stop("truth_from_panel() requires deterministic DGP columns missing from sim$data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  states <- sim$data %>%
    dplyr::filter(year %in% c(baseline_year, comparison_year)) %>%
    dplyr::select(code, sex, age_group, year, pop_mean, pi_mean, lambda_mean, kappa_mean)

  if (!all(c(baseline_year, comparison_year) %in% states$year)) {
    stop("baseline_year and comparison_year must both be present in the simulated panel.", call. = FALSE)
  }

  cell_keys <- c("code", "sex", "age_group")
  base <- states %>%
    dplyr::filter(year == baseline_year) %>%
    dplyr::select(dplyr::all_of(cell_keys),
                  N_base = pop_mean, PI_base = pi_mean,
                  LAMBDA_base = lambda_mean, KAPPA_base = kappa_mean)
  comp <- states %>%
    dplyr::filter(year == comparison_year) %>%
    dplyr::select(dplyr::all_of(cell_keys),
                  N_comp = pop_mean, PI_comp = pi_mean,
                  LAMBDA_comp = lambda_mean, KAPPA_comp = kappa_mean)

  state_tbl <- base %>% dplyr::inner_join(comp, by = cell_keys)
  if (nrow(state_tbl) == 0) stop("No matched cells found for truth_from_panel().", call. = FALSE)

  expected_total <- function(active) {
    use <- function(nm) nm %in% active
    N      <- if (use("N"))      state_tbl$N_comp      else state_tbl$N_base
    PI     <- if (use("PI"))     state_tbl$PI_comp     else state_tbl$PI_base
    LAMBDA <- if (use("LAMBDA")) state_tbl$LAMBDA_comp else state_tbl$LAMBDA_base
    KAPPA  <- if (use("KAPPA"))  state_tbl$KAPPA_comp  else state_tbl$KAPPA_base
    sum(N * PI * LAMBDA * KAPPA, na.rm = TRUE)
  }

  contrib_mat <- matrix(0, nrow = length(perms), ncol = 4, dimnames = list(NULL, comps))
  for (i in seq_along(perms)) {
    perm <- perms[[i]]
    active <- character(0)
    v_prev <- expected_total(active)
    for (component_name in perm) {
      active <- c(active, component_name)
      v_curr <- expected_total(active)
      contrib_mat[i, component_name] <- v_curr - v_prev
      v_prev <- v_curr
    }
  }

  sv <- colMeans(contrib_mat)
  V0 <- expected_total(character(0))
  V4 <- expected_total(comps)

  tibble::tibble(
    baseline_year = baseline_year, comparison_year = comparison_year,
    baseline_total = V0, comparison_total = V4, total_change = V4 - V0,
    N_effect = unname(sv["N"]), PI_effect = unname(sv["PI"]),
    LAMBDA_effect = unname(sv["LAMBDA"]), KAPPA_effect = unname(sv["KAPPA"]),
    residual = (V4 - V0) - sum(sv)
  )
}

#' Run the full simulation study: repeat (simulate -> fit -> decompose
#' with all 4 methods -> compare to truth) `reps` times, and summarise
#' bias and variance of the recovered component effects.
#'
#' @param model_cfg The `model:` config block (used to fit each replicate).
#' @param true_growth DGP parameters (see simulate_disease_panel()).
#' @param baseline_year,comparison_year Years to compare within each
#'   simulated panel's year range.
#' @param years Year range to simulate.
#' @param reps Number of independent simulated panels to draw.
#' @param seed Base RNG seed (replicate i uses seed + i for reproducibility).
#' @param methods Methods to evaluate, subset of
#'   c("naive","ofat","sequential","shapley").
#' @return list(raw = <all rep x method x component rows, with `truth`
#'              and `error = estimate - truth` attached>,
#'              summary = <bias, RMSE, SD of error per method x component>)
run_simulation_study <- function(model_cfg,
                                  true_growth = list(N = 0.01, PI = 0.03, LAMBDA = -0.01, KAPPA = 0.02),
                                  baseline_year = 2018, comparison_year = 2023,
                                  years = 2018:2023,
                                  n_units = 30,
                                  sex_levels = c("F", "M"),
                                  age_levels = c("0-17", "18-44", "45-64", "65+"),
                                  heterogeneous_slopes = TRUE,
                                  regional_random_slopes = TRUE,
                                  regional_slope_sd = list(N = 0.003, PI = 0.008, LAMBDA = 0.008, KAPPA = 0.006),
                                  pandemic_shock = FALSE,
                                  noise_scale = 1,
                                  composition_shift = FALSE,
                                  composition_shift_strength = 0.25,
                                  misspecification_stress = FALSE,
                                  misspecification_strength = 0.25,
                                  reps = 200, seed = 1,
                                  methods = c("naive", "ofat", "sequential", "shapley"),
                                  stop_on_fit_failure = TRUE,
                                  stop_on_method_failure = TRUE) {

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
      composition_shift = composition_shift,
      composition_shift_strength = composition_shift_strength,
      misspecification_stress = misspecification_stress,
      misspecification_strength = misspecification_strength,
      seed = seed + r
    )

    truth <- truth_from_panel(sim, baseline_year, comparison_year)

    fit_error <- NULL
    fit_r <- tryCatch(
      fit_model_chain_simulation_safe(sim$data, model_cfg),
      error = function(e) {
        fit_error <<- conditionMessage(e)
        NULL
      }
    )
    df_pair <- sim$data %>% filter(year %in% c(baseline_year, comparison_year))

    rows <- list()
    if ("naive" %in% methods) {
      rows$naive <- decompose_naive_2y(df_pair, NULL, baseline_year, comparison_year) %>% mutate(method = "naive")
    }

    model_based_methods <- intersect(methods, c("ofat", "sequential", "shapley"))
    if (length(model_based_methods) > 0) {
      if (is.null(fit_r)) {
        msg <- paste0(
          "Simulation replicate ", r, ": fit_model_chain_simulation_safe() failed, so model-based ",
          "methods (", paste(model_based_methods, collapse = ", "), ") cannot be computed. ",
          "Original error: ", fit_error
        )
        if (isTRUE(stop_on_fit_failure)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
      } else {
        run_method <- function(method_name, expr) {
          method_error <- NULL
          ans <- tryCatch(expr, error = function(e) {
            method_error <<- conditionMessage(e)
            NULL
          })
          if (is.null(ans)) {
            msg <- paste0(
              "Simulation replicate ", r, ": decomposition method '", method_name,
              "' failed. Original error: ", method_error
            )
            if (isTRUE(stop_on_method_failure)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
          }
          ans
        }

        if ("ofat" %in% methods) {
          rows$ofat <- run_method(
            "ofat",
            decompose_ofat_2y(df_pair, fit_r, baseline_year, comparison_year) %>% mutate(method = "ofat")
          )
        }
        if ("sequential" %in% methods) {
          rows$sequential <- run_method(
            "sequential",
            decompose_sequential_2y(df_pair, fit_r, baseline_year, comparison_year) %>% mutate(method = "sequential")
          )
        }
        if ("shapley" %in% methods) {
          rows$shapley <- run_method(
            "shapley",
            decompose_shapley_2y(df_pair, fit_r, baseline_year, comparison_year) %>% mutate(method = "shapley")
          )
        }
      }
    }

    rep_tbl <- bind_rows(rows)
    if (nrow(rep_tbl) > 0) {
      rep_tbl$rep_id      <- r
      rep_tbl$truth_N      <- truth$N_effect
      rep_tbl$truth_PI     <- truth$PI_effect
      rep_tbl$truth_LAMBDA <- truth$LAMBDA_effect
      rep_tbl$truth_KAPPA  <- truth$KAPPA_effect
    }

    out[[r]] <- rep_tbl
    setTxtProgressBar(pb, r)
  }
  close(pb)

  raw <- bind_rows(out) %>%
    tidyr::pivot_longer(
      cols = c(N_effect, PI_effect, LAMBDA_effect, KAPPA_effect),
      names_to = "component", values_to = "estimate"
    ) %>%
    mutate(
      component = sub("_effect$", "", component),
      truth = case_when(
        component == "N"      ~ truth_N,
        component == "PI"     ~ truth_PI,
        component == "LAMBDA" ~ truth_LAMBDA,
        component == "KAPPA"  ~ truth_KAPPA
      ),
      error = estimate - truth
    ) %>%
    dplyr::select(rep_id, method, component, estimate, truth, error)

  expected_methods <- intersect(methods, c("naive", "ofat", "sequential", "shapley"))
  observed_methods <- sort(unique(raw$method))
  missing_methods <- setdiff(expected_methods, observed_methods)
  if (length(missing_methods) > 0) {
    stop(
      "Simulation output is missing method(s): ", paste(missing_methods, collapse = ", "),
      ". This usually means the fitted model or a decomposition step failed in every replicate. ",
      "The recovery plot would be misleading, so the pipeline has been stopped.",
      call. = FALSE
    )
  }

  summary_tbl <- raw %>%
    group_by(method, component) %>%
    summarise(
      truth     = mean(truth, na.rm = TRUE),
      bias      = mean(error, na.rm = TRUE),
      rmse      = sqrt(mean(error^2, na.rm = TRUE)),
      sd_error  = sd(error, na.rm = TRUE),
      n_reps_ok = sum(!is.na(error)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      bias_pct_of_truth = dplyr::if_else(abs(truth) > 0, bias / abs(truth), NA_real_),
      rmse_pct_of_truth = dplyr::if_else(abs(truth) > 0, rmse / abs(truth), NA_real_)
    ) %>%
    dplyr::group_by(component) %>%
    dplyr::mutate(rmse_rank_within_component = dplyr::min_rank(rmse)) %>%
    dplyr::ungroup() %>%
    arrange(component, method)

  list(raw = raw, summary = summary_tbl)
}
