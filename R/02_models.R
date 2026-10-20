# ============================================================
# R/02_models.R
# ------------------------------------------------------------
# Fits the three offset-linked GAM/GLM models that generate
# the model-based counterfactual predictions:
#
#   patients   ~ sex * age_group + s(year, by = sex_age) + s(code, bs='re')
#                + offset(log(pop))
#   provision  ~ sex * age_group + s(year, by = sex_age) + s(code, bs='re')
#                + offset(log(n_patients))
#   value      ~ sex * age_group + s(year, by = sex_age) + s(code, bs='re')
#                + offset(log(n_provision))
#
# patients and provision use a negative binomial family (overdispersed
# counts). value uses, by default, a Tweedie family (continuous,
# right-skewed, with a point mass at zero) -- but this is configurable
# via config$model$value_family ("tweedie" or "gamma"). Use the
# diagnostics in R/13_distribution_diagnostics.R / scripts/run_08_*.R
# to check empirically which is better supported by your data: if
# value is confirmed strictly positive whenever n_provision > 0 (no
# genuine zero mass), Gamma fits at least as well with one fewer
# parameter, since Tweedie's zero-mass flexibility is then unused.
#
# This module operates on a single disease group's worth of data
# (already filtered by filter_and_aggregate()), spanning multiple years.
# ============================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
})

#' Build the right-hand side formula shared by all three stage models
#'
#' @param response Name of the response column.
#' @param offset_var Name of the column entering as log(offset).
#' @param use_spatial_smooth Whether to add a 2D spatial smooth s(x, y).
#' @param k_year Basis dimension for the year spline.
#' @param eps Floor applied inside log(pmax(offset_var, eps)) to avoid log(0).
build_formula <- function(response, offset_var, use_spatial_smooth = FALSE,
                           k_year = 5, eps = 1) {
  rhs <- c(
    "sex * age_group",
    paste0("s(year_num, by = sex_age, k = ", k_year, ")"),
    "s(code, bs = 're')"
  )
  if (use_spatial_smooth) rhs <- c(rhs, "s(x, y, bs = 'tp')")
  rhs <- c(rhs, paste0("offset(log(pmax(", offset_var, ", ", eps, ")))"))

  as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
}

#' Fit the three-stage counterfactual-generating model chain for one disease
#'
#' @param df Aggregated cells for a single disease, multiple years
#'   (output of filter_and_aggregate()).
#' @param model_cfg The `model` block of config.R.
#' @return A list with the three fitted mgcv::bam objects, deviance
#'   explained for each, and bookkeeping (year range, spline k used).
resolve_value_family_name <- function(model_cfg) {
  override <- NULL
  if (exists(".NFZ_VALUE_FAMILY_OVERRIDE", envir = .GlobalEnv, inherits = FALSE)) {
    override <- get(".NFZ_VALUE_FAMILY_OVERRIDE", envir = .GlobalEnv)
  }
  fam <- override %||% (model_cfg$value_family %||% "gamma")
  fam <- tolower(as.character(fam)[1])
  if (!fam %in% c("tweedie", "gamma")) {
    stop("value_family must be 'tweedie' or 'gamma', got: ", fam, call. = FALSE)
  }
  fam
}

make_value_family <- function(value_family_name) {
  switch(
    value_family_name,
    "tweedie" = tw(link = "log"),
    "gamma"   = Gamma(link = "log"),
    stop("Unknown value family: ", value_family_name, call. = FALSE)
  )
}

prepare_value_fit_data <- function(df_fit, value_family_name, stage_label = "value model") {
  if (identical(value_family_name, "gamma")) {
    keep <- is.finite(df_fit$value) & is.finite(df_fit$n_provision) & df_fit$n_provision > 0 & df_fit$value > 0
    n_excluded <- sum(!keep, na.rm = TRUE)
    n_zero_pos_prov <- sum(df_fit$n_provision > 0 & df_fit$value == 0, na.rm = TRUE)
    n_nonpos_value <- sum(df_fit$value <= 0, na.rm = TRUE)
    if (n_excluded > 0) {
      message(
        "Gamma ", stage_label, ": excluding ", n_excluded,
        " row(s) with non-positive value or non-positive/invalid service count ",
        "from the valuation-stage fit. Zero value with positive services: ", n_zero_pos_prov, "."
      )
    }
    d <- df_fit[keep, , drop = FALSE]
    if (nrow(d) < 10) {
      stop("Too few positive value/service observations for Gamma ", stage_label, ": ", nrow(d), call. = FALSE)
    }
    attr(d, "value_fit_n_excluded") <- n_excluded
    attr(d, "value_zero_positive_provision") <- n_zero_pos_prov
    attr(d, "value_nonpositive_value") <- n_nonpos_value
    return(d)
  }

  attr(df_fit, "value_fit_n_excluded") <- 0L
  attr(df_fit, "value_zero_positive_provision") <- sum(df_fit$n_provision > 0 & df_fit$value == 0, na.rm = TRUE)
  attr(df_fit, "value_nonpositive_value") <- sum(df_fit$value <= 0, na.rm = TRUE)
  df_fit
}

#' Fit the three-stage counterfactual-generating model chain for one disease
#'
#' @param df Aggregated cells for a single disease, multiple years
#'   (output of filter_and_aggregate()).
#' @param model_cfg The `model` block of config.R.
#' @return A list with the three fitted mgcv::bam objects, deviance
#'   explained for each, and bookkeeping (year range, spline k used).
fit_model_chain <- function(df, model_cfg) {
  yrs <- sort(unique(df$year))
  if (length(yrs) < 2) {
    stop("Need at least 2 distinct years to fit the model chain; found ", length(yrs), ".")
  }

  k_year <- min(model_cfg$k_year, length(yrs))
  if (k_year < 3) k_year <- length(yrs)

  df_fit <- df %>%
    mutate(
      year_num  = as.numeric(year),
      sex_age   = droplevels(interaction(sex, age_group, drop = TRUE)),
      code      = droplevels(code),
      sex       = droplevels(sex),
      age_group = droplevels(age_group)
    )

  f_pat  <- build_formula("n_patients",  "pop",         model_cfg$use_spatial_smooth, k_year, model_cfg$eps_offset)
  f_prov <- build_formula("n_provision", "n_patients",  model_cfg$use_spatial_smooth, k_year, model_cfg$eps_offset)
  f_val  <- build_formula("value",       "n_provision", model_cfg$use_spatial_smooth, k_year, model_cfg$eps_offset)

  ctrl <- gam.control(maxit = model_cfg$maxit_gam, trace = isTRUE(model_cfg$trace_gam))

  # Diagnostic update: for the NFZ neoplasm panel, zero settlement values
  # conditional on positive service counts were extremely rare (6/42,556;
  # 0.014%), and the Gamma log-link offset model had much better residual
  # behavior than the previous Tweedie specification. Gamma is therefore
  # supported as the main valuation-stage family, with Tweedie retained as
  # a sensitivity option. For Gamma we estimate the valuation stage only on
  # strictly positive value and service-count cells, as required by the family.
  value_family_name <- resolve_value_family_name(model_cfg)
  value_family <- make_value_family(value_family_name)
  df_val_fit <- prepare_value_fit_data(df_fit, value_family_name, stage_label = "main value model")

  m_pat <- bam(f_pat,  data = df_fit,     family = nb(),         discrete = TRUE, method = "fREML", control = ctrl)
  m_prov <- bam(f_prov, data = df_fit,    family = nb(),         discrete = TRUE, method = "fREML", control = ctrl)
  m_val <- bam(f_val,  data = df_val_fit, family = value_family, discrete = TRUE, method = "fREML", control = ctrl)

  list(
    m_pat = m_pat,
    m_prov = m_prov,
    m_val = m_val,
    dev_pat  = summary(m_pat)$dev.expl,
    dev_prov = summary(m_prov)$dev.expl,
    dev_val  = summary(m_val)$dev.expl,
    year_levels = yrs,
    k_year = k_year,
    value_family = value_family_name,
    value_fit_n_obs = nrow(df_val_fit),
    value_fit_n_excluded = attr(df_val_fit, "value_fit_n_excluded") %||% 0L,
    value_zero_positive_provision = attr(df_val_fit, "value_zero_positive_provision") %||% NA_integer_,
    value_nonpositive_value = attr(df_val_fit, "value_nonpositive_value") %||% NA_integer_
  )
}

#' One-row model-quality summary, useful for a paper's methods table
quality_summary <- function(fit, disease_label) {
  tibble::tibble(
    icd10    = disease_label,
    n_years  = length(fit$year_levels),
    k_year   = fit$k_year,
    value_family = fit$value_family %||% NA_character_,
    value_fit_n_obs = fit$value_fit_n_obs %||% NA_integer_,
    value_fit_n_excluded = fit$value_fit_n_excluded %||% NA_integer_,
    value_zero_positive_provision = fit$value_zero_positive_provision %||% NA_integer_,
    dev_pat  = fit$dev_pat,
    dev_prov = fit$dev_prov,
    dev_val  = fit$dev_val
  )
}

#' Build a simulation-safe model formula
#'
#' The main real-data formula uses sex-age-specific smooths of time. That is
#' appropriate for the NFZ panel, but in the simulation study some noisy
#' synthetic panels can be too sparse for that full specification. This simpler
#' formula is used ONLY as a fallback inside the simulation validation so that
#' OFAT / sequential / Shapley are evaluated for every replicate instead of
#' silently disappearing from the recovery plot.
build_formula_simulation_safe <- function(response, offset_var, k_year = 4, eps = 1) {
  rhs <- c(
    "sex * age_group",
    paste0("s(year_num, k = ", k_year, ")"),
    "s(code, bs = 're')",
    paste0("offset(log(pmax(", offset_var, ", ", eps, ")))" )
  )
  as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
}

#' Fit the model chain for simulation validation, with a robust fallback
#'
#' First tries the same model as the real-data analysis. If that fails on a
#' noisy synthetic replicate, falls back to a simpler but still three-stage
#' offset-linked GAM/GLM chain. This keeps the simulation study focused on the
#' full estimation + prediction + attribution pipeline while avoiding the old
#' failure mode where only the naive method was plotted.
fit_model_chain_simulation_safe <- function(df, model_cfg) {
  main_error <- NULL
  fit_main <- tryCatch(
    fit_model_chain(df, model_cfg),
    error = function(e) {
      main_error <<- conditionMessage(e)
      NULL
    }
  )
  if (!is.null(fit_main)) {
    fit_main$simulation_fit_spec <- "main"
    return(fit_main)
  }

  yrs <- sort(unique(df$year))
  if (length(yrs) < 2) {
    stop("Need at least 2 distinct years for simulation fallback; found ", length(yrs), ".")
  }

  # Safer basis for small synthetic panels. Keep k below the number of years.
  k_year <- min(4, length(yrs) - 1)
  if (k_year < 3) k_year <- min(3, length(yrs))

  df_fit <- df %>%
    mutate(
      year_num  = as.numeric(year),
      sex_age   = droplevels(interaction(sex, age_group, drop = TRUE)),
      code      = droplevels(code),
      sex       = droplevels(sex),
      age_group = droplevels(age_group)
    )

  f_pat  <- build_formula_simulation_safe("n_patients",  "pop",         k_year, model_cfg$eps_offset)
  f_prov <- build_formula_simulation_safe("n_provision", "n_patients",  k_year, model_cfg$eps_offset)
  f_val  <- build_formula_simulation_safe("value",       "n_provision", k_year, model_cfg$eps_offset)

  ctrl <- gam.control(maxit = model_cfg$maxit_gam %||% 200, trace = FALSE)

  value_family_name <- resolve_value_family_name(model_cfg)
  value_family <- make_value_family(value_family_name)
  df_val_fit <- prepare_value_fit_data(df_fit, value_family_name, stage_label = "simulation fallback value model")

  fallback_error <- NULL
  fit_fb <- tryCatch({
    m_pat <- bam(f_pat,  data = df_fit,     family = nb(), discrete = TRUE, method = "fREML", control = ctrl)
    m_prov <- bam(f_prov, data = df_fit,    family = nb(), discrete = TRUE, method = "fREML", control = ctrl)
    m_val <- bam(f_val,  data = df_val_fit, family = value_family, discrete = TRUE, method = "fREML", control = ctrl)

    list(
      m_pat = m_pat,
      m_prov = m_prov,
      m_val = m_val,
      dev_pat  = summary(m_pat)$dev.expl,
      dev_prov = summary(m_prov)$dev.expl,
      dev_val  = summary(m_val)$dev.expl,
      year_levels = yrs,
      k_year = k_year,
      value_family = value_family_name,
      value_fit_n_obs = nrow(df_val_fit),
      value_fit_n_excluded = attr(df_val_fit, "value_fit_n_excluded") %||% 0L,
      value_zero_positive_provision = attr(df_val_fit, "value_zero_positive_provision") %||% NA_integer_,
      value_nonpositive_value = attr(df_val_fit, "value_nonpositive_value") %||% NA_integer_,
      simulation_fit_spec = "fallback_simple_time_smooth",
      simulation_main_error = main_error
    )
  }, error = function(e) {
    fallback_error <<- conditionMessage(e)
    NULL
  })

  if (is.null(fit_fb)) {
    stop(
      "Both main and simulation-safe model fits failed. Main error: ", main_error,
      "; fallback error: ", fallback_error,
      call. = FALSE
    )
  }

  fit_fb
}
