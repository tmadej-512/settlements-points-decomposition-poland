# ============================================================
# scripts/run_14_family_selection_with_moran_by_distribution.R
# ------------------------------------------------------------
# Cel:
#   Pełna diagnostyka rodzin rozkładu + Moran's I dla KAŻDEGO
#   kandydata rodziny.
#
#   Stages:
#     patients  : Poisson, quasi-Poisson, Negative Binomial, Tweedie
#     provision : Poisson, quasi-Poisson, Negative Binomial, Tweedie
#     value     : Gamma, Tweedie, Inverse Gaussian
#
#   Moran's I:
#     - residuale Pearsona z każdego fitted candidate model
#     - agregacja residuali do administrative unit `code`
#     - k-nearest-neighbour weights, domyślnie k = 8
#     - Monte Carlo p-value, domyślnie 999 permutacji
#
# Wejście:
#   1) config.R, jeśli dostępny i zawiera:
#        load_config(), read_input_data(), normalize_schema()
#   albo
#   2) output/rds/aggregated_data.rds
#
# Wyjście:
#   output/tables/R14_family_selection_moran/
#     01_family_diagnostics_all_candidates.csv
#     02_moran_by_stage_family.csv
#     03_family_selection_recommendations_with_moran.csv
#     04_final_selected_model_quality_with_moran.csv
#     05_text_for_methods_limitations.csv
#
#   output/figures/R14_family_selection_moran/
#     moran_by_family.png
#     aic_by_family.png
#     pearson_dispersion_by_family.png
# ============================================================

options(scipen = 999)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(mgcv)
  library(purrr)
  library(tibble)
})

if (!requireNamespace("spdep", quietly = TRUE)) {
  stop(
    "Brakuje pakietu 'spdep'. Zainstaluj go przed uruchomieniem:\n",
    "install.packages('spdep')",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# 1. Ustawienia
# ------------------------------------------------------------

config_path <- "config.R"
fallback_rds <- "output/rds/aggregated_data.rds"

output_tables <- "output/tables/R14_family_selection_moran"
output_figures <- "output/figures/R14_family_selection_moran"
output_rds <- "output/rds/R14_family_selection_moran"

dir.create(output_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(output_rds, recursive = TRUE, showWarnings = FALSE)

# Główna choroba z artykułu.
# Ustaw NULL, jeśli wejściowy RDS jest już przefiltrowany.
disease_filter <- "Neoplasms"
icd_col_override <- NULL

# GAM.
k_year_default <- 5
include_code_re <- TRUE
use_discrete_bam <- TRUE

# Count stages: Tweedie jako sensitivity. Domyślnie nie wybiera Tweedie
# jako finalnej rodziny dla zmiennych licznikowych.
allow_tweedie_for_count_final <- FALSE

# Moran's I.
run_moran <- TRUE
moran_k <- 8
moran_nsim <- 999
moran_residual_type <- "pearson"

# Minimalne warunki.
min_rows_stage <- 50
min_positive_outcome <- 20

# ------------------------------------------------------------
# 2. Helpery ogólne
# ------------------------------------------------------------

safe_sum <- function(x) sum(x, na.rm = TRUE)

safe_div <- function(num, den) {
  num <- as.numeric(num)
  den <- as.numeric(den)
  n <- max(length(num), length(den))
  num <- rep(num, length.out = n)
  den <- rep(den, length.out = n)
  out <- rep(NA_real_, n)
  ok <- !is.na(num) & !is.na(den) & den != 0
  out[ok] <- num[ok] / den[ok]
  out
}

safe_pct <- function(num, den) {
  num <- as.numeric(num)
  den <- as.numeric(den)
  n <- max(length(num), length(den))
  num <- rep(num, length.out = n)
  den <- rep(den, length.out = n)
  out <- rep(NA_real_, n)
  ok <- !is.na(num) & !is.na(den) & den != 0
  out[ok] <- 100 * num[ok] / abs(den[ok])
  out
}

safe_aic <- function(model) {
  out <- tryCatch(
    AIC(model),
    error = function(e) NA_real_,
    warning = function(w) NA_real_
  )
  as.numeric(out)
}

safe_dev_expl <- function(model) {
  out <- tryCatch(summary(model)$dev.expl, error = function(e) NA_real_)
  as.numeric(out)
}

safe_pearson_dispersion <- function(model) {
  out <- tryCatch({
    rp <- residuals(model, type = "pearson")
    rdf <- stats::df.residual(model)
    if (is.na(rdf) || rdf <= 0) return(NA_real_)
    sum(rp^2, na.rm = TRUE) / rdf
  }, error = function(e) NA_real_)
  as.numeric(out)
}

safe_theta <- function(model) {
  out <- tryCatch({
    if (!is.null(model$family$getTheta)) model$family$getTheta(TRUE) else NA_real_
  }, error = function(e) NA_real_)
  as.numeric(out)
}

safe_rmse_response <- function(model, data, outcome) {
  out <- tryCatch({
    pred <- as.numeric(predict(model, newdata = data, type = "response"))
    obs <- data[[outcome]]
    sqrt(mean((obs - pred)^2, na.rm = TRUE))
  }, error = function(e) NA_real_)
  as.numeric(out)
}

safe_mae_response <- function(model, data, outcome) {
  out <- tryCatch({
    pred <- as.numeric(predict(model, newdata = data, type = "response"))
    obs <- data[[outcome]]
    mean(abs(obs - pred), na.rm = TRUE)
  }, error = function(e) NA_real_)
  as.numeric(out)
}

stop_if_missing_cols <- function(df, cols) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Brakuje wymaganych kolumn: ",
      paste(missing_cols, collapse = ", "),
      "\nDostępne kolumny:\n",
      paste(names(df), collapse = ", "),
      call. = FALSE
    )
  }
}

detect_icd_col <- function(df) {
  if (!is.null(icd_col_override)) {
    if (!icd_col_override %in% names(df)) {
      stop("icd_col_override nie istnieje w danych: ", icd_col_override, call. = FALSE)
    }
    return(icd_col_override)
  }

  candidates <- c("ICD", "icd", "icd10", "disease", "disease_group", "stan_chorobowy", "choroba")
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

detect_value_col <- function(df) {
  candidates <- c("value", "sum_points", "points", "liczba_punktow", "liczba_punktów")
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    stop(
      "Nie znaleziono kolumny value/sum_points. Dostępne kolumny:\n",
      paste(names(df), collapse = ", "),
      call. = FALSE
    )
  }
  hit[1]
}

# ------------------------------------------------------------
# 3. Wczytanie danych
# ------------------------------------------------------------

read_main_data <- function() {
  if (file.exists(config_path)) {
    source(config_path)
    has_config_funs <- all(c("load_config", "read_input_data", "normalize_schema") %in% ls(envir = .GlobalEnv))

    if (has_config_funs) {
      cfg <- load_config(config_path)
      message("Wczytuję dane z config.R: ", cfg$data$path)
      df_raw <- read_input_data(cfg$data$path)

      if (!is.null(cfg$data$value_var)) {
        return(normalize_schema(df_raw, value_var = cfg$data$value_var))
      }
      return(normalize_schema(df_raw))
    }
  }

  if (file.exists(fallback_rds)) {
    message("Wczytuję fallback RDS: ", fallback_rds)
    return(readRDS(fallback_rds))
  }

  stop("Nie znaleziono config.R ani ", fallback_rds, call. = FALSE)
}

dat_raw <- read_main_data()

message("Dane wejściowe: ", nrow(dat_raw), " wierszy")
print(names(dat_raw))

value_col <- detect_value_col(dat_raw)
icd_col <- detect_icd_col(dat_raw)

required_cols <- c("year", "code", "sex", "age_group", "pop", "n_patients", "n_provision", value_col)
stop_if_missing_cols(dat_raw, required_cols)

# x/y są wymagane tylko do Morana.
if (run_moran && !all(c("x", "y") %in% names(dat_raw))) {
  warning("Brakuje kolumn x/y. Moran's I zostanie pominięty.")
  run_moran <- FALSE
}

dat <- dat_raw %>%
  dplyr::mutate(
    ICD = if (!is.na(icd_col)) as.character(.data[[icd_col]]) else "ALL",
    value = as.numeric(.data[[value_col]]),
    year = as.integer(.data$year),
    code = as.character(.data$code),
    sex = as.character(.data$sex),
    age_group = as.character(.data$age_group),
    pop = as.numeric(.data$pop),
    n_patients = as.numeric(.data$n_patients),
    n_provision = as.numeric(.data$n_provision),
    x = if ("x" %in% names(dat_raw)) as.numeric(.data$x) else NA_real_,
    y = if ("y" %in% names(dat_raw)) as.numeric(.data$y) else NA_real_
  ) %>%
  dplyr::filter(
    !is.na(.data$year),
    !is.na(.data$code),
    !is.na(.data$sex),
    !is.na(.data$age_group),
    !is.na(.data$pop),
    !is.na(.data$n_patients),
    !is.na(.data$n_provision),
    !is.na(.data$value)
  )

if (!is.null(disease_filter) && !is.na(icd_col)) {
  dat <- dat %>% dplyr::filter(.data$ICD == disease_filter)
  if (nrow(dat) == 0) {
    stop("Po disease_filter = '", disease_filter, "' nie zostały żadne wiersze.", call. = FALSE)
  }
}

dat_model <- dat %>%
  dplyr::group_by(.data$year, .data$code, .data$sex, .data$age_group, .data$ICD) %>%
  dplyr::summarise(
    pop = safe_sum(.data$pop),
    n_patients = safe_sum(.data$n_patients),
    n_provision = safe_sum(.data$n_provision),
    value = safe_sum(.data$value),
    x = mean(.data$x, na.rm = TRUE),
    y = mean(.data$y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    year_num = as.numeric(.data$year),
    code = factor(.data$code),
    sex = factor(.data$sex),
    age_group = factor(.data$age_group),
    sex_age = interaction(.data$sex, .data$age_group, drop = TRUE)
  )

years_available <- sort(unique(dat_model$year))
k_year <- min(k_year_default, max(3, length(years_available) - 1))

message("Dane modelowe:")
print(dat_model %>%
  dplyr::summarise(
    ICD = paste(unique(.data$ICD), collapse = ", "),
    n_rows = dplyr::n(),
    n_years = dplyr::n_distinct(.data$year),
    n_codes = dplyr::n_distinct(.data$code),
    n_sex_age = dplyr::n_distinct(.data$sex_age),
    total_pop = sum(.data$pop, na.rm = TRUE),
    total_patients = sum(.data$n_patients, na.rm = TRUE),
    total_provision = sum(.data$n_provision, na.rm = TRUE),
    total_value = sum(.data$value, na.rm = TRUE)
  ))

saveRDS(dat_model, file.path(output_rds, "dat_model.rds"))

# ------------------------------------------------------------
# 4. Formuły i family candidates
# ------------------------------------------------------------

make_rhs <- function(include_code_re = TRUE) {
  rhs <- paste0("sex * age_group + s(year_num, by = sex_age, k = ", k_year, ")")
  if (include_code_re) rhs <- paste0(rhs, " + s(code, bs = 're')")
  rhs
}

rhs_main <- make_rhs(include_code_re = include_code_re)

stage_specs <- tibble::tribble(
  ~stage,      ~outcome,       ~offset_var,     ~family,
  "patients",  "n_patients",   "pop",          "Poisson",
  "patients",  "n_patients",   "pop",          "quasi-Poisson",
  "patients",  "n_patients",   "pop",          "Negative Binomial",
  "patients",  "n_patients",   "pop",          "Tweedie",
  "provision", "n_provision",  "n_patients",   "Poisson",
  "provision", "n_provision",  "n_patients",   "quasi-Poisson",
  "provision", "n_provision",  "n_patients",   "Negative Binomial",
  "provision", "n_provision",  "n_patients",   "Tweedie",
  "value",     "value",        "n_provision",  "Gamma",
  "value",     "value",        "n_provision",  "Tweedie",
  "value",     "value",        "n_provision",  "Inverse Gaussian"
)

make_formula <- function(outcome, offset_var) {
  stats::as.formula(paste0(outcome, " ~ ", rhs_main, " + offset(log(", offset_var, "))"))
}

make_family <- function(family_name) {
  switch(
    family_name,
    "Poisson" = poisson(link = "log"),
    "quasi-Poisson" = quasipoisson(link = "log"),
    "Negative Binomial" = mgcv::nb(link = "log"),
    "Tweedie" = mgcv::tw(link = "log"),
    "Gamma" = Gamma(link = "log"),
    "Inverse Gaussian" = inverse.gaussian(link = "log"),
    stop("Nieznana rodzina: ", family_name, call. = FALSE)
  )
}

prepare_stage_data <- function(df, outcome, offset_var, family_name) {
  out <- df %>%
    dplyr::filter(
      !is.na(.data[[outcome]]),
      !is.na(.data[[offset_var]]),
      .data[[offset_var]] > 0
    )

  if (family_name %in% c("Gamma", "Inverse Gaussian")) {
    out <- out %>% dplyr::filter(.data[[outcome]] > 0)
  } else {
    out <- out %>% dplyr::filter(.data[[outcome]] >= 0)
  }

  out
}

fit_candidate <- function(stage, outcome, offset_var, family_name) {
  data_use <- prepare_stage_data(dat_model, outcome, offset_var, family_name)
  formula_use <- make_formula(outcome, offset_var)

  message("------------------------------------------------------------")
  message("Stage: ", stage, " | Family: ", family_name, " | n = ", nrow(data_use))
  message("------------------------------------------------------------")

  if (nrow(data_use) < min_rows_stage || sum(data_use[[outcome]] > 0, na.rm = TRUE) < min_positive_outcome) {
    return(list(
      fit = NULL,
      data = data_use,
      formula = formula_use,
      error = "Za mało obserwacji lub dodatnich outcome."
    ))
  }

  fam <- make_family(family_name)
  warn_txt <- character(0)

  fit <- tryCatch(
    withCallingHandlers(
      mgcv::bam(
        formula = formula_use,
        data = data_use,
        family = fam,
        method = "fREML",
        discrete = use_discrete_bam
      ),
      warning = function(w) {
        warn_txt <<- c(warn_txt, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      message("Błąd: ", conditionMessage(e))
      NULL
    }
  )

  list(
    fit = fit,
    data = data_use,
    formula = formula_use,
    warning = paste(unique(warn_txt), collapse = " | "),
    error = if (is.null(fit)) "fit_failed" else NA_character_
  )
}

# ------------------------------------------------------------
# 5. Moran's I helper
# ------------------------------------------------------------

compute_moran_for_candidate <- function(fit, data_use, stage, family_name, outcome, offset_var) {
  base_empty <- tibble::tibble(
    stage = stage,
    outcome = outcome,
    offset = paste0("log(", offset_var, ")"),
    family = family_name,
    moran_i = NA_real_,
    moran_mc_p_value = NA_real_,
    moran_expectation = NA_real_,
    moran_variance = NA_real_,
    moran_k = moran_k,
    moran_nsim = moran_nsim,
    n_units_moran = NA_integer_,
    n_cells_residual = nrow(data_use),
    moran_status = NA_character_
  )

  if (!run_moran) {
    return(base_empty %>% dplyr::mutate(moran_status = "skipped_no_xy_or_disabled"))
  }

  if (is.null(fit)) {
    return(base_empty %>% dplyr::mutate(moran_status = "skipped_fit_failed"))
  }

  if (!all(c("code", "x", "y") %in% names(data_use))) {
    return(base_empty %>% dplyr::mutate(moran_status = "skipped_missing_code_xy"))
  }

  res <- tryCatch(
    as.numeric(residuals(fit, type = moran_residual_type)),
    error = function(e) rep(NA_real_, nrow(data_use))
  )

  if (length(res) != nrow(data_use)) {
    return(base_empty %>% dplyr::mutate(moran_status = "skipped_residual_length_mismatch"))
  }

  unit_resid <- data_use %>%
    dplyr::mutate(.resid = res) %>%
    dplyr::filter(
      is.finite(.data$.resid),
      is.finite(.data$x),
      is.finite(.data$y)
    ) %>%
    dplyr::group_by(.data$code) %>%
    dplyr::summarise(
      x = mean(.data$x, na.rm = TRUE),
      y = mean(.data$y, na.rm = TRUE),
      residual_mean = mean(.data$.resid, na.rm = TRUE),
      residual_sum = sum(.data$.resid, na.rm = TRUE),
      n_cells = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      is.finite(.data$x),
      is.finite(.data$y),
      is.finite(.data$residual_mean)
    )

  n_units <- nrow(unit_resid)
  if (n_units <= moran_k + 1) {
    return(base_empty %>%
      dplyr::mutate(
        n_units_moran = n_units,
        moran_status = "skipped_too_few_units"
      ))
  }

  if (stats::sd(unit_resid$residual_mean, na.rm = TRUE) == 0) {
    return(base_empty %>%
      dplyr::mutate(
        n_units_moran = n_units,
        moran_status = "skipped_constant_residual"
      ))
  }

  coords <- as.matrix(unit_resid[, c("x", "y")])

  out <- tryCatch({
    knn <- spdep::knearneigh(coords, k = moran_k)
    nb <- spdep::knn2nb(knn)
    lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

    mc <- spdep::moran.mc(
      x = unit_resid$residual_mean,
      listw = lw,
      nsim = moran_nsim,
      zero.policy = TRUE
    )

    mt <- spdep::moran.test(
      x = unit_resid$residual_mean,
      listw = lw,
      zero.policy = TRUE
    )

    tibble::tibble(
      stage = stage,
      outcome = outcome,
      offset = paste0("log(", offset_var, ")"),
      family = family_name,
      moran_i = as.numeric(mc$statistic),
      moran_mc_p_value = as.numeric(mc$p.value),
      moran_expectation = as.numeric(mt$estimate[["Expectation"]]),
      moran_variance = as.numeric(mt$estimate[["Variance"]]),
      moran_k = moran_k,
      moran_nsim = moran_nsim,
      n_units_moran = n_units,
      n_cells_residual = nrow(data_use),
      moran_status = "ok"
    )
  }, error = function(e) {
    base_empty %>%
      dplyr::mutate(
        n_units_moran = n_units,
        moran_status = paste0("moran_failed: ", conditionMessage(e))
      )
  })

  out
}

# ------------------------------------------------------------
# 6. Fit wszystkich kandydatów + diagnostyka + Moran
# ------------------------------------------------------------

candidate_objects <- vector("list", nrow(stage_specs))
quality_rows <- vector("list", nrow(stage_specs))
moran_rows <- vector("list", nrow(stage_specs))

for (i in seq_len(nrow(stage_specs))) {
  spec <- stage_specs[i, ]

  obj <- fit_candidate(
    stage = spec$stage,
    outcome = spec$outcome,
    offset_var = spec$offset_var,
    family_name = spec$family
  )

  candidate_objects[[i]] <- obj

  fit <- obj$fit
  data_use <- obj$data

  quality_rows[[i]] <- tibble::tibble(
    stage = spec$stage,
    outcome = spec$outcome,
    offset = paste0("log(", spec$offset_var, ")"),
    family = spec$family,
    n_rows = nrow(data_use),
    n_positive_outcome = sum(data_use[[spec$outcome]] > 0, na.rm = TRUE),
    n_zero_outcome = sum(data_use[[spec$outcome]] == 0, na.rm = TRUE),
    aic = if (!is.null(fit)) safe_aic(fit) else NA_real_,
    deviance_explained = if (!is.null(fit)) safe_dev_expl(fit) else NA_real_,
    pearson_dispersion = if (!is.null(fit)) safe_pearson_dispersion(fit) else NA_real_,
    rmse_response = if (!is.null(fit)) safe_rmse_response(fit, data_use, spec$outcome) else NA_real_,
    mae_response = if (!is.null(fit)) safe_mae_response(fit, data_use, spec$outcome) else NA_real_,
    theta_nb = if (!is.null(fit) && spec$family == "Negative Binomial") safe_theta(fit) else NA_real_,
    converged = !is.null(fit),
    warning = ifelse(is.null(obj$warning) || obj$warning == "", NA_character_, obj$warning),
    error = obj$error
  )

  moran_rows[[i]] <- compute_moran_for_candidate(
    fit = fit,
    data_use = data_use,
    stage = spec$stage,
    family_name = spec$family,
    outcome = spec$outcome,
    offset_var = spec$offset_var
  )
}

names(candidate_objects) <- paste(stage_specs$stage, stage_specs$family, sep = "__")

family_quality <- dplyr::bind_rows(quality_rows)
moran_results <- dplyr::bind_rows(moran_rows)

family_quality_with_moran <- family_quality %>%
  dplyr::left_join(
    moran_results %>%
      dplyr::select(.data$stage, .data$family, .data$moran_i, .data$moran_mc_p_value, .data$moran_status, .data$n_units_moran),
    by = c("stage", "family")
  ) %>%
  dplyr::arrange(
    factor(.data$stage, levels = c("patients", "provision", "value")),
    .data$aic
  )

readr::write_csv(
  family_quality_with_moran,
  file.path(output_tables, "01_family_diagnostics_all_candidates.csv")
)

readr::write_csv(
  moran_results,
  file.path(output_tables, "02_moran_by_stage_family.csv")
)

saveRDS(candidate_objects, file.path(output_rds, "candidate_model_objects.rds"))

# ------------------------------------------------------------
# 7. Rekomendacje rodzin z uwzględnieniem Morana
# ------------------------------------------------------------

select_family_for_stage <- function(stage_name) {
  df <- family_quality_with_moran %>%
    dplyr::filter(.data$stage == stage_name, .data$converged)

  if (nrow(df) == 0) {
    return(tibble::tibble(
      stage = stage_name,
      selected_family = NA_character_,
      best_aic_family_among_allowed = NA_character_,
      best_moran_family = NA_character_,
      rationale = "No converged candidate model."
    ))
  }

  if (stage_name %in% c("patients", "provision")) {
    allowed <- if (allow_tweedie_for_count_final) {
      c("Poisson", "Negative Binomial", "Tweedie")
    } else {
      c("Poisson", "Negative Binomial")
    }

    df_allowed <- df %>% dplyr::filter(.data$family %in% allowed, !is.na(.data$aic))
    best_aic <- df_allowed %>% dplyr::slice_min(.data$aic, n = 1, with_ties = FALSE)

    poisson <- df %>% dplyr::filter(.data$family == "Poisson")
    nb <- df %>% dplyr::filter(.data$family == "Negative Binomial")
    tweedie <- df %>% dplyr::filter(.data$family == "Tweedie")

    best_moran <- df %>%
      dplyr::filter(.data$moran_status == "ok", !is.na(.data$moran_i)) %>%
      dplyr::mutate(abs_moran = abs(.data$moran_i)) %>%
      dplyr::slice_min(.data$abs_moran, n = 1, with_ties = FALSE)

    selected <- if (nrow(best_aic) > 0) best_aic$family[1] else "Negative Binomial"

    rationale <- dplyr::case_when(
      selected == "Negative Binomial" & nrow(poisson) > 0 & nrow(nb) > 0 ~
        paste0(
          "Negative Binomial selected among count families: lower AIC than Poisson and reduced overdispersion. ",
          "Tweedie is reported as sensitivity unless allow_tweedie_for_count_final = TRUE. ",
          "Moran's I is reported as residual spatial diagnostic, not as the primary family-selection criterion."
        ),
      selected == "Poisson" ~
        "Poisson selected among allowed count families by AIC; check dispersion and Moran's I before using in final manuscript.",
      selected == "Tweedie" ~
        "Tweedie selected because allow_tweedie_for_count_final = TRUE and it had the best AIC among allowed count candidates.",
      TRUE ~
        "Family selected by AIC among allowed count candidates."
    )

    return(tibble::tibble(
      stage = stage_name,
      selected_family = selected,
      best_aic_family_among_allowed = if (nrow(best_aic) > 0) best_aic$family[1] else NA_character_,
      best_moran_family = if (nrow(best_moran) > 0) best_moran$family[1] else NA_character_,
      poisson_aic = if (nrow(poisson) > 0) poisson$aic[1] else NA_real_,
      nb_aic = if (nrow(nb) > 0) nb$aic[1] else NA_real_,
      tweedie_aic = if (nrow(tweedie) > 0) tweedie$aic[1] else NA_real_,
      delta_aic_nb_minus_poisson = if (nrow(poisson) > 0 && nrow(nb) > 0) nb$aic[1] - poisson$aic[1] else NA_real_,
      poisson_pearson_dispersion = if (nrow(poisson) > 0) poisson$pearson_dispersion[1] else NA_real_,
      nb_pearson_dispersion = if (nrow(nb) > 0) nb$pearson_dispersion[1] else NA_real_,
      tweedie_pearson_dispersion = if (nrow(tweedie) > 0) tweedie$pearson_dispersion[1] else NA_real_,
      selected_moran_i = df$moran_i[df$family == selected][1],
      selected_moran_p = df$moran_mc_p_value[df$family == selected][1],
      best_moran_i = if (nrow(best_moran) > 0) best_moran$moran_i[1] else NA_real_,
      rationale = rationale
    ))
  }

  if (stage_name == "value") {
    allowed <- c("Gamma", "Tweedie", "Inverse Gaussian")
    df_allowed <- df %>% dplyr::filter(.data$family %in% allowed, !is.na(.data$aic))
    best_aic <- df_allowed %>% dplyr::slice_min(.data$aic, n = 1, with_ties = FALSE)
    best_dev <- df %>% dplyr::filter(!is.na(.data$deviance_explained)) %>% dplyr::slice_max(.data$deviance_explained, n = 1, with_ties = FALSE)
    best_moran <- df %>%
      dplyr::filter(.data$moran_status == "ok", !is.na(.data$moran_i)) %>%
      dplyr::mutate(abs_moran = abs(.data$moran_i)) %>%
      dplyr::slice_min(.data$abs_moran, n = 1, with_ties = FALSE)

    selected <- if (nrow(best_aic) > 0) best_aic$family[1] else if (nrow(best_dev) > 0) best_dev$family[1] else NA_character_

    return(tibble::tibble(
      stage = stage_name,
      selected_family = selected,
      best_aic_family_among_allowed = if (nrow(best_aic) > 0) best_aic$family[1] else NA_character_,
      best_moran_family = if (nrow(best_moran) > 0) best_moran$family[1] else NA_character_,
      gamma_aic = df$aic[df$family == "Gamma"][1],
      tweedie_aic = df$aic[df$family == "Tweedie"][1],
      inverse_gaussian_aic = df$aic[df$family == "Inverse Gaussian"][1],
      gamma_deviance_explained = df$deviance_explained[df$family == "Gamma"][1],
      tweedie_deviance_explained = df$deviance_explained[df$family == "Tweedie"][1],
      inverse_gaussian_deviance_explained = df$deviance_explained[df$family == "Inverse Gaussian"][1],
      selected_moran_i = df$moran_i[df$family == selected][1],
      selected_moran_p = df$moran_mc_p_value[df$family == selected][1],
      best_moran_i = if (nrow(best_moran) > 0) best_moran$moran_i[1] else NA_real_,
      rationale = paste0(
        selected, " selected for value by AIC among positive-outcome continuous candidates. ",
        "Moran's I is reported as residual spatial diagnostic and can be used to describe remaining spatial autocorrelation."
      )
    ))
  }

  stop("Nieznany stage: ", stage_name, call. = FALSE)
}

recommendations <- dplyr::bind_rows(
  select_family_for_stage("patients"),
  select_family_for_stage("provision"),
  select_family_for_stage("value")
)

readr::write_csv(
  recommendations,
  file.path(output_tables, "03_family_selection_recommendations_with_moran.csv")
)

selected_quality <- family_quality_with_moran %>%
  dplyr::inner_join(
    recommendations %>% dplyr::select(.data$stage, .data$selected_family),
    by = "stage"
  ) %>%
  dplyr::filter(.data$family == .data$selected_family) %>%
  dplyr::select(-.data$selected_family)

readr::write_csv(
  selected_quality,
  file.path(output_tables, "04_final_selected_model_quality_with_moran.csv")
)

# ------------------------------------------------------------
# 8. Krótki tekst do Methods / Supplement
# ------------------------------------------------------------

text_for_methods <- tibble::tibble(
  location = c("Methods", "Supplement", "Limitations/diagnostics"),
  text = c(
    "For the two count stages, candidate generalized additive models were fitted using Poisson, quasi-Poisson, Negative Binomial, and Tweedie specifications with identical covariates and offsets. For the valuation stage, Gamma, Tweedie, and inverse-Gaussian specifications were compared.",
    "Residual spatial autocorrelation was assessed for every fitted candidate family by computing Moran's I on Pearson residuals aggregated to administrative unit, using k = 8 nearest-neighbour weights and 999 Monte Carlo permutations.",
    "Family selection was based primarily on distributional appropriateness, AIC where defined, explained deviance, and residual overdispersion; Moran's I was treated as a residual spatial diagnostic rather than as the primary distribution-selection criterion."
  )
)

readr::write_csv(
  text_for_methods,
  file.path(output_tables, "05_text_for_methods_limitations.csv")
)

# ------------------------------------------------------------
# 9. Wykresy
# ------------------------------------------------------------

theme_publication <- function(base_size = 10, legend_position = "bottom") {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      legend.box = "horizontal",
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = max(base_size - 1, 7)),
      strip.text = ggplot2::element_text(size = base_size, face = "bold"),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

label_pub_number <- function(accuracy = 0.01) {
  scales::label_number(accuracy = accuracy, big.mark = " ", decimal.mark = ".", trim = TRUE)
}

plot_moran <- family_quality_with_moran %>%
  dplyr::filter(.data$moran_status == "ok", !is.na(.data$moran_i))

if (nrow(plot_moran) > 0) {
  p_moran <- ggplot2::ggplot(
    plot_moran,
    ggplot2::aes(x = .data$family, y = .data$moran_i)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey40") +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::facet_wrap(~ stage, scales = "free_x") +
    ggplot2::scale_y_continuous(labels = label_pub_number(accuracy = 0.001)) +
    ggplot2::labs(x = "Candidate family", y = "Moran's I") +
    theme_publication(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  ggplot2::ggsave(
    file.path(output_figures, "moran_by_family.png"),
    p_moran,
    width = 9,
    height = 5,
    dpi = 600
  )
}

plot_aic <- family_quality_with_moran %>% dplyr::filter(!is.na(.data$aic))
if (nrow(plot_aic) > 0) {
  p_aic <- ggplot2::ggplot(plot_aic, ggplot2::aes(x = .data$family, y = .data$aic)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::facet_wrap(~ stage, scales = "free_x") +
    ggplot2::scale_y_continuous(labels = label_pub_number(accuracy = 1)) +
    ggplot2::labs(x = "Candidate family", y = "AIC") +
    theme_publication(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  ggplot2::ggsave(
    file.path(output_figures, "aic_by_family.png"),
    p_aic,
    width = 9,
    height = 5,
    dpi = 600
  )
}

plot_disp <- family_quality_with_moran %>% dplyr::filter(!is.na(.data$pearson_dispersion))
if (nrow(plot_disp) > 0) {
  p_disp <- ggplot2::ggplot(plot_disp, ggplot2::aes(x = .data$family, y = .data$pearson_dispersion)) +
    ggplot2::geom_hline(yintercept = 1, linewidth = 0.35, colour = "grey40") +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::facet_wrap(~ stage, scales = "free_x") +
    ggplot2::scale_y_continuous(labels = label_pub_number(accuracy = 0.01)) +
    ggplot2::labs(x = "Candidate family", y = "Pearson dispersion") +
    theme_publication(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  ggplot2::ggsave(
    file.path(output_figures, "pearson_dispersion_by_family.png"),
    p_disp,
    width = 9,
    height = 5,
    dpi = 600
  )
}

# ------------------------------------------------------------
# 10. Konsola
# ------------------------------------------------------------

message("============================================================")
message("GOTOWE — family selection + Moran's I for each distribution")
message("============================================================")
message("Tabele:")
message("  ", file.path(output_tables, "01_family_diagnostics_all_candidates.csv"))
message("  ", file.path(output_tables, "02_moran_by_stage_family.csv"))
message("  ", file.path(output_tables, "03_family_selection_recommendations_with_moran.csv"))
message("  ", file.path(output_tables, "04_final_selected_model_quality_with_moran.csv"))
message("  ", file.path(output_tables, "05_text_for_methods_limitations.csv"))
message("Figury:")
message("  ", file.path(output_figures, "moran_by_family.png"))
message("  ", file.path(output_figures, "aic_by_family.png"))
message("  ", file.path(output_figures, "pearson_dispersion_by_family.png"))
message("============================================================")

message("Wybrane rodziny:")
print(recommendations)

message("Moran's I by stage and family:")
print(moran_results)
