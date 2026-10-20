# ============================================================
# scripts/run_12_refresh_gamma_main_and_count_family_diagnostics.R
# ------------------------------------------------------------
# Cel:
#   Najnowsze przeliczenie glownego modelu:
#     patients  ~ Negative Binomial + offset(log(pop))
#     provision ~ Negative Binomial + offset(log(n_patients))
#     value     ~ Gamma + offset(log(n_provision))
#   oraz formalna diagnostyka rodzin rozkladu dla etapow licznikowych:
#     patients:   Poisson vs quasi-Poisson vs Negative Binomial
#     provision:  Poisson vs quasi-Poisson vs Negative Binomial
#   plus aktualny test Gamma vs Tweedie dla value / valuation.
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

# ------------------------------------------------------------
# 1. Ustawienia
# ------------------------------------------------------------

output_tables <- "output/tables/refresh_gamma_count"
output_figures <- "output/figures/refresh_gamma_count"
output_rds <- "output/rds/refresh_gamma_count"

dir.create(output_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(output_rds, recursive = TRUE, showWarnings = FALSE)

config_path <- "config.R"
fallback_rds <- "output/rds/aggregated_data.rds"

disease_filter <- "Neoplasms"
icd_col_override <- NULL

year0_override <- NULL
year1_override <- NULL

k_year_default <- 5
include_code_re <- TRUE
use_discrete_bam <- TRUE

min_rows_stage <- 50
min_positive_outcome <- 20

component_order <- c("N", "PI", "LAMBDA", "KAPPA")

component_labels <- c(
  N = "Population size (N)",
  PI = "Treated-patient reach (π)",
  LAMBDA = "Service intensity (λ)",
  KAPPA = "Point valuation (κ)"
)

component_labels_short <- c(
  N = "N: population",
  PI = "π: patients/pop",
  LAMBDA = "λ: services/patient",
  KAPPA = "κ: value/service"
)

# ------------------------------------------------------------
# 2. Formatowanie i helpery
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
  scales::label_number(
    accuracy = accuracy,
    big.mark = " ",
    decimal.mark = ".",
    trim = TRUE
  )
}

safe_sum <- function(x) sum(x, na.rm = TRUE)

safe_div <- function(num, den) {
  dplyr::if_else(!is.na(den) & den != 0, num / den, NA_real_)
}

safe_pct <- function(num, den) {
  dplyr::if_else(!is.na(den) & den != 0, 100 * num / abs(den), NA_real_)
}

safe_aic <- function(model) {
  out <- tryCatch(AIC(model), error = function(e) NA_real_, warning = function(w) NA_real_)
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

safe_theta <- function(model) {
  out <- tryCatch({
    if (!is.null(model$family$getTheta)) model$family$getTheta(TRUE) else NA_real_
  }, error = function(e) NA_real_)
  as.numeric(out)
}

stop_if_missing_cols <- function(df, cols) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Brakuje wymaganych kolumn: ", paste(missing_cols, collapse = ", "),
      "\nDostepne kolumny:\n", paste(names(df), collapse = ", "),
      call. = FALSE
    )
  }
}

detect_icd_col <- function(df) {
  if (!is.null(icd_col_override)) {
    if (!icd_col_override %in% names(df)) stop("icd_col_override nie istnieje w danych: ", icd_col_override, call. = FALSE)
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
    stop("Nie znaleziono kolumny value/sum_points.\nDostepne kolumny:\n", paste(names(df), collapse = ", "), call. = FALSE)
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
      message("Wczytuje dane z config.R: ", cfg$data$path)
      df_raw <- read_input_data(cfg$data$path)
      if (!is.null(cfg$data$value_var)) {
        dat0 <- normalize_schema(df_raw, value_var = cfg$data$value_var)
      } else {
        dat0 <- normalize_schema(df_raw)
      }
      return(dat0)
    }
  }
  if (file.exists(fallback_rds)) {
    message("Wczytuje fallback RDS: ", fallback_rds)
    return(readRDS(fallback_rds))
  }
  stop("Nie znaleziono danych. Oczekiwano config.R albo ", fallback_rds, call. = FALSE)
}

dat_raw <- read_main_data()

message("Dane wejsciowe:")
message("  wiersze: ", nrow(dat_raw))
message("  kolumny:")
print(names(dat_raw))

# ------------------------------------------------------------
# 4. Standaryzacja danych do modelu
# ------------------------------------------------------------

value_col <- detect_value_col(dat_raw)
icd_col <- detect_icd_col(dat_raw)

required_base <- c("year", "code", "sex", "age_group", "pop", "n_patients", "n_provision", value_col)
stop_if_missing_cols(dat_raw, required_base)

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
    n_provision = as.numeric(.data$n_provision)
  ) %>%
  dplyr::filter(
    !is.na(.data$year), !is.na(.data$code), !is.na(.data$sex), !is.na(.data$age_group),
    !is.na(.data$pop), !is.na(.data$n_patients), !is.na(.data$n_provision), !is.na(.data$value)
  )

if (!is.null(disease_filter) && !is.na(icd_col)) {
  dat <- dat %>% dplyr::filter(.data$ICD == disease_filter)
  if (nrow(dat) == 0) stop("Po disease_filter = '", disease_filter, "' nie zostaly zadne wiersze.", call. = FALSE)
}

dat_model <- dat %>%
  dplyr::group_by(.data$year, .data$code, .data$sex, .data$age_group, .data$ICD) %>%
  dplyr::summarise(
    pop = safe_sum(.data$pop),
    n_patients = safe_sum(.data$n_patients),
    n_provision = safe_sum(.data$n_provision),
    value = safe_sum(.data$value),
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
year0 <- if (is.null(year0_override)) min(years_available) else year0_override
year1 <- if (is.null(year1_override)) max(years_available) else year1_override
if (!all(c(year0, year1) %in% years_available)) stop("year0/year1 nie wystepuja w danych.", call. = FALSE)
k_year <- min(k_year_default, max(3, length(years_available) - 1))

message("Dane po standaryzacji:")
print(dat_model %>% dplyr::summarise(
  ICD = paste(unique(.data$ICD), collapse = ", "),
  n_rows = dplyr::n(), n_years = dplyr::n_distinct(.data$year),
  year0 = year0, year1 = year1,
  n_codes = dplyr::n_distinct(.data$code),
  n_sex_age = dplyr::n_distinct(.data$sex_age),
  total_pop = sum(.data$pop, na.rm = TRUE),
  total_patients = sum(.data$n_patients, na.rm = TRUE),
  total_provision = sum(.data$n_provision, na.rm = TRUE),
  total_value = sum(.data$value, na.rm = TRUE)
))

saveRDS(dat_model, file.path(output_rds, "dat_model_refresh_gamma_count.rds"))

# ------------------------------------------------------------
# 5. Formuly GAM
# ------------------------------------------------------------

make_rhs <- function(include_code_re = TRUE) {
  rhs <- paste0("sex * age_group + s(year_num, by = sex_age, k = ", k_year, ")")
  if (include_code_re) rhs <- paste0(rhs, " + s(code, bs = 're')")
  rhs
}

rhs_main <- make_rhs(include_code_re = include_code_re)
formula_pat <- stats::as.formula(paste0("n_patients ~ ", rhs_main, " + offset(log(pop))"))
formula_prov <- stats::as.formula(paste0("n_provision ~ ", rhs_main, " + offset(log(n_patients))"))
formula_val <- stats::as.formula(paste0("value ~ ", rhs_main, " + offset(log(n_provision))"))

# ------------------------------------------------------------
# 6. Fit glownego modelu Gamma/NB/NB
# ------------------------------------------------------------

fit_bam_safe <- function(formula, data, family, label) {
  message("Fit: ", label)
  print(formula)
  tryCatch(
    mgcv::bam(formula = formula, data = data, family = family, method = "fREML", discrete = use_discrete_bam),
    error = function(e) {
      message("Blad bam dla ", label, ": ", conditionMessage(e))
      message("Probuje gam() bez discrete.")
      mgcv::gam(formula = formula, data = data, family = family, method = "REML")
    }
  )
}

dat_pat <- dat_model %>% dplyr::filter(.data$pop > 0, .data$n_patients >= 0)
dat_prov <- dat_model %>% dplyr::filter(.data$n_patients > 0, .data$n_provision >= 0)
dat_val_gamma <- dat_model %>% dplyr::filter(.data$n_provision > 0, .data$value > 0)
n_value_zero_excluded <- dat_model %>% dplyr::filter(.data$n_provision > 0, .data$value <= 0) %>% nrow()

m_pat_nb <- fit_bam_safe(formula_pat, dat_pat, mgcv::nb(link = "log"), "patients Negative Binomial")
m_prov_nb <- fit_bam_safe(formula_prov, dat_prov, mgcv::nb(link = "log"), "provision Negative Binomial")
m_val_gamma <- fit_bam_safe(formula_val, dat_val_gamma, Gamma(link = "log"), "value Gamma")

main_model_quality <- tibble::tibble(
  stage = c("patients", "provision", "value"),
  outcome = c("n_patients", "n_provision", "value"),
  family = c("Negative Binomial", "Negative Binomial", "Gamma"),
  offset = c("log(pop)", "log(n_patients)", "log(n_provision)"),
  n_rows_fit = c(nrow(dat_pat), nrow(dat_prov), nrow(dat_val_gamma)),
  n_excluded_zero_value_positive_services = c(NA_integer_, NA_integer_, n_value_zero_excluded),
  aic = c(safe_aic(m_pat_nb), safe_aic(m_prov_nb), safe_aic(m_val_gamma)),
  deviance_explained = c(safe_dev_expl(m_pat_nb), safe_dev_expl(m_prov_nb), safe_dev_expl(m_val_gamma)),
  pearson_dispersion = c(safe_pearson_dispersion(m_pat_nb), safe_pearson_dispersion(m_prov_nb), safe_pearson_dispersion(m_val_gamma)),
  theta_nb = c(safe_theta(m_pat_nb), safe_theta(m_prov_nb), NA_real_)
)

readr::write_csv(main_model_quality, file.path(output_tables, "01_model_quality_updated.csv"))
saveRDS(list(m_pat = m_pat_nb, m_prov = m_prov_nb, m_val = m_val_gamma, dat_model = dat_model, year0 = year0, year1 = year1, k_year = k_year), file.path(output_rds, "main_gamma_model_chain.rds"))

# ------------------------------------------------------------
# 7. Diagnostyka rodzin dla patients/provision
# ------------------------------------------------------------

fit_one_count_family <- function(data, formula, family_name) {
  fam <- switch(family_name,
    "Poisson" = poisson(link = "log"),
    "quasi-Poisson" = quasipoisson(link = "log"),
    "Negative Binomial" = mgcv::nb(link = "log"),
    stop("Nieznana family_name: ", family_name, call. = FALSE)
  )
  tryCatch(
    mgcv::bam(formula = formula, data = data, family = fam, method = "fREML", discrete = use_discrete_bam),
    error = function(e) { message("Blad modelu ", family_name, ": ", conditionMessage(e)); NULL }
  )
}

diagnose_count_stage <- function(data, stage_name, outcome, offset_var, formula) {
  message("Diagnostyka rodzin count: ", stage_name)
  dat_stage <- data %>% dplyr::filter(!is.na(.data[[outcome]]), !is.na(.data[[offset_var]]), .data[[offset_var]] > 0, .data[[outcome]] >= 0)
  if (nrow(dat_stage) < min_rows_stage) { warning("Za malo obserwacji dla etapu: ", stage_name); return(tibble::tibble()) }
  if (sum(dat_stage[[outcome]] > 0, na.rm = TRUE) < min_positive_outcome) { warning("Za malo dodatnich obserwacji outcome dla etapu: ", stage_name); return(tibble::tibble()) }
  families <- c("Poisson", "quasi-Poisson", "Negative Binomial")
  fits <- lapply(families, function(fam) fit_one_count_family(dat_stage, formula, fam))
  names(fits) <- families
  purrr::imap_dfr(fits, function(fit, fam) {
    if (is.null(fit)) {
      return(tibble::tibble(stage = stage_name, outcome = outcome, offset = paste0("log(", offset_var, ")"), family = fam, n_rows = nrow(dat_stage), n_positive_outcome = sum(dat_stage[[outcome]] > 0, na.rm = TRUE), aic = NA_real_, deviance_explained = NA_real_, pearson_dispersion = NA_real_, rmse_response = NA_real_, mae_response = NA_real_, theta_nb = NA_real_, converged = FALSE))
    }
    tibble::tibble(stage = stage_name, outcome = outcome, offset = paste0("log(", offset_var, ")"), family = fam, n_rows = nrow(dat_stage), n_positive_outcome = sum(dat_stage[[outcome]] > 0, na.rm = TRUE), aic = safe_aic(fit), deviance_explained = safe_dev_expl(fit), pearson_dispersion = safe_pearson_dispersion(fit), rmse_response = safe_rmse_response(fit, dat_stage, outcome), mae_response = safe_mae_response(fit, dat_stage, outcome), theta_nb = if (fam == "Negative Binomial") safe_theta(fit) else NA_real_, converged = TRUE)
  })
}

recommend_count_family <- function(stage_df) {
  if (nrow(stage_df) == 0) return(tibble::tibble())
  stage_name <- unique(stage_df$stage); outcome <- unique(stage_df$outcome)
  poisson <- stage_df %>% dplyr::filter(.data$family == "Poisson")
  quasi <- stage_df %>% dplyr::filter(.data$family == "quasi-Poisson")
  nb <- stage_df %>% dplyr::filter(.data$family == "Negative Binomial")
  poisson_aic <- poisson$aic[1]; nb_aic <- nb$aic[1]
  delta_aic_nb_minus_poisson <- nb_aic - poisson_aic
  poisson_disp <- poisson$pearson_dispersion[1]; nb_disp <- nb$pearson_dispersion[1]
  recommendation <- dplyr::case_when(
    is.na(poisson_aic) | is.na(nb_aic) ~ "AIC comparison incomplete; use dispersion and deviance explained.",
    delta_aic_nb_minus_poisson < -10 & !is.na(poisson_disp) & poisson_disp > 1.5 ~ "Negative Binomial strongly preferred: lower AIC and clear Poisson overdispersion.",
    delta_aic_nb_minus_poisson < -2 ~ "Negative Binomial preferred: lower AIC than Poisson.",
    !is.na(poisson_disp) & poisson_disp > 2 & !is.na(nb_disp) & nb_disp < poisson_disp ~ "Negative Binomial preferred mainly due to reduced overdispersion.",
    abs(delta_aic_nb_minus_poisson) <= 2 & !is.na(poisson_disp) & poisson_disp <= 1.5 ~ "Poisson and Negative Binomial similar; Poisson may be adequate.",
    TRUE ~ "No strong family preference; retain Negative Binomial if consistency and overdispersion robustness are priorities."
  )
  tibble::tibble(stage = stage_name, outcome = outcome, poisson_aic = poisson_aic, nb_aic = nb_aic, delta_aic_nb_minus_poisson = delta_aic_nb_minus_poisson, poisson_pearson_dispersion = poisson_disp, quasi_poisson_pearson_dispersion = quasi$pearson_dispersion[1], nb_pearson_dispersion = nb_disp, poisson_deviance_explained = poisson$deviance_explained[1], quasi_poisson_deviance_explained = quasi$deviance_explained[1], nb_deviance_explained = nb$deviance_explained[1], recommendation = recommendation)
}

count_patients <- diagnose_count_stage(dat_model, "patients", "n_patients", "pop", formula_pat)
count_provision <- diagnose_count_stage(dat_model, "provision", "n_provision", "n_patients", formula_prov)

count_family_quality <- dplyr::bind_rows(count_patients, count_provision) %>%
  dplyr::mutate(family = factor(.data$family, levels = c("Poisson", "quasi-Poisson", "Negative Binomial"))) %>%
  dplyr::arrange(.data$stage, .data$family)
count_family_recommendation <- dplyr::bind_rows(recommend_count_family(count_patients), recommend_count_family(count_provision))

readr::write_csv(count_family_quality, file.path(output_tables, "02_count_family_diagnostics_model_quality.csv"))
readr::write_csv(count_family_recommendation, file.path(output_tables, "03_count_family_diagnostics_recommendation.csv"))

# ------------------------------------------------------------
# 8. Diagnostyka value: Gamma vs Tweedie
# ------------------------------------------------------------

fit_value_family <- function(family_name) {
  if (family_name == "Gamma") {
    data_use <- dat_model %>% dplyr::filter(.data$n_provision > 0, .data$value > 0)
    fam <- Gamma(link = "log")
  } else if (family_name == "Tweedie") {
    data_use <- dat_model %>% dplyr::filter(.data$n_provision > 0, .data$value >= 0)
    fam <- mgcv::tw(link = "log")
  } else stop("Nieznana rodzina value: ", family_name, call. = FALSE)
  fit <- tryCatch(mgcv::bam(formula = formula_val, data = data_use, family = fam, method = "fREML", discrete = use_discrete_bam), error = function(e) { message("Blad value family ", family_name, ": ", conditionMessage(e)); NULL })
  list(fit = fit, data = data_use)
}

value_gamma <- list(fit = m_val_gamma, data = dat_val_gamma)
value_tweedie <- fit_value_family("Tweedie")

value_family_quality <- dplyr::bind_rows(
  tibble::tibble(specification = "Gamma, log link, offset log(services)", family = "Gamma", n_fitted_cells = nrow(value_gamma$data), n_excluded_zero_value_positive_services = n_value_zero_excluded, aic = safe_aic(value_gamma$fit), deviance_explained = safe_dev_expl(value_gamma$fit), pearson_dispersion = safe_pearson_dispersion(value_gamma$fit)),
  tibble::tibble(specification = "Tweedie, log link, offset log(services)", family = "Tweedie", n_fitted_cells = nrow(value_tweedie$data), n_excluded_zero_value_positive_services = 0L, aic = if (!is.null(value_tweedie$fit)) safe_aic(value_tweedie$fit) else NA_real_, deviance_explained = if (!is.null(value_tweedie$fit)) safe_dev_expl(value_tweedie$fit) else NA_real_, pearson_dispersion = if (!is.null(value_tweedie$fit)) safe_pearson_dispersion(value_tweedie$fit) else NA_real_)
)
readr::write_csv(value_family_quality, file.path(output_tables, "04_value_family_gamma_vs_tweedie.csv"))

# ------------------------------------------------------------
# 9. Predykcja komponentow jako rates
# ------------------------------------------------------------

predict_components <- function(df, m_pat, m_prov, m_val) {
  nd_pi <- df; nd_pi$pop <- 1
  nd_lambda <- df; nd_lambda$n_patients <- 1
  nd_kappa <- df; nd_kappa$n_provision <- 1
  pi_hat <- as.numeric(predict(m_pat, newdata = nd_pi, type = "response"))
  lambda_hat <- as.numeric(predict(m_prov, newdata = nd_lambda, type = "response"))
  kappa_hat <- as.numeric(predict(m_val, newdata = nd_kappa, type = "response"))
  df %>% dplyr::mutate(N = .data$pop, PI = pmax(pi_hat, 0), LAMBDA = pmax(lambda_hat, 0), KAPPA = pmax(kappa_hat, 0), Y_model = .data$N * .data$PI * .data$LAMBDA * .data$KAPPA)
}

components_df <- predict_components(dat_model, m_pat_nb, m_prov_nb, m_val_gamma)
saveRDS(components_df, file.path(output_rds, "components_predicted_gamma_chain.rds"))

# ------------------------------------------------------------
# 10. Model-based Shapley
# ------------------------------------------------------------

all_permutations <- function(x) {
  if (length(x) == 1) return(list(x))
  out <- list()
  for (i in seq_along(x)) {
    rest <- x[-i]
    sub <- all_permutations(rest)
    out <- c(out, lapply(sub, function(s) c(x[i], s)))
  }
  out
}

scenario_value_from_components <- function(d, active_set) {
  N <- if ("N" %in% active_set) d$N_1 else d$N_0
  PI <- if ("PI" %in% active_set) d$PI_1 else d$PI_0
  LAMBDA <- if ("LAMBDA" %in% active_set) d$LAMBDA_1 else d$LAMBDA_0
  KAPPA <- if ("KAPPA" %in% active_set) d$KAPPA_1 else d$KAPPA_0
  sum(N * PI * LAMBDA * KAPPA, na.rm = TRUE)
}

compute_model_shapley <- function(component_df, year0, year1) {
  comps <- component_order
  id_vars <- intersect(c("code", "sex", "age_group", "sex_age", "ICD"), names(component_df))
  keep_vars <- c(id_vars, comps)
  d0 <- component_df %>% dplyr::filter(.data$year == year0) %>% dplyr::select(dplyr::all_of(keep_vars))
  d1 <- component_df %>% dplyr::filter(.data$year == year1) %>% dplyr::select(dplyr::all_of(keep_vars))
  d <- dplyr::full_join(d0, d1, by = id_vars, suffix = c("_0", "_1"))
  for (cc in comps) {
    d[[paste0(cc, "_0")]] <- dplyr::coalesce(d[[paste0(cc, "_0")]], 0)
    d[[paste0(cc, "_1")]] <- dplyr::coalesce(d[[paste0(cc, "_1")]], 0)
  }
  perms <- all_permutations(comps)
  contrib_mat <- matrix(0, nrow = length(perms), ncol = length(comps), dimnames = list(NULL, comps))
  for (i in seq_along(perms)) {
    perm <- perms[[i]]; active <- character(0)
    for (comp in perm) {
      v_prev <- scenario_value_from_components(d, active)
      active <- c(active, comp)
      v_curr <- scenario_value_from_components(d, active)
      contrib_mat[i, comp] <- v_curr - v_prev
    }
  }
  shapley <- colMeans(contrib_mat)
  v0 <- scenario_value_from_components(d, character(0))
  v1 <- scenario_value_from_components(d, comps)
  tibble::tibble(component = comps, component_label = component_labels[comps], component_label_short = component_labels_short[comps], contribution = as.numeric(shapley), contribution_bln = contribution / 1e9, model_y0 = v0, model_y1 = v1, model_y0_bln = v0 / 1e9, model_y1_bln = v1 / 1e9, model_change = v1 - v0, model_change_bln = (v1 - v0) / 1e9, contribution_share_pct = safe_pct(contribution, v1 - v0), residual = (v1 - v0) - sum(shapley), n_cells_joined = nrow(d), year0 = year0, year1 = year1)
}

shapley_updated <- compute_model_shapley(components_df, year0, year1)
readr::write_csv(shapley_updated, file.path(output_tables, "05_empirical_shapley_updated.csv"))

# ------------------------------------------------------------
# 11. Aggregate LMDI-style benchmark
# ------------------------------------------------------------

lmdi_weight <- function(y1, y0) {
  if (is.na(y1) || is.na(y0)) return(NA_real_)
  if (y1 > 0 && y0 > 0 && y1 != y0) return((y1 - y0) / (log(y1) - log(y0)))
  if (y1 > 0 && y0 > 0 && y1 == y0) return(y1)
  NA_real_
}

compute_observed_components <- function(df) {
  df %>% dplyr::mutate(N = .data$pop, PI = dplyr::if_else(.data$pop > 0, .data$n_patients / .data$pop, 0), LAMBDA = dplyr::if_else(.data$n_patients > 0, .data$n_provision / .data$n_patients, 0), KAPPA = dplyr::if_else(.data$n_provision > 0, .data$value / .data$n_provision, 0), Y_obs = .data$N * .data$PI * .data$LAMBDA * .data$KAPPA)
}

observed_components <- compute_observed_components(dat_model)
agg_obs <- observed_components %>%
  dplyr::filter(.data$year %in% c(year0, year1)) %>%
  dplyr::group_by(.data$year) %>%
  dplyr::summarise(Y = sum(.data$Y_obs, na.rm = TRUE), N = sum(.data$N, na.rm = TRUE), patients = sum(.data$n_patients, na.rm = TRUE), provision = sum(.data$n_provision, na.rm = TRUE), value = sum(.data$value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(PI = .data$patients / .data$N, LAMBDA = .data$provision / .data$patients, KAPPA = .data$value / .data$provision)

a0 <- agg_obs %>% dplyr::filter(.data$year == year0)
a1 <- agg_obs %>% dplyr::filter(.data$year == year1)
Y0 <- a0$Y[1]; Y1 <- a1$Y[1]
L_weight <- lmdi_weight(Y1, Y0)
lmdi_contrib <- c(N = L_weight * log(a1$N[1] / a0$N[1]), PI = L_weight * log(a1$PI[1] / a0$PI[1]), LAMBDA = L_weight * log(a1$LAMBDA[1] / a0$LAMBDA[1]), KAPPA = L_weight * log(a1$KAPPA[1] / a0$KAPPA[1]))

lmdi_updated <- tibble::tibble(component = component_order, component_label = component_labels[component_order], component_label_short = component_labels_short[component_order], contribution = as.numeric(lmdi_contrib[component_order]), contribution_bln = contribution / 1e9, observed_y0 = Y0, observed_y1 = Y1, observed_y0_bln = Y0 / 1e9, observed_y1_bln = Y1 / 1e9, observed_change = Y1 - Y0, observed_change_bln = (Y1 - Y0) / 1e9, contribution_share_pct = safe_pct(contribution, Y1 - Y0), residual = (Y1 - Y0) - sum(contribution), year0 = year0, year1 = year1)
readr::write_csv(lmdi_updated, file.path(output_tables, "06_empirical_lmdi_updated.csv"))

comparison_updated <- lmdi_updated %>%
  dplyr::select(.data$component, .data$component_label, lmdi_bln = .data$contribution_bln, lmdi_share_pct = .data$contribution_share_pct) %>%
  dplyr::left_join(shapley_updated %>% dplyr::select(.data$component, shapley_bln = .data$contribution_bln, shapley_share_pct = .data$contribution_share_pct), by = "component") %>%
  dplyr::arrange(factor(.data$component, levels = c("KAPPA", "LAMBDA", "N", "PI")))
readr::write_csv(comparison_updated, file.path(output_tables, "07_empirical_lmdi_vs_shapley_updated.csv"))

# ------------------------------------------------------------
# 12. Aktualizacja S1B i liczb do tekstu
# ------------------------------------------------------------

observed_y0_bln <- Y0 / 1e9; observed_y1_bln <- Y1 / 1e9; observed_change_bln <- (Y1 - Y0) / 1e9
model_y0 <- unique(shapley_updated$model_y0); model_y1 <- unique(shapley_updated$model_y1); model_change <- unique(shapley_updated$model_change)
model_y0_bln <- model_y0 / 1e9; model_y1_bln <- model_y1 / 1e9; model_change_bln <- model_change / 1e9
relative_difference_vs_observed_change_pct <- 100 * (model_change - (Y1 - Y0)) / abs(Y1 - Y0)

supp_s1b_updated <- tibble::tibble(
  quantity = c("ICD group", "Years available", "Patient component explained deviance", "Service component explained deviance", "Point-valuation component explained deviance", "Observed baseline total, billion points", "Observed comparison total, billion points", "Observed change, billion points", "Model-implied baseline total, billion points", "Model-implied comparison total, billion points", "Model-implied change, billion points", "Relative difference vs observed change"),
  value = c(paste(unique(dat_model$ICD), collapse = ", "), as.character(length(years_available)), sprintf("%.3f", main_model_quality$deviance_explained[main_model_quality$stage == "patients"]), sprintf("%.3f", main_model_quality$deviance_explained[main_model_quality$stage == "provision"]), sprintf("%.3f", main_model_quality$deviance_explained[main_model_quality$stage == "value"]), sprintf("%.2f", observed_y0_bln), sprintf("%.2f", observed_y1_bln), sprintf("%.2f", observed_change_bln), sprintf("%.2f", model_y0_bln), sprintf("%.2f", model_y1_bln), sprintf("%.2f", model_change_bln), paste0(sprintf("%.2f", relative_difference_vs_observed_change_pct), "%"))
)
readr::write_csv(supp_s1b_updated, file.path(output_tables, "08_supplement_s1b_values_updated.csv"))

get_comp_bln <- function(df, comp) df %>% dplyr::filter(.data$component == comp) %>% dplyr::pull(.data$contribution_bln)
get_comp_share <- function(df, comp) df %>% dplyr::filter(.data$component == comp) %>% dplyr::pull(.data$contribution_share_pct)

text_numbers <- tibble::tibble(
  item = c("observed_y0_bln", "observed_y1_bln", "observed_change_bln", "model_change_bln", "model_relative_difference_vs_observed_change_pct", "shapley_N_bln", "shapley_PI_bln", "shapley_LAMBDA_bln", "shapley_KAPPA_bln", "shapley_N_share_pct", "shapley_PI_share_pct", "shapley_LAMBDA_share_pct", "shapley_KAPPA_share_pct", "lmdi_N_bln", "lmdi_PI_bln", "lmdi_LAMBDA_bln", "lmdi_KAPPA_bln", "dev_patients", "dev_provision", "dev_value_gamma", "gamma_aic", "tweedie_aic", "gamma_dev_expl", "tweedie_dev_expl", "value_zero_cells_excluded", "value_zero_cells_excluded_pct"),
  value = c(observed_y0_bln, observed_y1_bln, observed_change_bln, model_change_bln, relative_difference_vs_observed_change_pct, get_comp_bln(shapley_updated, "N"), get_comp_bln(shapley_updated, "PI"), get_comp_bln(shapley_updated, "LAMBDA"), get_comp_bln(shapley_updated, "KAPPA"), get_comp_share(shapley_updated, "N"), get_comp_share(shapley_updated, "PI"), get_comp_share(shapley_updated, "LAMBDA"), get_comp_share(shapley_updated, "KAPPA"), get_comp_bln(lmdi_updated, "N"), get_comp_bln(lmdi_updated, "PI"), get_comp_bln(lmdi_updated, "LAMBDA"), get_comp_bln(lmdi_updated, "KAPPA"), main_model_quality$deviance_explained[main_model_quality$stage == "patients"], main_model_quality$deviance_explained[main_model_quality$stage == "provision"], main_model_quality$deviance_explained[main_model_quality$stage == "value"], value_family_quality$aic[value_family_quality$family == "Gamma"], value_family_quality$aic[value_family_quality$family == "Tweedie"], value_family_quality$deviance_explained[value_family_quality$family == "Gamma"], value_family_quality$deviance_explained[value_family_quality$family == "Tweedie"], n_value_zero_excluded, 100 * n_value_zero_excluded / nrow(dat_model))
)
readr::write_csv(text_numbers, file.path(output_tables, "09_text_numbers_updated.csv"))

# ------------------------------------------------------------
# 13. Wykresy
# ------------------------------------------------------------

plot_count_aic <- count_family_quality %>% dplyr::filter(.data$family %in% c("Poisson", "Negative Binomial"), !is.na(.data$aic))
if (nrow(plot_count_aic) > 0) {
  p_count_aic <- ggplot2::ggplot(plot_count_aic, ggplot2::aes(x = .data$stage, y = .data$aic, fill = .data$family)) + ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65) + ggplot2::scale_y_continuous(labels = label_pub_number(accuracy = 1)) + ggplot2::labs(x = "Model stage", y = "AIC") + theme_publication(base_size = 10)
  ggplot2::ggsave(file.path(output_figures, "count_family_aic.png"), p_count_aic, width = 7, height = 4.5, dpi = 600)
}

plot_count_disp <- count_family_quality %>% dplyr::filter(!is.na(.data$pearson_dispersion))
if (nrow(plot_count_disp) > 0) {
  p_count_disp <- ggplot2::ggplot(plot_count_disp, ggplot2::aes(x = .data$stage, y = .data$pearson_dispersion, fill = .data$family)) + ggplot2::geom_hline(yintercept = 1, linewidth = 0.35, colour = "grey40") + ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65) + ggplot2::scale_y_continuous(labels = label_pub_number(accuracy = 0.01)) + ggplot2::labs(x = "Model stage", y = "Pearson dispersion") + theme_publication(base_size = 10)
  ggplot2::ggsave(file.path(output_figures, "count_family_dispersion.png"), p_count_disp, width = 7, height = 4.5, dpi = 600)
}

plot_value_aic <- value_family_quality %>% dplyr::filter(!is.na(.data$aic))
if (nrow(plot_value_aic) > 0) {
  p_value_aic <- ggplot2::ggplot(plot_value_aic, ggplot2::aes(x = .data$family, y = .data$aic)) + ggplot2::geom_col(width = 0.65) + ggplot2::scale_y_continuous(labels = label_pub_number(accuracy = 1)) + ggplot2::labs(x = "Value-stage family", y = "AIC") + theme_publication(base_size = 10)
  ggplot2::ggsave(file.path(output_figures, "value_family_aic.png"), p_value_aic, width = 6, height = 4.5, dpi = 600)
}

plot_comp <- dplyr::bind_rows(lmdi_updated %>% dplyr::mutate(method = "Aggregate LMDI-style"), shapley_updated %>% dplyr::mutate(method = "Model-based Shapley")) %>% dplyr::mutate(component_label = factor(.data$component_label, levels = component_labels[c("KAPPA", "LAMBDA", "N", "PI")]), method = factor(.data$method, levels = c("Aggregate LMDI-style", "Model-based Shapley")))
p_comp <- ggplot2::ggplot(plot_comp, ggplot2::aes(x = .data$component_label, y = .data$contribution_bln, fill = .data$method)) + ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65) + ggplot2::coord_flip() + ggplot2::scale_y_continuous(labels = label_pub_number(accuracy = 0.01)) + ggplot2::labs(x = NULL, y = "Contribution, billion points") + theme_publication(base_size = 10)
ggplot2::ggsave(file.path(output_figures, "empirical_lmdi_vs_shapley_updated.png"), p_comp, width = 8, height = 5, dpi = 600)

# ------------------------------------------------------------
# 14. Podsumowanie
# ------------------------------------------------------------

message("============================================================")
message("GOTOWE — najnowsze przeliczenie Gamma + diagnostyka rodzin")
message("============================================================")
message("Tabele: ", output_tables)
message("Figury: ", output_figures)
message("RDS: ", output_rds)
message("============================================================")
print(main_model_quality)
print(count_family_recommendation)
print(value_family_quality)
print(supp_s1b_updated)
print(comparison_updated)
