# ============================================================
# scripts/run_13_full_family_selection_and_refresh_results.R
# ------------------------------------------------------------
# Pełna diagnostyka rodzin rozkładu + dobór finalnego modelu.
#
# Etapy:
#   patients  : Poisson, quasi-Poisson, Negative Binomial,
#               Tweedie sensitivity
#   provision : Poisson, quasi-Poisson, Negative Binomial,
#               Tweedie sensitivity
#   value     : Gamma, Tweedie, Inverse Gaussian
#
# Finalny model chain:
#   wybierany z diagnostyki, domyślnie:
#   - count stages: tylko Poisson albo Negative Binomial jako finalne
#   - value stage : Gamma, Tweedie albo Inverse Gaussian
#
# Wyjście:
#   output/tables/R13_family_selection/
#   output/figures/R13_family_selection/
#   output/rds/R13_family_selection/
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

output_tables <- "output/tables/R13_family_selection"
output_figures <- "output/figures/R13_family_selection"
output_rds <- "output/rds/R13_family_selection"

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

# Dla liczników Tweedie jest raportowany jako sensitivity, ale domyślnie
# nie będzie wybrany jako finalna rodzina. Zmień na TRUE, jeśli chcesz
# pozwolić, by Tweedie wygrał też dla patients/provision.
allow_tweedie_for_count_final <- FALSE

# Reguły selection dla count stages.
delta_aic_strong <- 10
poisson_overdispersion_threshold <- 1.5

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
# 2. Helpery
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

# Wektorowa wersja: działa zarówno dla skalara w mianowniku
# i wektora w liczniku, np. 4 komponenty / 1 total_change.
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
      "\nDostępne kolumny:\n", paste(names(df), collapse = ", "),
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
  if (length(hit) == 0) NA_character_ else hit[1]
}

detect_value_col <- function(df) {
  candidates <- c("value", "sum_points", "points", "liczba_punktow", "liczba_punktów")
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    stop("Nie znaleziono kolumny value/sum_points.\nDostępne kolumny:\n", paste(names(df), collapse = ", "), call. = FALSE)
  }
  hit[1]
}

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
      } else {
        return(normalize_schema(df_raw))
      }
    }
  }
  if (file.exists(fallback_rds)) {
    message("Wczytuję fallback RDS: ", fallback_rds)
    return(readRDS(fallback_rds))
  }
  stop("Nie znaleziono danych. Oczekiwano config.R albo ", fallback_rds, call. = FALSE)
}

# ------------------------------------------------------------
# 3. Dane
# ------------------------------------------------------------

dat_raw <- read_main_data()
message("Dane wejściowe: ", nrow(dat_raw), " wierszy")
print(names(dat_raw))

value_col <- detect_value_col(dat_raw)
icd_col <- detect_icd_col(dat_raw)

required_base <- c("year", "code", "sex", "age_group", "pop", "n_patients", "n_provision", value_col)
stop_if_missing_cols(dat_raw, required_base)

dat <- dat_raw %>%
  mutate(
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
  filter(!is.na(year), !is.na(code), !is.na(sex), !is.na(age_group),
         !is.na(pop), !is.na(n_patients), !is.na(n_provision), !is.na(value))

if (!is.null(disease_filter) && !is.na(icd_col)) {
  dat <- dat %>% filter(.data$ICD == disease_filter)
  if (nrow(dat) == 0) {
    stop("Po disease_filter = '", disease_filter, "' nie zostały żadne wiersze.", call. = FALSE)
  }
}

dat_model <- dat %>%
  group_by(year, code, sex, age_group, ICD) %>%
  summarise(
    pop = safe_sum(pop),
    n_patients = safe_sum(n_patients),
    n_provision = safe_sum(n_provision),
    value = safe_sum(value),
    .groups = "drop"
  ) %>%
  mutate(
    year_num = as.numeric(year),
    code = factor(code),
    sex = factor(sex),
    age_group = factor(age_group),
    sex_age = interaction(sex, age_group, drop = TRUE)
  )

years_available <- sort(unique(dat_model$year))
year0 <- if (is.null(year0_override)) min(years_available) else year0_override
year1 <- if (is.null(year1_override)) max(years_available) else year1_override
if (!all(c(year0, year1) %in% years_available)) stop("year0/year1 nie występują w danych.", call. = FALSE)

k_year <- min(k_year_default, max(3, length(years_available) - 1))

saveRDS(dat_model, file.path(output_rds, "dat_model.rds"))

message("Dane modelowe:")
print(dat_model %>% summarise(
  ICD = paste(unique(ICD), collapse = ", "),
  n_rows = n(), n_years = n_distinct(year), year0 = year0, year1 = year1,
  n_codes = n_distinct(code), n_sex_age = n_distinct(sex_age),
  total_pop = sum(pop), total_patients = sum(n_patients),
  total_provision = sum(n_provision), total_value = sum(value)
))

# ------------------------------------------------------------
# 4. Kandydaci rodzin
# ------------------------------------------------------------

make_rhs <- function(include_code_re = TRUE) {
  rhs <- paste0("sex * age_group + s(year_num, by = sex_age, k = ", k_year, ")")
  if (include_code_re) rhs <- paste0(rhs, " + s(code, bs = 're')")
  rhs
}

rhs_main <- make_rhs(include_code_re = include_code_re)

make_formula <- function(outcome, offset_var) {
  as.formula(paste0(outcome, " ~ ", rhs_main, " + offset(log(", offset_var, "))"))
}

formula_pat <- make_formula("n_patients", "pop")
formula_prov <- make_formula("n_provision", "n_patients")
formula_val <- make_formula("value", "n_provision")

family_object <- function(family_name) {
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

filter_for_family <- function(data, outcome, offset_var, family_name) {
  base <- data %>%
    filter(!is.na(.data[[outcome]]), !is.na(.data[[offset_var]]), .data[[offset_var]] > 0)
  if (family_name %in% c("Gamma", "Inverse Gaussian")) {
    base %>% filter(.data[[outcome]] > 0)
  } else {
    base %>% filter(.data[[outcome]] >= 0)
  }
}

fit_candidate <- function(stage, outcome, offset_var, formula, family_name) {
  data_use <- filter_for_family(dat_model, outcome, offset_var, family_name)
  n_excluded <- nrow(dat_model %>% filter(!is.na(.data[[offset_var]]), .data[[offset_var]] > 0)) - nrow(data_use)

  message("------------------------------------------------------------")
  message("Stage: ", stage, " | Family: ", family_name, " | n = ", nrow(data_use))
  message("------------------------------------------------------------")

  fit <- tryCatch(
    mgcv::bam(
      formula = formula,
      data = data_use,
      family = family_object(family_name),
      method = "fREML",
      discrete = use_discrete_bam
    ),
    error = function(e1) {
      message("bam() failed: ", conditionMessage(e1))
      tryCatch(
        mgcv::gam(
          formula = formula,
          data = data_use,
          family = family_object(family_name),
          method = "REML"
        ),
        error = function(e2) {
          message("gam() failed: ", conditionMessage(e2))
          NULL
        }
      )
    }
  )

  if (is.null(fit)) {
    return(list(
      fit = NULL,
      data = data_use,
      metrics = tibble(
        stage = stage, outcome = outcome, offset = paste0("log(", offset_var, ")"),
        family = family_name, n_fitted = nrow(data_use), n_excluded = n_excluded,
        aic = NA_real_, deviance_explained = NA_real_, pearson_dispersion = NA_real_,
        rmse_response = NA_real_, mae_response = NA_real_, theta_nb = NA_real_, converged = FALSE
      )
    ))
  }

  list(
    fit = fit,
    data = data_use,
    metrics = tibble(
      stage = stage,
      outcome = outcome,
      offset = paste0("log(", offset_var, ")"),
      family = family_name,
      n_fitted = nrow(data_use),
      n_excluded = n_excluded,
      aic = safe_aic(fit),
      deviance_explained = safe_dev_expl(fit),
      pearson_dispersion = safe_pearson_dispersion(fit),
      rmse_response = safe_rmse_response(fit, data_use, outcome),
      mae_response = safe_mae_response(fit, data_use, outcome),
      theta_nb = if (family_name == "Negative Binomial") safe_theta(fit) else NA_real_,
      converged = TRUE
    )
  )
}

count_candidates <- c("Poisson", "quasi-Poisson", "Negative Binomial", "Tweedie")
value_candidates <- c("Gamma", "Tweedie", "Inverse Gaussian")

fits_pat <- setNames(lapply(count_candidates, function(fam) {
  fit_candidate("patients", "n_patients", "pop", formula_pat, fam)
}), count_candidates)

fits_prov <- setNames(lapply(count_candidates, function(fam) {
  fit_candidate("provision", "n_provision", "n_patients", formula_prov, fam)
}), count_candidates)

fits_val <- setNames(lapply(value_candidates, function(fam) {
  fit_candidate("value", "value", "n_provision", formula_val, fam)
}), value_candidates)

family_diagnostics <- bind_rows(
  map_dfr(fits_pat, "metrics"),
  map_dfr(fits_prov, "metrics"),
  map_dfr(fits_val, "metrics")
) %>% arrange(stage, family)

write_csv(family_diagnostics, file.path(output_tables, "01_family_diagnostics_all_candidates.csv"))

# ------------------------------------------------------------
# 5. Wybór rodziny
# ------------------------------------------------------------

select_count_family <- function(diag_stage) {
  poisson <- diag_stage %>% filter(family == "Poisson")
  nb <- diag_stage %>% filter(family == "Negative Binomial")
  tw <- diag_stage %>% filter(family == "Tweedie")

  eligible <- diag_stage %>%
    filter(family %in% c("Poisson", "Negative Binomial", if (allow_tweedie_for_count_final) "Tweedie" else character(0))) %>%
    filter(converged, !is.na(aic))

  best_aic_family <- if (nrow(eligible) > 0) eligible$family[which.min(eligible$aic)] else NA_character_

  poisson_disp <- poisson$pearson_dispersion[1]
  poisson_aic <- poisson$aic[1]
  nb_aic <- nb$aic[1]
  delta_nb_minus_poisson <- nb_aic - poisson_aic

  selected <- dplyr::case_when(
    !is.na(best_aic_family) & best_aic_family == "Tweedie" & allow_tweedie_for_count_final ~ "Tweedie",
    !is.na(delta_nb_minus_poisson) & delta_nb_minus_poisson < -2 ~ "Negative Binomial",
    !is.na(poisson_disp) & poisson_disp > poisson_overdispersion_threshold & !is.na(nb_aic) ~ "Negative Binomial",
    !is.na(best_aic_family) ~ best_aic_family,
    TRUE ~ "Negative Binomial"
  )

  rationale <- dplyr::case_when(
    selected == "Negative Binomial" & !is.na(delta_nb_minus_poisson) & delta_nb_minus_poisson < -delta_aic_strong ~
      "Negative Binomial selected: substantially lower AIC than Poisson.",
    selected == "Negative Binomial" & !is.na(poisson_disp) & poisson_disp > poisson_overdispersion_threshold ~
      "Negative Binomial selected: Poisson showed overdispersion and NB is count-appropriate.",
    selected == "Poisson" ~
      "Poisson selected: no strong AIC or dispersion evidence requiring NB.",
    selected == "Tweedie" ~
      "Tweedie selected because allow_tweedie_for_count_final=TRUE and it had best AIC.",
    TRUE ~ "Defaulted to Negative Binomial for count-stage robustness."
  )

  tibble(
    stage = unique(diag_stage$stage),
    selected_family = selected,
    best_aic_family_among_allowed = best_aic_family,
    poisson_aic = poisson_aic,
    nb_aic = nb_aic,
    tweedie_aic = tw$aic[1],
    delta_aic_nb_minus_poisson = delta_nb_minus_poisson,
    poisson_pearson_dispersion = poisson_disp,
    nb_pearson_dispersion = nb$pearson_dispersion[1],
    tweedie_pearson_dispersion = tw$pearson_dispersion[1],
    rationale = rationale
  )
}

select_value_family <- function(diag_stage) {
  eligible <- diag_stage %>% filter(converged, !is.na(aic))
  selected <- if (nrow(eligible) > 0) eligible$family[which.min(eligible$aic)] else "Gamma"
  best_dev <- if (nrow(eligible) > 0) eligible$family[which.max(eligible$deviance_explained)] else NA_character_
  tibble(
    stage = unique(diag_stage$stage),
    selected_family = selected,
    best_aic_family_among_allowed = selected,
    gamma_aic = diag_stage$aic[diag_stage$family == "Gamma"][1],
    tweedie_aic = diag_stage$aic[diag_stage$family == "Tweedie"][1],
    inverse_gaussian_aic = diag_stage$aic[diag_stage$family == "Inverse Gaussian"][1],
    gamma_deviance_explained = diag_stage$deviance_explained[diag_stage$family == "Gamma"][1],
    tweedie_deviance_explained = diag_stage$deviance_explained[diag_stage$family == "Tweedie"][1],
    inverse_gaussian_deviance_explained = diag_stage$deviance_explained[diag_stage$family == "Inverse Gaussian"][1],
    best_deviance_family = best_dev,
    rationale = paste0(selected, " selected by minimum AIC among positive/nonnegative value-stage families.")
  )
}

selection_recommendations <- bind_rows(
  select_count_family(family_diagnostics %>% filter(stage == "patients")),
  select_count_family(family_diagnostics %>% filter(stage == "provision")),
  select_value_family(family_diagnostics %>% filter(stage == "value"))
)

write_csv(selection_recommendations, file.path(output_tables, "02_family_selection_recommendations.csv"))

selected_pat <- selection_recommendations$selected_family[selection_recommendations$stage == "patients"]
selected_prov <- selection_recommendations$selected_family[selection_recommendations$stage == "provision"]
selected_val <- selection_recommendations$selected_family[selection_recommendations$stage == "value"]

message("Wybrane rodziny:")
print(selection_recommendations)

get_fit <- function(fits, family_name) fits[[family_name]]$fit

m_pat <- get_fit(fits_pat, selected_pat)
m_prov <- get_fit(fits_prov, selected_prov)
m_val <- get_fit(fits_val, selected_val)

if (is.null(m_pat) || is.null(m_prov) || is.null(m_val)) stop("Co najmniej jeden finalny model jest NULL.", call. = FALSE)

final_model_quality <- family_diagnostics %>%
  semi_join(selection_recommendations %>% select(stage, selected_family), by = c("stage" = "stage", "family" = "selected_family")) %>%
  arrange(stage)

write_csv(final_model_quality, file.path(output_tables, "03_final_model_quality.csv"))

saveRDS(list(m_pat = m_pat, m_prov = m_prov, m_val = m_val,
             selected_pat = selected_pat, selected_prov = selected_prov, selected_val = selected_val,
             dat_model = dat_model, year0 = year0, year1 = year1, k_year = k_year),
        file.path(output_rds, "final_selected_model_chain.rds"))

# ------------------------------------------------------------
# 6. Predykcja komponentów i Shapley
# ------------------------------------------------------------

predict_components <- function(df, m_pat, m_prov, m_val) {
  nd_pi <- df; nd_pi$pop <- 1
  nd_lambda <- df; nd_lambda$n_patients <- 1
  nd_kappa <- df; nd_kappa$n_provision <- 1

  pi_hat <- as.numeric(predict(m_pat, newdata = nd_pi, type = "response"))
  lambda_hat <- as.numeric(predict(m_prov, newdata = nd_lambda, type = "response"))
  kappa_hat <- as.numeric(predict(m_val, newdata = nd_kappa, type = "response"))

  df %>% mutate(
    N = pop,
    PI = pmax(pi_hat, 0),
    LAMBDA = pmax(lambda_hat, 0),
    KAPPA = pmax(kappa_hat, 0),
    Y_model = N * PI * LAMBDA * KAPPA
  )
}

components_df <- predict_components(dat_model, m_pat, m_prov, m_val)
saveRDS(components_df, file.path(output_rds, "components_predicted_selected_chain.rds"))

all_permutations <- function(x) {
  if (length(x) == 1) return(list(x))
  out <- list()
  for (i in seq_along(x)) {
    rest <- x[-i]
    out <- c(out, lapply(all_permutations(rest), function(s) c(x[i], s)))
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

  d0 <- component_df %>% filter(year == year0) %>% select(all_of(keep_vars))
  d1 <- component_df %>% filter(year == year1) %>% select(all_of(keep_vars))

  d <- full_join(d0, d1, by = id_vars, suffix = c("_0", "_1"))
  for (cc in comps) {
    d[[paste0(cc, "_0")]] <- coalesce(d[[paste0(cc, "_0")]], 0)
    d[[paste0(cc, "_1")]] <- coalesce(d[[paste0(cc, "_1")]], 0)
  }

  perms <- all_permutations(comps)
  contrib_mat <- matrix(0, nrow = length(perms), ncol = length(comps), dimnames = list(NULL, comps))

  for (i in seq_along(perms)) {
    active <- character(0)
    for (comp in perms[[i]]) {
      v_prev <- scenario_value_from_components(d, active)
      active <- c(active, comp)
      v_curr <- scenario_value_from_components(d, active)
      contrib_mat[i, comp] <- v_curr - v_prev
    }
  }

  shapley <- colMeans(contrib_mat)
  v0 <- scenario_value_from_components(d, character(0))
  v1 <- scenario_value_from_components(d, comps)

  tibble(
    component = comps,
    component_label = component_labels[comps],
    component_label_short = component_labels_short[comps],
    contribution = as.numeric(shapley),
    contribution_bln = contribution / 1e9,
    model_y0 = v0,
    model_y1 = v1,
    model_y0_bln = v0 / 1e9,
    model_y1_bln = v1 / 1e9,
    model_change = v1 - v0,
    model_change_bln = (v1 - v0) / 1e9,
    contribution_share_pct = safe_pct(contribution, v1 - v0),
    residual = (v1 - v0) - sum(shapley),
    n_cells_joined = nrow(d),
    year0 = year0,
    year1 = year1
  )
}

shapley_updated <- compute_model_shapley(components_df, year0, year1)
write_csv(shapley_updated, file.path(output_tables, "04_empirical_shapley_updated.csv"))

# ------------------------------------------------------------
# 7. LMDI benchmark
# ------------------------------------------------------------

lmdi_weight <- function(y1, y0) {
  if (is.na(y1) || is.na(y0)) return(NA_real_)
  if (y1 > 0 && y0 > 0 && y1 != y0) return((y1 - y0) / (log(y1) - log(y0)))
  if (y1 > 0 && y0 > 0 && y1 == y0) return(y1)
  NA_real_
}

observed_components <- dat_model %>% mutate(
  N = pop,
  PI = if_else(pop > 0, n_patients / pop, 0),
  LAMBDA = if_else(n_patients > 0, n_provision / n_patients, 0),
  KAPPA = if_else(n_provision > 0, value / n_provision, 0),
  Y_obs = N * PI * LAMBDA * KAPPA
)

agg_obs <- observed_components %>%
  filter(year %in% c(year0, year1)) %>%
  group_by(year) %>%
  summarise(
    Y = sum(Y_obs, na.rm = TRUE),
    N = sum(N, na.rm = TRUE),
    patients = sum(n_patients, na.rm = TRUE),
    provision = sum(n_provision, na.rm = TRUE),
    value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(PI = patients / N, LAMBDA = provision / patients, KAPPA = value / provision)

a0 <- agg_obs %>% filter(year == year0)
a1 <- agg_obs %>% filter(year == year1)
Y0 <- a0$Y[1]
Y1 <- a1$Y[1]
L_weight <- lmdi_weight(Y1, Y0)

lmdi_contrib <- c(
  N = L_weight * log(a1$N[1] / a0$N[1]),
  PI = L_weight * log(a1$PI[1] / a0$PI[1]),
  LAMBDA = L_weight * log(a1$LAMBDA[1] / a0$LAMBDA[1]),
  KAPPA = L_weight * log(a1$KAPPA[1] / a0$KAPPA[1])
)

lmdi_updated <- tibble(
  component = component_order,
  component_label = component_labels[component_order],
  component_label_short = component_labels_short[component_order],
  contribution = as.numeric(lmdi_contrib[component_order]),
  contribution_bln = contribution / 1e9,
  observed_y0 = Y0,
  observed_y1 = Y1,
  observed_y0_bln = Y0 / 1e9,
  observed_y1_bln = Y1 / 1e9,
  observed_change = Y1 - Y0,
  observed_change_bln = (Y1 - Y0) / 1e9,
  contribution_share_pct = safe_pct(contribution, Y1 - Y0),
  residual = (Y1 - Y0) - sum(contribution),
  year0 = year0,
  year1 = year1
)

write_csv(lmdi_updated, file.path(output_tables, "05_empirical_lmdi_updated.csv"))

comparison_updated <- lmdi_updated %>%
  select(component, component_label, lmdi_bln = contribution_bln, lmdi_share_pct = contribution_share_pct) %>%
  left_join(shapley_updated %>% select(component, shapley_bln = contribution_bln, shapley_share_pct = contribution_share_pct), by = "component") %>%
  arrange(factor(component, levels = c("KAPPA", "LAMBDA", "N", "PI")))

write_csv(comparison_updated, file.path(output_tables, "06_empirical_lmdi_vs_shapley_updated.csv"))

# ------------------------------------------------------------
# 8. Aktualne liczby do supplementu i tekstu
# ------------------------------------------------------------

observed_y0_bln <- Y0 / 1e9
observed_y1_bln <- Y1 / 1e9
observed_change_bln <- (Y1 - Y0) / 1e9
model_y0 <- unique(shapley_updated$model_y0)
model_y1 <- unique(shapley_updated$model_y1)
model_change <- unique(shapley_updated$model_change)
model_change_bln <- model_change / 1e9
relative_difference_vs_observed_change_pct <- 100 * (model_change - (Y1 - Y0)) / abs(Y1 - Y0)

get_dev <- function(stage) final_model_quality$deviance_explained[final_model_quality$stage == stage][1]

supp_s1b_updated <- tibble(
  quantity = c(
    "ICD group", "Years available",
    "Patient component family", "Service component family", "Point-valuation component family",
    "Patient component explained deviance", "Service component explained deviance", "Point-valuation component explained deviance",
    "Observed baseline total, billion points", "Observed comparison total, billion points", "Observed change, billion points",
    "Model-implied baseline total, billion points", "Model-implied comparison total, billion points", "Model-implied change, billion points",
    "Relative difference vs observed change"
  ),
  value = c(
    paste(unique(dat_model$ICD), collapse = ", "), as.character(length(years_available)),
    selected_pat, selected_prov, selected_val,
    sprintf("%.3f", get_dev("patients")), sprintf("%.3f", get_dev("provision")), sprintf("%.3f", get_dev("value")),
    sprintf("%.2f", observed_y0_bln), sprintf("%.2f", observed_y1_bln), sprintf("%.2f", observed_change_bln),
    sprintf("%.2f", model_y0 / 1e9), sprintf("%.2f", model_y1 / 1e9), sprintf("%.2f", model_change_bln),
    paste0(sprintf("%.2f", relative_difference_vs_observed_change_pct), "%")
  )
)

write_csv(supp_s1b_updated, file.path(output_tables, "07_supplement_s1b_values_updated.csv"))

get_comp_bln <- function(df, comp) df %>% filter(component == comp) %>% pull(contribution_bln)
get_comp_share <- function(df, comp) df %>% filter(component == comp) %>% pull(contribution_share_pct)

text_numbers <- tibble(
  item = c(
    "selected_patients_family", "selected_provision_family", "selected_value_family",
    "observed_y0_bln", "observed_y1_bln", "observed_change_bln",
    "model_change_bln", "model_relative_difference_vs_observed_change_pct",
    "shapley_N_bln", "shapley_PI_bln", "shapley_LAMBDA_bln", "shapley_KAPPA_bln",
    "shapley_N_share_pct", "shapley_PI_share_pct", "shapley_LAMBDA_share_pct", "shapley_KAPPA_share_pct",
    "lmdi_N_bln", "lmdi_PI_bln", "lmdi_LAMBDA_bln", "lmdi_KAPPA_bln",
    "dev_patients", "dev_provision", "dev_value"
  ),
  value = c(
    selected_pat, selected_prov, selected_val,
    observed_y0_bln, observed_y1_bln, observed_change_bln,
    model_change_bln, relative_difference_vs_observed_change_pct,
    get_comp_bln(shapley_updated, "N"), get_comp_bln(shapley_updated, "PI"), get_comp_bln(shapley_updated, "LAMBDA"), get_comp_bln(shapley_updated, "KAPPA"),
    get_comp_share(shapley_updated, "N"), get_comp_share(shapley_updated, "PI"), get_comp_share(shapley_updated, "LAMBDA"), get_comp_share(shapley_updated, "KAPPA"),
    get_comp_bln(lmdi_updated, "N"), get_comp_bln(lmdi_updated, "PI"), get_comp_bln(lmdi_updated, "LAMBDA"), get_comp_bln(lmdi_updated, "KAPPA"),
    get_dev("patients"), get_dev("provision"), get_dev("value")
  )
)

write_csv(text_numbers, file.path(output_tables, "08_text_numbers_updated.csv"))

# ------------------------------------------------------------
# 9. Figury diagnostyczne
# ------------------------------------------------------------

p_aic <- family_diagnostics %>%
  filter(!is.na(aic)) %>%
  ggplot(aes(x = family, y = aic)) +
  geom_col(width = 0.65) +
  facet_wrap(~ stage, scales = "free_y") +
  scale_y_continuous(labels = label_pub_number(accuracy = 1)) +
  labs(x = NULL, y = "AIC") +
  theme_publication(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(output_figures, "family_diagnostics_aic_all_candidates.png"), p_aic, width = 9, height = 5, dpi = 600)

p_disp <- family_diagnostics %>%
  filter(!is.na(pearson_dispersion)) %>%
  ggplot(aes(x = family, y = pearson_dispersion)) +
  geom_hline(yintercept = 1, linewidth = 0.35, colour = "grey40") +
  geom_col(width = 0.65) +
  facet_wrap(~ stage, scales = "free_y") +
  scale_y_continuous(labels = label_pub_number(accuracy = 0.01)) +
  labs(x = NULL, y = "Pearson dispersion") +
  theme_publication(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(output_figures, "family_diagnostics_pearson_dispersion_all_candidates.png"), p_disp, width = 9, height = 5, dpi = 600)

plot_comp <- bind_rows(
  lmdi_updated %>% mutate(method = "Aggregate LMDI-style"),
  shapley_updated %>% mutate(method = "Model-based Shapley")
) %>%
  mutate(
    component_label = factor(component_label, levels = component_labels[c("KAPPA", "LAMBDA", "N", "PI")]),
    method = factor(method, levels = c("Aggregate LMDI-style", "Model-based Shapley"))
  )

p_comp <- ggplot(plot_comp, aes(x = component_label, y = contribution_bln, fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  coord_flip() +
  scale_y_continuous(labels = label_pub_number(accuracy = 0.01)) +
  labs(x = NULL, y = "Contribution, billion points") +
  theme_publication(base_size = 10)

ggsave(file.path(output_figures, "empirical_lmdi_vs_shapley_updated.png"), p_comp, width = 8, height = 5, dpi = 600)

# ------------------------------------------------------------
# 10. Konsola
# ------------------------------------------------------------

message("============================================================")
message("GOTOWE — R13 family selection + refreshed results")
message("============================================================")
message("Tabele:")
message("  ", file.path(output_tables, "01_family_diagnostics_all_candidates.csv"))
message("  ", file.path(output_tables, "02_family_selection_recommendations.csv"))
message("  ", file.path(output_tables, "03_final_model_quality.csv"))
message("  ", file.path(output_tables, "04_empirical_shapley_updated.csv"))
message("  ", file.path(output_tables, "05_empirical_lmdi_updated.csv"))
message("  ", file.path(output_tables, "06_empirical_lmdi_vs_shapley_updated.csv"))
message("  ", file.path(output_tables, "07_supplement_s1b_values_updated.csv"))
message("  ", file.path(output_tables, "08_text_numbers_updated.csv"))
message("")
message("Figury:")
message("  ", file.path(output_figures, "family_diagnostics_aic_all_candidates.png"))
message("  ", file.path(output_figures, "family_diagnostics_pearson_dispersion_all_candidates.png"))
message("  ", file.path(output_figures, "empirical_lmdi_vs_shapley_updated.png"))
message("============================================================")

print(selection_recommendations)
print(final_model_quality)
print(supp_s1b_updated)
print(comparison_updated)
