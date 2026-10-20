# ============================================================
# R/01_data_schema.R
# ------------------------------------------------------------
# Normalises a raw NFZ administrative extract into the fixed
# internal schema used by every downstream step, then
# aggregates it to (icd10 x code x sex x age_group x year) cells.
#
# The decomposition itself operates on a SINGLE disease group
# (e.g. "Cardiovascular diseases"), selected via
# config$data$disease_filter. The disease groups in the source
# extract are already pre-aggregated clinical categories
# ("Jednostka chorobowa"), not raw ICD-10 codes, so no further
# grouping of codes is performed here.
#
# Internal schema after normalize_schema():
#   code        chr   administrative unit identifier (e.g. TERYT)
#   powiat      chr   administrative unit display name
#   x, y        dbl   coordinates (optional; NA if absent)
#   year        int
#   icd10       chr   disease group label ("Jednostka chorobowa")
#   age_group   chr
#   sex         chr   "F" / "M"
#   pop         dbl   population at risk (offset for the patient model)
#   n_patients  dbl   number of patients ("Liczba osob (bez <5)")
#   n_provision dbl   number of services/claims
#   value       dbl   value_var column (sum_points or val_provision)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(data.table)
})

#' Normalise a raw data frame to the internal 4-component schema
#'
#' Column matching is tolerant of several naming conventions seen
#' in NFZ extracts, including the Polish "Data Preview" export
#' format (TERYT, NAZWA, "Jednostka chorobowa", "Wiek pacjenta",
#' "Plec pacjenta", Rok, "Liczba osob (bez <5)").
#'
#' @param x Raw data frame.
#' @param value_var Name of the outcome column to use as `value`
#'   ("sum_points" or "val_provision").
#' @return A tibble in the internal schema, with invalid rows removed.
normalize_schema <- function(x, value_var) {
  nm <- names(x)

  c_code   <- pick_first(nm, c("code", "id", "TERYT", "teryt", "region", "powiat_code"))
  c_powiat <- pick_first(nm, c("powiat", "county", "powiat_name", "NAZWA", "nazwa"))
  c_x      <- pick_first(nm, c("x", "lon", "lng", "longitude"))
  c_y      <- pick_first(nm, c("y", "lat", "latitude"))
  c_year   <- pick_first(nm, c("year", "rok", "Rok", "YEAR"))
  c_icd10  <- pick_first(nm, c("icd10", "disease", "Jednostka chorobowa", "code_icd10"))
  c_age    <- pick_first(nm, c("age_group", "age", "Wiek pacjenta", "grupa_wiek"))
  c_sex    <- pick_first(nm, c("sex", "plec", "Plec pacjenta", "Płeć pacjenta"))
  c_pop    <- pick_first(nm, c("pop", "populacja", "population", "ludnosc", "N"))
  c_pat    <- pick_first(nm, c(
    "n_patients", "cases", "Liczba osob", "Liczba osob (bez <5)",
    "Liczba osób", "Liczba  osób (bez <5)"
  ))
  c_prov   <- pick_first(nm, c("n_provision", "provision_count", "n_services", "liczba_swiadczen"))
  c_val    <- pick_first(nm, c(value_var))

  required <- c(year = c_year, icd10 = c_icd10, age = c_age, sex = c_sex,
                pop = c_pop, n_patients = c_pat, n_provision = c_prov, value = c_val)
  missing <- names(required)[is.na(required)]
  if (length(missing) > 0) {
    stop(
      "Could not find required column(s): ", paste(missing, collapse = ", "), "\n",
      "Available columns: ", paste(nm, collapse = ", "), "\n",
      "Minimum required fields: year, icd10, age_group, sex, pop, ",
      "n_patients, n_provision, and the chosen value_var ('", value_var, "')."
    )
  }

  x %>%
    mutate(
      .code        = if (!is.na(c_code))   .data[[c_code]]   else "ALL",
      .powiat      = if (!is.na(c_powiat)) .data[[c_powiat]] else as.character(.code),
      .x           = if (!is.na(c_x))      .data[[c_x]]      else NA,
      .y           = if (!is.na(c_y))      .data[[c_y]]      else NA
    ) %>%
    transmute(
      code        = as.character(.code),
      powiat      = as.character(.powiat),
      x           = safe_numeric(.x),
      y           = safe_numeric(.y),
      year        = readr::parse_number(as.character(.data[[c_year]])),
      icd10       = as.character(.data[[c_icd10]]),
      age_group   = as.character(.data[[c_age]]),
      sex         = str_to_upper(as.character(.data[[c_sex]])),
      pop         = safe_numeric(.data[[c_pop]]),
      n_patients  = safe_numeric(.data[[c_pat]]),
      n_provision = safe_numeric(.data[[c_prov]]),
      value       = safe_numeric(.data[[c_val]])
    ) %>%
    mutate(
      sex = recode(sex, "K" = "F", "KOBIETA" = "F", "MEZCZYZNA" = "M",
                   "MĘŻCZYZNA" = "M", .default = sex)
    ) %>%
    filter(
      !is.na(year), !is.na(icd10), !is.na(age_group), !is.na(sex),
      !is.na(pop), !is.na(n_patients), !is.na(n_provision), !is.na(value),
      pop >= 0, n_patients >= 0, n_provision >= 0, value >= 0
    )
}

#' Filter to a single disease group and aggregate to analysis cells
#'
#' @param dat Output of normalize_schema().
#' @param disease_filter Value of `icd10` to keep (e.g. "Cardiovascular diseases").
#' @return A tibble aggregated to (code, powiat, x, y, sex, age_group, year),
#'   with `icd10` constant and equal to `disease_filter`.
filter_and_aggregate <- function(dat, disease_filter = NULL) {
  available <- sort(unique(dat$icd10))

  if (is.null(disease_filter) || length(disease_filter) == 0 || is.na(disease_filter) || identical(disease_filter, "")) {
    if (length(available) == 1) {
      disease_filter <- available[[1]]
      message("No disease_filter set; using the only disease group available: '", disease_filter, "'")
    } else {
      stop(
        "config$data$disease_filter is NULL, but the input file contains multiple disease groups.\n",
        "Set config$data$disease_filter to one exact value from the icd10 column.\n",
        "Available values: ", paste(available, collapse = " | ")
      )
    }
  }

  if (!disease_filter %in% available) {
    stop(
      "disease_filter = '", disease_filter, "' not found in the icd10 column.\n",
      "Available values: ", paste(available, collapse = " | ")
    )
  }

  out <- dat %>%
    filter(icd10 == disease_filter) %>%
    group_by(icd10, code, powiat, x, y, sex, age_group, year) %>%
    summarise(
      pop         = sum(pop, na.rm = TRUE),
      n_patients  = sum(n_patients, na.rm = TRUE),
      n_provision = sum(n_provision, na.rm = TRUE),
      value       = sum(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      code      = factor(code),
      powiat    = as.character(powiat),
      sex       = factor(sex),
      age_group = factor(age_group),
      year      = as.integer(year)
    ) %>%
    arrange(code, sex, age_group, year)

  attr(out, "disease_filter") <- disease_filter
  out
}

#' Read a data file, dispatching on extension (.csv, .rds, .xlsx, .xls)
read_input_data <- function(path) {
  path_expanded <- path.expand(path)
  if (!file.exists(path_expanded)) {
    stop("Input data file not found: ", path_expanded)
  }
  ext <- tolower(tools::file_ext(path_expanded))

  if (ext == "csv")              return(data.table::fread(path_expanded, data.table = FALSE))
  if (ext == "rds")               return(readRDS(path_expanded))
  if (ext %in% c("xlsx", "xls"))  return(readxl::read_excel(path_expanded))

  stop("Unsupported data format: .", ext, ". Use .csv, .rds, .xlsx, or .xls.")
}

#' Raw descriptive identity table for the supplement
#'
#' This table makes the accounting factorisation transparent before any
#' model-based counterfactuals are used. It should be reported as a
#' descriptive/supplementary check, not as the primary decomposition.
#'
#' At the aggregate level we define:
#'   value = population * (patients / population) *
#'           (services / patients) * (value / services)
#' using period totals. This is an accounting identity for the observed
#' aggregate; it is not a claim that the payer assigns unit point values
#' as a function of population size.
#'
#' @param df Aggregated analysis data from filter_and_aggregate().
#' @param years Optional vector of years to keep; NULL = all years.
#' @return One row per year with raw totals and aggregate ratios.
make_raw_identity_table <- function(df, years = NULL) {
  if (!is.null(years)) df <- df %>% dplyr::filter(year %in% years)

  df %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      total_population = sum(pop, na.rm = TRUE),
      total_patients   = sum(n_patients, na.rm = TRUE),
      total_services   = sum(n_provision, na.rm = TRUE),
      total_value      = sum(value, na.rm = TRUE),
      n_cells          = dplyr::n(),
      n_admin_units    = dplyr::n_distinct(code),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      pi_patients_per_population = dplyr::if_else(total_population > 0, total_patients / total_population, NA_real_),
      lambda_services_per_patient = dplyr::if_else(total_patients > 0, total_services / total_patients, NA_real_),
      kappa_value_per_service     = dplyr::if_else(total_services > 0, total_value / total_services, NA_real_),
      identity_reconstructed_value = total_population * pi_patients_per_population *
        lambda_services_per_patient * kappa_value_per_service,
      identity_gap = total_value - identity_reconstructed_value,
      identity_gap_relative = dplyr::if_else(abs(total_value) > 0, identity_gap / total_value, NA_real_)
    ) %>%
    dplyr::arrange(year)
}

#' Compare raw observed totals with the totals used by each decomposition
#'
#' For naive log-share decomposition, baseline_total/comparison_total are
#' raw aggregate totals. For model-based OFAT/sequential/Shapley, they are
#' model-implied counterfactual totals V(empty) and V(full). This table is
#' useful in the paper to avoid conflating the observed raw change with the
#' model-implied aggregate change being decomposed.
make_model_vs_observed_table <- function(method_results, raw_identity) {
  obs <- raw_identity %>%
    dplyr::filter(year %in% c(min(year), max(year))) %>%
    dplyr::arrange(year)

  if (nrow(obs) != 2) {
    stop("raw_identity must contain exactly the two comparison years for make_model_vs_observed_table().")
  }

  observed_baseline_total   <- obs$total_value[1]
  observed_comparison_total <- obs$total_value[2]
  observed_total_change     <- observed_comparison_total - observed_baseline_total

  method_results %>%
    dplyr::transmute(
      method,
      baseline_year,
      comparison_year,
      decomposition_total_type = dplyr::if_else(
        as.character(method) == "naive",
        "observed aggregate log-share",
        "model-implied counterfactual"
      ),
      decomposition_baseline_total = baseline_total,
      decomposition_comparison_total = comparison_total,
      decomposition_total_change = total_change,
      observed_baseline_total = observed_baseline_total,
      observed_comparison_total = observed_comparison_total,
      observed_total_change = observed_total_change,
      change_difference_vs_observed = total_change - observed_total_change,
      change_difference_pct_observed = dplyr::if_else(
        abs(observed_total_change) > 0,
        (total_change - observed_total_change) / observed_total_change,
        NA_real_
      )
    )
}
