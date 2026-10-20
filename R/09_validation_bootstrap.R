# ============================================================
# R/09_validation_bootstrap.R
# ------------------------------------------------------------
# Bootstrap uncertainty for the decomposition effects.
#
# Approach: cluster (block) bootstrap over administrative units
# (`code`), refitting the full model chain and recomputing all four
# decomposition methods on each resample. Resampling whole
# administrative units (rather than individual rows) respects the
# within-unit dependence induced by the random effect s(code, bs='re')
# and by repeated observations of the same unit across sex/age_group
# cells.
#
# This produces percentile confidence intervals for each component's
# effect under each method, used in the paper to assess whether the
# differences between e.g. sequential and Shapley estimates are
# distinguishable from resampling noise, and to report uncertainty
# around the headline N/pi/lambda/kappa contributions.
#
# COST WARNING: each bootstrap replicate refits THREE GAM models.
# This is the most expensive step in the pipeline. Start with a small
# `reps` (e.g. 20) to check timing before committing to
# config$validation$bootstrap_reps (default 200).
# ============================================================

suppressPackageStartupMessages(library(dplyr))

#' One bootstrap replicate: resample administrative units with
#' replacement, refit the model chain, and run all requested
#' decomposition methods for one (baseline, comparison) year pair.
#'
#' @param df Full aggregated single-disease data (all years).
#' @param model_cfg The `model:` config block.
#' @param baseline_year,comparison_year Years to compare.
#' @param methods Character vector of methods to run, subset of
#'   c("naive","ofat","sequential","shapley").
#' @return A tibble with one row per requested method, or NULL if the
#'   refit fails on this resample (logged, not fatal -- some resamples
#'   may not have enough distinct units/years for identifiability).
bootstrap_replicate <- function(df, model_cfg, baseline_year, comparison_year,
                                 methods = c("naive", "ofat", "sequential", "shapley")) {
  units <- unique(df$code)
  resampled_units <- sample(units, length(units), replace = TRUE)

  # Rebuild a resampled data frame; duplicated units get a unique
  # synthetic `code` suffix so the random effect treats each draw as
  # a distinct level (standard cluster-bootstrap practice for mixed/
  # random-effect models).
  df_boot <- purrr::map_dfr(seq_along(resampled_units), function(i) {
    df %>%
      filter(code == resampled_units[i]) %>%
      mutate(code = paste0(code, "_b", i))
  }) %>%
    mutate(code = factor(code))

  fit_boot <- tryCatch(
    fit_model_chain(df_boot, model_cfg),
    error = function(e) {
      message("  [bootstrap] refit failed on this resample: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit_boot)) return(NULL)

  df_pair <- df_boot %>% filter(year %in% c(baseline_year, comparison_year))
  if (length(unique(df_pair$year)) < 2) return(NULL)

  out <- list()
  if ("naive" %in% methods) {
    out$naive <- decompose_naive_2y(df_pair, NULL, baseline_year, comparison_year) %>%
      mutate(method = "naive")
  }
  if ("ofat" %in% methods) {
    out$ofat <- tryCatch(
      decompose_ofat_2y(df_pair, fit_boot, baseline_year, comparison_year) %>% mutate(method = "ofat"),
      error = function(e) NULL
    )
  }
  if ("sequential" %in% methods) {
    out$sequential <- tryCatch(
      decompose_sequential_2y(df_pair, fit_boot, baseline_year, comparison_year) %>% mutate(method = "sequential"),
      error = function(e) NULL
    )
  }
  if ("shapley" %in% methods) {
    out$shapley <- tryCatch(
      decompose_shapley_2y(df_pair, fit_boot, baseline_year, comparison_year) %>% mutate(method = "shapley"),
      error = function(e) NULL
    )
  }

  bind_rows(out)
}

#' Run the full bootstrap and summarise as percentile confidence intervals
#'
#' @param df Full aggregated single-disease data.
#' @param model_cfg The `model:` config block.
#' @param baseline_year,comparison_year Years to compare.
#' @param reps Number of bootstrap replicates.
#' @param seed RNG seed for reproducibility.
#' @param methods Methods to bootstrap.
#' @param ci_level Confidence level for the percentile interval (default 0.95).
#' @return list(raw = <all replicate rows, with `rep_id`>,
#'              summary = <one row per method x component, with point
#'                         estimate from the FULL (non-resampled) fit,
#'                         bootstrap mean/SD, and percentile CI>)
run_bootstrap <- function(df, model_cfg, baseline_year, comparison_year,
                           reps = 10, seed = 1,
                           methods = c("naive", "ofat", "sequential", "shapley"),
                           ci_level = 0.95) {
  set.seed(seed)

  message("Bootstrap: ", reps, " replicates over administrative units (this refits ",
          length(methods[methods != "naive"]), " model-based method(s) x ", reps, " resamples)...")

  pb <- txtProgressBar(min = 0, max = reps, style = 3)
  raw <- vector("list", reps)
  for (b in seq_len(reps)) {
    raw[[b]] <- bootstrap_replicate(df, model_cfg, baseline_year, comparison_year, methods)
    if (!is.null(raw[[b]])) raw[[b]]$rep_id <- b
    setTxtProgressBar(pb, b)
  }
  close(pb)

  raw_df <- bind_rows(raw)
  n_failed <- reps - length(unique(raw_df$rep_id))
  if (n_failed > 0) {
    message("Note: ", n_failed, " / ", reps, " bootstrap resamples failed to refit and were dropped.")
  }

  alpha <- 1 - ci_level
  comp_cols <- c("N_effect", "PI_effect", "LAMBDA_effect", "KAPPA_effect")

  summary_tbl <- raw_df %>%
    tidyr::pivot_longer(all_of(comp_cols), names_to = "component", values_to = "effect") %>%
    mutate(component = sub("_effect$", "", component)) %>%
    group_by(method, component) %>%
    summarise(
      boot_mean   = mean(effect, na.rm = TRUE),
      boot_sd     = sd(effect, na.rm = TRUE),
      ci_lower    = quantile(effect, alpha / 2, na.rm = TRUE),
      ci_upper    = quantile(effect, 1 - alpha / 2, na.rm = TRUE),
      n_reps_used = sum(!is.na(effect)),
      .groups = "drop"
    )

  list(raw = raw_df, summary = summary_tbl)
}
