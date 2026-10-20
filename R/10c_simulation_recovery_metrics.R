# ============================================================
# R/10c_simulation_recovery_metrics.R
# ============================================================

compute_simulation_recovery_metrics <- function(sim_raw) {
  required <- c("scenario", "rep_id", "method", "component", "estimate", "truth")
  missing <- setdiff(required, names(sim_raw))
  if (length(missing) > 0) {
    stop("compute_simulation_recovery_metrics(): missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  sim_raw %>%
    dplyr::mutate(
      error = .data$estimate - .data$truth,
      abs_error = abs(.data$error),
      sq_error = .data$error^2,
      rel_error_pct = dplyr::if_else(abs(.data$truth) > 1e-9, 100 * .data$error / abs(.data$truth), NA_real_),
      abs_rel_error_pct = abs(.data$rel_error_pct)
    ) %>%
    dplyr::group_by(.data$scenario, .data$method, .data$component) %>%
    dplyr::summarise(
      n_reps_ok = sum(!is.na(.data$error)),
      truth = mean(.data$truth, na.rm = TRUE),
      estimate_mean = mean(.data$estimate, na.rm = TRUE),
      bias = mean(.data$error, na.rm = TRUE),
      median_error = median(.data$error, na.rm = TRUE),
      mae = mean(.data$abs_error, na.rm = TRUE),
      rmse = sqrt(mean(.data$sq_error, na.rm = TRUE)),
      sd_error = stats::sd(.data$error, na.rm = TRUE),
      mean_abs_relative_error_pct = mean(.data$abs_rel_error_pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      bias_pct_of_truth = dplyr::if_else(abs(.data$truth) > 1e-9, .data$bias / abs(.data$truth), NA_real_),
      rmse_pct_of_truth = dplyr::if_else(abs(.data$truth) > 1e-9, .data$rmse / abs(.data$truth), NA_real_)
    ) %>%
    dplyr::group_by(.data$scenario, .data$component) %>%
    dplyr::mutate(rmse_rank_within_component = dplyr::min_rank(.data$rmse)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$scenario, .data$component, .data$method)
}

compute_simulation_recovery_metrics_overall <- function(metrics_by_component) {
  metrics_by_component %>%
    dplyr::group_by(.data$scenario, .data$method) %>%
    dplyr::summarise(
      mean_abs_bias_pct_of_truth = mean(abs(.data$bias_pct_of_truth), na.rm = TRUE),
      mean_rmse_pct_of_truth = mean(.data$rmse_pct_of_truth, na.rm = TRUE),
      mean_mae = mean(.data$mae, na.rm = TRUE),
      mean_rmse = mean(.data$rmse, na.rm = TRUE),
      max_abs_bias_pct_of_truth = max(abs(.data$bias_pct_of_truth), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$scenario, .data$mean_rmse_pct_of_truth)
}
