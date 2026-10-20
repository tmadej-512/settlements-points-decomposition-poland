# ============================================================
# R/10d_simulation_recovery_plots.R
# ------------------------------------------------------------
# Scenario-aware simulation plots.
# ============================================================

compute_simulation_plot_limits <- function(sim_raw) {
  list(
    truth = range(c(sim_raw$truth, sim_raw$estimate), na.rm = TRUE),
    error = range(sim_raw$estimate - sim_raw$truth, na.rm = TRUE)
  )
}

plot_simulation_recovery_vs_truth_by_scenario <- function(sim_raw,
                                                          scenario_id,
                                                          expected_methods = c("naive", "ofat", "sequential", "shapley"),
                                                          y_limits = NULL) {
  expected_methods <- intersect(expected_methods, names(METHOD_LABELS))
  plot_df <- sim_raw %>%
    dplyr::filter(.data$scenario == scenario_id) %>%
    dplyr::mutate(
      method = factor(.data$method, levels = expected_methods),
      component = factor(.data$component, levels = c("N", "PI", "LAMBDA", "KAPPA"),
                              labels = COMPONENT_STRIP_LABELS[c("N", "PI", "LAMBDA", "KAPPA")])
    )

  if (nrow(plot_df) == 0) stop("No simulation rows for scenario: ", scenario_id, call. = FALSE)

  truth_df <- plot_df %>%
    dplyr::group_by(.data$component) %>%
    dplyr::summarise(truth = mean(.data$truth, na.rm = TRUE), .groups = "drop")

  title <- dplyr::case_when(
    scenario_id == "oracle_decomposition" ~ "Oracle decomposition vs known ground truth",
    scenario_id == "clean_recovery" ~ "Simulation recovery vs known ground truth: clean model-based recovery",
    scenario_id == "heterogeneous_stress" ~ "Simulation recovery vs known ground truth: heterogeneous stress test",
    scenario_id == "cell_composition_shift" ~ "Simulation recovery vs known ground truth: cell-composition shift",
    scenario_id == "model_misspecification" ~ "Simulation recovery vs known ground truth: model-misspecification stress",
    TRUE ~ paste0("Simulation recovery vs known ground truth: ", scenario_id)
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$method, y = .data$estimate, fill = .data$method)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.25) +
    ggplot2::geom_hline(
      data = truth_df,
      ggplot2::aes(yintercept = .data$truth),
      linewidth = 0.45,
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~ component, scales = "fixed") +
    ggplot2::scale_fill_manual(values = METHOD_COLORS[expected_methods], labels = METHOD_LABELS[expected_methods], guide = "none") +
    ggplot2::scale_x_discrete(labels = METHOD_AXIS_LABELS[expected_methods], drop = FALSE) +
    ggplot2::scale_y_continuous(labels = label_number_pl(), limits = y_limits) +
    ggplot2::labs(
      title = title,
      subtitle = "Boxplots = estimated component contributions; dashed line = known ground truth",
      x = NULL,
      y = "Estimated component contribution"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1, size = 7.8))
}

plot_simulation_error_by_scenario <- function(sim_raw,
                                              scenario_id,
                                              expected_methods = c("naive", "ofat", "sequential", "shapley"),
                                              y_limits = NULL) {
  expected_methods <- intersect(expected_methods, names(METHOD_LABELS))
  plot_df <- sim_raw %>%
    dplyr::filter(.data$scenario == scenario_id) %>%
    dplyr::mutate(
      method = factor(.data$method, levels = expected_methods),
      component = factor(.data$component, levels = c("N", "PI", "LAMBDA", "KAPPA"),
                              labels = COMPONENT_STRIP_LABELS[c("N", "PI", "LAMBDA", "KAPPA")]),
      error = .data$estimate - .data$truth
    )

  if (nrow(plot_df) == 0) stop("No simulation rows for scenario: ", scenario_id, call. = FALSE)

  title <- dplyr::case_when(
    scenario_id == "oracle_decomposition" ~ "Oracle decomposition error relative to known ground truth",
    scenario_id == "clean_recovery" ~ "Simulation error relative to known ground truth: clean model-based recovery",
    scenario_id == "heterogeneous_stress" ~ "Simulation error relative to known ground truth: heterogeneous stress test",
    scenario_id == "cell_composition_shift" ~ "Simulation error relative to known ground truth: cell-composition shift",
    scenario_id == "model_misspecification" ~ "Simulation error relative to known ground truth: model-misspecification stress",
    TRUE ~ paste0("Simulation error relative to known ground truth: ", scenario_id)
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$method, y = .data$error, fill = .data$method)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.45, linetype = "dashed") +
    ggplot2::geom_boxplot(outlier.alpha = 0.25) +
    ggplot2::facet_wrap(~ component, scales = "fixed") +
    ggplot2::scale_fill_manual(values = METHOD_COLORS[expected_methods], labels = METHOD_LABELS[expected_methods], guide = "none") +
    ggplot2::scale_x_discrete(labels = METHOD_AXIS_LABELS[expected_methods], drop = FALSE) +
    ggplot2::scale_y_continuous(labels = label_number_pl(), limits = y_limits) +
    ggplot2::labs(
      title = title,
      subtitle = "Error = estimate − truth; dashed line = perfect recovery",
      x = NULL,
      y = "Estimation error"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1, size = 7.8))
}
