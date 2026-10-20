# ============================================================
# R/12_plots.R
# ------------------------------------------------------------
# Plots comparing the four decomposition methods (aggregate LMDI-style, OFAT,
# sequential, Shapley) for a single disease group, plus diagnostic
# plots used in the validation sections (order sensitivity,
# simulation recovery, bootstrap intervals).
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

METHOD_COLORS <- c(
  "naive"      = "#7A7A7A",  # neutral grey for aggregate LMDI-style benchmark
  "ofat"       = "#00468B",  # Lancet-style deep blue
  "sequential" = "#0099B4",  # blue-teal
  "shapley"    = "#D55E00"   # burnt orange for the preferred method
)

METHOD_LABELS <- c(
  naive      = "Aggregate LMDI-style",
  ofat       = "OFAT",
  sequential = "Sequential (1 order)",
  shapley    = "Shapley (avg. of 24 orders)"
)

METHOD_AXIS_LABELS <- c(
  naive      = "LMDI-style",
  ofat       = "OFAT",
  sequential = "Sequential",
  shapley    = "Shapley"
)

COMPONENT_STRIP_LABELS <- c(
  N = "N",
  PI = "\u03c0",
  LAMBDA = "\u03bb",
  KAPPA = "\u03ba"
)

COMP_LABELS <- c(
  N_effect      = "N: demography",
  PI_effect     = "\u03c0: patients/pop",
  LAMBDA_effect = "\u03bb: services/patient",
  KAPPA_effect  = "\u03ba: value/service"
)

#' Stack the four methods' results into one tidy data frame for plotting
#'
#' @param naive,ofat,sequential,shapley One-row (or multi-row, for
#'   trend-type results) tibbles from each decompose_*_2y() call,
#'   each already tagged with its own `method` column OR untagged
#'   (this function adds the tag if missing).
stack_methods <- function(naive = NULL, ofat = NULL, sequential = NULL, shapley = NULL) {
  add_tag <- function(df, tag) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    if (!"method" %in% names(df)) df$method <- tag
    df
  }
  bind_rows(
    add_tag(naive, "naive"),
    add_tag(ofat, "ofat"),
    add_tag(sequential, "sequential"),
    add_tag(shapley, "shapley")
  ) %>%
    mutate(method = factor(method, levels = names(METHOD_COLORS)))
}

#' Grouped bar chart: each component's effect, one bar cluster per method
#'
#' The diagnostic plot for "how much do the four methods disagree, for
#' this disease and this year pair".
plot_method_comparison_bars <- function(stacked) {
  long <- stacked %>%
    tidyr::pivot_longer(
      cols = c(N_effect, PI_effect, LAMBDA_effect, KAPPA_effect),
      names_to = "component", values_to = "effect"
    ) %>%
    mutate(component = factor(COMP_LABELS[component], levels = COMP_LABELS))

  ggplot(long, aes(x = component, y = effect, fill = method)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = METHOD_COLORS, labels = METHOD_LABELS, name = NULL) +
    scale_y_continuous(labels = label_number_pl()) +
    labs(
      title = "Component effects by decomposition method",
      subtitle = paste0(unique(stacked$baseline_year), " vs ", unique(stacked$comparison_year)),
      x = NULL, y = "Effect on total settlement points"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

#' Residual comparison: shows that sequential and Shapley are
#' numerically additive (residual ~ 0) while aggregate LMDI-style and OFAT may not be
#' (LMDI-style is additive by its log-ratio construction; OFAT is not).
plot_residual_comparison <- function(stacked) {
  ggplot(stacked, aes(x = method, y = residual, fill = method)) +
    geom_col(width = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = METHOD_COLORS, labels = METHOD_LABELS, guide = "none") +
    scale_x_discrete(labels = METHOD_AXIS_LABELS) +
    scale_y_continuous(labels = label_number_pl()) +
    labs(
      title = "Additivity check: residual by method",
      subtitle = "Sequential and Shapley are additive by construction (residual \u2248 0); OFAT need not be",
      x = NULL, y = "Residual (total_change \u2212 sum of effects)"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7.8))
}

#' Order-sensitivity plot: range of the sequential effect across all 24
#' orderings, with the Shapley value overlaid as a point.
plot_order_sensitivity <- function(range_summary) {
  range_summary <- range_summary %>%
    mutate(component = factor(COMP_LABELS[paste0(component, "_effect")],
                               levels = COMP_LABELS))

  ggplot(range_summary, aes(x = component)) +
    geom_linerange(aes(ymin = seq_min, ymax = seq_max), linewidth = 4, color = "#4E79A7", alpha = 0.4) +
    geom_point(aes(y = shapley_value), color = "#F28E2B", size = 3) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_y_continuous(labels = label_number_pl()) +
    labs(
      title = "Order sensitivity of the sequential decomposition",
      subtitle = "Bar = range across all 24 substitution orders; dot = Shapley average",
      x = NULL, y = "Effect on total settlement points"
    ) +
    theme_minimal(base_size = 12)
}

#' Simulation recovery plot: estimate vs truth, one panel per component,
#' one boxplot per method, across all simulation replicates.
plot_simulation_recovery <- function(sim_raw,
                                     expected_methods = c("naive", "ofat", "sequential", "shapley")) {
  expected_methods <- intersect(expected_methods, names(METHOD_LABELS))
  observed_methods <- sort(unique(as.character(sim_raw$method)))
  missing_methods <- setdiff(expected_methods, observed_methods)
  if (length(missing_methods) > 0) {
    stop(
      "simulation_raw is missing method(s): ", paste(missing_methods, collapse = ", "),
      ". The plot would show only a subset of methods. Re-run scripts/run_05_simulation_study.R ",
      "after fixing the model/decomposition error reported there.",
      call. = FALSE
    )
  }

  sim_raw <- sim_raw %>%
    mutate(
      method = factor(method, levels = expected_methods),
      component = factor(component, levels = c("N", "PI", "LAMBDA", "KAPPA"),
                         labels = COMPONENT_STRIP_LABELS[c("N", "PI", "LAMBDA", "KAPPA")])
    )

  ggplot(sim_raw, aes(x = method, y = estimate - truth, fill = method)) +
    geom_boxplot(outlier.alpha = 0.3) +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
    facet_wrap(~component, scales = "fixed") +
    scale_fill_manual(values = METHOD_COLORS[expected_methods], labels = METHOD_LABELS[expected_methods], guide = "none") +
    scale_x_discrete(labels = METHOD_AXIS_LABELS[expected_methods], drop = FALSE) +
    scale_y_continuous(labels = label_number_pl()) +
    labs(
      title = "Simulation study: estimation error vs. known ground truth",
      subtitle = "estimate \u2212 truth, by component and method (0 = perfect recovery)",
      x = NULL, y = "Error"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7.8))
}

#' Bootstrap interval plot: point estimate plus percentile CI, by
#' method and component.
plot_bootstrap_intervals <- function(boot_summary, point_estimates) {
  point_long <- point_estimates %>%
    tidyr::pivot_longer(
      cols = c(N_effect, PI_effect, LAMBDA_effect, KAPPA_effect),
      names_to = "component", values_to = "point_estimate"
    ) %>%
    mutate(component = sub("_effect$", "", component))

  merged <- boot_summary %>%
    left_join(point_long, by = c("method", "component")) %>%
    mutate(
      component = factor(COMP_LABELS[paste0(component, "_effect")], levels = COMP_LABELS),
      method = factor(method, levels = names(METHOD_COLORS))
    )

  ggplot(merged, aes(x = component, y = point_estimate, color = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 2.5) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  position = position_dodge(width = 0.5), width = 0.2) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_color_manual(values = METHOD_COLORS, labels = METHOD_LABELS, name = NULL) +
    scale_y_continuous(labels = label_number_pl()) +
    labs(
      title = "Bootstrap confidence intervals by component and method",
      x = NULL, y = "Effect on total settlement points"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

#' Simulation recovery plot: estimated contributions with known ground truth line
#'
#' Easier to read than the error plot: each facet shows estimated component
#' contributions across simulation replicates, and the dashed line is the known
#' true contribution for that component.
plot_simulation_recovery_vs_truth <- function(sim_raw,
                                              expected_methods = c("naive", "ofat", "sequential", "shapley")) {
  expected_methods <- intersect(expected_methods, names(METHOD_LABELS))
  observed_methods <- sort(unique(as.character(sim_raw$method)))
  missing_methods <- setdiff(expected_methods, observed_methods)
  if (length(missing_methods) > 0) {
    stop(
      "simulation_raw is missing method(s): ", paste(missing_methods, collapse = ", "),
      ". Re-run scripts/run_05_simulation_study.R after fixing the simulation step.",
      call. = FALSE
    )
  }

  plot_df <- sim_raw %>%
    mutate(
      method = factor(method, levels = expected_methods),
      component = factor(component, levels = c("N", "PI", "LAMBDA", "KAPPA"),
                         labels = COMPONENT_STRIP_LABELS[c("N", "PI", "LAMBDA", "KAPPA")])
    )

  truth_df <- plot_df %>%
    group_by(component) %>%
    summarise(truth = mean(truth, na.rm = TRUE), .groups = "drop")

  ggplot(plot_df, aes(x = method, y = estimate, fill = method)) +
    geom_boxplot(outlier.alpha = 0.25) +
    geom_hline(data = truth_df, aes(yintercept = truth),
               linewidth = 0.45, linetype = "dashed", inherit.aes = FALSE) +
    facet_wrap(~component, scales = "fixed") +
    scale_fill_manual(values = METHOD_COLORS[expected_methods], labels = METHOD_LABELS[expected_methods], guide = "none") +
    scale_x_discrete(labels = METHOD_AXIS_LABELS[expected_methods], drop = FALSE) +
    scale_y_continuous(labels = label_number_pl()) +
    labs(
      title = "Simulation recovery against known ground truth",
      subtitle = "Boxplots = estimated contribution; dashed line = mean known true contribution across replicates",
      x = NULL, y = "Estimated component contribution"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7.8))
}

#' Diagnostic plots for the synthetic DGP used in simulation study
#'
#' Saves separate trajectories for N, PI, LAMBDA, and KAPPA, plus a combined
#' 2x2 figure when {patchwork} is available. Uses deterministic DGP means,
#' not observed noisy draws, so slopes/curvature are visible.
plot_simulation_dgp_diagnostics <- function(seed = 123,
                                            out_dir = "output/figures/simulation_dgp",
                                            years = 2018:2023,
                                            n_units = 30,
                                            sex_levels = c("F", "M"),
                                            age_levels = c("0-17", "18-44", "45-64", "65+"),
                                            true_growth = list(N = 0.01, PI = 0.03, LAMBDA = -0.01, KAPPA = 0.02),
                                            heterogeneous_slopes = TRUE,
                                            regional_random_slopes = TRUE,
                                            regional_slope_sd = list(N = 0.003, PI = 0.008, LAMBDA = 0.008, KAPPA = 0.006),
                                            width = 11,
                                            height = 7,
                                            dpi = 600) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  sim <- simulate_disease_panel(
    n_units = n_units,
    years = years,
    sex_levels = sex_levels,
    age_levels = age_levels,
    true_growth = true_growth,
    heterogeneous_slopes = heterogeneous_slopes,
    regional_random_slopes = regional_random_slopes,
    regional_slope_sd = regional_slope_sd,
    seed = seed
  )

  d <- sim$data %>%
    mutate(
      pi_obs = if_else(pop > 0, n_patients / pop, NA_real_),
      lambda_obs = if_else(n_patients > 0, n_provision / n_patients, NA_real_),
      kappa_obs = if_else(n_provision > 0, value / n_provision, NA_real_)
    )

  dgp_long <- d %>%
    group_by(year, sex, age_group) %>%
    summarise(
      N = mean(pop_mean, na.rm = TRUE),
      PI = weighted.mean(pi_mean, w = pop_mean, na.rm = TRUE),
      LAMBDA = weighted.mean(lambda_mean, w = patients_mean, na.rm = TRUE),
      KAPPA = weighted.mean(kappa_mean, w = provision_mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(N, PI, LAMBDA, KAPPA), names_to = "component", values_to = "value") %>%
    mutate(
      component = factor(component, levels = c("N", "PI", "LAMBDA", "KAPPA")),
      sex_age = interaction(sex, age_group, sep = " / ", drop = TRUE)
    )

  make_panel <- function(comp, ylab, title, percent = FALSE) {
    p <- dgp_long %>% filter(component == comp) %>%
      ggplot(aes(x = year, y = value, group = sex_age, linetype = sex)) +
      geom_line(linewidth = 0.85) +
      geom_point(size = 1.4) +
      facet_wrap(~ age_group, scales = "fixed") +
      labs(
        title = title,
        subtitle = if (isTRUE(heterogeneous_slopes) && isTRUE(regional_random_slopes))
          "Synthetic DGP with age-sex-specific slopes, regional random slopes, and nonlinear trajectories" else if (isTRUE(heterogeneous_slopes))
          "Synthetic DGP with age-sex-specific slopes and nonlinear trajectories" else
          "Synthetic DGP with common component slopes",
        x = "Year", y = ylab, linetype = "Sex"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank())
    if (percent) {
      p <- p + scale_y_continuous(labels = scales::label_percent(accuracy = 0.1))
    } else {
      p <- p + scale_y_continuous(labels = label_number_pl())
    }
    p
  }

  p_N <- make_panel("N", "Mean population", "Synthetic DGP: population trajectories")
  p_PI <- make_panel("PI", "Patients / population", "Synthetic DGP: patient-share trajectories", percent = TRUE)
  p_LA <- make_panel("LAMBDA", "Services / patient", "Synthetic DGP: service-intensity trajectories")
  p_KA <- make_panel("KAPPA", "Value / service", "Synthetic DGP: unit-value trajectories")

  safe_ggsave(p_N, file.path(out_dir, "DGP_population_trajectories.png"), width, height, dpi)
  safe_ggsave(p_PI, file.path(out_dir, "DGP_pi_patient_share.png"), width, height, dpi)
  safe_ggsave(p_LA, file.path(out_dir, "DGP_lambda_services_per_patient.png"), width, height, dpi)
  safe_ggsave(p_KA, file.path(out_dir, "DGP_kappa_value_per_service.png"), width, height, dpi)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- (p_N + p_PI) / (p_LA + p_KA)
    ggplot2::ggsave(
      file.path(out_dir, "DGP_combined_trajectories.png"), combined,
      width = width + 5, height = (height + 1) * 2, dpi = dpi
    )
  }

  unit_slopes <- d %>%
    dplyr::distinct(code, N_unit_slope, PI_unit_slope, LAMBDA_unit_slope, KAPPA_unit_slope) %>%
    tidyr::pivot_longer(
      cols = c(N_unit_slope, PI_unit_slope, LAMBDA_unit_slope, KAPPA_unit_slope),
      names_to = "component", values_to = "unit_slope"
    ) %>%
    dplyr::mutate(
      component = dplyr::recode(component,
        N_unit_slope = "N", PI_unit_slope = "PI",
        LAMBDA_unit_slope = "LAMBDA", KAPPA_unit_slope = "KAPPA"
      ),
      component = factor(component, levels = c("N", "PI", "LAMBDA", "KAPPA"))
    )

  p_slopes <- unit_slopes %>%
    ggplot(aes(x = component, y = unit_slope)) +
    geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed") +
    geom_boxplot(outlier.alpha = 0.45) +
    labs(
      title = "Synthetic DGP: regional random slopes",
      subtitle = "Unit-level deviations from the age-sex component trend",
      x = "Component", y = "Unit-level annual log-slope deviation"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  safe_ggsave(p_slopes, file.path(out_dir, "DGP_regional_random_slopes.png"), width, height, dpi)
  readr::write_csv(unit_slopes, file.path(out_dir, "DGP_regional_random_slopes.csv"))
  readr::write_csv(dgp_long, file.path(out_dir, "DGP_component_trajectories.csv"))

  invisible(list(data = d, trajectories = dgp_long, unit_slopes = unit_slopes,
                 plots = list(N = p_N, PI = p_PI, LAMBDA = p_LA, KAPPA = p_KA,
                              regional_slopes = p_slopes)))
}

#' All 24 sequential orders: component-specific attribution dispersion
#'
#' This plot makes the practical order-dependence of the sequential
#' decomposition visible. Each point is one of the 24 possible substitution
#' orders. The red/orange Shapley reference is added by the separate
#' range plot; this plot focuses on how wide the sequential cloud is.
plot_sequential_all_orders <- function(all_orders) {
  long <- all_orders %>%
    dplyr::mutate(order_id = dplyr::row_number()) %>%
    tidyr::pivot_longer(
      cols = c(N_effect, PI_effect, LAMBDA_effect, KAPPA_effect),
      names_to = "component", values_to = "effect"
    ) %>%
    dplyr::mutate(
      component = sub("_effect$", "", component),
      component = factor(COMPONENT_STRIP_LABELS[component],
                         levels = COMPONENT_STRIP_LABELS[c("N", "PI", "LAMBDA", "KAPPA")]),
      sequential_order = factor(sequential_order, levels = unique(sequential_order))
    )

  ggplot(long, aes(x = sequential_order, y = effect)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_point(size = 1.8, alpha = 0.75) +
    facet_wrap(~component, scales = "fixed") +
    scale_y_continuous(labels = label_number_pl()) +
    labs(
      title = "Sequential decomposition across all 24 substitution orders",
      subtitle = "Each point is one possible component ordering; dispersion shows order dependence",
      x = "Sequential substitution order",
      y = "Attributed component effect"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      panel.grid.minor = element_blank()
    )
}
