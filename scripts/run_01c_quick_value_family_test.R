#!/usr/bin/env Rscript
# ============================================================
# scripts/run_01c_quick_value_family_test.R
# ------------------------------------------------------------
# Quick distribution-family test for the value / kappa model.
#
# Compares:
#   1. Tweedie log-link with offset log(n_provision)
#   2. Gamma log-link with offset log(n_provision), positive values only
#   3. Gaussian on log(value / n_provision), positive values only
#
# Run after:
#   source("scripts/run_01_fit_model.R")
#
# Output:
#   output/tables/quick_value_family_test.csv
#   output/figures/model_diagnostics/value_family_quick_test.png
# ============================================================

source("R/99_load_all.R")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(mgcv)
})

cfg <- load_config("config.R")
ensure_dir(cfg$output$tables_dir)
ensure_dir(cfg$output$figures_dir)
ensure_dir(file.path(cfg$output$figures_dir, "model_diagnostics"))

dat_path <- "output/rds/aggregated_data.rds"

if (!file.exists(dat_path)) {
  stop(
    "Missing output/rds/aggregated_data.rds. Run source('scripts/run_01_fit_model.R') first.",
    call. = FALSE
  )
}

dat <- readRDS(dat_path)

value_col <- if ("value" %in% names(dat)) "value" else if ("sum_points" %in% names(dat)) "sum_points" else NA_character_
if (is.na(value_col)) stop("Could not find value column: expected 'value' or 'sum_points'.", call. = FALSE)

# Standardize to `value` for formulas.
df <- dat %>%
  mutate(
    value     = .data[[value_col]],
    year_num  = as.numeric(year),
    sex_age   = droplevels(interaction(sex, age_group, drop = TRUE)),
    code      = droplevels(code),
    sex       = droplevels(sex),
    age_group = droplevels(age_group)
  )

eps_offset <- cfg$model$eps_offset %||% 1
k_year <- min(cfg$model$k_year %||% 5, length(unique(df$year)))
if (k_year < 3) k_year <- length(unique(df$year))

# ------------------------------------------------------------
# Basic zero-mass diagnostics
# ------------------------------------------------------------

zero_diag <- tibble::tibble(
  n_rows = nrow(df),
  n_positive_provision = sum(df$n_provision > 0, na.rm = TRUE),
  n_zero_value_given_positive_provision = sum(df$n_provision > 0 & df$value == 0, na.rm = TRUE),
  prop_zero_value_given_positive_provision =
    mean(df$value == 0 & df$n_provision > 0, na.rm = TRUE),
  n_negative_value = sum(df$value < 0, na.rm = TRUE),
  value_min = suppressWarnings(min(df$value, na.rm = TRUE)),
  value_p50 = suppressWarnings(as.numeric(stats::quantile(df$value, 0.50, na.rm = TRUE))),
  value_p99 = suppressWarnings(as.numeric(stats::quantile(df$value, 0.99, na.rm = TRUE))),
  value_max = suppressWarnings(max(df$value, na.rm = TRUE))
)

message("== Zero / tail diagnostics for value ==")
print(zero_diag)

# ------------------------------------------------------------
# Shared formula pieces
# ------------------------------------------------------------

rhs_common <- paste(
  "sex * age_group",
  paste0("s(year_num, by = sex_age, k = ", k_year, ")"),
  "s(code, bs = 're')",
  sep = " + "
)

f_tweedie <- as.formula(
  paste(
    "value ~",
    rhs_common,
    "+ offset(log(pmax(n_provision,", eps_offset, ")))"
  )
)

f_gamma <- f_tweedie

df_pos <- df %>%
  filter(n_provision > 0, value > 0) %>%
  mutate(
    kappa_obs = value / n_provision,
    log_kappa_obs = log(kappa_obs)
  )

f_lognormal <- as.formula(
  paste(
    "log_kappa_obs ~",
    rhs_common
  )
)

ctrl <- gam.control(
  maxit = cfg$model$maxit_gam %||% 200,
  trace = FALSE
)

fit_one <- function(label, formula, data, family, discrete = TRUE) {
  message("\nFitting: ", label)

  t0 <- Sys.time()

  fit <- tryCatch(
    mgcv::bam(
      formula = formula,
      data = data,
      family = family,
      discrete = discrete,
      method = "fREML",
      control = ctrl
    ),
    error = function(e) e
  )

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (inherits(fit, "error")) {
    warning("Model failed: ", label, " | ", conditionMessage(fit))

    return(list(
      label = label,
      fit = NULL,
      summary = tibble::tibble(
        family_tested = label,
        model_ok = FALSE,
        error = conditionMessage(fit),
        n_obs = nrow(data),
        aic = NA_real_,
        dev_expl = NA_real_,
        adj_r_sq = NA_real_,
        mean_pearson = NA_real_,
        sd_pearson = NA_real_,
        mean_abs_pearson = NA_real_,
        pct_abs_pearson_gt_2 = NA_real_,
        cor_abs_pearson_fitted = NA_real_,
        elapsed_seconds = elapsed
      )
    ))
  }

  r_pearson <- tryCatch(
    residuals(fit, type = "pearson"),
    error = function(e) rep(NA_real_, nrow(data))
  )

  fitted_values <- tryCatch(as.numeric(fitted(fit)), error = function(e) rep(NA_real_, nrow(data)))

  list(
    label = label,
    fit = fit,
    summary = tibble::tibble(
      family_tested = label,
      model_ok = TRUE,
      error = NA_character_,
      n_obs = nrow(data),
      aic = tryCatch(AIC(fit), error = function(e) NA_real_),
      dev_expl = summary(fit)$dev.expl %||% NA_real_,
      adj_r_sq = summary(fit)$r.sq %||% NA_real_,
      mean_pearson = mean(r_pearson, na.rm = TRUE),
      sd_pearson = sd(r_pearson, na.rm = TRUE),
      mean_abs_pearson = mean(abs(r_pearson), na.rm = TRUE),
      pct_abs_pearson_gt_2 = mean(abs(r_pearson) > 2, na.rm = TRUE),
      cor_abs_pearson_fitted = suppressWarnings(stats::cor(abs(r_pearson), fitted_values, use = "complete.obs")),
      elapsed_seconds = elapsed
    )
  )
}

# ------------------------------------------------------------
# Fit candidate families
# ------------------------------------------------------------

fits <- list(
  tweedie = fit_one(
    label = "Tweedie_log_offset_services",
    formula = f_tweedie,
    data = df,
    family = mgcv::tw(link = "log"),
    discrete = TRUE
  )
)

if (nrow(df_pos) > 0) {
  fits$gamma <- fit_one(
    label = "Gamma_log_offset_services_positive_only",
    formula = f_gamma,
    data = df_pos,
    family = Gamma(link = "log"),
    discrete = TRUE
  )

  fits$lognormal <- fit_one(
    label = "Gaussian_log_kappa_positive_only",
    formula = f_lognormal,
    data = df_pos,
    family = gaussian(),
    discrete = TRUE
  )
}

family_summary <- bind_rows(lapply(fits, `[[`, "summary")) %>%
  mutate(
    delta_aic = ifelse(all(is.na(aic)), NA_real_, aic - min(aic, na.rm = TRUE)),
    rank_aic = ifelse(is.na(aic), NA_real_, min_rank(aic))
  ) %>%
  arrange(rank_aic)

out_table <- bind_cols(
  zero_diag[rep(1, nrow(family_summary)), ],
  family_summary
)

readr::write_csv(
  out_table,
  file.path(cfg$output$tables_dir, "quick_value_family_test.csv")
)

message("\n== Quick value-family comparison ==")
print(out_table)

# ------------------------------------------------------------
# Quick residual plot for fitted candidates
# ------------------------------------------------------------

plot_rows <- lapply(fits, function(x) {
  if (is.null(x$fit)) return(NULL)

  fit <- x$fit

  tibble::tibble(
    family_tested = x$label,
    fitted = as.numeric(fitted(fit)),
    pearson_residual = as.numeric(residuals(fit, type = "pearson"))
  )
})

plot_dt <- bind_rows(plot_rows)

if (nrow(plot_dt) > 0) {
  p <- ggplot(
    plot_dt,
    aes(x = fitted, y = pearson_residual)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed") +
    geom_point(alpha = 0.12, size = 0.6, color = "#00468B") +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x, linewidth = 0.7, color = "#D55E00") +
    facet_wrap(~ family_tested, scales = "free_x") +
    scale_x_continuous(labels = label_number_pl()) +
    scale_y_continuous(labels = label_number_pl()) +
    labs(
      title = "Quick value-family diagnostic",
      subtitle = "Pearson residuals versus fitted values",
      x = "Fitted value",
      y = "Pearson residual"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  safe_ggsave(
    p,
    file.path(cfg$output$figures_dir, "model_diagnostics", "value_family_quick_test.png"),
    cfg$output$fig_width + 3,
    cfg$output$fig_height,
    cfg$output$fig_dpi
  )

  print(p)
}

message("\nSaved:")
message("  ", file.path(cfg$output$tables_dir, "quick_value_family_test.csv"))
message("  ", file.path(cfg$output$figures_dir, "model_diagnostics", "value_family_quick_test.png"))
