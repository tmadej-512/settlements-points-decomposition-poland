# ============================================================
# R/15_model_residual_diagnostics.R
# ------------------------------------------------------------
# Residual diagnostics for the three main fitted models:
#   patients, provision, value
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

extract_main_model_residuals <- function(fit, dat) {
  stopifnot(is.list(fit), is.data.frame(dat))

  model_map <- list(
    patients  = fit$m_pat,
    provision = fit$m_prov,
    value     = fit$m_val
  )

  take_or_na <- function(mf, candidates, n) {
    for (cc in candidates) {
      if (cc %in% names(mf)) return(mf[[cc]])
    }
    rep(NA, n)
  }

  out <- lapply(names(model_map), function(nm) {
    m <- model_map[[nm]]
    mf <- stats::model.frame(m)

    fitted_vec <- as.numeric(stats::fitted(m))
    n <- length(fitted_vec)

    # Critical for Gamma value model:
    # the valuation fit may exclude value <= 0 or invalid service-count rows,
    # so residual vectors can be shorter than the original data panel.
    # Therefore all row-level metadata must come from model.frame(m), not from dat.
    year_vec <- take_or_na(mf, c("year", "year_num"), n)
    code_vec <- take_or_na(mf, c("code"), n)
    sex_vec <- take_or_na(mf, c("sex"), n)
    age_vec <- take_or_na(mf, c("age_group"), n)

    tibble::tibble(
      model = rep(nm, n),
      fitted = fitted_vec,
      pearson_residual = as.numeric(stats::residuals(m, type = "pearson")),
      deviance_residual = as.numeric(stats::residuals(m, type = "deviance")),
      response = as.numeric(stats::model.response(mf)),
      year = year_vec,
      code = code_vec,
      sex = sex_vec,
      age_group = age_vec
    )
  })

  dplyr::bind_rows(out)
}

summarise_main_model_diagnostics <- function(fit, resid_df) {
  model_map <- list(
    patients  = fit$m_pat,
    provision = fit$m_prov,
    value     = fit$m_val
  )

  dplyr::bind_rows(lapply(names(model_map), function(nm) {
    m <- model_map[[nm]]
    d <- resid_df %>% dplyr::filter(.data$model == nm)
    kc <- tryCatch(mgcv::k.check(m), error = function(e) NULL)
    k_index_min <- NA_real_
    k_p_min <- NA_real_
    if (!is.null(kc) && is.matrix(kc) && nrow(kc) > 0) {
      cn <- colnames(kc)
      if ("k-index" %in% cn) k_index_min <- suppressWarnings(min(kc[, "k-index"], na.rm = TRUE))
      if ("p-value" %in% cn) k_p_min <- suppressWarnings(min(kc[, "p-value"], na.rm = TRUE))
      if (!is.finite(k_index_min)) k_index_min <- NA_real_
      if (!is.finite(k_p_min)) k_p_min <- NA_real_
    }

    tibble::tibble(
      model = nm,
      family = family(m)$family,
      link = family(m)$link,
      n_obs = nrow(d),
      dev_expl = summary(m)$dev.expl %||% NA_real_,
      adj_r_sq = summary(m)$r.sq %||% NA_real_,
      aic = tryCatch(AIC(m), error = function(e) NA_real_),
      mean_pearson = mean(d$pearson_residual, na.rm = TRUE),
      sd_pearson = stats::sd(d$pearson_residual, na.rm = TRUE),
      mean_abs_pearson = mean(abs(d$pearson_residual), na.rm = TRUE),
      pct_abs_pearson_gt_2 = mean(abs(d$pearson_residual) > 2, na.rm = TRUE),
      cor_abs_pearson_fitted = suppressWarnings(stats::cor(abs(d$pearson_residual), d$fitted, use = "complete.obs")),
      k_index_min = k_index_min,
      k_check_p_min = k_p_min
    )
  }))
}

plot_main_model_residuals_vs_fitted <- function(resid_df) {
  ggplot(resid_df, aes(x = .data$fitted, y = .data$pearson_residual)) +
    geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed") +
    geom_point(alpha = 0.18, size = 0.7, color = "#00468B") +
    geom_smooth(se = FALSE, color = "#D55E00", linewidth = 0.7, method = "loess", formula = y ~ x) +
    facet_wrap(~ model, scales = "free_x") +
    scale_y_continuous(labels = label_number_pl()) +
    scale_x_continuous(labels = label_number_pl()) +
    labs(x = "Fitted value", y = "Pearson residual") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

plot_main_model_residual_histogram <- function(resid_df) {
  ggplot(resid_df, aes(x = .data$pearson_residual)) +
    geom_histogram(bins = 40, fill = "#0099B4", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed") +
    facet_wrap(~ model, scales = "free_y") +
    scale_x_continuous(labels = label_number_pl()) +
    labs(x = "Pearson residual", y = "Count") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

plot_main_model_residual_qq <- function(resid_df) {
  ggplot(resid_df, aes(sample = .data$deviance_residual)) +
    stat_qq(alpha = 0.18, size = 0.8, color = "#00468B") +
    stat_qq_line(color = "#D55E00", linewidth = 0.6) +
    facet_wrap(~ model, scales = "free") +
    labs(x = "Theoretical quantiles", y = "Deviance residual quantiles") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

plot_main_model_residuals_by_year <- function(resid_df) {
  ggplot(resid_df, aes(x = factor(.data$year), y = .data$pearson_residual)) +
    geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed") +
    geom_boxplot(outlier.alpha = 0.15, fill = "#7A7A7A") +
    facet_wrap(~ model, scales = "free_y") +
    scale_y_continuous(labels = label_number_pl()) +
    labs(x = "Year", y = "Pearson residual") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}


make_dharma_diagnostics <- function(fit, dat, out_dir, nsim = 250, width = 1800, height = 1400, res = 200) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    warning("Package {DHARMa} is not installed. Install it with install.packages('DHARMa') to generate simulation-based residual diagnostics.")
    return(NULL)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  model_map <- list(
    patients  = fit$m_pat,
    provision = fit$m_prov,
    value     = fit$m_val
  )

  test_rows <- lapply(names(model_map), function(nm) {
    m <- model_map[[nm]]
    message("Running DHARMa diagnostics for model: ", nm)

    simres <- tryCatch(
      DHARMa::simulateResiduals(fittedModel = m, n = nsim, plot = FALSE),
      error = function(e) e
    )

    if (inherits(simres, "error")) {
      warning("DHARMa failed for model '", nm, "': ", conditionMessage(simres))
      return(tibble::tibble(
        model = nm,
        dharma_ok = FALSE,
        dharma_message = conditionMessage(simres),
        uniformity_p = NA_real_,
        dispersion_p = NA_real_,
        outlier_p = NA_real_,
        zero_inflation_p = NA_real_
      ))
    }

    # overview plot
    grDevices::png(file.path(out_dir, paste0("dharma_overview_", nm, ".png")), width = width, height = height, res = res)
    plot(simres)
    grDevices::dev.off()

    # QQ / uniformity plot
    grDevices::png(file.path(out_dir, paste0("dharma_qq_", nm, ".png")), width = width, height = height, res = res)
    DHARMa::plotQQunif(simres)
    grDevices::dev.off()

    # residual vs predicted
    grDevices::png(file.path(out_dir, paste0("dharma_residual_vs_predicted_", nm, ".png")), width = width, height = height, res = res)
    tryCatch(
      DHARMa::plotResiduals(simres, form = simres$predictedResponse,
                            xlab = "Predicted value", ylab = "Scaled DHARMa residual"),
      error = function(e) warning("DHARMa residual-vs-predicted plot failed for ", nm, ": ", conditionMessage(e))
    )
    grDevices::dev.off()

    # residual vs year
    year_vec <- tryCatch(model.frame(m)$year, error = function(e) NULL)
    if (is.null(year_vec) || length(year_vec) != length(simres$scaledResiduals)) {
      year_vec <- dat$year
      if (length(year_vec) != length(simres$scaledResiduals)) year_vec <- seq_along(simres$scaledResiduals)
    }
    grDevices::png(file.path(out_dir, paste0("dharma_residual_vs_year_", nm, ".png")), width = width, height = height, res = res)
    tryCatch(
      DHARMa::plotResiduals(simres, form = year_vec,
                            xlab = "Year", ylab = "Scaled DHARMa residual"),
      error = function(e) warning("DHARMa residual-vs-year plot failed for ", nm, ": ", conditionMessage(e))
    )
    grDevices::dev.off()

    test_uniformity <- tryCatch(DHARMa::testUniformity(simres), error = function(e) NULL)
    test_dispersion <- tryCatch(DHARMa::testDispersion(simres), error = function(e) NULL)
    test_outliers <- tryCatch(DHARMa::testOutliers(simres, type = "bootstrap"), error = function(e) NULL)
    test_zi <- tryCatch(DHARMa::testZeroInflation(simres), error = function(e) NULL)

    tibble::tibble(
      model = nm,
      dharma_ok = TRUE,
      dharma_message = NA_character_,
      uniformity_p = if (!is.null(test_uniformity)) unname(test_uniformity$p.value) else NA_real_,
      dispersion_p = if (!is.null(test_dispersion)) unname(test_dispersion$p.value) else NA_real_,
      outlier_p = if (!is.null(test_outliers)) unname(test_outliers$p.value) else NA_real_,
      zero_inflation_p = if (!is.null(test_zi)) unname(test_zi$p.value) else NA_real_
    )
  })

  dplyr::bind_rows(test_rows)
}
