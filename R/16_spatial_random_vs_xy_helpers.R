# ============================================================
# R/16_spatial_random_vs_xy_helpers.R
# ============================================================
#
# Diagnostics for spatial structure in GAM residuals.
#
# Purpose:
#   Check whether an administrative-unit random effect is sufficient,
#   or whether an additional smooth spatial surface s(x, y) is needed.
#
# Model variants tested:
#   1. base
#      no administrative-unit random effect, no spatial smooth
#
#   2. random_region
#      + s(code, bs = "re") or equivalent region column
#
#   3. spatial_xy
#      + s(x, y, bs = "tp")
#
#   4. random_region_plus_spatial_xy
#      + s(code, bs = "re") + s(x, y, bs = "tp")
#
# Moran's I is calculated on residuals aggregated to administrative units
# by default. This avoids pseudo-replication caused by many age/sex/year
# rows sharing the same coordinates.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(spdep)
})

# ------------------------------------------------------------
# Add terms to a formula
# ------------------------------------------------------------

formula_add_terms <- function(formula_base, add_terms = NULL) {
  if (is.null(add_terms) || length(add_terms) == 0L) {
    return(formula_base)
  }

  formula_txt <- paste(deparse(formula_base), collapse = " ")
  add_txt <- paste(add_terms, collapse = " + ")

  as.formula(paste(formula_txt, "+", add_txt))
}

# ------------------------------------------------------------
# k-nearest-neighbour spatial weights
# ------------------------------------------------------------

make_knn_listw <- function(dt, x_col = "x", y_col = "y", k = 8L) {
  d <- as.data.table(dt)

  if (!x_col %in% names(d)) stop("Missing x column: ", x_col)
  if (!y_col %in% names(d)) stop("Missing y column: ", y_col)

  coords <- as.matrix(d[, .SD, .SDcols = c(x_col, y_col)])

  ok <- is.finite(coords[, 1]) & is.finite(coords[, 2])
  coords_ok <- coords[ok, , drop = FALSE]

  if (nrow(coords_ok) <= k + 2L) {
    stop("Too few observations for kNN weights after coordinate filtering.")
  }

  nb <- spdep::knn2nb(
    spdep::knearneigh(coords_ok, k = k)
  )

  lw <- spdep::nb2listw(
    nb,
    style = "W",
    zero.policy = TRUE
  )

  list(
    listw = lw,
    ok = ok
  )
}

# ------------------------------------------------------------
# Moran's I for residuals
# ------------------------------------------------------------

moran_for_model <- function(
  fit,
  data,
  x_col = "x",
  y_col = "y",
  group_cols = "code",
  k = 8L,
  residual_type = "pearson",
  nsim = 999L
) {
  d <- as.data.table(copy(data))

  r <- residuals(fit, type = residual_type)

  if (length(r) != nrow(d)) {
    stop(
      "Residual length does not match data rows. ",
      "Pass the exact data used in model fitting."
    )
  }

  d[, residual := as.numeric(r)]

  d <- d[
    is.finite(residual) &
      is.finite(get(x_col)) &
      is.finite(get(y_col))
  ]

  # Aggregate residuals to spatial units by default. This is usually the
  # right diagnostic object because the analysis data contain repeated
  # rows per spatial unit: age x sex x year x code.
  if (!is.null(group_cols) && length(group_cols) > 0L) {
    missing_group <- setdiff(group_cols, names(d))
    if (length(missing_group) > 0L) {
      stop("Missing group_cols in data: ", paste(missing_group, collapse = ", "))
    }

    d <- d[, c(
      .(residual = mean(residual, na.rm = TRUE)),
      lapply(.SD, function(z) z[which(is.finite(z))[1]])
    ), by = group_cols, .SDcols = c(x_col, y_col)]

    setnames(d, c(x_col, y_col), c(".x_moran", ".y_moran"))
    x_col_moran <- ".x_moran"
    y_col_moran <- ".y_moran"
  } else {
    x_col_moran <- x_col
    y_col_moran <- y_col
  }

  d <- d[
    is.finite(residual) &
      is.finite(get(x_col_moran)) &
      is.finite(get(y_col_moran))
  ]

  if (nrow(d) <= k + 2L) {
    return(data.table(
      n_moran = nrow(d),
      moran_i = NA_real_,
      expected_i = NA_real_,
      variance_i = NA_real_,
      p_value_asymptotic = NA_real_,
      p_value_mc = NA_real_
    ))
  }

  w <- make_knn_listw(
    d,
    x_col = x_col_moran,
    y_col = y_col_moran,
    k = k
  )

  r_ok <- d$residual[w$ok]

  mt <- spdep::moran.test(
    r_ok,
    w$listw,
    zero.policy = TRUE
  )

  mc <- spdep::moran.mc(
    r_ok,
    w$listw,
    nsim = nsim,
    zero.policy = TRUE
  )

  data.table(
    n_moran = length(r_ok),
    moran_i = unname(mt$estimate[["Moran I statistic"]]),
    expected_i = unname(mt$estimate[["Expectation"]]),
    variance_i = unname(mt$estimate[["Variance"]]),
    p_value_asymptotic = mt$p.value,
    p_value_mc = mc$p.value
  )
}

# ------------------------------------------------------------
# Extract smooth information
# ------------------------------------------------------------

extract_smooth_info <- function(fit, region_col = "code") {
  st <- tryCatch(
    as.data.table(summary(fit)$s.table, keep.rownames = "smooth"),
    error = function(e) NULL
  )

  if (is.null(st) || nrow(st) == 0L) {
    return(data.table(
      spatial_edf = NA_real_,
      spatial_p_value = NA_real_,
      region_re_edf = NA_real_,
      region_re_p_value = NA_real_
    ))
  }

  p_col <- intersect(c("p-value", "p.value", "Pr(>F)", "Pr(>Chi)"), names(st))
  if (length(p_col) == 0L) {
    st[, pval_tmp := NA_real_]
  } else {
    setnames(st, p_col[1], "pval_tmp")
  }

  edf_col <- intersect(c("edf", "EDF"), names(st))
  if (length(edf_col) == 0L) {
    st[, edf_tmp := NA_real_]
  } else {
    setnames(st, edf_col[1], "edf_tmp")
  }

  spatial_row <- st[grepl("x.*y|y.*x", smooth)]
  region_pat <- paste0("s\\(", region_col, "\\)")
  region_row <- st[grepl(region_pat, smooth)]

  data.table(
    spatial_edf = if (nrow(spatial_row) > 0) spatial_row$edf_tmp[1] else NA_real_,
    spatial_p_value = if (nrow(spatial_row) > 0) spatial_row$pval_tmp[1] else NA_real_,
    region_re_edf = if (nrow(region_row) > 0) region_row$edf_tmp[1] else NA_real_,
    region_re_p_value = if (nrow(region_row) > 0) region_row$pval_tmp[1] else NA_real_
  )
}

# ------------------------------------------------------------
# Fit metrics
# ------------------------------------------------------------

model_fit_metrics <- function(fit) {
  sm <- summary(fit)

  data.table(
    dev_expl = unname(sm$dev.expl),
    adj_r_sq = if (!is.null(sm$r.sq)) unname(sm$r.sq) else NA_real_,
    reml_score = if (!is.null(fit$gcv.ubre)) unname(fit$gcv.ubre) else NA_real_,
    aic = tryCatch(unname(AIC(fit)), error = function(e) NA_real_),
    n_obs_model = length(fit$fitted.values)
  )
}

# ------------------------------------------------------------
# Fit one variant
# ------------------------------------------------------------

fit_gam_variant <- function(
  formula_base,
  data,
  family,
  variant = c(
    "base",
    "random_region",
    "spatial_xy",
    "random_region_plus_spatial_xy"
  ),
  region_col = "code",
  x_col = "x",
  y_col = "y",
  spatial_k = 80L,
  method = "fREML",
  discrete = TRUE
) {
  variant <- match.arg(variant)

  add_terms <- character(0)

  if (variant %in% c("random_region", "random_region_plus_spatial_xy")) {
    add_terms <- c(
      add_terms,
      paste0("s(", region_col, ", bs = 're')")
    )
  }

  if (variant %in% c("spatial_xy", "random_region_plus_spatial_xy")) {
    add_terms <- c(
      add_terms,
      paste0("s(", x_col, ", ", y_col, ", bs = 'tp', k = ", spatial_k, ")")
    )
  }

  f <- formula_add_terms(formula_base, add_terms)

  mgcv::bam(
    formula = f,
    data = data,
    family = family,
    method = method,
    discrete = discrete,
    select = TRUE
  )
}

# ------------------------------------------------------------
# Full test for one component
# ------------------------------------------------------------

test_spatial_need_one_component <- function(
  component_name,
  formula_base,
  data,
  family,
  region_col = "code",
  x_col = "x",
  y_col = "y",
  spatial_k = 80L,
  knn_k = 8L,
  moran_group_cols = "code",
  nsim = 999L,
  residual_type = "pearson"
) {
  variants <- c(
    "base",
    "random_region",
    "spatial_xy",
    "random_region_plus_spatial_xy"
  )

  out <- lapply(variants, function(vv) {
    message("Fitting ", component_name, " | ", vv)

    fit <- tryCatch(
      fit_gam_variant(
        formula_base = formula_base,
        data = data,
        family = family,
        variant = vv,
        region_col = region_col,
        x_col = x_col,
        y_col = y_col,
        spatial_k = spatial_k
      ),
      error = function(e) {
        message("Model failed: ", component_name, " | ", vv, " | ", e$message)
        NULL
      }
    )

    if (is.null(fit)) {
      return(data.table(
        component = component_name,
        variant = vv,
        model_ok = FALSE
      ))
    }

    moran <- tryCatch(
      moran_for_model(
        fit = fit,
        data = data,
        x_col = x_col,
        y_col = y_col,
        group_cols = moran_group_cols,
        k = knn_k,
        residual_type = residual_type,
        nsim = nsim
      ),
      error = function(e) {
        message("Moran failed: ", component_name, " | ", vv, " | ", e$message)

        data.table(
          n_moran = NA_integer_,
          moran_i = NA_real_,
          expected_i = NA_real_,
          variance_i = NA_real_,
          p_value_asymptotic = NA_real_,
          p_value_mc = NA_real_
        )
      }
    )

    metrics <- model_fit_metrics(fit)
    smooths <- extract_smooth_info(fit, region_col = region_col)

    cbind(
      data.table(
        component = component_name,
        variant = vv,
        model_ok = TRUE
      ),
      metrics,
      smooths,
      moran
    )
  })

  rbindlist(out, fill = TRUE)
}

# ------------------------------------------------------------
# Decision table
# ------------------------------------------------------------

decide_spatial_need <- function(result_dt, alpha = 0.05) {
  d <- as.data.table(result_dt)

  out <- dcast(
    d,
    component ~ variant,
    value.var = c(
      "moran_i",
      "p_value_mc",
      "dev_expl",
      "reml_score",
      "aic",
      "spatial_edf",
      "spatial_p_value",
      "region_re_edf",
      "region_re_p_value"
    )
  )

  out[, random_sufficient := (
    !is.na(p_value_mc_random_region) &
      p_value_mc_random_region >= alpha
  )]

  out[, spatial_needed_after_random := (
    !random_sufficient &
      !is.na(p_value_mc_random_region_plus_spatial_xy) &
      p_value_mc_random_region_plus_spatial_xy >= alpha
  )]

  out[, spatial_still_problem := (
    !is.na(p_value_mc_random_region_plus_spatial_xy) &
      p_value_mc_random_region_plus_spatial_xy < alpha
  )]

  out[, recommendation := fifelse(
    random_sufficient,
    "random_region_wystarczy",
    fifelse(
      spatial_needed_after_random,
      "dodaj_s_xy",
      fifelse(
        spatial_still_problem,
        "s_xy_nie_wystarcza_rozwaz_MRF_lub_inny_spatial",
        "niejednoznaczne_sprawdz_diagnostyke"
      )
    )
  )]

  out[]
}
