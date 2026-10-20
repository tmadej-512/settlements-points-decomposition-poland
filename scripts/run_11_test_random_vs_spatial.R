#!/usr/bin/env Rscript
# ============================================================
# scripts/run_11_test_random_vs_spatial.R
# ============================================================
#
# Purpose:
#   Test whether the production models need an explicit spatial smooth
#   s(x, y), or whether the administrative-unit random effect alone is
#   sufficient.
#
# It fits four variants for each component:
#   1. base
#   2. random_region                 + s(code, bs = "re")
#   3. spatial_xy                    + s(x, y)
#   4. random_region_plus_spatial_xy + s(code, bs = "re") + s(x, y)
#
# Then it calculates Moran's I on residuals aggregated to administrative
# units. This avoids pseudo-replication from repeated age/sex/year cells
# at the same coordinates.
#
# Outputs:
#   output/tables/TABLE_random_vs_spatial_moran_long.csv
#   output/tables/TABLE_random_vs_spatial_decision.csv
#   output/figures/FIGURE_random_vs_spatial_moran.png
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(ggplot2)
  library(scales)
})

source("R/99_load_all.R")
source("R/16_spatial_random_vs_xy_helpers.R")

if (!exists("fwrite_safe", mode = "function")) {
  fwrite_safe <- function(x, file, ...) {
    data.table::fwrite(make_fwrite_safe(x), file, ...)
  }
}

cfg <- load_config("config.R")
ensure_dir(cfg$output$tables_dir)
ensure_dir(cfg$output$figures_dir)

# ------------------------------------------------------------
# Options
# ------------------------------------------------------------

SPATIAL_K <- get0(".NFZ_SPATIAL_K", ifnotfound = 80L)
KNN_K <- get0(".NFZ_MORAN_K", ifnotfound = 8L)
MORAN_NSIM <- get0(".NFZ_MORAN_NSIM", ifnotfound = if (isTRUE(get0(".NFZ_QUICK", ifnotfound = FALSE))) 199L else 999L)

# ------------------------------------------------------------
# Load model data
# ------------------------------------------------------------

if (file.exists("output/rds/aggregated_data.rds")) {
  dt <- readRDS("output/rds/aggregated_data.rds")
  message("Loaded: output/rds/aggregated_data.rds")
} else {
  message("output/rds/aggregated_data.rds not found; rebuilding aggregated data from input file.")
  df_raw <- read_input_data(cfg$data$path)
  dat0 <- normalize_schema(df_raw, value_var = cfg$data$value_var)
  dt <- filter_and_aggregate(dat0, disease_filter = cfg$data$disease_filter)
}

dt <- as.data.table(dt)

# ------------------------------------------------------------
# Prepare variables used in formulas
# ------------------------------------------------------------

required_cols <- c(
  "code", "x", "y", "year", "sex", "age_group",
  "pop", "n_patients", "n_provision", "value"
)
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0L) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (all(is.na(dt$x)) || all(is.na(dt$y))) {
  stop(
    "Columns x and y are missing or all NA. Spatial diagnostics require coordinates.\n",
    "Add projected or consistent coordinates to the source data before running this script."
  )
}

dt <- dt[
  is.finite(x) & is.finite(y) &
    is.finite(pop) & is.finite(n_patients) &
    is.finite(n_provision) & is.finite(value)
]

if (nrow(dt) == 0L) {
  stop("No rows left after spatial-coordinate filtering.")
}

yrs <- sort(unique(dt$year))
k_year <- min(cfg$model$k_year, length(yrs))
if (k_year < 3) k_year <- length(yrs)

dt[, `:=`(
  year_num = as.numeric(year),
  sex_age = droplevels(interaction(sex, age_group, drop = TRUE)),
  code = droplevels(factor(code)),
  sex = droplevels(factor(sex)),
  age_group = droplevels(factor(age_group))
)]

# ------------------------------------------------------------
# Families
# ------------------------------------------------------------

value_family <- switch(
  cfg$model$value_family %||% "tweedie",
  "tweedie" = tw(link = "log"),
  "gamma" = Gamma(link = "log"),
  stop("config$model$value_family must be 'tweedie' or 'gamma'.")
)

# ------------------------------------------------------------
# Base formulas deliberately omit s(code, bs='re') and s(x,y).
# The four variants add them explicitly.
# ------------------------------------------------------------

base_rhs <- c(
  "sex * age_group",
  paste0("s(year_num, by = sex_age, k = ", k_year, ")")
)

make_diag_formula <- function(response, offset_var, eps = cfg$model$eps_offset) {
  rhs <- c(
    base_rhs,
    paste0("offset(log(pmax(", offset_var, ", ", eps, ")))")
  )
  as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
}

component_specs <- list(
  patients = list(
    formula = make_diag_formula("n_patients", "pop"),
    family = nb()
  ),
  provision = list(
    formula = make_diag_formula("n_provision", "n_patients"),
    family = nb()
  ),
  value = list(
    formula = make_diag_formula("value", "n_provision"),
    family = value_family
  )
)

# ------------------------------------------------------------
# Run diagnostics
# ------------------------------------------------------------

all_results <- lapply(names(component_specs), function(cc) {
  spec <- component_specs[[cc]]

  test_spatial_need_one_component(
    component_name = cc,
    formula_base = spec$formula,
    data = dt,
    family = spec$family,
    region_col = "code",
    x_col = "x",
    y_col = "y",
    spatial_k = SPATIAL_K,
    knn_k = KNN_K,
    moran_group_cols = "code",
    nsim = MORAN_NSIM,
    residual_type = "pearson"
  )
})

random_vs_spatial_long <- rbindlist(all_results, fill = TRUE)

random_vs_spatial_long[, `:=`(
  spatial_k = SPATIAL_K,
  knn_k = KNN_K,
  moran_nsim = MORAN_NSIM,
  moran_residual_aggregation = "code"
)]

fwrite_safe(
  random_vs_spatial_long,
  file.path(cfg$output$tables_dir, "TABLE_random_vs_spatial_moran_long.csv")
)

cat("\n============================================================\n")
cat("RANDOM VS SPATIAL: LONG TABLE\n")
cat("============================================================\n")
print(random_vs_spatial_long)

# ------------------------------------------------------------
# Decision table
# ------------------------------------------------------------

random_vs_spatial_decision <- decide_spatial_need(
  random_vs_spatial_long,
  alpha = 0.05
)

fwrite_safe(
  random_vs_spatial_decision,
  file.path(cfg$output$tables_dir, "TABLE_random_vs_spatial_decision.csv")
)

cat("\n============================================================\n")
cat("RANDOM VS SPATIAL: DECISION TABLE\n")
cat("============================================================\n")
print(random_vs_spatial_decision)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

plot_dt <- random_vs_spatial_long[
  model_ok == TRUE &
    is.finite(moran_i)
]

plot_dt[, variant_label := factor(
  variant,
  levels = c(
    "base",
    "random_region",
    "spatial_xy",
    "random_region_plus_spatial_xy"
  ),
  labels = c(
    "Base",
    "Region random",
    "Spatial s(x,y)",
    "Region random + s(x,y)"
  )
)]

plot_dt[, component_label := factor(
  component,
  levels = c("patients", "provision", "value"),
  labels = c("Patients", "Services", "Settlement value")
)]

p_random_vs_spatial_moran <- ggplot(
  plot_dt,
  aes(x = variant_label, y = moran_i)
) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  facet_wrap(~ component_label, scales = "fixed") +
  labs(
    title = "Residual spatial autocorrelation by model variant",
    subtitle = "Moran's I calculated on residuals aggregated to administrative units",
    x = NULL,
    y = "Moran's I"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1, size = 7.8)
  )

ggsave(
  file.path(cfg$output$figures_dir, "FIGURE_random_vs_spatial_moran.png"),
  p_random_vs_spatial_moran,
  width = 8.0,
  height = 4.8,
  dpi = cfg$output$fig_dpi %||% 600
)

print(p_random_vs_spatial_moran)

cat("\n============================================================\n")
cat("RANDOM VS SPATIAL DIAGNOSTICS FINISHED\n")
cat("Saved:\n")
cat("  ", file.path(cfg$output$tables_dir, "TABLE_random_vs_spatial_moran_long.csv"), "\n", sep = "")
cat("  ", file.path(cfg$output$tables_dir, "TABLE_random_vs_spatial_decision.csv"), "\n", sep = "")
cat("  ", file.path(cfg$output$figures_dir, "FIGURE_random_vs_spatial_moran.png"), "\n", sep = "")
cat("============================================================\n")
