# ============================================================
# config.R
# ------------------------------------------------------------
# Single source of truth for all run parameters -- plain R, no YAML.
# Edit the values below, not the scripts in R/ or scripts/, to
# change a study setting.
#
# This file defines one variable, `config`, a nested list with
# exactly the same structure used throughout the codebase as `cfg`
# (e.g. cfg$data$path, cfg$model$k_year). Loading it is just
# source("config.R") -- see load_config() in R/00_utils.R.
# ============================================================

config <- list(

  data = list(
    path            = "data/dane_2018_2024.csv",  # .csv, .rds, .xlsx, or .xls
    disease_filter  = "Neoplasms",
                              # if the file contains multiple groups, set exact icd10 label
    value_var       = "sum_points",  # "sum_points" or "val_provision" -- both columns are
                                      # present in this extract; sum_points is the NFZ
                                      # settlement-point total (the paper's primary outcome),
                                      # val_provision is its monetary equivalent
    baseline_year   = 2018
  ),

  model = list(
    k_year             = 5,       # max spline basis dimension for s(year_num)
    use_spatial_smooth = FALSE,   # adds s(x, y, bs = "tp") if TRUE and x/y present
    maxit_gam          = 200,     # mgcv default is much higher than the old maxit=10;
                                   # raise further (e.g. 500) if you see "algorithm did
                                   # not converge" warnings on the real data
    trace_gam          = FALSE,
    eps_offset         = 1,       # floor for log(offset) to avoid log(0)
    value_family       = "gamma"  # "tweedie" (default, safe if value can be 0 when
                                     # n_provision > 0) or "gamma" (preferred on parsimony
                                     # grounds if value is confirmed strictly positive
                                     # whenever n_provision > 0 -- check with
                                     # scripts/run_08_distribution_diagnostics.R first)
  ),

  decomposition = list(
    components       = c("N", "PI", "LAMBDA", "KAPPA"),
    run_naive        = TRUE,
    run_ofat         = TRUE,
    run_sequential   = TRUE,
    run_shapley      = TRUE,
    sequential_order = c("N", "PI", "LAMBDA", "KAPPA")  # order used for the single-order
                                                          # sequential comparator
  ),

  validation = list(
    run_axiom_checks      = TRUE,

    run_bootstrap         = TRUE,
    bootstrap_reps        = 200,
    bootstrap_seed        = 20260620,
    bootstrap_use_cache   = TRUE,
    bootstrap_force_rerun = FALSE,
    bootstrap_cache_file  = "output/cache/bootstrap_results.rds",

    run_simulation_study    = TRUE,
    simulation_reps         = 200,
    simulation_seed         = 1,
    simulation_use_cache    = TRUE,
    simulation_force_rerun  = FALSE,
    simulation_cache_file   = "output/cache/simulation_results.rds",
    simulation_two_scenarios_cache_file = "output/cache/simulation_two_scenarios_results.rds",

    # Synthetic DGP settings. These control the simulation validation only,
    # not the empirical NFZ model or bootstrap.
    simulation_n_units      = 60,
    simulation_years        = 2018:2024,
    simulation_comparison_year = 2024,
    simulation_sex_levels   = c("F", "M"),
    simulation_age_levels   = c("0-17", "18-44", "45-64", "65+"),
    simulation_heterogeneous_slopes = TRUE,
    simulation_regional_random_slopes = TRUE,
    simulation_regional_slope_sd = list(N = 0.003, PI = 0.008, LAMBDA = 0.008, KAPPA = 0.006),
    simulation_true_growth  = list(N = 0.005, PI = 0.025, LAMBDA = 0.035, KAPPA = -0.005),
    simulation_make_dgp_plots = TRUE,
    simulation_dgp_seed     = 123
  ),

  output = list(
    tables_dir  = "output/tables",
    figures_dir = "output/figures",
    fig_width   = 10,
    fig_height  = 6,
    fig_dpi     = 600
  )
)
