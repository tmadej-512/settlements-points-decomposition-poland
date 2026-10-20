#!/usr/bin/env Rscript
# ============================================================
# scripts/run_06_bootstrap.R
# ------------------------------------------------------------
# Cluster bootstrap with RDS cache.
#
# Expensive bootstrap results are saved to
# cfg$validation$bootstrap_cache_file. Re-running this script will
# reuse the cache when bootstrap_use_cache = TRUE and
# bootstrap_force_rerun = FALSE, provided the cache matches the
# requested number of reps/seed.
# ============================================================

source("R/99_load_all.R")

cfg <- load_config("config.R")
make_plots <- !(exists(".NFZ_SKIP_PLOTS") && isTRUE(.NFZ_SKIP_PLOTS)) && !("--no-plots" %in% commandArgs(trailingOnly = TRUE))
dat <- readRDS("output/rds/aggregated_data.rds")

reps <- cfg$validation$bootstrap_reps
args <- commandArgs(trailingOnly = TRUE)
if (exists(".NFZ_REPS_OVERRIDE")) {
  reps <- .NFZ_REPS_OVERRIDE
} else if ("--reps" %in% args) {
  reps <- as.integer(args[which(args == "--reps") + 1])
}

boot_cache_file <- cfg$validation$bootstrap_cache_file %||% "output/cache/bootstrap_results.rds"
boot_use_cache <- cfg$validation$bootstrap_use_cache %||% TRUE
boot_force_rerun <- cfg$validation$bootstrap_force_rerun %||% FALSE
if ("--force" %in% args || "--force-bootstrap" %in% args ||
    (exists(".NFZ_FORCE_RERUN") && isTRUE(.NFZ_FORCE_RERUN)) ||
    (exists(".NFZ_BOOTSTRAP_FORCE_RERUN") && isTRUE(.NFZ_BOOTSTRAP_FORCE_RERUN))) {
  boot_force_rerun <- TRUE
}

baseline_year   <- cfg$data$baseline_year
comparison_year <- max(dat$year)

ensure_dir(cfg$output$tables_dir)
ensure_dir(cfg$output$figures_dir)
ensure_dir_for_file(boot_cache_file)

cache_ok <- FALSE
if (should_use_cache(boot_use_cache, boot_force_rerun, boot_cache_file)) {
  message("Checking cached bootstrap results: ", boot_cache_file)
  cached <- tryCatch(cache_read(boot_cache_file), error = function(e) NULL)
  cache_ok <- !is.null(cached) &&
    cache_metadata_matches(cached, expected_reps = reps, expected_seed = cfg$validation$bootstrap_seed) &&
    !is.null(cached$raw) && !is.null(cached$summary) &&
    is.data.frame(cached$raw) && is.data.frame(cached$summary) &&
    nrow(cached$raw) > 0 && nrow(cached$summary) > 0

  if (cache_ok) {
    message("Using cached bootstrap results: ", boot_cache_file)
    boot <- cached
  } else {
    message("Bootstrap cache exists but is stale/incomplete or uses different reps/seed; recomputing.")
  }
}

if (!cache_ok) {
  message("Bootstrap (", reps, " reps) for ", baseline_year, " vs ", comparison_year, " ...")

  boot <- run_bootstrap(
    df = dat, model_cfg = cfg$model,
    baseline_year = baseline_year, comparison_year = comparison_year,
    reps = reps, seed = cfg$validation$bootstrap_seed,
    methods = c("naive", "ofat", "sequential", "shapley")
  )

  boot$metadata <- list(
    step = "bootstrap",
    reps = as.integer(reps),
    seed = as.integer(cfg$validation$bootstrap_seed),
    baseline_year = as.integer(baseline_year),
    comparison_year = as.integer(comparison_year),
    methods = c("naive", "ofat", "sequential", "shapley"),
    created_at = as.character(Sys.time()),
    cache_version = "bootstrap-cache-v1"
  )

  cache_write(boot, boot_cache_file)
  message("Saved bootstrap cache: ", boot_cache_file)
}

message("\n== Bootstrap summary (percentile 95% CIs) ==")
print(boot$summary)

readr::write_csv(boot$summary, file.path(cfg$output$tables_dir, "bootstrap_summary.csv"))
readr::write_csv(boot$raw, file.path(cfg$output$tables_dir, "bootstrap_raw.csv"))

if (make_plots) {
  method_comparison_path <- file.path(cfg$output$tables_dir, "method_comparison.csv")
  if (file.exists(method_comparison_path)) {
    point_estimates <- readr::read_csv(method_comparison_path, show_col_types = FALSE)
    p <- plot_bootstrap_intervals(boot$summary, point_estimates)
    safe_ggsave(p, file.path(cfg$output$figures_dir, "bootstrap_intervals.png"),
                cfg$output$fig_width, cfg$output$fig_height, cfg$output$fig_dpi)
    message("Saved: ", file.path(cfg$output$figures_dir, "bootstrap_intervals.png"))
  } else {
    message("Skipping bootstrap interval plot: run scripts/run_02_compare_methods.R first ",
            "to get point estimates (", method_comparison_path, " not found).")
  }
} else {
  message("Skipping bootstrap interval plot because .NFZ_SKIP_PLOTS=TRUE or --no-plots was used.")
}

message("\nSaved: ", file.path(cfg$output$tables_dir, "bootstrap_summary.csv"))
message("Saved: ", file.path(cfg$output$tables_dir, "bootstrap_raw.csv"))
message("Done.")
