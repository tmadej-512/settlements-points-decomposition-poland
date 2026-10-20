# ============================================================
# R/00_utils.R
# ------------------------------------------------------------
# Small, generic helpers used throughout the pipeline:
# config loading, safe coercions, column-name matching,
# and plot label formatters.
#
# No project-specific logic lives here.
# ============================================================

suppressPackageStartupMessages({
  library(scales)
})

#' Load the project configuration
#'
#' Sources config.R, a plain-R file defining a single nested list
#' variable named `config` (e.g. config$data$path, config$model$k_year).
#' This replaces an earlier YAML-based config.yml -- the structure
#' returned is identical (a nested list with the same names), so every
#' other script in this project that uses `cfg$...` works unchanged
#' regardless of which config format was used to produce it.
#'
#' @param path Path to config.R (default: "config.R" in the working directory).
#' @return A nested list with the parsed configuration.
load_config <- function(path = "config.R") {
  if (!file.exists(path)) {
    stop("Config file not found at: ", path,
         "\nRun this script from the project root, or pass an explicit path.")
  }
  config_env <- new.env()
  source(path, local = config_env)
  if (!exists("config", envir = config_env, inherits = FALSE)) {
    stop("config.R was found but did not define a variable named `config`. ",
         "Check that config.R contains `config <- list(...)`.")
  }
  get("config", envir = config_env)
}

#' `%||%`: return `a` unless it is NULL / length-0 / all-NA, else `b`
`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
}

#' Coerce to numeric, suppressing the "NAs introduced by coercion" warning
safe_numeric <- function(x) suppressWarnings(as.numeric(x))

#' Pick the first candidate name that exists in `nm`
#'
#' Used to make the schema normaliser tolerant of different
#' column-naming conventions across data extracts (English, Polish,
#' TERYT-style administrative exports, etc.).
#'
#' @param nm Character vector of available column names.
#' @param cands Character vector of candidate names, in priority order.
#' @return The first match, or NA_character_ if none match.
pick_first <- function(nm, cands) {
  hit <- intersect(cands, nm)
  if (length(hit)) hit[1] else NA_character_
}

#' Polish-locale number formatter for plot axes (space thousands sep, comma decimal)
label_number_pl <- function(accuracy = 1, big.mark = " ") {
  scales::label_number(accuracy = accuracy, big.mark = big.mark, decimal.mark = ",")
}

#' Polish-locale percent formatter for plot axes
label_percent_pl <- function(accuracy = 1) {
  scales::label_percent(accuracy = accuracy, decimal.mark = ",")
}

#' Ensure a directory exists, creating it (recursively) if needed
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

#' Save a ggplot only if it is non-NULL (skip gracefully on failed plots)
safe_ggsave <- function(plot, filename, width, height, dpi) {
  if (is.null(plot)) {
    message("  [skip] plot not generated: ", basename(filename))
    return(invisible(NULL))
  }
  ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi)
  invisible(filename)
}


# ============================================================
# CACHE HELPERS
# ============================================================

#' Ensure the parent directory for a file path exists
ensure_dir_for_file <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

#' Read an RDS cache object
cache_read <- function(path) {
  readRDS(path)
}

#' Save an object to an RDS cache file
cache_write <- function(object, path) {
  ensure_dir_for_file(path)
  saveRDS(object, path)
  invisible(path)
}

#' Basic cache gate
should_use_cache <- function(use_cache, force_rerun, cache_file) {
  isTRUE(use_cache) && !isTRUE(force_rerun) && file.exists(cache_file)
}

#' Get a nested config value with a fallback
cfg_get <- function(cfg, path, default = NULL) {
  x <- cfg
  for (nm in path) {
    if (is.null(x) || is.null(x[[nm]])) return(default)
    x <- x[[nm]]
  }
  x
}

#' Validate cached object metadata for expensive steps
cache_metadata_matches <- function(obj, expected_reps = NULL, expected_seed = NULL) {
  md <- obj$metadata
  if (is.null(md)) return(FALSE)
  if (!is.null(expected_reps) && !identical(as.integer(md$reps), as.integer(expected_reps))) return(FALSE)
  if (!is.null(expected_seed) && !identical(as.integer(md$seed), as.integer(expected_seed))) return(FALSE)
  TRUE
}

#' Check whether a simulation raw table contains all required method/component cells
simulation_has_required_methods <- function(sim_raw,
                                            required_methods = c("naive", "ofat", "sequential", "shapley"),
                                            required_components = c("N", "PI", "LAMBDA", "KAPPA")) {
  if (is.null(sim_raw) || !is.data.frame(sim_raw) || nrow(sim_raw) == 0) return(FALSE)
  if (!all(c("method", "component") %in% names(sim_raw))) return(FALSE)
  present <- unique(paste(as.character(sim_raw$method), as.character(sim_raw$component), sep = "::"))
  required <- as.vector(outer(required_methods, required_components, paste, sep = "::"))
  all(required %in% present)
}
