#!/usr/bin/env Rscript
# Optional helper: permanently change config.R model$value_family from tweedie to gamma.
# The full recompute script can also run without this, using .NFZ_VALUE_FAMILY_OVERRIDE.

if (!file.exists("config.R")) stop("config.R not found. Run from project root.", call. = FALSE)

x <- readLines("config.R", warn = FALSE)
old <- x
x <- gsub('value_family\s*=\s*"tweedie"', 'value_family = "gamma"', x)
x <- gsub("value_family\s*=\s*'tweedie'", 'value_family = "gamma"', x)

if (identical(x, old)) {
  message("No explicit value_family = 'tweedie' found. Please check config.R manually.")
} else {
  backup <- paste0("config.R.backup_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  file.copy("config.R", backup, overwrite = FALSE)
  writeLines(x, "config.R")
  message("Updated config.R: value_family = 'gamma'. Backup: ", backup)
}
