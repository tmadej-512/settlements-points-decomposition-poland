# ============================================================
# R/00_safe_fwrite.R
# ------------------------------------------------------------
# Helper for writing data frames/data.tables that may contain list-columns.
# CSV is flat, so list-cells are converted to character/JSON before fwrite().
# ============================================================

safe_list_cell_to_chr <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (length(x) == 0) return(NA_character_)
  if (is.atomic(x)) return(paste(as.character(x), collapse = "; "))
  if (is.list(x)) {
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      return(as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null")))
    }
    return(paste(unlist(x), collapse = "; "))
  }
  as.character(x)
}

make_fwrite_safe <- function(x) {
  x <- data.table::as.data.table(x)
  list_cols <- names(x)[vapply(x, is.list, logical(1))]
  if (length(list_cols) > 0) {
    message("Converting list-columns before fwrite: ", paste(list_cols, collapse = ", "))
    x[, (list_cols) := lapply(
      .SD,
      function(col) vapply(col, safe_list_cell_to_chr, character(1))
    ), .SDcols = list_cols]
  }
  x
}
