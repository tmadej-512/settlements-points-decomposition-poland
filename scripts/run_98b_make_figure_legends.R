#!/usr/bin/env Rscript
# ============================================================
# scripts/run_98b_make_figure_legends.R
# ============================================================
# Creates separate Figure Legends files for Value in Health submission.
# This all-figures version writes legends for main and supplementary figures.
#
# Outputs:
#   submission/Figure_Legends_VIH.txt
#   submission/Figure_Legends_VIH.rtf
#
# The RTF uses double line spacing. Main manuscript files should cite figures
# in text only; image files are uploaded separately.
# ============================================================

suppressPackageStartupMessages({ library(data.table) })

dir.create("submission", recursive = TRUE, showWarnings = FALSE)

manifest_path <- "submission/figures/FIGURE_MANIFEST.csv"
manifest <- if (file.exists(manifest_path)) data.table::fread(manifest_path) else NULL

legend_blocks <- list(
  list(
    figure = "Figure 1",
    title = "Empirical method comparison.",
    text = "Component contributions to the change in total public-payer settlement points between 2018 and 2024. Components are population size, treated-patient reach, service intensity, and point valuation. Values are shown in billion settlement points. LMDI denotes logarithmic mean Divisia index; OFAT, one-factor-at-a-time substitution."
  ),
  list(
    figure = "Figure 2",
    title = "Multi-year decomposition.",
    text = "Component contributions to changes in public-payer settlement points from the 2018 baseline to each comparison year. Lines compare the aggregate LMDI-style benchmark and model-based Shapley decomposition. Values are shown in billion settlement points."
  ),
  list(
    figure = "Figure 3",
    title = "Simulation performance.",
    text = "Root mean squared error of estimated component contributions across decomposition methods and validation scenarios. Scenarios include oracle decomposition, clean model-based recovery, heterogeneous stress, cell-composition shift, and model misspecification. Values are shown in million settlement points. Lower values indicate better recovery of known component contributions. RMSE denotes root mean squared error."
  ),
  list(
    figure = "Supplementary Figure S1",
    title = "Order sensitivity.",
    text = "Distribution of component contributions across all sequential component-replacement orders. The Shapley contribution is the average marginal contribution across replacement orders. Values are shown in billion settlement points."
  ),
  list(
    figure = "Supplementary Figure S2",
    title = "Bootstrap intervals.",
    text = "Bootstrap uncertainty intervals for empirical component contributions. Points are original estimates and intervals are bootstrap percentile intervals. Values are shown in billion settlement points."
  ),
  list(
    figure = "Supplementary Figure S3",
    title = "Additivity residual.",
    text = "Residual equal to total change minus the sum of component effects by decomposition method. Sequential substitution and Shapley satisfy additivity by construction; OFAT is not a full decomposition. Values are shown in billion settlement points. OFAT denotes one-factor-at-a-time substitution."
  ),
  list(
    figure = "Supplementary Figure S4",
    title = "Validation panel.",
    text = "Four-panel validation summary. Panel A compares empirical component effects by method; Panel B shows order sensitivity across sequential substitution orders; Panel C compares estimated versus known true component contributions in simulation; Panel D shows bootstrap uncertainty intervals. Values are shown in billion or million settlement points as indicated by the axis labels."
  ),
  list(
    figure = "Supplementary Figure S5",
    title = "Shapley versus LMDI-style differences.",
    text = "Component-level differences between model-based Shapley attribution and the aggregate LMDI-style benchmark. Positive values indicate larger Shapley contributions; negative values indicate larger LMDI-style contributions. LMDI denotes logarithmic mean Divisia index. Values are shown in million settlement points."
  ),
  list(
    figure = "Supplementary Figure S6",
    title = "Residual spatial autocorrelation.",
    text = "Moran's I for residuals under alternative generalized additive model specifications: no unit adjustment, administrative-unit random effect, two-dimensional spatial smooth, and combined model. Spatial weights are based on nearest-neighbour relationships among administrative units. This is a diagnostic sensitivity figure rather than a main display item."
  )
)

if (!is.null(manifest) && nrow(manifest) > 0L) {
  exported <- manifest$figure
  keep <- vapply(legend_blocks, function(z) z$figure %in% exported, logical(1))
  legend_blocks <- legend_blocks[keep]
}

txt_lines <- c("Figure Legends", "")
for (b in legend_blocks) {
  txt_lines <- c(txt_lines, paste0(b$figure, ". ", b$title), b$text, "")
}
writeLines(txt_lines, "submission/Figure_Legends_VIH.txt")

escape_rtf <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\\{", "\\\\{", x)
  x <- gsub("\\}", "\\\\}", x)
  x
}

rtf_lines <- c(
  "{\\rtf1\\ansi\\deff0",
  "{\\fonttbl{\\f0 Times New Roman;}}",
  "\\paperw12240\\paperh15840\\margl1440\\margr1440\\margt1440\\margb1440",
  "\\f0\\fs24\\sl480\\slmult1",
  "\\b Figure Legends\\b0\\par\\par"
)
for (b in legend_blocks) {
  rtf_lines <- c(
    rtf_lines,
    paste0("\\b ", escape_rtf(paste0(b$figure, ". ", b$title)), "\\b0\\par"),
    paste0(escape_rtf(b$text), "\\par\\par")
  )
}
rtf_lines <- c(rtf_lines, "}")
writeLines(rtf_lines, "submission/Figure_Legends_VIH.rtf", useBytes = TRUE)

cat("Saved:\n")
cat("  submission/Figure_Legends_VIH.txt\n")
cat("  submission/Figure_Legends_VIH.rtf\n")
