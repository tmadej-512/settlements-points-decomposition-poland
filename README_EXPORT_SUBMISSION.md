# Value in Health submission export: clean, versioned and audited

Place `scripts/export_submission.R` inside the refactored project. Start RStudio in the project root. **Use the newest ZIP source** as agreed; this script never modifies original `output/`, `config.R`, models or manuscript.

```r
system2('Rscript', c('scripts/export_submission.R', '--strict'))
```

This creates a unique `submission_exports/VIH_YYYYmmdd_HHMMSS/` with `tables/`, `figures/`, and `QA/`. It runs existing `run_98_export_vih_figures.R` in an isolated workspace, exporting TIFFs (600 DPI per config). It copies a curated list of paper tables and records hashes and input dates. `QA/issues.csv` lists missing tables, old split-format main comparison, and discrepancies between R13 and the main export. All generated figure names are sourced from the original exporter manifest. The script does not copy old `output/figures/` or old `submission/` graphics.

**First test, without figures:**

```r
system2('Rscript', c('scripts/export_submission.R', '--tables-only', '--strict'))
```

For each re-run, choose a fresh timestamp or `--run-id=RUN42`. Nothing overwrites previous exports. `--output-dir=other_folder` changes export root. On errors the export folder remains for diagnosis.

**IMPORTANT:** This is an exporter, not a re-run of model estimation. It cannot determine whether old `output/tables/` rows came from current model or stale caches. To certify scientific reproducibility, re-run `scripts/run_pipeline.R --full` and compare its numeric outputs. In `QA/issues.csv`, ERROR means do not submit until resolved, WARNING needs human review. The split-format `TABLE_empirical_shapley_vs_lmdi.csv` is *not* silently replaced; a paired `TABLE_empirical_comparison_RECONSTRUCTED_QA_ONLY.csv` is generated for checking. No raw NFZ data are bundled by design. Never upload `QA/build_workspace` if a run failed until confirming it contains no confidential material.

No historical execution logs are required. A new `QA/export.log` and exporter stdout/stderr are produced by this script.


## Curated supplementary figures (October 2026)

Submission export generates exactly three supplementary TIFFs:

- S1: order sensitivity (former S1)
- S2: bootstrap intervals (former S2)
- S3: Shapley–LMDI differences across simulation scenarios (former S5)

Former figures S3 (additivity residual), S4 (validation panel), and S6
(residual diagnostics) are deliberately excluded from submission figure export.
No GAM, bootstrap or simulation is rerun by `scripts/export_submission.R`.
Old exports remain untouched; use only the newest `submission_exports/VIH_*` directory.

**Manuscript/supplement synchronization:** remove the former S3/S4/S6 figure
blocks, rename former Figure S5 to Figure S3, and update the figure citations
in the manuscript and supplement. Supplementary **Table S3** is unchanged.
The script modifies code exports, not DOCX manuscripts.

Run: `Rscript scripts/export_submission.R --strict` from project root.
Expected QA: three main and three supplementary TIFFs.
