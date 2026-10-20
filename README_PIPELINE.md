# NFZ / Value in Health — canonical reproducible run

**Canonical entry point:** `scripts/run_pipeline.R`. The 8-stage legacy main driver and
`run_12`–`run_14` are retained for audit and backward compatibility; do not use them
as alternative sources of published results. Statistical functions in `R/` were not rewritten.

## Requirements

1. Install R and dependencies listed in the source scripts (including `mgcv`, `dplyr`,
   `readr`, `ggplot2`, `tidyr`, `purrr`, `stringr`, `tibble`, `scales`; `spdep` for Moran diagnostics).
2. Put your licensed NFZ extract at `data/dane_2018_2024.csv` or edit
   `config.R` -> `config$data$path`. Input data are not shipped.
3. Keep `config.R` in the project root (supplied by the author); `config.example.R` is a template.
4. Run from the project root.

```sh
Rscript scripts/run_pipeline.R --preflight
Rscript scripts/run_pipeline.R --quick
Rscript scripts/run_pipeline.R --full --diagnostics --figures
```

The default run **forces bootstrap and simulation recomputation**, because the prior
caches verified only seeds/replicate counts, not hashes of all input data and model settings.
Use `--reuse-cache` only when deliberately reproducing a previously verified run.
Quick mode uses 20 bootstrap/simulation replicates and is NOT for publication.
`--diagnostics` runs R14 *after* the main chain; its model-selection recommendations
are diagnostic, and do not silently overwrite the manuscript's chosen models.
`--figures` additionally produces figures and figure legends.

## Outputs and audits

- `output/pipeline_execution_log.csv`: completed stages and failures.
- `output/sessionInfo.txt`: R environment versions.
- `output/pipeline_source_manifest.csv`: MD5 hashes of configuration and executed entry scripts.
- Original outputs of each retained script are under `output/`, as before.

## Important limitations

- No empirical rerun has been completed in this package: the user did not provide
  the NFZ input dataset and the build environment had no R runtime.
- Script R14 is diagnostic and independently refits candidates. R12/R13 are
  historical alternate analysis paths and are not part of the canonical pipeline.
- The three-stage GAM models use the setting `model$value_family = "gamma"`,
  `model$use_spatial_smooth = FALSE` in the supplied config. Changing these
  settings changes the analysis and requires reviewing the manuscript numbers.
- Some legacy outputs may be sensitive to scripts' global state. The canonical
  execution order has not yet been empirically integration-tested.
- Do not publish the NFZ dataset or any identifying proprietary extracts.
