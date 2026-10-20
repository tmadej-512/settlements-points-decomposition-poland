# Refactor change log

- Integrated user-provided `config.R` (not guessed).
- Added canonical `scripts/run_pipeline.R` with deterministic stage order and flags.
- Kept scientific estimators and original scripts unchanged for audit comparison.
- Explicit preflight for data, config and packages; no fallback to stale aggregated RDS.
- Default to forced recompute of simulation and bootstrap; cache opt-in.
- Added execution log, sessionInfo and manifest; documentation and .gitignore.
- Isolated R14 as diagnostic stage; did not automatically swap family/specification.
- Supplied non-publishing quick smoke-test option.

NOT VERIFIED: R syntax/integration/runtime outputs, because R executable and NFZ
raw data are unavailable. Do not interpret this package as an independently
validated replication of manuscript figures.
