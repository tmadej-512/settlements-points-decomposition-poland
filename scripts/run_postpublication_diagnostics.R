#!/usr/bin/env Rscript
# Optional analyses of saved GAM objects. Does not refit the primary pipeline.
# Usage: Rscript scripts/run_postpublication_diagnostics.R --help
args <- commandArgs(trailingOnly=TRUE)
allowed <- c('--help','--moran-re','--moran-block','--sensitivity','--audit','--all')
if (any(!args %in% allowed)) stop('Unknown options: ', paste(setdiff(args,allowed),collapse=', '))
if (!length(args) || '--help' %in% args) {
 cat('Run optional checks with EXISTING local RDS models, no refitting:\n',
 '  --moran-re   RE cell slices\n',
 '  --moran-block  RE/MRF/RE+MRF block Moran\n',
 '  --sensitivity  16-state fast Shapley for saved patient variants\n',
 '  --audit  independent original-helper audit of 16 RE states\n',
 '  --all  all the above (requires all cached models)\n',sep='')
 quit(status=0)
}
requested <- if ('--all' %in% args) setdiff(allowed,c('--all','--help')) else args
steps <- list('--moran-re'='scripts/run_moran_re_only.R',
              '--moran-block'='scripts/run_moran_cell_block_pilot.R',
              '--sensitivity'='scripts/run_fast_shapley_model_sensitivity.R',
              '--audit'='scripts/run_audit_fast_shapley_RE.R')
if (!file.exists('config.R')) stop('Run from project root; copy config.example.R to config.R first')
for (name in names(steps)) if (name %in% requested) {
 cat('Starting ',name,'\n')
 source(steps[[name]],local=FALSE)
}
