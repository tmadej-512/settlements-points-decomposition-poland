#!/usr/bin/env Rscript
# Canonical Value in Health pipeline runner (refactor of orchestration only).
# Run from project root: Rscript scripts/run_pipeline.R --full
# Optional: --quick (20 simulation/bootstrap reps), --diagnostics, --figures,
#           --reuse-cache, --preflight
# No statistical estimator has been changed by this orchestrator.
args <- commandArgs(trailingOnly = TRUE)
allowed <- c('--full','--quick','--diagnostics','--figures','--reuse-cache','--preflight','--help')
unknown <- setdiff(args, allowed)
if (length(unknown)) stop('Unknown arguments: ', paste(unknown, collapse=', '))
if ('--help' %in% args) {
  cat('Rscript scripts/run_pipeline.R [--full|--quick] [--diagnostics] [--figures] [--reuse-cache] [--preflight]\n')
  quit(status=0)
}
if (!file.exists('config.R') || !file.exists('R/99_load_all.R'))
  stop('Run from project root (config.R and R/ required).')
source('R/00_utils.R', local=TRUE)
cfg <- load_config('config.R')
if (!is.list(cfg$data) || !is.list(cfg$model) || !is.list(cfg$validation))
  stop('Malformed config.R: expected data, model, validation sections.')
required <- c('data','model','validation','output','decomposition')
if (!all(required %in% names(cfg))) stop('Missing config sections: ',paste(setdiff(required,names(cfg)),collapse=', '))
if (!identical(tolower(cfg$model$value_family),'gamma'))
  warning('Manuscript specifies Gamma for valuation, but config.R uses ',cfg$model$value_family)
if (!identical(cfg$data$value_var,'sum_points'))
  warning('Manuscript uses settlement points; config$data$value_var is ',cfg$data$value_var)
if (!identical(cfg$data$disease_filter,'Neoplasms'))
  warning('Manuscript illustration uses Neoplasms; configuration differs.')
if (!file.exists(cfg$data$path))
  stop('Missing NFZ input: ',cfg$data$path,'\nPlace licensed data at configured path; data are not distributed.')
required_packages <- c('dplyr','tidyr','readr','stringr','ggplot2','mgcv','purrr','tibble','scales','data.table','MASS')
if ('--diagnostics' %in% args) required_packages <- c(required_packages,'spdep')
missing_packages <- required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing_packages)) stop('Install R packages: ',paste(missing_packages,collapse=', '))
if ('--preflight' %in% args) {message('Preflight passed.');quit(status=0)}
quick <- '--quick' %in% args
reuse <- '--reuse-cache' %in% args
if (quick) .NFZ_REPS_OVERRIDE <- 20L
# Force expensive recomputation by default: existing legacy caches check only
# seeds and replicate counts, not the complete model/data configuration.
if (!reuse) .NFZ_FORCE_RERUN <- TRUE
.NFZ_SKIP_PLOTS <- !('--figures' %in% args)
ensure_dir('output')
run_log <- file.path('output','pipeline_execution_log.csv')
stages <- list(
 c('01_fit_model','scripts/run_01_fit_model.R'),
 c('02_compare_methods','scripts/run_02_compare_methods.R'),
 c('03_order_sensitivity','scripts/run_03_order_sensitivity.R'),
 c('04_axioms','scripts/run_04_validate_axioms.R'),
 c('05_simulation_all_five_scenarios','scripts/run_05b_simulation_two_scenarios.R'),
 c('06_bootstrap','scripts/run_06_bootstrap.R'),
 c('07_paper_tables','scripts/run_07_make_paper_outputs.R'),
 c('08_vih_tables','scripts/run_10_value_in_health_outputs.R')
)
if ('--diagnostics' %in% args) stages <- c(stages,list(
 c('09_residual_diagnostics','scripts/run_01b_residual_diagnostics.R'),
 c('10_distribution_diagnostics','scripts/run_08_distribution_diagnostics.R'),
 c('11_family_and_moran_diagnostics','scripts/run_14_family_selection_with_moran_by_distribution.R')
))
if ('--figures' %in% args) stages <- c(stages,list(
 c('12_all_plots','scripts/run_99_make_all_plots.R'),
 c('13_validation_plots','scripts/run_97_make_diagnostics_validation_plots.R'),
 c('14_submission_figures','scripts/run_98_export_vih_figures.R'),
 c('15_figure_legends','scripts/run_98b_make_figure_legends.R')
))
if (file.exists(run_log)) unlink(run_log)
write.csv(data.frame(stage=character(),status=character(),duration_sec=numeric(),completed=character()),run_log,row.names=FALSE)
start_time <- Sys.time()
for (stage in stages) {
  label <- stage[[1]]; path <- stage[[2]]
  message('\n===== ',label,' [',path,'] =====')
  began <- Sys.time()
  tryCatch({
    # local = FALSE preserves the existing scripts' shared-session conventions.
    source(path,local=FALSE)
    entry <- data.frame(stage=label,status='OK',duration_sec=as.numeric(difftime(Sys.time(),began,units='secs')),completed=as.character(Sys.time()))
    write.table(entry,run_log,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)
  },error=function(e){
    entry <- data.frame(stage=label,status='FAILED',duration_sec=as.numeric(difftime(Sys.time(),began,units='secs')),completed=as.character(Sys.time()))
    write.table(entry,run_log,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)
    stop('Pipeline stopped at ',label,': ',conditionMessage(e),call.=FALSE)
  })
}
writeLines(capture.output(sessionInfo()),'output/sessionInfo.txt')
source_files <- c('config.R', vapply(stages, function(stage) stage[[2]], character(1)))
manifest <- data.frame(file=source_files, md5=unname(tools::md5sum(source_files)))
write.csv(manifest,'output/pipeline_source_manifest.csv',row.names=FALSE)
message('Pipeline complete in ',round(as.numeric(difftime(Sys.time(),start_time,units='mins')),2),' minutes.')
