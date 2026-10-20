#!/usr/bin/env Rscript
# Value in Health submission export, non-destructive, with provenance and QA.
# Invoke from project root: Rscript scripts/export_submission.R
# Flags: --tables-only, --strict, --output-dir=PATH, --run-id=NAME
options(stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)
opt <- function(prefix, default=NULL) { x <- grep(paste0('^',prefix,'='),args,value=TRUE); if(length(x)) substring(tail(x,1),nchar(prefix)+2) else default }
tables_only <- '--tables-only' %in% args
strict <- '--strict' %in% args
run_id <- opt('--run-id',format(Sys.time(),'%Y%m%d_%H%M%S'))
if(!grepl('^[A-Za-z0-9_-]+$',run_id)) stop('Invalid run ID')
root <- normalizePath('.',winslash='/',mustWork=TRUE)
if(!file.exists('config.R') || !file.exists('scripts/run_98_export_vih_figures.R')) stop('Run from refactored project root containing config.R and scripts/.')
source('config.R')
tdir <- config$output$tables_dir
if(!dir.exists(tdir)) stop('Source tables do not exist: ',tdir,'; run the analysis first.')
base <- opt('--output-dir','submission_exports')
dir.create(base,recursive=TRUE,showWarnings=FALSE)
out <- file.path(base,paste0('VIH_',run_id))
if(file.exists(out)) stop('Export directory already exists: ',out, '; choose a new --run-id')
dir.create(out,recursive=TRUE)
dir.create(file.path(out,'tables'))
dir.create(file.path(out,'QA'))
logfile <- file.path(out,'QA','export.log')
log_line <- function(level,msg) { line <- paste(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),level,msg); cat(line,'\n'); cat(line,'\n',file=logfile,append=TRUE) }
issues <- data.frame(severity=character(),code=character(),detail=character())
issue <- function(severity,code,detail) { issues <<- rbind(issues,data.frame(severity=severity,code=code,detail=detail)); log_line(severity,paste(code,detail)) }
# Explicit source mapping: no automatic copying from old output folders.
required <- c('TABLE_empirical_shapley_vs_lmdi.csv','TABLE_empirical_multiyear_lmdi_shapley.csv',
              'TABLE_simulation_performance_overall.csv','TABLE0_model_quality.csv',
              'TABLE1_method_comparison_with_CI.csv','TABLE2_order_sensitivity.csv',
              'TABLE3_simulation_recovery.csv','TABLE4_axiom_validation.csv',
              'TABLE_model_vs_observed_change.csv','TABLE_method_decision_guide.csv')
optional <- c('bootstrap_summary.csv','order_sensitivity_summary.csv','sequential_all_24_orders.csv',
              'TABLE_methods_feature_comparison.csv','TABLE_shapley_vs_lmdi_simulation.csv',
              'TABLE_random_vs_spatial_moran_long.csv')
provenance <- data.frame(source=character(),destination=character(),md5=character(),size_bytes=numeric(),modified=character())
copy_checked <- function(rel, dest_rel=rel, is_required=TRUE) {
 src <- file.path(tdir,rel); dst <- file.path(out,'tables',dest_rel)
 if(!file.exists(src)) {issue(if(is_required) 'ERROR' else 'WARNING','MISSING_TABLE',rel);return(invisible(FALSE))}
 if(file.info(src)$size<=0) {issue('ERROR','EMPTY_TABLE',rel);return(invisible(FALSE))}
 dir.create(dirname(dst),recursive=TRUE,showWarnings=FALSE)
 if(!file.copy(src,dst,overwrite=FALSE)) stop('Cannot copy ',src)
 dat <- tryCatch(read.csv(src,check.names=FALSE),error=function(e)e)
 if(inherits(dat,'error') || !is.data.frame(dat) || nrow(dat)==0) issue('ERROR','INVALID_CSV',rel)
 provenance <<- rbind(provenance,data.frame(source=src,destination=dst,md5=unname(tools::md5sum(dst)),size_bytes=file.info(dst)$size,modified=as.character(file.info(src)$mtime)))
 TRUE
}
log_line('INFO',paste('Started; project',root,'; target',out))
for(f in required) copy_checked(f)
for(f in optional) copy_checked(f,is_required=FALSE)
# Materialize correct 4-row comparison from method_comparison.csv without refitting models.
# Prefer the authoritative raw method comparison, not the previously broken derived table.
source_cmp <- file.path(tdir, 'method_comparison.csv')
if (file.exists(source_cmp)) {
  d <- tryCatch(read.csv(source_cmp, check.names=FALSE), error=function(e) NULL)
  comps <- c('N','PI','LAMBDA','KAPPA')
  if (!is.null(d) && all(c('method','baseline_year','comparison_year',paste0(comps,'_effect')) %in% names(d))) {
    d <- d[d$method %in% c('naive','shapley'),,drop=FALSE]
    valid <- nrow(d)==2 && setequal(d$method,c('naive','shapley')) &&
      length(unique(d$baseline_year))==1 && length(unique(d$comparison_year))==1
    if (valid) {
      raw_a <- d[d$method=='naive',,drop=FALSE]
      raw_b <- d[d$method=='shapley',,drop=FALSE]
      x <- data.frame(baseline_year=raw_a$baseline_year,comparison_year=raw_a$comparison_year,
                      component=comps,lmdi_style=as.numeric(raw_a[1,paste0(comps,'_effect')]),
                      shapley=as.numeric(raw_b[1,paste0(comps,'_effect')]))
      x$difference_shapley_minus_lmdi <- x$shapley-x$lmdi_style
      x$abs_difference <- abs(x$difference_shapley_minus_lmdi)
      if (all(is.finite(x$lmdi_style)) && all(is.finite(x$shapley))) {
        write.csv(x,file.path(out,'tables','TABLE_empirical_shapley_vs_lmdi.csv'),row.names=FALSE)
        log_line('INFO','Created paired empirical comparison from method_comparison.csv (no model refit)')
      } else issue('ERROR','COMPARISON_NONFINITE','Comparison source contains non-finite estimates')
    } else issue('ERROR','COMPARISON_SOURCE_AMBIGUOUS','Expected one naive and one shapley row for one year pair')
  } else issue('ERROR','COMPARISON_SOURCE_SCHEMA','method_comparison.csv absent or unexpected schema')
} else issue('ERROR','MISSING_COMPARISON_SOURCE','Missing method_comparison.csv; cannot guarantee 4-row comparison')
paired_path <- file.path(out,'tables','TABLE_empirical_shapley_vs_lmdi.csv')
if (file.exists(paired_path)) {
 d <- tryCatch(read.csv(paired_path), error=function(e)NULL)
 if (is.null(d) || nrow(d)!=4 || !all(c('component','lmdi_style','shapley') %in% names(d)) ||
     anyNA(d$lmdi_style) || anyNA(d$shapley)) issue('ERROR','UNPAIRED_MAIN_COMPARISON','Submission comparison does not contain 4 complete pairs')
}
# Newer R13 estimates may be different: report numeric discrepancies, do not silently overwrite main analysis.
r13 <- file.path(tdir,'R13_family_selection','06_empirical_lmdi_vs_shapley_updated.csv')
if(file.exists(r13)) {
 a <- tryCatch(read.csv(file.path(out,'tables','TABLE_empirical_shapley_vs_lmdi.csv')),error=function(e)NULL)
 b <- tryCatch(read.csv(r13),error=function(e)NULL)
 if(!is.null(a)&&!is.null(b)&&all(c('component','shapley_bln')%in%names(b))) {
  chk <- merge(a,b[,c('component','shapley_bln')],by='component')
  chk$absolute_difference_points <- abs(chk$shapley-chk$shapley_bln*1e9)
  write.csv(chk,file.path(out,'QA','R13_vs_main_comparison.csv'),row.names=FALSE)
  maxdiff <- max(chk$absolute_difference_points,na.rm=TRUE)
  if(is.finite(maxdiff)&&maxdiff>1) issue('WARNING','R13_DIFF',paste('R13 and main differ: max absolute Shapley difference',format(maxdiff,scientific=FALSE),'points'))
 }
}
# Compare reported simulation 1.81 million with sources: descriptive, no forced correction.
sim <- file.path(tdir,'TABLE_simulation_performance_overall.csv')
if(file.exists(sim)) {
 d <- tryCatch(read.csv(sim),error=function(e)NULL)
 if(!is.null(d) && all(c('scenario','method','rmse')%in%names(d))) {
  x <- d[d$scenario=='cell_composition_shift' & d$method=='shapley',,drop=FALSE]
  if(nrow(x)==1) log_line('INFO',paste('Cell composition shift Shapley RMSE =',format(x$rmse,scientific=FALSE),'points; million =',round(x$rmse/1e6,2)))
 } else issue('WARNING','SIM_SCHEMA','Cannot check simulation RMSE automatically')
}
# Figures: existing exporter is used in a clean isolated project snapshot.
if(!tables_only) {
 if(!requireNamespace('data.table',quietly=TRUE)||!requireNamespace('ggplot2',quietly=TRUE)||!requireNamespace('scales',quietly=TRUE)) {
  issue('ERROR','MISSING_R_PACKAGES','Need data.table, ggplot2, scales to render submission figures')
 } else {
  stage <- file.path(out,'QA','build_workspace');dir.create(stage,recursive=TRUE)
  file.copy('config.R',file.path(stage,'config.R'))
  file.copy('R',stage,recursive=TRUE)
  file.copy('scripts',stage,recursive=TRUE)
  dir.create(file.path(stage,'output'),recursive=TRUE)
  # Copy only selected complete source tables for figure exporter; no stale figures carried over.
  file.copy(tdir,file.path(stage,'output'),recursive=TRUE)
  # Note: cfg$output$tables_dir typically output/tables; retain relative paths from config in stage.
  staged_tdir <- file.path(stage,tdir)
  if(!dir.exists(staged_tdir)) dir.create(staged_tdir,recursive=TRUE)
  source_files <- list.files(tdir,recursive=TRUE,full.names=TRUE)
  source_files <- source_files[!dir.exists(source_files)]
  for(src in source_files){rel<-substring(src,nchar(tdir)+2);dst<-file.path(staged_tdir,rel);dir.create(dirname(dst),recursive=TRUE,showWarnings=FALSE);file.copy(src,dst,overwrite=TRUE)}
  old <- getwd(); setwd(stage)
  exit_code <- tryCatch({
    cmd <- file.path(R.home('bin'),'Rscript')
    system2(cmd,c('scripts/run_98_export_vih_figures.R','--tiff','--no-scenario-specific'),stdout=file.path(old,out,'QA','figures_stdout.log'),stderr=file.path(old,out,'QA','figures_stderr.log'))
  },error=function(e){issue('ERROR','FIGURE_EXPORT_EXCEPTION',conditionMessage(e));1},finally={setwd(old)})
  if(!identical(as.integer(exit_code),0L)) issue('ERROR','FIGURE_EXPORT_FAILED',paste('exit status',exit_code))
  produced <- file.path(stage,'submission','figures')
  if(dir.exists(produced)) {
   dir.create(file.path(out,'figures'),showWarnings=FALSE)
   files <- list.files(produced,full.names=TRUE,recursive=TRUE)
   for(src in files) {if(dir.exists(src)) next;rel<-substring(src,nchar(produced)+2);dst<-file.path(out,'figures',rel);dir.create(dirname(dst),recursive=TRUE,showWarnings=FALSE);file.copy(src,dst)}
  }
  if(!file.exists(file.path(out,'figures','FIGURE_MANIFEST.csv'))) issue('ERROR','NO_FIGURE_MANIFEST','Figure exporter did not produce a manifest')
  figfiles <- list.files(file.path(out,'figures'),full.names=TRUE,recursive=TRUE,pattern='\\.(tiff|tif|png|jpg)$',ignore.case=TRUE)
  for(f in figfiles) if(file.info(f)$size==0) issue('ERROR','EMPTY_FIGURE',basename(f))
  log_line('INFO',paste('Exported image files:',length(figfiles)))
  # Only remove temporary workspace if complete figure creation; keep on failure for debugging.
  if(identical(as.integer(exit_code),0L)) unlink(stage,recursive=TRUE)
 }
}
# Curated set: main Figure 1-3 and supplementary Figure S1-S3 only.
# No S6 or Moran figure is expected in the submission package.
if (!tables_only) {
  fdir <- file.path(out, 'figures')
  figs <- list.files(fdir, pattern='\\.(tiff|tif)$', recursive=TRUE, full.names=TRUE, ignore.case=TRUE)
  figs <- figs[!grepl('scenario_specific', figs, fixed=TRUE)]
  expected <- c(sprintf('Figure_%02d_', 1:3), sprintf('Supplementary_Figure_%02d_', 1:3))
  for (prefix in expected) {
    n <- sum(startsWith(basename(figs), prefix))
    if (n != 1L) issue('ERROR', 'FIGURE_SET', paste('Expected one figure with prefix', prefix, 'got', n))
  }
  if (length(figs) != 6L) issue('ERROR', 'FIGURE_COUNT', paste('Expected six TIFFs (3 main + 3 supplement); found', length(figs)))
}
# provenance metadata for outputs including exported figures
result_files <- list.files(out,recursive=TRUE,full.names=TRUE)
result_files <- result_files[!dir.exists(result_files)]
result_files <- result_files[!grepl('build_workspace',result_files,fixed=TRUE)]
all_manifest <- data.frame(file=substring(result_files,nchar(out)+2),size_bytes=file.info(result_files)$size,md5=unname(tools::md5sum(result_files)))
write.csv(provenance,file.path(out,'QA','source_tables_manifest.csv'),row.names=FALSE)
write.csv(all_manifest,file.path(out,'QA','files_manifest.csv'),row.names=FALSE)
write.csv(issues,file.path(out,'QA','issues.csv'),row.names=FALSE)
writeLines(c(paste('Timestamp:',as.character(Sys.time())),paste('Project root:',root),paste('Config md5:',unname(tools::md5sum('config.R'))),paste('Tables source:',tdir),paste('Figure exporter:', 'scripts/run_98_export_vih_figures.R'),paste('Figure exporter md5:',unname(tools::md5sum('scripts/run_98_export_vih_figures.R'))),paste('R version:',R.version.string)),file.path(out,'QA','provenance.txt'))
log_line('INFO',paste('Finished: errors',sum(issues$severity=='ERROR'),'warnings',sum(issues$severity=='WARNING')))
cat('\nSUBMISSION EXPORT: ',normalizePath(out,winslash='/'),'\n',sep='')
cat('QA issues: ',file.path(out,'QA','issues.csv'),'\n',sep='')
if(strict && any(issues$severity=='ERROR')) quit(status=2)
