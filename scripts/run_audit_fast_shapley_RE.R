#!/usr/bin/env Rscript
# Validate each of 16 cached Fast-Shapley RE states against original NFZ helpers.
# No refit, no bootstrap. Run in project root.
suppressPackageStartupMessages({library(data.table);library(dplyr);library(mgcv)})
source('R/99_load_all.R')
source('R/03_predict_helpers.R')
cfg <- load_config('config.R')
OUT <- 'output/pilot_fast_shapley_sensitivity/audit_RE_16'
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
fast_path <- 'output/pilot_fast_shapley_sensitivity/cache_scenarios_RE.rds'
chain_path <- 'output/rds/model_fit.rds'
data_path <- 'output/rds/aggregated_data.rds'
for(p in c(fast_path,chain_path,data_path))if(!file.exists(p))stop('Missing: ',p)
fast <- readRDS(fast_path)
if(!is.numeric(fast)||is.null(names(fast)))stop('Invalid named numeric fast cache')
fit <- readRDS(chain_path)
if(!all(c('m_pat','m_prov','m_val')%in%names(fit)))stop('Invalid original model chain')
d <- as.data.frame(readRDS(data_path))
d$year <- as.integer(as.character(d$year))
base <- as.integer(cfg$data$baseline_year)
last <- max(d$year)
if(base==last)stop('Baseline and target coincide')
d2 <- d[d$year %in% c(base,last),,drop=FALSE]
if(anyDuplicated(d2[,c('code','year','sex','age_group')]))stop('Duplicate analytical cells')
p <- prep_pair(d2,base,last)
comps <- c('N','PI','LAMBDA','KAPPA')
subsets <- unlist(lapply(0:4,function(k) combn(comps,k,simplify=FALSE)),recursive=FALSE)
key <- function(v)if(!length(v))'EMPTY' else paste(comps[comps%in%v],collapse='_')
expected <- vapply(subsets,key,character(1))
if(!setequal(names(fast),expected)||length(fast)!=16L)stop('Expected exactly 16 fast keys; found: ',paste(names(fast),collapse=','))
# Original fit is intentionally used, not the RE pilot model.
# This lets the audit detect mismatched RE model objects as well as differences in scenario construction.
message('Evaluating original 16 states (no fitting).')
bar <- utils::txtProgressBar(min=0,max=16,style=3,width=35)
rows <- vector('list',16)
for(i in seq_along(subsets)){
  active <- subsets[[i]]
  k <- key(active)
  reference <- scenario_from_set(p$target,p$base,fit,base,last,active)
  if(length(reference)!=1L||!is.finite(reference))stop('Non-finite original state ',k)
  candidate <- unname(fast[[k]])
  rows[[i]] <- data.table(state=k,active_components=if(length(active))paste(active,collapse='+') else '(none)',
                         original=reference,fast=candidate,fast_minus_original=candidate-reference,
                         abs_diff=abs(candidate-reference),rel_diff=abs(candidate-reference)/max(abs(reference),1))
  utils::setTxtProgressBar(bar,i)
}
close(bar)
tab <- rbindlist(rows)
fwrite(tab,file.path(OUT,'TABLE_16_scenario_RE_audit.csv'))
get <- function(values,which)values[match(which,tab$state)]
refv <- tab$original; fastv <- tab$fast
component <- function(v){
  perms <- combinat_permutations(comps)
  contrib <- matrix(0,length(perms),length(comps),dimnames=list(NULL,comps))
  for(j in seq_along(perms)){
    act <- character(0);prev <- get(v,key(act))
    for(c in perms[[j]]){
      act <- c(act,c);nxt <- get(v,key(act));contrib[j,c]<-nxt-prev;prev<-nxt
    }
  }
  colMeans(contrib)
}
oref <- component(refv);ofast <- component(fastv)
endpoints <- data.table(source=c('ORIGINAL_03_predict_helpers','FAST_CACHE'),
                        value_2018=c(get(refv,'EMPTY'),get(fastv,'EMPTY')),
                        value_2024=c(get(refv,key(comps)),get(fastv,key(comps))),
                        change=c(get(refv,key(comps))-get(refv,'EMPTY'),get(fastv,key(comps))-get(fastv,'EMPTY')))
fwrite(endpoints,file.path(OUT,'TABLE_endpoints_RE_audit.csv'))
contrib <- data.table(component=comps,original=as.numeric(oref[comps]),fast=as.numeric(ofast[comps]))
contrib[,diff_fast_minus_original:=fast-original]
fwrite(contrib,file.path(OUT,'TABLE_shapley_RE_audit.csv'))
# Strict/loose status thresholds in settlement POINTS, not model p-values.
max_abs <- max(tab$abs_diff)
status <- if(max_abs <= 1) 'MATCH_WITHIN_1_POINT' else if(max_abs <= 1e5) 'SMALL_DIFFERENCE_CHECK_MODEL_OBJECT' else 'MISMATCH_INVESTIGATE'
meta <- data.table(status=status, max_abs_point_difference=max_abs,
                   mean_abs_point_difference=mean(tab$abs_diff),
                   n_states=16L,baseline=base,target=last,
                   original_patients_AIC=AIC(fit$m_pat),
                   original_patients_n=nobs(fit$m_pat))
# Optional: compare RE pilot model AIC and coefficients as diagnostic for source mismatch.
pilot_path <- 'output/pilot_mrf_patients/patient_models_RE_MRF.rds'
if(file.exists(pilot_path)){
  tmp <- readRDS(pilot_path)
  if(!is.null(tmp$RE)){
    meta[,pilot_RE_AIC:=AIC(tmp$RE)]
    meta[,original_vs_pilot_AIC:=AIC(fit$m_pat)-AIC(tmp$RE)]
  }
}
fwrite(meta,file.path(OUT,'TABLE_audit_status.csv'))
cat('\nAUDIT:',status,'; maximum absolute state difference:',format(max_abs,digits=10),'points\n')
cat('Files:',OUT,'\n')
print(endpoints)
print(contrib)
if(status!='MATCH_WITHIN_1_POINT')message('Important: original model or scenario computations differ. Do not assert equivalence from additive closure alone.')
