#!/usr/bin/env Rscript
# Patient NB GAM: RE vs RE+AGE vs RE+SLOPE vs RE+AGE+SLOPE.
# Original cell-level Pearson residuals, blocked spatial Moran, kNN 4/8/12.
# No bootstrap, simulation or Shapley.
suppressPackageStartupMessages({library(data.table);library(mgcv);library(spdep)})
source('R/99_load_all.R')
cfg <- load_config('config.R')
out <- 'output/pilot_patients_random_slopes'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
ref <- 'output/pilot_mrf_patients/patient_models_RE_MRF.rds'
if (!file.exists(ref)) stop('Missing RE reference RDS: ', ref)
old <- readRDS(ref)
if (is.null(old$RE)) stop('Reference RDS lacks RE model')
if (file.exists('output/rds/aggregated_data.rds')) {
 d <- as.data.table(readRDS('output/rds/aggregated_data.rds'))
} else {
 raw <- read_input_data(cfg$data$path)
 d <- as.data.table(filter_and_aggregate(normalize_schema(raw,value_var=cfg$data$value_var),disease_filter=cfg$data$disease_filter))
}
req <- c('code','x','y','year','sex','age_group','pop','n_patients')
if(length(setdiff(req,names(d)))) stop('Missing data columns: ',paste(setdiff(req,names(d)),collapse=', '))
d <- d[is.finite(x)&is.finite(y)&is.finite(pop)&is.finite(n_patients)&pop>0&n_patients>=0]
d[,`:=`(code=as.character(code),sex=factor(sex),age_group=factor(age_group),year_num=as.numeric(year))]
coords <- d[,.(x=median(x),y=median(y),ncoords=uniqueN(paste(x,y))),by=code][order(code)]
if(any(coords$ncoords!=1L)||anyDuplicated(coords[,.(x,y)])) stop('Invalid / duplicate coordinates')
ids <- coords$code
d[,code:=factor(code,levels=ids)]
d[,sex_age:=interaction(sex,age_group,drop=TRUE)]
d[,code_age:=interaction(code,age_group,drop=TRUE,lex.order=TRUE)]
d[,year_c:=as.numeric(year_num - mean(range(year_num)))]
keycols <- c('year','sex','age_group','code')
if(anyDuplicated(d[,..keycols])) stop('Duplicate original cells; no aggregation permitted')
eps <- as.numeric(cfg$model$eps_offset)
if(length(eps)!=1L||!is.finite(eps)||eps<=0) stop('Invalid eps_offset')
d[,log_pop_offset:=base::log(pmax(as.numeric(pop),eps))]
if(any(!is.finite(d$log_pop_offset))) stop('Invalid population offset')
compatible <- function(m){
 inherits(m,'gam') && nobs(m)==nrow(d) && length(m$y)==nrow(d) &&
 isTRUE(all.equal(as.numeric(m$y),as.numeric(d$n_patients),tolerance=0,check.attributes=FALSE))
}
if(!compatible(old$RE)) stop('Stored RE incompatible with input row ordering / outcomes')
# Reuse AGE only if it has the correct formula and matches data.
age_paths <- c('output/pilot_patients_age_vs_none/pilot_new_models.rds',
 'output/pilot_patients_code_age_re/patients_RE_code_age.rds')
age <- NULL
for (p in age_paths) if(file.exists(p)) {
 obj <- readRDS(p)
 candidate <- if(is.list(obj) && !inherits(obj,'gam')) obj$RE_code_age else obj
 if(!is.null(candidate) && compatible(candidate) &&
    grepl('code_age',paste(deparse(formula(candidate)),collapse=' '),fixed=TRUE)) {
  age <- candidate; break
 }
}
if(is.null(age)) stop('Missing compatible RE+AGE fit. First run scripts/run_pilot_patients_age_vs_none.R')
# Preserve baseline model terms (and settings) and add random slopes, no altered fixed effects.
common <- paste0('n_patients ~ sex * age_group + s(year_num, by=sex_age, k=',
 min(cfg$model$k_year,uniqueN(d$year)),')')
random_slope <- 's(code, by=year_c, bs="re")'
fo_slope <- as.formula(paste(common,'+ s(code,bs="re") +',random_slope,'+ offset(log_pop_offset)'))
fo_age_slope <- as.formula(paste(common,'+ s(code,bs="re") + s(code_age,bs="re") +',random_slope,'+ offset(log_pop_offset)'))
cat('Patients: RE vs RE+AGE vs RE+SLOPE vs RE+AGE+SLOPE (cell-level).\n')
cat('Fitting only two new GAMs. Progress bar updates after completed stages.\n')
pb <- utils::txtProgressBar(min=0,max=8,style=3,width=35)
fit_one <- function(f) mgcv::bam(f,data=d,family=mgcv::nb(),method='fREML',discrete=TRUE,
 control=mgcv::gam.control(maxit=cfg$model$maxit_gam,trace=FALSE))
fit_slope <- tryCatch(fit_one(fo_slope),error=function(e){close(pb);stop(e)})
utils::setTxtProgressBar(pb,1)
fit_age_slope <- tryCatch(fit_one(fo_age_slope),error=function(e){close(pb);stop(e)})
utils::setTxtProgressBar(pb,2)
models <- list(RE=old$RE,RE_AGE=age,RE_SLOPE=fit_slope,RE_AGE_SLOPE=fit_age_slope)
if(!all(vapply(models,compatible,logical(1)))) stop('Comparison model data mismatch')
metrics <- rbindlist(lapply(names(models),function(nm){m<-models[[nm]]
 data.table(model=nm,n=nobs(m),aic=as.numeric(AIC(m)),deviance_explained=summary(m)$dev.expl,
 edf_total=sum(m$edf),theta=m$family$getTheta(TRUE),total_predicted=sum(as.numeric(predict(m,newdata=d,type='response'))))}))
fwrite(metrics,file.path(out,'TABLE_model_comparison.csv'))
by_year <- rbindlist(lapply(names(models),function(nm){
 data.table(year=d$year,obs=d$n_patients,pred=as.numeric(predict(models[[nm]],newdata=d,type='response')))[,
 .(observed=sum(obs),predicted=sum(pred)),by=year][,model:=nm]}))
fwrite(by_year,file.path(out,'TABLE_prediction_by_year.csv'))
utils::setTxtProgressBar(pb,3)
d[,block:=interaction(year,sex,age_group,drop=TRUE,lex.order=TRUE)]
blocks <- split(seq_len(nrow(d)),d$block)
if(length(blocks)!=112L) warning('Expected 112 strata, found ',length(blocks))
resids <- lapply(models,function(m){r<-as.numeric(residuals(m,type='pearson'))
 if(length(r)!=nrow(d)||any(!is.finite(r))) stop('Bad residuals');r})
make_graph <- function(k){
 nb <- spdep::make.sym.nb(spdep::knn2nb(spdep::knearneigh(as.matrix(coords[,.(x,y)]),k=k)))
 lapply(blocks,function(ix){
  idx <- match(as.character(d$code[ix]),ids)
  lookup <- integer(length(ids));lookup[idx]<-seq_along(idx)
  ii <- integer();jj <- integer();ww <- numeric()
  for(u in seq_along(idx)){
   near <- lookup[nb[[idx[u]]]];near<-near[near>0L]
   if(length(near)){ii<-c(ii,rep.int(u,length(near)));jj<-c(jj,near);ww<-c(ww,rep.int(1/length(near),length(near)))}
  }
  list(ix=ix,i=ii,j=jj,w=ww,n=length(ix),s0=sum(ww))
 })
}
calc_I <- function(r,gr,slices=FALSE){
 n<-0;den<-0;num<-0;s0<-0
 ans<-if(slices)vector('list',length(gr))else NULL
 for(b in seq_along(gr)){
  g<-gr[[b]];z<-r[g$ix]-mean(r[g$ix]);d2<-sum(z*z)
  nt<-sum(g$w*z[g$i]*z[g$j]);num<-num+nt;den<-den+d2;n<-n+g$n;s0<-s0+g$s0
  if(slices)ans[[b]]<-data.table(block=names(blocks)[b],n_cells=g$n,
   moran_i=if(d2>0&&g$s0>0)(g$n/g$s0)*nt/d2 else NA_real_)
 }
 list(I=if(den>0&&s0>0)(n/s0)*num/den else NA_real_,slices=if(slices)rbindlist(ans)else NULL)
}
set.seed(20261019)
perm <- as.integer(Sys.getenv('VIH_MORAN_PERMUTATIONS','199'))
if(is.na(perm)||perm<99L) stop('Invalid VIH_MORAN_PERMUTATIONS: >=99 required')
global <- list();all_slices<-list()
labels <- d[,.(year=year[1],sex=as.character(sex[1]),age_group=as.character(age_group[1])),by=block]
for(k in c(4L,8L,12L)){
 gr<-make_graph(k)
 for(nm in names(resids)){
  r<-resids[[nm]];obs<-calc_I(r,gr,TRUE)
  st<-merge(copy(obs$slices),labels,by='block',sort=FALSE);st[,`:=`(model=nm,k=k)]
  all_slices[[length(all_slices)+1L]]<-st
  null<-numeric(perm)
  for(p in seq_len(perm)){
   rp<-r
   for(g in gr)rp[g$ix]<-sample(r[g$ix],length(g$ix),replace=FALSE)
   null[p]<-calc_I(rp,gr,FALSE)$I
  }
  mu<-mean(null)
  global[[length(global)+1L]]<-data.table(model=nm,k=k,n_cells=nrow(d),n_blocks=length(blocks),
   moran_global_block=obs$I,median_slice_i=median(st$moran_i,na.rm=TRUE),
   q25_slice_i=as.numeric(quantile(st$moran_i,.25,na.rm=TRUE)),
   q75_slice_i=as.numeric(quantile(st$moran_i,.75,na.rm=TRUE)),
   p_within_block=(1+sum(abs(null-mu)>=abs(obs$I-mu)))/(perm+1),n_permutations=perm)
 }
 utils::setTxtProgressBar(pb,3+match(k,c(4L,8L,12L)))
}
ss<-rbindlist(all_slices);gg<-rbindlist(global)
fwrite(gg,file.path(out,'TABLE_global_block_moran.csv'))
fwrite(ss,file.path(out,'TABLE_moran_112_cell_slices.csv'))
fwrite(ss[,.(n_slices=.N,median_i=median(moran_i),q25_i=as.numeric(quantile(moran_i,.25)),
 q75_i=as.numeric(quantile(moran_i,.75))),by=.(model,k,age_group)],file.path(out,'TABLE_moran_by_age.csv'))
wide<-dcast(ss,year+sex+age_group+k~model,value.var='moran_i')
wide[,`:=`(slope_minus_RE=RE_SLOPE-RE,age_slope_minus_age=RE_AGE_SLOPE-RE_AGE,
 age_slope_minus_RE=RE_AGE_SLOPE-RE)]
fwrite(wide,file.path(out,'TABLE_paired_slices.csv'))
fwrite(data.table(metric=c('model_cells','unique_cell_keys','strata','codes','age_groups','years','sexes','year_center'),
 value=c(nrow(d),uniqueN(d,by=keycols),length(blocks),uniqueN(d$code),uniqueN(d$age_group),uniqueN(d$year),uniqueN(d$sex),mean(range(d$year_num)))),
 file.path(out,'TABLE_input_integrity.csv'))
# Add same-cell accuracy; in-sample only, not external validation.
calib<-rbindlist(lapply(names(models),function(nm){m<-models[[nm]];p<-as.numeric(predict(m,newdata=d,type='response'))
 data.table(model=nm,mae=mean(abs(d$n_patients-p)),rmse=sqrt(mean((d$n_patients-p)^2)),
 mean_observed=mean(d$n_patients),mean_predicted=mean(p))}))
fwrite(calib,file.path(out,'TABLE_in_sample_accuracy.csv'))
writeLines(c('One model observation = code x year x sex x age_group; no residual aggregation.',
 'Same input cells, NB GAM family, demographic terms, offsets and global sex-age year smooth across all fits.',
 'RE and RE+AGE are reused if compatible; only two slope models are fitted.',
 'year_c = year_num minus midpoint of years; slopes are code-specific LINEAR deviations in year.',
 'RE+SLOPE adds s(code,by=year_c,bs=re); RE+AGE+SLOPE adds same to RE+AGE.',
 'Moran is Pearson residuals centered within each year x sex x age block, on symmetric kNN-4/8/12 (proxy, not borders).',
 'Global block Moran is variance-weighted; permutation of residuals is within blocks, exploratory conditional on fits.',
 'The provided AIC and in-sample accuracy are NOT out-of-sample validation.',
 'No Shapley, bootstrap or simulation calculations are included.'),file.path(out,'README_method.txt'))
writeLines(capture.output(summary(fit_slope)),file.path(out,'diagnostics_RE_SLOPE.txt'))
writeLines(capture.output(summary(fit_age_slope)),file.path(out,'diagnostics_RE_AGE_SLOPE.txt'))
saveRDS(list(RE_SLOPE=fit_slope,RE_AGE_SLOPE=fit_age_slope),file.path(out,'models_random_slopes.rds'))
utils::setTxtProgressBar(pb,8)
close(pb)
cat('\nResults saved: ',out,'\n',sep='')
print(metrics)
print(gg[,.(model,k,moran_global_block,median_slice_i,p_within_block)])
