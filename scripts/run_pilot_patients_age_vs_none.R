#!/usr/bin/env Rscript
# One-pipeline comparison of no-region, RE and RE + (region x age) NB GAMs.
# Uses original model cells and saved reference RE model. No Shapley or bootstrap.
suppressPackageStartupMessages({library(data.table);library(mgcv);library(spdep)})
source('R/99_load_all.R')
cfg <- load_config('config.R')
out <- 'output/pilot_patients_age_vs_none'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
ref <- 'output/pilot_mrf_patients/patient_models_RE_MRF.rds'
if (!file.exists(ref)) stop('Missing RE reference: ',ref)
old <- readRDS(ref)
if (is.null(old$RE)) stop('Reference RDS lacks RE model')
if (file.exists('output/rds/aggregated_data.rds')) {
 d <- as.data.table(readRDS('output/rds/aggregated_data.rds'))
} else {
 raw <- read_input_data(cfg$data$path)
 d <- as.data.table(filter_and_aggregate(normalize_schema(raw,value_var=cfg$data$value_var), disease_filter=cfg$data$disease_filter))
}
needed<-c('code','x','y','year','sex','age_group','pop','n_patients')
if(length(setdiff(needed,names(d)))) stop('Missing columns: ',paste(setdiff(needed,names(d)),collapse=', '))
d<-d[is.finite(x)&is.finite(y)&is.finite(pop)&is.finite(n_patients)&pop>0&n_patients>=0]
d[,`:=`(code=as.character(code),sex=factor(sex),age_group=factor(age_group),year_num=as.numeric(year))]
coords<-d[,.(x=median(x),y=median(y),ncoords=uniqueN(paste(x,y))),by=code][order(code)]
if(any(coords$ncoords!=1L)||anyDuplicated(coords[,.(x,y)])) stop('Invalid/duplicated code coordinates')
ids<-coords$code
d[,code:=factor(code,levels=ids)]
d[,sex_age:=interaction(sex,age_group,drop=TRUE)]
keycols<-c('year','sex','age_group','code')
if(anyDuplicated(d[,..keycols])) stop('Duplicate model cells: no aggregation is allowed')
eps<-as.numeric(cfg$model$eps_offset)
if(length(eps)!=1L||!is.finite(eps)||eps<=0) stop('Invalid eps_offset')
d[,log_pop_offset:=base::log(pmax(as.numeric(pop),eps))]
if(any(!is.finite(d$log_pop_offset))) stop('Non-finite offset')
if(nobs(old$RE)!=nrow(d) || length(old$RE$y)!=nrow(d) ||
 !isTRUE(all.equal(as.numeric(old$RE$y),as.numeric(d$n_patients),tolerance=0,check.attributes=FALSE)))
 stop('Saved RE is incompatible with model data/row order')
# No spatial effects means no administrative-unit term, not removal of geographic records.
common<-paste0('n_patients ~ sex * age_group + s(year_num, by=sex_age, k=',
 min(cfg$model$k_year,uniqueN(d$year)),') + offset(log_pop_offset)')
f_none<-as.formula(common)
f_age<-as.formula(sub(' + offset(log_pop_offset)',' + s(code, bs="re") + s(code_age, bs="re") + offset(log_pop_offset)',common,fixed=TRUE))
d[,code_age:=interaction(code,age_group,drop=TRUE,lex.order=TRUE)]
# Saved pilot can be reused. No uncontrolled change to baseline specification.
prior_age <- 'output/pilot_patients_code_age_re/patients_RE_code_age.rds'
if(file.exists(prior_age)) {
 age_existing<-readRDS(prior_age)
 if(nobs(age_existing)!=nrow(d)||length(age_existing$y)!=nrow(d)||
  !isTRUE(all.equal(as.numeric(age_existing$y),as.numeric(d$n_patients),tolerance=0,check.attributes=FALSE)))
   stop('Saved code-age model incompatible with input data')
} else age_existing<-NULL
cat('Patients: NO spatial effects vs RE vs RE+code_age (cell-level analysis).\n')
cat('Progress updates per fitting/diagnostic stage. GAM fitting may take considerable time.\n')
pb<-utils::txtProgressBar(min=0,max=6,style=3,width=35)
fail<-function(e){try(close(pb),silent=TRUE);stop(e)}
fit_none<-tryCatch(mgcv::bam(f_none,data=d,family=mgcv::nb(),method='fREML',discrete=TRUE,
 control=mgcv::gam.control(maxit=cfg$model$maxit_gam,trace=FALSE)),error=fail)
utils::setTxtProgressBar(pb,1)
fit_age<-if(!is.null(age_existing)) age_existing else tryCatch(mgcv::bam(f_age,data=d,family=mgcv::nb(),method='fREML',discrete=TRUE,
 control=mgcv::gam.control(maxit=cfg$model$maxit_gam,trace=FALSE)),error=fail)
utils::setTxtProgressBar(pb,2)
models<-list(NONE=fit_none,RE=old$RE,RE_code_age=fit_age)
metrics<-rbindlist(lapply(names(models),function(nm){m<-models[[nm]]
 data.table(model=nm,n=nobs(m),aic=as.numeric(AIC(m)),deviance_explained=summary(m)$dev.expl,
 edf_total=sum(m$edf),theta=m$family$getTheta(TRUE),
 total_predicted=sum(as.numeric(predict(m,newdata=d,type='response'))))}))
fwrite(metrics,file.path(out,'TABLE_model_comparison.csv'))
pred_year<-rbindlist(lapply(names(models),function(nm){
 data.table(year=d$year,obs=d$n_patients,pred=as.numeric(predict(models[[nm]],newdata=d,type='response')))[,
 .(observed=sum(obs),predicted=sum(pred)),by=year][,model:=nm]}))
fwrite(pred_year,file.path(out,'TABLE_prediction_by_year.csv'))
utils::setTxtProgressBar(pb,3)
# Build a block-diagonal graph over UNIQUE original cell observations. Never average residuals.
d[,block:=interaction(year,sex,age_group,drop=TRUE,lex.order=TRUE)]
blocks<-split(seq_len(nrow(d)),d$block)
if(length(blocks)!=112L) warning('Expected 112 year/sex/age strata; got ',length(blocks))
resids<-lapply(models,function(m) {r<-as.numeric(residuals(m,type='pearson'));
 if(length(r)!=nrow(d)||any(!is.finite(r)))stop('Residual mismatch or nonfinite residuals')
 r})
make_graph<-function(k){
 nb<-spdep::make.sym.nb(spdep::knn2nb(spdep::knearneigh(as.matrix(coords[,.(x,y)]),k=k)))
 lapply(blocks,function(ix){
  idx<-match(as.character(d$code[ix]),ids)
  lookup<-integer(length(ids));lookup[idx]<-seq_along(idx)
  ii<-integer();jj<-integer();ww<-numeric()
  for(u in seq_along(idx)){
   near<-lookup[nb[[idx[u]]]];near<-near[near>0]
   if(length(near)){
    ii<-c(ii,rep.int(u,length(near)));jj<-c(jj,near)
    ww<-c(ww,rep.int(1/length(near),length(near)))
   }
  }
  list(ix=ix,i=ii,j=jj,w=ww,n=length(ix),s0=sum(ww))
 })
}
calc_I<-function(r,gr,slices=FALSE){
 n<-0;den<-0;num<-0;s0<-0
 ans<-if(slices)vector('list',length(gr))else NULL
 for(b in seq_along(gr)){
  g<-gr[[b]];z<-r[g$ix]-mean(r[g$ix]);d2<-sum(z*z)
  nt<-sum(g$w*z[g$i]*z[g$j]);num<-num+nt;den<-den+d2;n<-n+g$n;s0<-s0+g$s0
  if(slices)ans[[b]]<-data.table(block=names(blocks)[b],n_cells=g$n,
    moran_i=if(d2>0&&g$s0>0)(g$n/g$s0)*nt/d2 else NA_real_)
 }
 list(I=if(den>0&&s0>0)(n/s0)*num/den else NA_real_,
 slices=if(slices)rbindlist(ans)else NULL)
}
# Permutation null keeps labels within year/sex/age; p-values descriptive after GAM estimation.
set.seed(20261019)
perm<-as.integer(Sys.getenv('VIH_MORAN_PERMUTATIONS','199'))
if(is.na(perm)||perm<99)stop('Invalid VIH_MORAN_PERMUTATIONS; expected >= 99')
global<-list();all_slices<-list()
for(k in c(4L,8L,12L)){
 gr<-make_graph(k)
 for(nm in names(resids)){
  r<-resids[[nm]]; obs<-calc_I(r,gr,TRUE)
  st<-copy(obs$slices)
  lab<-d[,.(year=year[1],sex=as.character(sex[1]),age_group=as.character(age_group[1])),by=block]
  st<-merge(st,lab,by='block',sort=FALSE);st[,`:=`(model=nm,k=k)]
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
wide[,`:=`(age_minus_re=RE_code_age-RE,none_minus_re=NONE-RE)]
fwrite(wide,file.path(out,'TABLE_paired_slices.csv'))
fwrite(data.table(metric=c('model_cells','unique_cell_keys','strata','codes','age_groups','years','sexes'),
 value=c(nrow(d),uniqueN(d,by=keycols),length(blocks),uniqueN(d$code),uniqueN(d$age_group),uniqueN(d$year),uniqueN(d$sex))),
 file.path(out,'TABLE_input_integrity.csv'))
writeLines(c('All fits are Negative Binomial GAM on identical input cells, offset and year/sex/age terms.',
 'NONE: no code or spatial term. RE: published s(code,bs=re). RE_code_age: RE plus s(code_age,bs=re).',
 'Same saved reference RE is used; age pilot is reused if compatible.',
 'Moran: unaggregated Pearson residuals, block by year x sex x age, centered inside blocks.',
 'Graphs symmetric centroid kNN 4/8/12, NOT geographic polygon borders.',
 'Global Moran is variance-weighted, not the simple mean of slice Moran values.',
 'Permutation p values conditional on fitted residuals; exploratory, not formal model comparison.',
 'Please examine AIC, absolute calibration, residuals, and predictive performance jointly.'),
 file.path(out,'README_method.txt'))
writeLines(capture.output(summary(fit_age)),file.path(out,'age_model_diagnostics.txt'))
saveRDS(list(NONE=fit_none,RE_code_age=fit_age),file.path(out,'pilot_new_models.rds'))
close(pb)
cat('\nResults saved in: ',out,'\n',sep='')
print(metrics)
print(gg[,.(model,k,moran_global_block,median_slice_i,p_within_block)])
