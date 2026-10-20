#!/usr/bin/env Rscript
# Patient NB GAM: inferential specification tests for RE+AGE+SLOPE baseline.
# Original cell-level Pearson residuals, blocked spatial Moran, kNN 4/8/12.
# No bootstrap, simulation or Shapley.
suppressPackageStartupMessages({library(data.table);library(mgcv);library(spdep)})
source('R/99_load_all.R')
cfg <- load_config('config.R')
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
# Explicitly keep the publication-compatible baseline.
ref_slope_path <- 'output/pilot_patients_random_slopes/models_random_slopes.rds'
if (!file.exists(ref_slope_path)) stop('Missing baseline RE+AGE+SLOPE model: ', ref_slope_path)
ref_slope <- readRDS(ref_slope_path)
base_fit <- ref_slope$RE_AGE_SLOPE
if (is.null(base_fit) || !compatible(base_fit)) stop('Stored baseline is incompatible with the current rows')

# k must be no larger than the 7 available distinct time points.
k_ref <- as.integer(min(cfg$model$k_year, uniqueN(d$year_num)))
if (k_ref != 5L) warning('Baseline k_year is ', k_ref, ', not 5; check reference consistency.')
ref_formula_text <- paste(deparse(formula(base_fit)), collapse=' ')
if (!grepl('code_age',ref_formula_text,fixed=TRUE) || !grepl('year_c',ref_formula_text,fixed=TRUE))
 stop('Cached baseline is not RE+AGE+SLOPE')

# Center population-size covariate *within* year x sex x age to isolate differences across units.
# Keep population offset unchanged. Linear and smooth versions are ALTERNATIVES.
d[,log_pop:=base::log(pop)]
d[,log_pop_c:=log_pop-mean(log_pop),by=.(year,sex,age_group)]
if(any(!is.finite(d$log_pop_c)))stop('Nonfinite centered population predictor')
d[,code_sex:=interaction(code,sex,drop=TRUE,lex.order=TRUE)]
d[,pandemic_c:=as.numeric(year_num %in% c(2020,2021)) - 2/uniqueN(year_num)]
# Slightly safer for unbalanced designs: per-unit centering of the indicator.
d[,pandemic_c:=pandemic_c-mean(pandemic_c),by=code]

core <- function(k) paste0(
 'n_patients ~ sex * age_group + s(year_num, by=sex_age, k=',k,')',
 ' + s(code,bs="re") + s(code_age,bs="re") + s(code,by=year_c,bs="re")')
form <- function(k=k_ref,sex=FALSE,pandemic=FALSE,pop='none') {
 extra <- character()
 if(sex) extra <- c(extra,'s(code_sex,bs="re")')
 if(pandemic) extra <- c(extra,'s(code,by=pandemic_c,bs="re")')
 if(pop=='linear') extra <- c(extra,'log_pop_c')
 if(pop=='smooth') extra <- c(extra,'s(log_pop_c,k=6)')
 as.formula(paste(core(k), if(length(extra))paste0(' + ',paste(extra,collapse=' + '))else'',
   '+ offset(log_pop_offset)'))
}
# Predeclared hypotheses; toggle individual variants using environment variable.
specs <- list(
 SEX=list(k=k_ref,sex=TRUE),
 PANDEMIC=list(k=k_ref,pandemic=TRUE),
 POP_LINEAR=list(k=k_ref,pop='linear'),
 POP_SMOOTH=list(k=k_ref,pop='smooth'),
 SEX_PANDEMIC=list(k=k_ref,sex=TRUE,pandemic=TRUE),
 ALL=list(k=k_ref,sex=TRUE,pandemic=TRUE,pop='smooth'),
 K6=list(k=6L),
 K7=list(k=7L),
 ALL_K7=list(k=7L,sex=TRUE,pandemic=TRUE,pop='smooth')
)
selected <- Sys.getenv('VIH_VARIANTS','')
if(nzchar(selected)) {
 selected <- trimws(strsplit(selected,',',fixed=TRUE)[[1]])
 if(length(setdiff(selected,names(specs))))stop('Unknown variants: ',paste(setdiff(selected,names(specs)),collapse=', '))
 specs <- specs[selected]
}
if(uniqueN(d$year_num)!=7L)stop('Expected exactly 7 distinct years, got ',uniqueN(d$year_num))
if(k_ref>6L)stop('Expected baseline k <=6 to compare k=6 and 7')

out <- 'output/pilot_patients_missing_terms'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Resume completed experimental fits when present, to avoid expensive refits.
cache_dir <- file.path(out,'model_cache');dir.create(cache_dir,showWarnings=FALSE)
cat('Patient inference GAM: additional sex/pandemic/population terms and year k 5/6/7.\n')
cat('Model fits can take time; progress advances after each fit or diagnostic stage.\n')
pb <- utils::txtProgressBar(min=0,max=length(specs)+3,style=3,width=35)
fit_one <- function(f) mgcv::bam(f,data=d,family=mgcv::nb(),method='fREML',discrete=TRUE,
 control=mgcv::gam.control(maxit=cfg$model$maxit_gam,trace=FALSE))
models <- list(BASE=base_fit)
fit_log <- list()
for(j in seq_along(specs)){
 nm <- names(specs)[j]
 fo <- do.call(form,specs[[j]])
 p <- file.path(cache_dir,paste0('model_',nm,'.rds'))
 t0 <- Sys.time()
 reused <- FALSE
 if(file.exists(p)){
  z <- tryCatch(readRDS(p),error=function(e)NULL)
  if(!is.null(z)&&compatible(z) &&
   identical(paste(deparse(formula(z)),collapse=' '),paste(deparse(fo),collapse=' '))) {
   m <- z;reused <- TRUE
  }
 }
 if(!reused){
  m <- tryCatch(fit_one(fo),error=function(e){close(pb);stop('Fit ',nm,' failed: ',conditionMessage(e))})
  saveRDS(m,p)
 }
 models[[nm]] <- m
 fit_log[[j]]<-data.table(model=nm,reused=reused,elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units='secs')),
 formula=paste(deparse(fo),collapse=' '))
 utils::setTxtProgressBar(pb,j)
}
fwrite(rbindlist(fit_log),file.path(out,'TABLE_fitting_log.csv'))

metrics <- rbindlist(lapply(names(models),function(nm){m<-models[[nm]]
 pred <- as.numeric(predict(m,newdata=d,type='response'))
 data.table(model=nm,n=nobs(m),aic=as.numeric(AIC(m)),delta_aic_vs_base=as.numeric(AIC(m)-AIC(base_fit)),
 deviance_explained=summary(m)$dev.expl,edf_total=sum(m$edf),theta=m$family$getTheta(TRUE),
 total_predicted=sum(pred),pred_pct_vs_base=100*(sum(pred)/sum(as.numeric(predict(base_fit,newdata=d,type='response')))-1),
 mae=mean(abs(d$n_patients-pred)),rmse=sqrt(mean((d$n_patients-pred)^2)))
}))
fwrite(metrics,file.path(out,'TABLE_model_comparison.csv'))
by_year <- rbindlist(lapply(names(models),function(nm){
 data.table(year=d$year,obs=d$n_patients,pred=as.numeric(predict(models[[nm]],newdata=d,type='response')))[,
 .(observed=sum(obs),predicted=sum(pred)),by=year][,model:=nm]}))
fwrite(by_year,file.path(out,'TABLE_prediction_by_year.csv'))
utils::setTxtProgressBar(pb,length(specs)+1)

# Use original cells and fixed neighboring units; never average raw residuals.
d[,block:=interaction(year,sex,age_group,drop=TRUE,lex.order=TRUE)]
blocks <- split(seq_len(nrow(d)),d$block)
if(length(blocks)!=112L)warning('Unexpected number of age/sex/year blocks: ',length(blocks))
resids <- lapply(models,function(m){r<-as.numeric(residuals(m,type='pearson'))
 if(length(r)!=nrow(d)||any(!is.finite(r)))stop('Invalid residuals');r})
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
all_slices<-list();global<-list()
labels <- d[,.(year=year[1],sex=as.character(sex[1]),age_group=as.character(age_group[1])),by=block]
for(k in c(4L,8L,12L)){
 gr<-make_graph(k)
 for(nm in names(resids)){
  z<-calc_I(resids[[nm]],gr,TRUE)
  st<-merge(copy(z$slices),labels,by='block',sort=FALSE);st[,`:=`(model=nm,k=k)]
  all_slices[[length(all_slices)+1L]]<-st
  global[[length(global)+1L]]<-data.table(model=nm,k=k,n_cells=nrow(d),n_blocks=length(blocks),
   moran_global_block=z$I,median_slice_i=median(st$moran_i,na.rm=TRUE),
   q25_slice_i=as.numeric(quantile(st$moran_i,.25,na.rm=TRUE)),
   q75_slice_i=as.numeric(quantile(st$moran_i,.75,na.rm=TRUE)))
 }
}
ss<-rbindlist(all_slices);gg<-rbindlist(global)
fwrite(gg,file.path(out,'TABLE_global_block_moran.csv'))
fwrite(ss,file.path(out,'TABLE_moran_112_cell_slices.csv'))
fwrite(ss[,.(median_i=median(moran_i),n_slices=.N),by=.(model,k,age_group)],file.path(out,'TABLE_moran_by_age.csv'))
fwrite(ss[,.(median_i=median(moran_i),n_slices=.N),by=.(model,k,year)],file.path(out,'TABLE_moran_by_year.csv'))
fwrite(ss[,.(median_i=median(moran_i),n_slices=.N),by=.(model,k,sex)],file.path(out,'TABLE_moran_by_sex.csv'))
utils::setTxtProgressBar(pb,length(specs)+2)
# Compare per-cell predictions and age-by-year trends against baseline.
ref_pred <- as.numeric(predict(base_fit,newdata=d,type='response'))
pred_changes <- rbindlist(lapply(names(models),function(nm){p<-as.numeric(predict(models[[nm]],newdata=d,type='response'))
 data.table(model=nm,year=d$year,sex=d$sex,age_group=d$age_group,pred=p,ref=ref_pred)[,
 .(predicted=sum(pred),base_predicted=sum(ref),diff=sum(pred)-sum(ref)),by=.(model,year,sex,age_group)]}))
fwrite(pred_changes,file.path(out,'TABLE_prediction_changes_demographic.csv'))
fwrite(data.table(n=nrow(d),distinct_key=uniqueN(d,by=keycols),units=uniqueN(d$code),
 years=uniqueN(d$year),sexes=uniqueN(d$sex),age_groups=uniqueN(d$age_group),
 n_blocks=length(blocks),baseline_k=k_ref),file.path(out,'TABLE_input_integrity.csv'))
writeLines(c('NB GAM, same 42556 original cells; no residual aggregation or GAM refit of BASE.',
 'All variants include sex*age + sex-age-specific year smooth + RE(code) + RE(code-age) + RE(code)*year_c.',
 'Sex extension: RE(code-sex); pandemic: RE(code) multiplied by centered indicator of 2020/2021.',
 'Population predictor uses within-(year,sex,age) centered log(pop), offset log(pop) remains present.',
 'Population linear and smooth are alternative tests, not the same term entered twice.',
 'k_year=5/6/7 tests flexibility with only seven observed years: boundary k=7 is exploratory.',
 'AIC / deviance are in-sample comparisons, not predictive validation.',
 'Global block Moran: centered Pearson residuals within year x sex x age; symmetric kNN 4/8/12 proxy graph, not borders.',
 'No p-values printed: diagnostic focus on effect sizes. No Shapley / bootstrap / simulations.',
 'Cache each new model to allow resuming if fitting interrupted.'),file.path(out,'README_method.txt'))
utils::setTxtProgressBar(pb,length(specs)+3);close(pb)
cat('\nResults saved: ',out,'\n',sep='')
print(metrics[,.(model,aic,delta_aic_vs_base,deviance_explained,edf_total,total_predicted)])
print(gg[k==8,.(model,moran_global_block,median_slice_i)])
