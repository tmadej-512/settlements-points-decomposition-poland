#!/usr/bin/env Rscript
# One additional fit: patients GAM with RE + graph-MRF.
# Existing RE/MRF fits are READ ONLY. No changes to production analyses.
suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(spdep)
})
source('R/99_load_all.R')
cfg <- load_config('config.R')
ref_dir <- 'output/pilot_mrf_patients'
out_dir <- 'output/pilot_mrf_patients_re_plus_mrf'
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
ref_file <- file.path(ref_dir, 'patient_models_RE_MRF.rds')
if (!file.exists(ref_file)) stop('Missing ', ref_file, '; copy original pilot output into the project.')
old <- readRDS(ref_file)
if (!all(c('RE','MRF') %in% names(old))) stop('Saved pilot must contain models RE and MRF.')

if (file.exists('output/rds/aggregated_data.rds')) {
  d <- as.data.table(readRDS('output/rds/aggregated_data.rds'))
} else {
  raw <- read_input_data(cfg$data$path)
  d <- as.data.table(filter_and_aggregate(
    normalize_schema(raw, value_var=cfg$data$value_var),
    disease_filter=cfg$data$disease_filter))
}
needed <- c('code','x','y','year','sex','age_group','pop','n_patients')
if(length(setdiff(needed,names(d)))) stop('Missing columns: ',paste(setdiff(needed,names(d)),collapse=', '))
d <- d[is.finite(x)&is.finite(y)&is.finite(pop)&is.finite(n_patients)&pop>0&n_patients>=0]
d[,`:=`(code=as.character(code),sex=factor(sex),age_group=factor(age_group),year_num=as.numeric(year))]
coords <- d[,.(x=median(x),y=median(y),n_coords=uniqueN(paste(x,y))),by=code][order(code)]
if (any(coords$n_coords>1)) warning('Varying coordinates within unit; medians used.')
if (anyDuplicated(coords[,.(x,y)])) stop('Duplicate coordinates across units; cannot build kNN graph.')
ids <- as.character(coords$code)
kn <- spdep::knn2nb(spdep::knearneigh(as.matrix(coords[,.(x,y)]),k=8L))
nb <- spdep::make.sym.nb(kn)
names(nb) <- ids
nb <- lapply(nb,as.integer)
names(nb) <- ids
stopifnot(all(vapply(seq_along(nb),function(i) all(vapply(nb[[i]],function(j) i %in% nb[[j]],logical(1))),logical(1))))
d[,code:=factor(code,levels=ids)]
d[,sex_age:=interaction(sex,age_group,drop=TRUE)]
eps <- as.numeric(cfg$model$eps_offset)
if(length(eps)!=1L||!is.finite(eps)||eps<=0) stop('Invalid eps_offset')
d[,log_pop_offset:=base::log(pmax(as.numeric(pop),eps))]
if(any(!is.finite(d$log_pop_offset))) stop('Invalid offset')

# Verify exactly the same fitting sample and the same baseline GAM specification.
if(nrow(d)!=nobs(old$RE)||nrow(d)!=nobs(old$MRF)) stop('Observation count differs from saved models; do not compare.')
expected_re <- as.formula(paste('n_patients ~ sex * age_group +',
  sprintf('s(year_num, by=sex_age, k=%d)',min(cfg$model$k_year,uniqueN(d$year))),
  "+ s(code, bs='re') + offset(log_pop_offset)"))
if(!identical(gsub('[[:space:]]','',deparse(formula(old$RE))),
              gsub('[[:space:]]','',deparse(expected_re)))) {
  warning('Saved RE formula differs in formatting/specification; inspect formulas before interpreting comparison.')
}
# Compare response/offset alignment where available. Refit from same input cache if uncertain.
old_response <- tryCatch(as.numeric(old$RE$y),error=function(e) numeric(0))
if(length(old_response)==nrow(d) && !isTRUE(all.equal(old_response,as.numeric(d$n_patients))))
 stop('Row order or response differs from saved RE pilot; results are not comparable.')

form <- as.formula(paste('n_patients ~ sex * age_group +',
  sprintf('s(year_num, by=sex_age, k=%d)',min(cfg$model$k_year,uniqueN(d$year))),
  "+ s(code, bs='re') + s(code, bs='mrf', xt=list(nb=nb), k=100) + offset(log_pop_offset)"))
# mgcv smooths use the same unit factor: concurvity/identifiability requires review.
ctrl <- mgcv::gam.control(maxit=cfg$model$maxit_gam,trace=FALSE)
cat('RE + MRF patients: fitting 1 model (42k cells). Progress updates after fit only.\n')
pb <- utils::txtProgressBar(min=0,max=2,style=3,width=40)
fit <- tryCatch(mgcv::bam(form,data=d,family=mgcv::nb(),method='fREML',discrete=TRUE,control=ctrl),
                error=function(e){close(pb);stop(e)})
utils::setTxtProgressBar(pb,1)

models <- list(RE=old$RE,MRF=old$MRF,RE_plus_MRF=fit)
fitstats <- rbindlist(lapply(names(models),function(nm){m <- models[[nm]];
 data.table(model=nm,n=nobs(m),AIC=AIC(m),deviance_explained=summary(m)$dev.expl,
            edf_total=sum(m$edf),theta=m$family$getTheta(TRUE),
            total_predicted_patients=sum(predict(m,newdata=d,type='response')))}))
fwrite(fitstats,file.path(out_dir,'comparison_three_models.csv'))

# Slice-wise Moran for same graph across the SAME cell/administrative-unit support.
lwfull <- spdep::nb2listw(spdep::make.sym.nb(kn),style='W',zero.policy=TRUE)
calc_slices <- function(m,model_name){
 rd <- data.table(code=as.character(d$code),year=d$year,sex=as.character(d$sex),
   age_group=as.character(d$age_group),r=as.numeric(residuals(m,type='pearson')))
 rd <- rd[is.finite(r),.(r=mean(r)),by=.(year,sex,age_group,code)]
 rd[,{
   ids_here <- ids[ids %in% code]
   vec <- r[match(ids_here,code)]
   if(length(vec)<20||sd(vec)==0) list(n_units=length(vec),moran_i=NA_real_,p_asymptotic=NA_real_)
   else {
     lw <- spdep::subset.listw(lwfull,subset=ids %in% ids_here,zero.policy=TRUE)
     tst <- spdep::moran.test(vec,lw,zero.policy=TRUE)
     list(n_units=length(vec),moran_i=unname(tst$estimate[1]),p_asymptotic=tst$p.value)
   }
 },by=.(year,sex,age_group)][,model:=model_name][]
}
slices <- rbindlist(lapply(names(models),function(nm) calc_slices(models[[nm]],nm)))
fwrite(slices,file.path(out_dir,'moran_slices_three_models.csv'))
summary_tab <- slices[,.(n_slices=.N,median_i=median(moran_i,na.rm=TRUE),
  q25_i=quantile(moran_i,0.25,na.rm=TRUE),q75_i=quantile(moran_i,0.75,na.rm=TRUE),
  median_abs_i=median(abs(moran_i),na.rm=TRUE),
  frac_p_below_005=mean(p_asymptotic<0.05,na.rm=TRUE)),by=model]
fwrite(summary_tab,file.path(out_dir,'moran_summary_three_models.csv'))
# Explicit same-slice comparison RE against RE+MRF.
ref <- slices[model=='RE',.(year,sex,age_group,i_re=moran_i)]
new <- slices[model=='RE_plus_MRF',.(year,sex,age_group,i_combined=moran_i)]
by_slice <- merge(ref,new,by=c('year','sex','age_group'))
by_slice[,`:=`(delta_i=i_combined-i_re,improves_abs=abs(i_combined)<abs(i_re))]
fwrite(by_slice,file.path(out_dir,'moran_slices_RE_vs_combined.csv'))
# Report duplicate/spatial smooth diagnostics; avoid claiming fully separate effects are identifiable.
con <- tryCatch(mgcv::concurvity(fit,full=TRUE),error=function(e)NULL)
if(!is.null(con)) write.csv(con,file.path(out_dir,'concurvity.csv'))
writeLines(c('Graph: symmetric kNN-8, NOT shared-border adjacency.',
  'Single new RE+MRF fit. RE and MRF from existing saved pilot.',
  'Moran p-values are asymptotic and exploratory across dependent slices.',
  'RE and MRF smooths on the same code may suffer concurvity; inspect edf, summary and concurvity.',
  capture.output(summary(fit))), file.path(out_dir,'model_diagnostics.txt'))
writeLines(capture.output(sessionInfo()),file.path(out_dir,'environment.txt'))
saveRDS(fit,file.path(out_dir,'patients_RE_plus_MRF.rds'))
utils::setTxtProgressBar(pb,2);close(pb)
cat('\nSaved outputs to: ',out_dir,'\n',sep='')
