#!/usr/bin/env Rscript
# Pilot only. Fits one additional patient GAM; does not overwrite the paper's outputs.
suppressPackageStartupMessages({library(data.table);library(mgcv);library(spdep)})
source('R/99_load_all.R')
cfg <- load_config('config.R')
out <- 'output/pilot_patients_code_age_re'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
ref <- 'output/pilot_mrf_patients/patient_models_RE_MRF.rds'
if (!file.exists(ref)) stop('Missing reference file: ', ref)
old <- readRDS(ref)
if (is.null(old$RE)) stop('Saved reference model RE is missing')
if (file.exists('output/rds/aggregated_data.rds')) {
  d <- as.data.table(readRDS('output/rds/aggregated_data.rds'))
} else {
  raw <- read_input_data(cfg$data$path)
  d <- as.data.table(filter_and_aggregate(normalize_schema(raw, value_var=cfg$data$value_var),
                                           disease_filter=cfg$data$disease_filter))
}
required <- c('code','x','y','year','sex','age_group','pop','n_patients')
if (length(setdiff(required, names(d)))) stop('Missing fields: ', paste(setdiff(required,names(d)),collapse=', '))
d <- d[is.finite(x)&is.finite(y)&is.finite(pop)&is.finite(n_patients)&pop>0&n_patients>=0]
d[, `:=`(code=as.character(code), sex=factor(sex), age_group=factor(age_group),
          year_num=as.numeric(year))]
coords <- d[,.(x=median(x),y=median(y),n_coords=uniqueN(paste(x,y))),by=code][order(code)]
if(anyDuplicated(coords[,.(x,y)])) stop('Duplicate unit coordinates; cannot build kNN graph')
ids <- coords$code
d[,code:=factor(code, levels=ids)]
d[,sex_age:=interaction(sex,age_group,drop=TRUE)]
d[,code_age:=interaction(code, age_group, drop=TRUE, lex.order=TRUE)]
eps <- as.numeric(cfg$model$eps_offset)
if(length(eps)!=1L||!is.finite(eps)||eps<=0) stop('Invalid eps_offset')
d[,log_pop_offset:=base::log(pmax(as.numeric(pop),eps))]
if(any(!is.finite(d$log_pop_offset))) stop('Invalid log offset')
if(nrow(d)!=nobs(old$RE)) stop('Fitting sample has changed relative to saved RE')
if(length(old$RE$y)==nrow(d) && !isTRUE(all.equal(as.numeric(old$RE$y),as.numeric(d$n_patients))))
  stop('Data row order or patient response has changed; cannot reuse saved RE')
kn <- spdep::knn2nb(spdep::knearneigh(as.matrix(coords[,.(x,y)]),k=8L))
lwfull <- spdep::nb2listw(spdep::make.sym.nb(kn),style='W',zero.policy=TRUE)

# Input-only descriptive stratification of ALREADY computed reference results, when available.
previous <- 'output/pilot_mrf_patients/moran_slices_re.csv'
if (file.exists(previous)) {
  previous_dt <- fread(previous)
  if(all(c('year','sex','age_group','moran_i') %in% names(previous_dt))) {
    fwrite(previous_dt[,.(median_i=median(moran_i,na.rm=TRUE),
                           min_i=min(moran_i,na.rm=TRUE),max_i=max(moran_i,na.rm=TRUE)),by=age_group],
           file.path(out,'existing_RE_by_age.csv'))
    fwrite(previous_dt[,.(median_i=median(moran_i,na.rm=TRUE)),by=year],
           file.path(out,'existing_RE_by_year.csv'))
  }
}

formula_new <- as.formula(paste0('n_patients ~ sex * age_group + ',
  's(year_num, by=sex_age, k=',min(cfg$model$k_year,uniqueN(d$year)),') + ',
  's(code, bs="re") + s(code_age, bs="re") + offset(log_pop_offset)'))
cat(sprintf('Patients pilot: %d rows, %d units, %d code-age levels\n',
            nrow(d),uniqueN(d$code),uniqueN(d$code_age)))
cat('Fitting ONE extra GAM (code-age RE). No Shapley/bootstrap/simulations.\n')
cat('Progress advances after GAM fit; fitting may take time and substantial RAM.\n')
pb <- utils::txtProgressBar(min=0,max=3,style=3,width=35)
on.exit <- NULL
fit <- tryCatch(mgcv::bam(formula_new,data=d,family=mgcv::nb(),method='fREML',
                          discrete=TRUE,control=mgcv::gam.control(maxit=cfg$model$maxit_gam,trace=FALSE)),
                error=function(e){close(pb);stop(e)})
utils::setTxtProgressBar(pb,1)
models <- list(RE=old$RE,RE_code_age=fit)
metrics <- rbindlist(lapply(names(models),function(nm){m<-models[[nm]];
  data.table(model=nm,n=nobs(m),AIC=AIC(m),deviance_explained=summary(m)$dev.expl,
             edf_total=sum(m$edf),total_predicted_patients=sum(predict(m,newdata=d,type='response')))}))
fwrite(metrics,file.path(out,'comparison.csv'))
utils::setTxtProgressBar(pb,2)
calc_slices <- function(m,nm) {
  rd <- data.table(code=as.character(d$code),year=d$year,sex=as.character(d$sex),
                   age_group=as.character(d$age_group),r=as.numeric(residuals(m,type='pearson')))
  rd <- rd[is.finite(r),.(r=mean(r)),by=.(year,sex,age_group,code)]
  res <- rd[, {
    present <- ids[ids %in% code]
    v <- r[match(present,code)]
    if(length(v)<20L||!is.finite(sd(v))||sd(v)==0) list(n_units=length(v),moran_i=NA_real_,p_asymptotic=NA_real_)
    else {
      lw <- spdep::subset.listw(lwfull,subset=ids %in% present,zero.policy=TRUE)
      mt <- spdep::moran.test(v,lw,zero.policy=TRUE)
      list(n_units=length(v),moran_i=unname(mt$estimate[1]),p_asymptotic=mt$p.value)
    }
  },by=.(year,sex,age_group)]
  res[,model:=nm]
  res
}
slices <- rbindlist(lapply(names(models),function(nm) calc_slices(models[[nm]],nm)))
fwrite(slices,file.path(out,'moran_slices.csv'))
summary_tab <- slices[,.(n_slices=.N,valid=sum(is.finite(moran_i)),median_i=median(moran_i,na.rm=TRUE),
   q25_i=as.numeric(quantile(moran_i,.25,na.rm=TRUE)),q75_i=as.numeric(quantile(moran_i,.75,na.rm=TRUE)),
   fraction_p_lt_005=mean(p_asymptotic<.05,na.rm=TRUE)),by=model]
fwrite(summary_tab,file.path(out,'moran_summary.csv'))
for (strata in c('age_group','year','sex')) {
  fwrite(slices[,.(median_i=median(moran_i,na.rm=TRUE),n=.N),by=c('model',strata)],
         file.path(out,paste0('moran_by_',strata,'.csv')))
}
base <- slices[model=='RE',.(year,sex,age_group,i_re=moran_i)]
extra <- slices[model=='RE_code_age',.(year,sex,age_group,i_new=moran_i)]
diff <- merge(base,extra,by=c('year','sex','age_group'))
diff[,`:=`(delta_i=i_new-i_re,improves_abs=abs(i_new)<abs(i_re))]
fwrite(diff,file.path(out,'moran_paired_difference.csv'))
byyear <- rbindlist(lapply(names(models),function(nm){data.table(year=d$year,pred=predict(models[[nm]],newdata=d,type='response'))[,.(sum_pred=sum(pred)),by=year][,model:=nm]}))
fwrite(byyear,file.path(out,'prediction_by_year.csv'))
con <- tryCatch(mgcv::concurvity(fit,full=TRUE),error=function(e)NULL)
if(!is.null(con)) write.csv(con,file.path(out,'concurvity.csv'))
writeLines(c('Design: RE versus RE + code-by-age RE; same NB formula, rows and kNN-8 residual graph.',
  'Moran tested per year x sex x age group. P-values asymptotic and exploratory.',
  'This is a pilot: NOT production results. High-dimensional random effects may be expensive.',
  capture.output(summary(fit))),file.path(out,'model_diagnostics.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'environment.txt'))
saveRDS(fit,file.path(out,'patients_RE_code_age.rds'))
utils::setTxtProgressBar(pb,3);close(pb)
cat('\nSaved results to ',out,'\n',sep='')
print(metrics)
print(summary_tab)
