#!/usr/bin/env Rscript
# Pilot: PATIENTS only. Fit RE vs kNN-MRF on identical cell rows.
# Does not run Shapley/bootstrap/simulation or modify production outputs.
suppressPackageStartupMessages({library(data.table); library(mgcv); library(spdep)})
source('R/99_load_all.R')
cfg <- load_config('config.R')
out <- 'output/pilot_mrf_patients'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
logfile <- file.path(out,'pilot.log')
log_msg <- function(...) { msg <- paste0(...,collapse=''); cat(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),msg,'\n'); cat(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),msg,'\n',file=logfile,append=TRUE) }
tryCatch({
 if(file.exists('output/rds/aggregated_data.rds')) {
   d <- as.data.table(readRDS('output/rds/aggregated_data.rds'))
 } else {
   raw <- read_input_data(cfg$data$path)
   d <- as.data.table(filter_and_aggregate(normalize_schema(raw,value_var=cfg$data$value_var), disease_filter=cfg$data$disease_filter))
 }
 needed <- c('code','x','y','year','sex','age_group','pop','n_patients')
 stopifnot(length(setdiff(needed,names(d)))==0)
 d <- d[is.finite(x)&is.finite(y)&is.finite(pop)&is.finite(n_patients)&pop>0 & n_patients>=0]
 d[, `:=`(code=as.character(code),sex=factor(sex),age_group=factor(age_group),year_num=as.numeric(year))]
 # One coordinate pair per administrative unit; detect inconsistencies.
 coords <- d[,.(x=median(x),y=median(y),n_coords=uniqueN(paste(x,y))),by=code][order(code)]
 if(any(coords$n_coords>1)) log_msg('WARNING: varying coordinates within code; using medians')
 if(anyDuplicated(coords[,.(x,y)])>0) stop('Multiple codes have identical coordinates. Supply a real adjacency map first.')
 if(nrow(coords)<20) stop('Insufficient spatial units')
 k <- 8L
 coordsm <- as.matrix(coords[,.(x,y)])
 kn <- spdep::knn2nb(spdep::knearneigh(coordsm,k=k))
 # Symmetrize graph: MRF requires reciprocal neighbours.
 nb <- spdep::make.sym.nb(kn)
 ids <- as.character(coords$code)
 names(nb) <- ids
 nb <- lapply(nb, as.integer)
 names(nb) <- ids
 stopifnot(all(vapply(seq_along(nb), function(i) all(vapply(nb[[i]],function(j) i %in% nb[[j]],logical(1))),logical(1))))
 graph <- data.table(code=ids,n_neighbors=lengths(nb))
 fwrite(graph,file.path(out,'graph_degree.csv'))
 d[, code:=factor(code,levels=ids)]
 d[, sex_age:=interaction(sex,age_group,drop=TRUE)]
 yrs <- sort(unique(d$year))
 ky <- min(cfg$model$k_year,length(yrs)); if(ky<3) stop('At least 3 years required')
 rhs <- c('sex * age_group',sprintf('s(year_num, by=sex_age, k=%d)',ky))
 # Precompute offset: mgcv::bam() model.frame may fail to evaluate
 # nested pmax() inside offset() when reconstructing its model frame.
 eps <- as.numeric(cfg$model$eps_offset)
 if(length(eps)!=1L || !is.finite(eps) || eps<=0) stop('Invalid eps_offset')
 if (!'pop' %in% names(d) || !is.numeric(d[['pop']])) stop('Missing numeric population column pop')
 d[, log_pop_offset := base::log(pmax(as.numeric(pop), eps))]
 if (length(d[['log_pop_offset']]) != nrow(d) || any(!is.finite(d[['log_pop_offset']]))) stop('Invalid population offset; check pop values')
 offs <- 'offset(log_pop_offset)'
 fre <- as.formula(paste('n_patients ~',paste(c(rhs,"s(code, bs='re')",offs),collapse=' + ')))
 fmrf <- as.formula(paste('n_patients ~',paste(c(rhs,"s(code, bs='mrf', xt=list(nb=nb), k=100)",offs),collapse=' + ')))
 # Keep k at or below available rank, and well below 380; code-level MRF graph is kNN proxy.
 # 100 basis functions, not a fully saturated spatial effect; may be tuned separately.
 log_msg('N rows=',nrow(d),'; units=',nrow(coords),'; years=',length(yrs),'; graph=symmetric kNN-8 (not border adjacency)')
 ctrl <- mgcv::gam.control(maxit=cfg$model$maxit_gam,trace=isTRUE(cfg$model$trace_gam))
 log_msg('Fitting RE patients model...')
 mre <- mgcv::bam(fre,data=d,family=mgcv::nb(),method='fREML',discrete=TRUE,control=ctrl)
 log_msg('Fitting MRF patients model...')
 mmrf <- mgcv::bam(fmrf,data=d,family=mgcv::nb(),method='fREML',discrete=TRUE,control=ctrl)
 log_msg('Fits complete')
 fitstats <- rbindlist(lapply(list(RE=mre,MRF=mmrf),function(m) data.table(n=nobs(m),AIC=AIC(m),deviance_explained=summary(m)$dev.expl,edf_total=sum(m$edf),theta=m$family$getTheta(TRUE))),idcol='model')
 fwrite(fitstats,file.path(out,'model_comparison.csv'))
 # Per-slice Moran's I. Each year x sex x age group has one spatial value per code;
 # no pooling across repeats as if independent spatial observations.
 for(name in c('RE','MRF')) {
  fit <- if(name=='RE') mre else mmrf
  rd <- data.table(code=as.character(d$code),year=d$year,sex=as.character(d$sex),age_group=as.character(d$age_group),r=as.numeric(residuals(fit,type='pearson')))
  rd <- rd[is.finite(r),.(r=mean(r)),by=.(year,sex,age_group,code)]
  setkey(rd,code)
  lwfull <- spdep::nb2listw(spdep::make.sym.nb(kn),style='W',zero.policy=TRUE)
  slices <- rd[,{
    ids_here <- ids[ids %in% code]
    vec <- r[match(ids_here,code)]
    if(length(ids_here)<20 || sd(vec)==0) .(n_units=length(ids_here),moran_i=NA_real_,p_asymptotic=NA_real_) else {
     subset_idx <- match(ids_here,ids)
     lw <- spdep::subset.listw(lwfull,subset=seq_along(ids) %in% subset_idx,zero.policy=TRUE)
     tst <- spdep::moran.test(vec,lw,zero.policy=TRUE)
     .(n_units=length(vec),moran_i=unname(tst$estimate[1]),p_asymptotic=tst$p.value)
    }
  },by=.(year,sex,age_group)]
  slices[,model:=name]
  fwrite(slices,file.path(out,paste0('moran_slices_',tolower(name),'.csv')))
 }
 a<-fread(file.path(out,'moran_slices_re.csv'))
 b<-fread(file.path(out,'moran_slices_mrf.csv'))
 cmp<-merge(a,b,by=c('year','sex','age_group'),suffixes=c('_re','_mrf'))
 cmp[,`:=`(delta_i=moran_i_mrf-moran_i_re,abs_re=abs(moran_i_re),abs_mrf=abs(moran_i_mrf))]
 fwrite(cmp,file.path(out,'moran_slice_comparison.csv'))
 summary_tab<-cmp[,.(n_slices=.N,median_re=median(moran_i_re,na.rm=TRUE),median_mrf=median(moran_i_mrf,na.rm=TRUE),median_abs_re=median(abs_re,na.rm=TRUE),median_abs_mrf=median(abs_mrf,na.rm=TRUE),fraction_mrf_lower_abs=mean(abs_mrf<abs_re,na.rm=TRUE),fraction_significant_re=mean(p_asymptotic_re<.05,na.rm=TRUE),fraction_significant_mrf=mean(p_asymptotic_mrf<.05,na.rm=TRUE))]
 fwrite(summary_tab,file.path(out,'moran_summary.csv'))
 # Prediction check in original response scale, keeping data identical.
 pred <- data.table(model=c('RE','MRF'),total_predicted_patients=c(sum(predict(mre,newdata=d,type='response')),sum(predict(mmrf,newdata=d,type='response'))))
 fwrite(pred,file.path(out,'prediction_totals.csv'))
 saveRDS(list(RE=mre,MRF=mmrf),file.path(out,'patient_models_RE_MRF.rds'))
 writeLines(c(capture.output(sessionInfo()),'',paste('GRAPH: symmetric 8-nearest-neighbour proximity proxy, NOT geographic boundary adjacency'),paste('Run at:',Sys.time())),file.path(out,'environment.txt'))
 log_msg('Finished. Outputs at: ',out)
 print(fitstats);print(summary_tab);print(pred)
},error=function(e){log_msg('ERROR: ',conditionMessage(e));stop(e)})
