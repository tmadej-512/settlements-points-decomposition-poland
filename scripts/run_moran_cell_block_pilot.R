#!/usr/bin/env Rscript
# No model refits. Raw cell-level Pearson residuals, no averaging or collapsing.
suppressPackageStartupMessages({library(data.table);library(spdep);library(mgcv)})
source('R/99_load_all.R')
cfg <- load_config('config.R')
out <- 'output/pilot_moran_cell_block'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ref_file <- 'output/pilot_mrf_patients/patient_models_RE_MRF.rds'
comb_file <- 'output/pilot_mrf_patients_re_plus_mrf/patients_RE_plus_MRF.rds'
if (!file.exists(ref_file) || !file.exists(comb_file)) stop('Missing saved RDS models: ',ref_file,' or ',comb_file)
old <- readRDS(ref_file)
if(!all(c('RE','MRF') %in% names(old))) stop('RDS missing RE/MRF.')
models <- list(RE=old$RE, MRF=old$MRF, RE_plus_MRF=readRDS(comb_file))
if(file.exists('output/rds/aggregated_data.rds')) {
 d <- as.data.table(readRDS('output/rds/aggregated_data.rds'))
} else {
 raw <- read_input_data(cfg$data$path)
 d <- as.data.table(filter_and_aggregate(normalize_schema(raw,value_var=cfg$data$value_var),disease_filter=cfg$data$disease_filter))
}
needed <- c('code','x','y','year','sex','age_group','pop','n_patients')
if(length(setdiff(needed,names(d)))) stop('Missing columns: ',paste(setdiff(needed,names(d)),collapse=', '))
d <- d[is.finite(x)&is.finite(y)&is.finite(pop)&is.finite(n_patients)&pop>0&n_patients>=0]
d[,`:=`(code=as.character(code),sex=factor(sex),age_group=factor(age_group),year_num=as.numeric(year))]
coords <- d[,.(x=median(x),y=median(y),ncoords=uniqueN(paste(x,y))),by=code][order(code)]
if(any(coords$ncoords>1)) stop('Coordinates vary within code. Resolve before spatial inference.')
if(anyDuplicated(coords[,.(x,y)])) stop('Duplicated centroid coordinates across units.')
ids <- coords$code
d[,code:=factor(code,levels=ids)]
d[,sex_age:=interaction(sex,age_group,drop=TRUE)]
eps <- as.numeric(cfg$model$eps_offset)
d[,log_pop_offset:=base::log(pmax(as.numeric(pop),eps))]
if(any(!is.finite(d$log_pop_offset))) stop('Invalid population offset.')
keycols <- c('year','sex','age_group','code')
if(anyDuplicated(d[,..keycols])) stop('Duplicate cells detected. NO aggregation permitted: inspect input cell key.')
for(nm in names(models)) {
 m <- models[[nm]]
 if(nobs(m)!=nrow(d)) stop('Model ',nm,' nobs != data rows; refuse mismatched residuals.')
 yy <- tryCatch(as.numeric(m$y),error=function(e) numeric(0))
 if(length(yy)!=nrow(d) || !isTRUE(all.equal(yy,as.numeric(d$n_patients),tolerance=0,check.attributes=FALSE)))
   stop('Response/row order cannot be verified for model ',nm,'. Refuse unsafe residual alignment.')
}
# One block per year x sex x age, each observation = one individual spatial unit within the block.
d[,block:=interaction(year,sex,age_group,drop=TRUE,lex.order=TRUE)]
blocks <- split(seq_len(nrow(d)),d$block)
if(length(blocks)!=112L) warning('Number of available strata differs from 112: ',length(blocks))
if(any(lengths(blocks)<20L)) stop('Slice has fewer than 20 cells; inspect input.')
resids <- lapply(models,function(m) as.numeric(residuals(m,type='pearson')))
if(any(vapply(resids,function(z) any(!is.finite(z)),logical(1)))) stop('Non-finite Pearson residuals.')
set.seed(20261019)
nperm <- suppressWarnings(as.integer(Sys.getenv('VIH_MORAN_PERMUTATIONS','199')))
if(is.na(nperm)||nperm<99L) stop('VIH_MORAN_PERMUTATIONS must be integer >=99')
# Build symmetrized centroid kNN for each k; then restrict neighbours to the SAME demographic-year slice.
# The resulting raw-cell graph has no across-block links, ever.
make_graph <- function(k) {
 kn <- spdep::knearneigh(as.matrix(coords[,.(x,y)]),k=k)
 nb <- spdep::make.sym.nb(spdep::knn2nb(kn))
 grp <- vector('list',length(blocks))
 for(b in seq_along(blocks)) {
  ix <- blocks[[b]]
  idx <- match(as.character(d$code[ix]),ids)
  if(anyNA(idx)) stop('Unknown administrative code.')
  local <- integer(length(ids));local[idx] <- seq_along(idx)
  ii <- integer();jj <- integer();ww <- numeric()
  for(u in seq_along(idx)) {
   neigh <- local[nb[[idx[u]]]]
   neigh <- neigh[neigh>0L]
   if(length(neigh)) {
    ii <- c(ii,rep.int(u,length(neigh)))
    jj <- c(jj,neigh)
    ww <- c(ww,rep.int(1/length(neigh),length(neigh)))
   }
  }
  grp[[b]] <- list(ix=ix,i=ii,j=jj,w=ww,S0=sum(ww),n=length(ix),isolates=length(ix)-length(unique(ii)))
 }
 grp
}
# Pooled Moran on within-stratum-centered original residuals. This uses all cells without aggregating.
# Each directed edge has row-standardized weight within its stratum.
calc_stats <- function(r,graph,do_slices=TRUE) {
 numerator <- 0;denominator <- 0;S0 <- 0; n <- 0
 per <- if(do_slices) vector('list',length(graph)) else NULL
 for(b in seq_along(graph)) {
  g <- graph[[b]];v <- r[g$ix];z <- v-mean(v)
  num <- sum(g$w*z[g$i]*z[g$j]);den <- sum(z*z)
  numerator<-numerator+num; denominator<-denominator+den;S0<-S0+g$S0;n<-n+g$n
  if(do_slices) per[[b]] <- data.table(block=names(blocks)[b],n_cells=g$n,n_isolates=g$isolates,
       moran_i=if(den>0&&g$S0>0) (g$n/g$S0)*num/den else NA_real_,
       res_var=stats::var(v))
 }
 list(I=if(denominator>0&&S0>0) (n/S0)*numerator/denominator else NA_real_,
      n_cells=n,S0=S0,slices=if(do_slices)rbindlist(per)else NULL)
}
# Progress bar has stages: graph construction and 3 x 3 models+graphs with bounded block permutations.
cat('Cell-level Moran (RE / MRF / RE+MRF), no GAM fitting.\n')
pb <- utils::txtProgressBar(min=0,max=12,style=3,width=35);step<-0L
all_slices<-list(); all_global<-list(); all_detail<-list()
for(k in c(4L,8L,12L)) {
 graph <- make_graph(k);step<-step+1L;utils::setTxtProgressBar(pb,step)
 for(nm in names(resids)) {
  r <- resids[[nm]];obs <- calc_stats(r,graph,TRUE)
  slices <- copy(obs$slices)
  # Bring back slice labels without parsing the interaction string.
  labels <- d[,.(year=year[1],sex=as.character(sex[1]),age_group=as.character(age_group[1])),by=block]
  slices <- merge(slices,labels,by='block',sort=FALSE)
  slices[,`:=`(model=nm,k=k)]
  all_slices[[length(all_slices)+1L]] <- slices
  perm_i <- numeric(nperm)
  for(p in seq_len(nperm)) {
    rp <- r
    for(g in graph) rp[g$ix] <- sample(r[g$ix],length(g$ix),replace=FALSE)
    perm_i[p] <- calc_stats(rp,graph,FALSE)$I
  }
  # Two-sided p-value based on absolute deviations around the permutation null mean.
  null_mean<-mean(perm_i)
  pval<-(1+sum(abs(perm_i-null_mean)>=abs(obs$I-null_mean)))/(nperm+1)
  all_global[[length(all_global)+1L]]<-data.table(model=nm,k=k,n_cells=obs$n_cells,
       n_blocks=length(blocks),moran_global_block=obs$I,perm_null_mean=null_mean,
       p_two_sided_within_block=pval,n_permutations=nperm,
       median_slice_i=median(slices$moran_i,na.rm=TRUE),
       q25_slice_i=as.numeric(quantile(slices$moran_i,.25,na.rm=TRUE)),
       q75_slice_i=as.numeric(quantile(slices$moran_i,.75,na.rm=TRUE)),
       fraction_positive=mean(slices$moran_i>0,na.rm=TRUE))
  all_detail[[length(all_detail)+1L]]<-data.table(model=nm,k=k,perm=seq_len(nperm),moran_i=perm_i)
  step<-step+1L;utils::setTxtProgressBar(pb,step)
 }
}
close(pb)
sg<-rbindlist(all_global);ss<-rbindlist(all_slices);pd<-rbindlist(all_detail)
fwrite(sg,file.path(out,'TABLE_global_block_moran.csv'))
fwrite(ss,file.path(out,'TABLE_moran_112_cell_slices.csv'))
fwrite(pd,file.path(out,'TABLE_permutation_null_distribution.csv'))
# Compare baseline RE vs two alternatives, matched strata; no mean residuals at any stage.
wide<-dcast(ss,year+sex+age_group+k~model,value.var='moran_i')
wide[,`:=`(MRF_minus_RE=MRF-RE,RE_plus_MRF_minus_RE=RE_plus_MRF-RE)]
fwrite(wide,file.path(out,'TABLE_stratum_model_comparison.csv'))
by_age <- ss[,.(n_slices=.N,median_i=median(moran_i),q25_i=as.numeric(quantile(moran_i,.25)),
                q75_i=as.numeric(quantile(moran_i,.75))),by=.(model,k,age_group)]
fwrite(by_age,file.path(out,'TABLE_age_summary.csv'))
checks <- data.table(metric=c('observations','blocks','units','unique_cell_keys','min_stratum_size','max_stratum_size'),
 value=as.character(c(nrow(d),length(blocks),length(ids),uniqueN(d,by=keycols),min(lengths(blocks)),max(lengths(blocks)))))
fwrite(checks,file.path(out,'TABLE_input_integrity.csv'))
writeLines(c('RAW MODEL CELL = code x year x sex x age_group. No residual aggregation or replication.',
  'Within each cell stratum (year x sex x age), Moran compares administrative units.',
  'Pooled block I: within-block-centered Pearson residuals; row standardized graph restricted to blocks.',
  'Pooled statistic is variance-weighted; not an arithmetic mean of slice I.',
  'Null: residual labels permuted within each stratum independently. p-values exploratory for fitted model residuals.',
  'Graphs: symmetric k-nearest centroids, NOT shared-border adjacency.',
  'No refits; original models are read-only. All model response/row positions checked.',
  paste('R:',R.version.string),paste('mgcv:',as.character(packageVersion('mgcv'))),
  paste('spdep:',as.character(packageVersion('spdep')))),file.path(out,'README_method.txt'))
cat('\nOutputs: ',out,'\n',sep='')
print(sg[,.(model,k,moran_global_block,p_two_sided_within_block,median_slice_i)])
