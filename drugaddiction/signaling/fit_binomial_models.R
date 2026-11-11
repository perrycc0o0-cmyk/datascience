#!/usr/bin/env Rscript
# Binomial mixed models on new accepted IEG-high count summaries.
args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Usage: Rscript fit_binomial_models.R NEURON_METADATA.csv RUN_ROOT")
metadata_path <- normalizePath(args[1L], mustWork = TRUE)
root <- normalizePath(args[2L], mustWork = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])))
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({library(lme4); library(jsonlite); library(digest)})
out <- file.path(root,"04_proportions_GLMM")
for (d in c("tables","source_data","qa","qa/checkpoints","logs")) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
stopifnot(!file.exists(file.path(out,"qa/GLMM_complete.json")))
mdfile <- metadata_path
md <- read.csv(mdfile,stringsAsFactors=FALSE)
stopifnot(all(c("Cell","pool_id","condition","neuron_type","subtype","IEGhigh") %in% names(md)), !anyDuplicated(md$Cell), is.logical(md$IEGhigh), !anyNA(md$IEGhigh), any(md$IEGhigh))
orders <- list(Fig3L=c("GLU_L2_L3","GLU_L4_5","GLU_L5","GLU_L5_6","GLU_L6_CLA","GLU_ADGM"),
  FigS9B=c("GABA_Amb","GABA_Etv1","GABA_Lamp5","GABA_MSN","GABA_Sst","GABA_Hs3st2","GABA_Pvalb","GABA_Vip"))
pools <- c(paste0("Ext",1:3),paste0("Ret-Ext",1:3))
rows <- list(); datasets <- list()
for (panel in names(orders)) for (subtype in orders[[panel]]) {
  type <- if(panel=="Fig3L") "Excitatory" else "Inhibitory"
  m <- md[md$neuron_type==type,,drop=FALSE]
  for (pool in pools) {
    x <- m[m$pool_id==pool,,drop=FALSE]
    num <- sum(x$IEGhigh & x$subtype==subtype)
    den <- if(panel=="Fig3L") sum(x$IEGhigh) else sum(x$subtype==subtype)
    rows[[length(rows)+1L]] <- data.frame(panel,neuron_type=type,subtype,pool_id=pool,
      condition=if(startsWith(pool,"Ext")) "Ext" else "Ret-Ext",numerator=num,denominator=den,
      total_class_IEGhigh=sum(x$IEGhigh),total_subtype=sum(x$subtype==subtype),
      value=if(den>0) num/den else NA_real_,
      definition=if(panel=="Fig3L") "high_composition" else "activation_fraction")
  }
}
counts <- do.call(rbind,rows)
stopifnot(all(counts$numerator>=0),all(counts$numerator<=counts$denominator),all(counts$denominator>=0))
write.csv(counts,file.path(out,"source_data/pool_counts_original.csv"),row.names=FALSE)
for (panel in names(orders)) for (subtype in orders[[panel]]) {
  d <- counts[counts$panel==panel & counts$subtype==subtype & counts$denominator>0,,drop=FALSE]
  d$condition <- factor(d$condition,levels=c("Ext","Ret-Ext"));d$pool_id <- factor(d$pool_id)
  datasets[[paste(panel,subtype,sep="|")]] <- d
}
B <- 1999L; AGQ <- 15L
# Import only audited model helper definitions, never the historical runner or data.
helper <- file.path(script_dir,"binomial_model_helpers.R")
allowed <- c("control","form_full","form_null","fit_model","fit_pair","true_integral_ll")
for (e in parse(helper)) if (is.call(e) && is.symbol(e[[1]]) && as.character(e[[1]]) %in% c("<-","=") &&
    is.symbol(e[[2]]) && as.character(e[[2]]) %in% allowed) eval(e)
stopifnot(is.function(fit_pair),is.function(true_integral_ll))
save_csv <- function(x,p) write.csv(x,file.path(out,p),row.names=FALSE,na="NA")
run_one <- function(index) {
  d <- datasets[[index]]; key <- names(datasets)[index]; pp <- strsplit(key,"|",fixed=TRUE)[[1]]
  panel <- pp[1]; subtype <- pp[2]; seed <- 42L+index
  log <- file.path(out,"qa",paste0("progress_",panel,"_",subtype,".txt"))
  progress <- function(...) cat(format(Sys.time(),tz="UTC",usetz=TRUE),...,"\n",file=log,append=TRUE)
  progress("START",nrow(d),"nonempty pools")
  base <- data.frame(panel,subtype,n_Ext=sum(d$condition=="Ext"),n_RetExt=sum(d$condition=="Ret-Ext"),
    bootstrap_seed=seed,bootstrap_replicates=B,bootstrap_P=NA_real_,bootstrap_BH=NA_real_,
    Wald_P_asymptotic=NA_real_,Wald_BH_asymptotic=NA_real_,LR_chisq_P_asymptotic=NA_real_,
    beta_log_OR=NA_real_,conditional_OR=NA_real_,Wald_SE=NA_real_,LR=NA_real_,
    bootstrap_failures=NA_integer_,bootstrap_MCSE=NA_real_,full_singular=NA,random_SD=NA_real_,
    observed_valid=FALSE,pearson_dispersion=NA_real_,
    Ext_mean_percent=mean(d$value[d$condition=="Ext"])*100,
    RetExt_mean_percent=mean(d$value[d$condition=="Ret-Ext"])*100,
    status="not estimable: fewer than two nonempty pools per condition")
  if(any(table(d$condition)<2L)) return(base)
  pair <- fit_pair(d); highq <- fit_pair(d,25L)
  integral_lr <- if(!is.null(pair$full$fit) && !is.null(pair$null$fit)) tryCatch(
    2*(true_integral_ll(pair$full$fit,d)-true_integral_ll(pair$null$fit,d)),error=function(e) NA_real_) else NA_real_
  beta <- if(!is.null(pair$full$fit)) unname(fixef(pair$full$fit)[2]) else NA_real_
  beta25 <- if(!is.null(highq$full$fit)) unname(fixef(highq$full$fit)[2]) else NA_real_
  observed <- pair$valid && highq$valid && is.finite(integral_lr) &&
    abs(pair$lr-highq$lr)<1e-3 && abs(beta-beta25)<1e-3 && abs(pair$lr-integral_lr)<1e-3
  diagnostics <- do.call(rbind,lapply(c("full","null"),function(name) {
    f<-pair[[name]]; data.frame(panel,subtype,model=name,valid=f$valid,logLik=f$ll,random_SD=f$sd,
      singular=f$singular,min_hessian=f$hessian_min,max_gradient=f$max_gradient,
      boundary_glm_verified=f$boundary_glm_verified,notes=f$notes,LR15=pair$lr,LR25=highq$lr,
      independent_integral_LR=integral_lr,quadrature_valid=observed)
  }))
  save_csv(diagnostics,paste0("qa/fits_",panel,"_",subtype,".csv"))
  saveRDS(list(full=pair$full$fit,null=pair$null$fit),file.path(out,"qa",paste0("models_",panel,"_",subtype,".rds")))
  base$beta_log_OR<-beta;base$conditional_OR<-exp(beta);base$LR<-pair$lr;base$observed_valid<-observed
  base$full_singular<-pair$full$singular;base$random_SD<-pair$full$sd
  if(!is.null(pair$full$fit)) {
    cf<-coef(summary(pair$full$fit));base$Wald_SE<-unname(cf[2,"Std. Error"])
    base$Wald_P_asymptotic<-unname(cf[2,"Pr(>|z|)"])
    base$pearson_dispersion<-sum(residuals(pair$full$fit,type="pearson")^2)/max(1,nrow(d)-3L)
  }
  if(pair$valid) base$LR_chisq_P_asymptotic<-pchisq(pair$lr,1,lower.tail=FALSE)
  if(!observed) {base$status<-"P withheld: observed model/integration diagnostic failure";return(base)}
  sim<-data.frame(replicate=seq_len(B),LR=NA_real_,failed=FALSE,full_singular=NA,null_singular=NA)
  set.seed(seed); intercept<-unname(fixef(pair$null$fit)[1]); sd0<-pair$null$sd
  for(b in seq_len(B)) {
    ds<-d;ds$numerator<-rbinom(nrow(d),d$denominator,plogis(intercept+rnorm(nrow(d),0,sd0)))
    bp<-fit_pair(ds);sim$LR[b]<-bp$lr;sim$failed[b]<-!bp$valid
    sim$full_singular[b]<-bp$full$singular;sim$null_singular[b]<-bp$null$singular
    if(b%%100L==0L) {
      progress("bootstrap",b,"/",B,"failures",sum(sim$failed[seq_len(b)]))
      saveRDS(list(index=index,completed=b,sim=sim[seq_len(b),],RNG_state=.Random.seed),
        file.path(out,"qa/checkpoints",paste0(panel,"_",subtype,".rds")))
    }
  }
  save_csv(sim,paste0("qa/bootstrap_",panel,"_",subtype,".csv"))
  base$bootstrap_failures<-sum(sim$failed)
  if(mean(sim$failed)<=.01) {
    base$bootstrap_P<-(1+sum(sim$failed | (!is.na(sim$LR)&sim$LR>=pair$lr-1e-10)))/(B+1)
    base$bootstrap_MCSE<-sqrt(base$bootstrap_P*(1-base$bootstrap_P)/(B+1))
    base$status<-if(pair$full$singular||pair$null$singular) "boundary variance; interpret cautiously" else "interior model fit"
  } else base$status<-"P withheld: more than 1% bootstrap numerical failures"
  save_csv(base,paste0("tables/result_",panel,"_",subtype,".csv"))
  progress("COMPLETE",base$bootstrap_P,base$status);base
}
cluster <- parallel::makePSOCKcluster(as.integer(Sys.getenv("GLMM_WORKERS", "1")),outfile=file.path(out,"qa/worker_console.txt"))
parallel::clusterExport(cluster,c("out","B","AGQ","datasets","save_csv","control","form_full","form_null","fit_model","fit_pair","true_integral_ll","run_one"),envir=environment())
parallel::clusterEvalQ(cluster,{Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1");suppressPackageStartupMessages(library(lme4))})
results <- tryCatch(do.call(rbind,parallel::parLapplyLB(cluster,seq_along(datasets),run_one)),finally=parallel::stopCluster(cluster))
for(panel in names(orders)) {
  idx<-which(results$panel==panel);m<-length(orders[[panel]])
  for(cols in list(c("bootstrap_P","bootstrap_BH"),c("Wald_P_asymptotic","Wald_BH_asymptotic"))) {
    good<-idx[is.finite(results[[cols[1]]][idx])]
    results[[cols[2]]][good]<-p.adjust(results[[cols[1]]][good],"BH",n=m)
  }
  save_csv(results[idx,],paste0("tables/",panel,"_GLMM_statistics.csv"))
}
save_csv(results,"tables/All_subtypes_GLMM_statistics.csv")
capture.output(sessionInfo(),file=file.path(out,"qa/sessionInfo.txt"))
write_json(list(status="complete",method="binomial GLMM, condition fixed and pooled-sample random intercept; null-parametric-bootstrap likelihood ratio",B=B,nAGQ=AGQ,
  primary_P="bootstrap_P",primary_adjusted_P="bootstrap_BH",BH_family=list(Fig3L=6L,FigS9B=8L),
  n_models=nrow(results),n_valid=sum(results$observed_valid),n_primary_P_available=sum(is.finite(results$bootstrap_P)),
  count_MD5=unname(tools::md5sum(file.path(out,"source_data/pool_counts_original.csv"))),
  source_metadata_SHA256=digest(mdfile,file=TRUE,algo="sha256"),source_helpers_MD5=unname(tools::md5sum(helper)),
  partial_sample_variance_adjustment=NULL,no_method_search=TRUE),file.path(out,"qa/GLMM_complete.json"),auto_unbox=TRUE,pretty=TRUE)
cat("GLMM_COMPLETE\n");print(results[,c("panel","subtype","bootstrap_P","bootstrap_BH","status")],row.names=FALSE)
