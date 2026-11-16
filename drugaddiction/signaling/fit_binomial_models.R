# fit binomial models

# 1. Packages
library(lme4)
library(jsonlite)
library(digest)

# 2. Paths
work_dir <- "/path/to/work"
neuron_metadata <- file.path(work_dir, "stage3/01_inputs/neuron_metadata.csv")
result_dir <- file.path(work_dir, "stage3")
bootstrap_workers <- 1L
figure_dir <- file.path(result_dir, "figures")
script_path <- "/path/to/fit_binomial_models.R"

# 3. Samples and conditions
sample_cols <- c("WD1", "WD2", "WD3", "Ext1", "Ext2", "Ext3", "Ret-Ext1", "Ret-Ext2", "Ret-Ext3")
metadata <- data.frame(
  pool_id = sample_cols,
  source_sample_id = c("Con1", "Con2", "Con3", "Ext1", "Ext2", "Ext3", "Ret_ext1", "Ret_ext2", "Ret_ext3"),
  model = rep("Rat IL", 9L),
  condition = factor(rep(c("WD", "Ext", "Ret-Ext"), each = 3L), levels = c("WD", "Ext", "Ret-Ext")),
  rats_per_pool = rep(5L, 9L), stringsAsFactors = FALSE)

# 6. Analysis and output
# Binomial mixed models on new accepted IEG-high count summaries.
metadata_path <- normalizePath(neuron_metadata, mustWork = TRUE)
root <- normalizePath(result_dir, mustWork = TRUE)
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
# Standard pooled-sample binomial model fitting and quadrature checks.
control <- glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000),
  check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore",check.rankX="stop.deficient")
form_full <- cbind(numerator,denominator-numerator) ~ condition + (1|pool_id)
form_null <- cbind(numerator,denominator-numerator) ~ 1 + (1|pool_id)

fit_model <- function(d,full=TRUE,agq=AGQ) {
  warnings <- messages <- character()
  fit <- tryCatch(withCallingHandlers(glmer(if(full) form_full else form_null,
    data=d,family=binomial(link="logit"),nAGQ=agq,control=control),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")},
    message=function(m){messages<<-c(messages,conditionMessage(m));invokeRestart("muffleMessage")}),
    error=function(e){warnings<<-c(warnings,conditionMessage(e));NULL})
  if(is.null(fit)) return(list(fit=NULL,valid=FALSE,ll=NA_real_,sd=NA_real_,singular=NA,
    hessian_min=NA_real_,max_gradient=NA_real_,boundary_glm_verified=FALSE,
    scaled_fixed_score=NA_real_,notes=paste(warnings,collapse=" | ")))
  oi <- fit@optinfo
  conv_messages <- oi$conv$lme4$messages
  serious <- conv_messages[!grepl("boundary.*singular",conv_messages,ignore.case=TRUE)]
  optimizer_ok <- all(oi$conv$opt==0)
  hess <- oi$derivs$Hessian
  he <- if(is.null(hess)) NA_real_ else tryCatch(min(eigen(hess,symmetric=TRUE,only.values=TRUE)$values),error=function(e) NA_real_)
  grad <- if(is.null(oi$derivs$gradient)) NA_real_ else max(abs(oi$derivs$gradient))
  ll <- as.numeric(logLik(fit))
  sd <- as.numeric(getME(fit,"theta")); singular <- isSingular(fit,tol=1e-4)
  boundary_glm_verified <- FALSE; scaled_fixed_score <- NA_real_
  # lme4 2.0 may omit derivatives at a legitimate zero-variance boundary.
  # Verify its fixed-effect score and exact binomial likelihood, instead of
  # misclassifying a missing random-variance Hessian as optimizer failure.
  if(singular&&sd<=1e-4&&optimizer_ok&&length(serious)==0&&is.finite(ll)) {
    boundary_glm_verified <- tryCatch({
      X <- getME(fit,"X"); mu <- plogis(drop(X%*%fixef(fit)))
      score <- drop(crossprod(X,d$numerator-d$denominator*mu))
      information <- crossprod(X,X*as.vector(d$denominator*mu*(1-mu)))
      scaled_fixed_score <- max(abs(backsolve(chol(information),score,transpose=TRUE)))
      relative_ll <- sum(dbinom(d$numerator,d$denominator,mu,log=TRUE))-
        sum(dbinom(d$numerator,d$denominator,d$numerator/d$denominator,log=TRUE))
      scaled_fixed_score<1e-3&&abs(relative_ll-ll)<1e-3
    },error=function(e) FALSE)
  }
  curvature_ok <- (is.finite(he)&&he>1e-6)||boundary_glm_verified
  valid <- optimizer_ok&&length(serious)==0&&is.finite(ll)&&all(is.finite(fixef(fit)))&&curvature_ok
  list(fit=fit,valid=valid,ll=ll,sd=sd,singular=singular,
       hessian_min=he,max_gradient=grad,boundary_glm_verified=boundary_glm_verified,
       scaled_fixed_score=scaled_fixed_score,
       notes=paste(unique(c(warnings,messages,conv_messages)),collapse=" | "))
}
fit_pair <- function(d,agq=AGQ) {
  null <- fit_model(d,FALSE,agq); full <- fit_model(d,TRUE,agq)
  valid <- null$valid&&full$valid&&full$ll>=null$ll-1e-6
  list(null=null,full=full,valid=valid,lr=if(valid) max(0,2*(full$ll-null$ll)) else NA_real_)
}
true_integral_ll <- function(fit,d) {
  X <- model.matrix(delete.response(terms(fit,fixed.only=TRUE)),d)
  eta <- drop(X%*%fixef(fit)); sd <- as.numeric(getME(fit,"theta"))
  sum(vapply(seq_len(nrow(d)),function(j) {
    if(sd<1e-8) return(dbinom(d$numerator[j],d$denominator[j],plogis(eta[j]),log=TRUE))
    log_integrand <- function(z) dbinom(d$numerator[j],d$denominator[j],plogis(eta[j]+sd*z),log=TRUE)+dnorm(z,log=TRUE)
    opt <- optimize(function(z) -log_integrand(z),interval=c(-12,12))
    peak <- log_integrand(opt$minimum)
    val <- integrate(function(z) exp(log_integrand(z)-peak),lower=-Inf,upper=Inf,
      rel.tol=1e-9,abs.tol=1e-11,subdivisions=300L,stop.on.error=TRUE)$value
    peak+log(val)
  },numeric(1)))
}

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
cluster <- parallel::makePSOCKcluster(as.integer(bootstrap_workers),outfile=file.path(out,"qa/worker_console.txt"))
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
  source_metadata_SHA256=digest(mdfile,file=TRUE,algo="sha256"),source_helpers_MD5=unname(tools::md5sum(script_path)),
  partial_sample_variance_adjustment=NULL,no_method_search=TRUE),file.path(out,"qa/GLMM_complete.json"),auto_unbox=TRUE,pretty=TRUE)
cat("GLMM_COMPLETE\n");print(results[,c("panel","subtype","bootstrap_P","bootstrap_BH","status")],row.names=FALSE)
