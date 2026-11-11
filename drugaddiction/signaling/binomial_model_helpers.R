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
