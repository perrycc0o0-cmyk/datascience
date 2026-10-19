#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(lme4)
})
source('R/common.R')
ctx <- pipeline_context()
out <- new_stage(ctx, '08_GLMM')
samples <- read_samples(ctx)
samples <- samples[samples$condition %in% c('Ext', 'Ret-Ext'), , drop = FALSE]
settings <- ctx$cfg$glmm
if (!identical(settings$family, 'binomial') || !identical(settings$link, 'logit') ||
    !identical(settings$fixed_effect, 'condition') ||
    !identical(settings$random_intercept, 'sample') ||
    !identical(settings$p_method, 'Wald_z') || !identical(settings$adjust, 'BH')) {
  stop('GLMM settings must specify binomial logit, condition, sample, Wald_z and BH')
}

subtypes <- list(
  Excitatory = c('GLU_L6_CLA', 'GLU_L5_6', 'GLU_L2_L3',
                 'GLU_L5', 'GLU_L4_5', 'GLU_ADGM'),
  Inhibitory = c('GABA_Amb', 'GABA_Etv1', 'GABA_Lamp5', 'GABA_MSN',
                 'GABA_Sst', 'GABA_Hs3st2', 'GABA_Pvalb', 'GABA_Vip')
)

read_neuron_metadata <- function(neuron_type) {
  obj <- readRDS(stage_path(ctx, '04_ieg', paste0(neuron_type, '.rds')))
  require_neuron_metadata(obj)
  m <- obj[[]]
  if (!'IEGhigh' %in% names(m)) stop('IEGhigh metadata is required: ', neuron_type)
  high <- m$IEGhigh
  if (is.factor(high)) high <- as.character(high)
  if (!is.logical(high)) {
    if (anyNA(high) || !all(high %in% c('TRUE', 'FALSE', 0, 1))) {
      stop('IEGhigh must contain logical values: ', neuron_type)
    }
    high <- if (is.character(high) && all(high %in% c('TRUE', 'FALSE'))) {
      high == 'TRUE'
    } else as.logical(as.integer(high))
  }
  if (anyNA(high)) stop('IEGhigh contains missing values: ', neuron_type)
  m$IEGhigh <- high
  if (any(!as.character(m$subtype) %in% subtypes[[neuron_type]])) {
    stop('Unexpected subtype annotation: ', neuron_type)
  }
  sample_conditions <- setNames(samples$condition, samples$sample)
  m <- m[m$condition %in% c('Ext', 'Ret-Ext'), , drop = FALSE]
  if (any(!m$sample %in% samples$sample) ||
      any(as.character(m$condition) != sample_conditions[as.character(m$sample)])) {
    stop('Metadata sample and condition disagree with sample table: ', neuron_type)
  }
  m
}

make_pool_counts <- function(m, neuron_type) {
  ans <- list()
  analysis <- if (neuron_type == 'Excitatory') {
    'excitatory_IEGhigh_composition'
  } else 'inhibitory_IEGhigh_fraction'
  for (subtype in subtypes[[neuron_type]]) {
    for (j in seq_len(nrow(samples))) {
      same_pool <- as.character(m$sample) == samples$sample[j]
      same_subtype <- as.character(m$subtype) == subtype
      success <- sum(same_pool & same_subtype & m$IEGhigh)
      total <- if (neuron_type == 'Excitatory') {
        sum(same_pool & m$IEGhigh)
      } else {
        sum(same_pool & same_subtype)
      }
      proportion <- if (total > 0L) success / total else NA_real_
      ans[[length(ans) + 1L]] <- data.frame(
        analysis = analysis, neuron_type = neuron_type, subtype = subtype,
        sample = samples$sample[j], condition = samples$condition[j],
        success = success, total = total, proportion = proportion,
        percent = 100 * proportion, excluded = total == 0L,
        exclusion_reason = if (total == 0L) 'zero_denominator' else '',
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, ans)
}

counts <- do.call(rbind, lapply(names(subtypes), function(neuron_type) {
  make_pool_counts(read_neuron_metadata(neuron_type), neuron_type)
}))
write_table(counts, file.path(out, 'pool_counts.csv'))

models <- list()
fit_one <- function(d) {
  result <- data.frame(
    analysis = d$analysis[1], neuron_type = d$neuron_type[1], subtype = d$subtype[1],
    contrast = 'Ret-Ext vs Ext', denominator = if (d$neuron_type[1] == 'Excitatory') {
      'All IEG-high excitatory nuclei within the same pooled sample'
    } else 'All nuclei of the same inhibitory subtype within the same pooled sample',
    planned_tests = length(subtypes[[d$neuron_type[1]]]),
    pools_Ext = sum(d$condition == 'Ext' & !d$excluded),
    pools_RetExt = sum(d$condition == 'Ret-Ext' & !d$excluded),
    coefficient = NA_real_, SE = NA_real_, z = NA_real_, pvalue = NA_real_,
    padj = NA_real_, odds_ratio = NA_real_, OR_CI_low = NA_real_,
    OR_CI_high = NA_real_, random_intercept_SD = NA_real_,
    singular = NA, converged = NA, status = 'not_fitted',
    warning = '', error = '', stringsAsFactors = FALSE
  )
  d <- d[!d$excluded, , drop = FALSE]
  if (any(d$success < 0 | d$success > d$total) || any(d$total < 0)) {
    result$status <- 'invalid_counts'
    result$error <- 'Success counts must be between zero and the denominator'
    return(result)
  }
  if (nrow(d) < 4L || any(c(result$pools_Ext, result$pools_RetExt) < 2L)) {
    result$status <- 'insufficient_pools'
    return(result)
  }
  if (sum(d$success) == 0L || sum(d$total - d$success) == 0L) {
    result$status <- 'no_binary_outcome_variation'
    return(result)
  }
  d$condition <- factor(d$condition, levels = c('Ext', 'Ret-Ext'))
  d$sample <- factor(d$sample)
  messages <- character()
  fit <- tryCatch(withCallingHandlers(
    lme4::glmer(
      cbind(success, total - success) ~ condition + (1 | sample),
      data = d, family = stats::binomial(link = 'logit'),
      nAGQ = as.integer(settings$nAGQ),
      # Each grouped row represents many Bernoulli nuclei from one pooled sample.
      control = lme4::glmerControl(optimizer = settings$optimizer,
                                   check.nobs.vs.nlev = 'ignore',
                                   check.nobs.vs.nRE = 'ignore')
    ),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart('muffleWarning')
    }
  ), error = function(e) e)
  result$warning <- paste(unique(messages), collapse = ' | ')
  if (inherits(fit, 'error')) {
    result$status <- 'fit_error'
    result$error <- conditionMessage(fit)
    return(result)
  }
  models[[paste(result$analysis, result$subtype, sep = '__')]] <<- fit
  convergence_messages <- fit@optinfo$conv$lme4$messages
  optimizer_code <- fit@optinfo$conv$opt
  numerical_messages <- convergence_messages[
    !grepl('boundary \\(singular\\)', convergence_messages)
  ]
  result$converged <- (is.null(optimizer_code) || all(optimizer_code == 0)) &&
    !length(numerical_messages)
  result$singular <- lme4::isSingular(fit, tol = 1e-4)
  if (length(convergence_messages)) {
    result$warning <- paste(c(result$warning, convergence_messages), collapse = ' | ')
  }
  result$random_intercept_SD <- attr(lme4::VarCorr(fit)$sample, 'stddev')[1]
  coefficients <- summary(fit)$coefficients
  term <- setdiff(rownames(coefficients), '(Intercept)')
  if (length(term) != 1L || !is.finite(coefficients[term, 'Std. Error'])) {
    result$status <- 'coefficient_not_estimable'
    return(result)
  }
  result$coefficient <- coefficients[term, 'Estimate']
  result$SE <- coefficients[term, 'Std. Error']
  result$z <- coefficients[term, 'z value']
  result$odds_ratio <- exp(result$coefficient)
  result$OR_CI_low <- exp(result$coefficient - stats::qnorm(0.975) * result$SE)
  result$OR_CI_high <- exp(result$coefficient + stats::qnorm(0.975) * result$SE)
  if (!isTRUE(result$converged)) {
    result$status <- 'nonconverged'
  } else {
    result$pvalue <- coefficients[term, 'Pr(>|z|)']
    result$status <- if (result$singular) 'fitted_singular' else 'fitted'
  }
  result
}

statistics <- do.call(rbind, lapply(split(counts, interaction(counts$analysis, counts$subtype,
                                                            drop = TRUE)), fit_one))
for (analysis in unique(statistics$analysis)) {
  i <- statistics$analysis == analysis
  statistics$padj[i] <- stats::p.adjust(statistics$pvalue[i], method = 'BH',
                                       n = unique(statistics$planned_tests[i]))
}
write_table(statistics, file.path(out, 'model_statistics.csv'))
save_object(models, file.path(out, 'models.rds'))
write_table(statistics[, c('analysis', 'subtype', 'pools_Ext', 'pools_RetExt',
                           'status', 'converged', 'singular', 'warning', 'error')],
            file.path(out, 'model_diagnostics.csv'))
write_session(file.path(out, 'sessionInfo.txt'))
