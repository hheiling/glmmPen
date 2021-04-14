#' @name lambdaControl
#' @aliases selectControl
#' 
#' @title Control of Penalization Parameters and Selection Criteria
#' 
#' @description Constructs control structures for penalized mixed model fitting.
#' 
#' @details If unspecified, the \code{lambda0_seq} and \code{lambda1_seq} numeric sequences are 
#' automatically calculated. For an initial
#' grid search using \code{\link{glmmPen}}, the range will be calculated in the same manner as 
#' \code{ncvreg} calculates the range: max penalizes all fixed and random effects to 0, min is a 
#' small portion of max (\code{lambda.min}*(lambda max)), range is composed of 
#' \code{nlambda} values spread evenly on the log scale. Unlike \code{ncvreg}, the order of penalty
#' values used in the algorithm must run from the min lambda to the max lambda (as opposed to 
#' running from max lambda to min lambda).
#' For secondary rounds run using \code{\link{glmmPen_FineSearch}}, the optimal lambda combination 
#' from the initial round (based on the specified \code{BIC_option}) will be used to calculate a 
#' finer grid search. The new max and min lambda values are the lambda values \code{idx_lambda} 
#' positions away from the best lambda in the original lambda sequences for the fixed and random
#' effects.
#' 
#' The \code{lambda0} and \code{lambda1} arguments allow for a user to fit a model with a single 
#' non-zero penalty parameter combination. However, this is generally not recommended.  
#' 
#' @param lambda0 a non-negative numeric penalty parameter for the fixed effects parameters
#' @param lambda1 a non-negative numeric penalty parameter for the (grouped) random effects
#' covariance parameters
#' @param lambda0_seq,lambda1_seq a range of non-negative numeric penalty parameters for the fixed 
#' and random effect parameters, respectively. If \code{NULL}, then a range will be automatically 
#' calculated. See 'Details' section for more details on these default calculations.
#' @param nlambda positive integer specifying number of lambda with which to fit a model.
#' Ignored if \code{lambda0_seq} and \code{lambda1_seq} are specified by the user.
#' @param BIC_option character string specifing the selection criteria used to select the 'best' model.
#' Default "BICq" option specifies the BIC-ICQ criterion, which requires a fit of 
#' a full model; a small penalty of 0.001*(lambda max) is used for the fixed and/or random effects if the number of 
#' effects is 10 or greater. The "BICh" option utilizes the hybrid BIC value described in Delattre, Lavielle, and Poursat (2014).
#' The regular "BIC" option penalty term uses (total non-zero coefficients)*(length(y) = total number
#' observations). The "BICNgrp" option penalty term uses (total non-zero coefficients)*(nlevels(group) = number
#' groups).
#' @param logLik_calc logical value specifying if the log likelihood (and log-likelihood based 
#' calculations BIC, BICh, and BICNgrp) should be calculated for all of the models in the selection procedure. 
#' If BIC-ICQ is used for selection, the log-likelihood is not needed. However, if users are interested
#' in comparing the best models from BIC-ICQ and other BIC-type selection criteria, setting
#' \code{logLik_calc} to \code{TRUE} will calculate these other options for all of the models.
#' @param pre_screen logical value indicating if pre-screening of random effects should be performed
#' during selection when the number random effect predictors is 10 or more. Default \code{TRUE}. 
#' Pre-screening not performed when the number of random effect predictors is 9 or less.
#' @param lambda.min numeric fraction between 0 and 1. The sequence of the lambda penalty parameters
#' ranges from the maximum lambda where all fixed and random effects are penalized to 0 and 
#' a minimum lambda value, which equals a small fraction of the maximum lambda. The parameter 
#' \code{lambda.min} specifiex this fraction. Default value is set to \code{NULL}, which
#' automatically selects \code{lambda.min} to equal 0.01 when n < p and 0.05 when p >= n.
#' 
#' @return The *Control functions return a list (inheriting from class "\code{pglmmControl}") 
#' containing penalization parameter values, presented either as an individual set or as a range of
#' possible values.
#' 
#' @export
lambdaControl = function(lambda0 = 0, lambda1 = 0){
  if(!is.numeric(c(lambda0,lambda1))){
    stop("lambda0 and lambda1 must be numeric values")
  }
  if((lambda0 < 0) | (lambda1 < 0)){
    stop("lambda0 and lambda1 must be postive numeric values")
  }
  structure(list(lambda0 = lambda0, lambda1 = lambda1), 
            class = c("lambdaControl","pglmmControl"))
}

#' @rdname lambdaControl
#' @export
selectControl = function(lambda0_seq = NULL, lambda1_seq = NULL, nlambda = 10,
                         BIC_option = c("BICq","BICh","BIC","BICNgrp"), 
                         logLik_calc = switch(BIC_option[1], BICq = F, T), 
                         lambda.min = NULL, pre_screen = T){
  
  # Perform input checks
  
  if(length(BIC_option) > 1){
    BIC_option = BIC_option[1]
  }
  if(!(BIC_option %in% c("BICh","BIC","BICq","BICNgrp"))){
    stop("BIC_option must be 'BICq', 'BICh', 'BIC', or 'BICNgrp'")
  }
  
  if(!is.null(lambda0_seq)){
    if(!is.numeric(lambda0_seq)){
      stop("lambda0_seq must be numeric")
    }
    if(any(lambda0_seq) < 0){
      stop("lambda0_seq cannot be negative")
    }
    for(i in 2:length(lambda0_seq)){
      if(lambda0_seq[i-1] > lambda0_seq[i]){
        stop("lambda0_seq must be a sequence of ascending order (min value to max value)")
      }
    }
  }
  if(!is.null(lambda1_seq)){
    if(!is.numeric(lambda1_seq)){
      stop("lambda1_seq must be numeric")
    }
    if(any(lambda1_seq < 0)){
      stop("lambda1_seq cannot be negative")
    }
    for(i in 2:length(lambda1_seq)){
      if(lambda1_seq[i-1] > lambda1_seq[i]){
        stop("lambda1_seq must be a sequence of ascending order (min value to max value)")
      }
    }
  }
  
  if(!(floor(nlambda) == nlambda)){
    stop("nlambda must be an integer")
  }
  if(nlambda < 2){
    stop("nlambda must be at least 2")
  }
  
  if (!is.logical(logLik_calc)) {
    stop("'logLik_calc' must be a logical value (T or F)")
  }
  
  if (!is.logical(pre_screen)) {
    stop("'pre_screen' must be a logical value (T or F)")
  }
  
  if(!is.null(lambda.min)){
    if(!is.numeric(lambda.min)){stop("lambda.min must be numeric")}
    if((lambda.min >= 1) | (lambda.min <=0 )){stop("lambda.min must be a fraction between 0 and 1")}
  }
  
  structure(list(lambda0_seq = lambda0_seq, lambda1_seq = lambda1_seq,
                 nlambda = nlambda, BIC_option = BIC_option, logLik_calc = logLik_calc,
                 lambda.min = lambda.min, pre_screen = pre_screen),
            class = c("selectControl", "pglmmControl"))
}

#' @title Control of Metropolis-within-Gibbs Adaptive Random Walk Sampling Procedure
#' 
#' Controls the adaptive random walk Metropolis-within-Gibbs sampling procedure.
#' 
#' @param batch_length positive integer specifying the number of posterior draws to collect
#' before the proposal variance is adjusted based on the acceptance rate of the last 
#' \code{batch_length} accepted posterior draws. Default is set to 100. Batch length restricted
#' to be no less than 50.
#' @param offset non-negative integer specifying an offset value for the increment of the proposal
#' variance adjustment. Optionally used to ensure the required diminishing adaptation condition. 
#' Default set to 0.
#' \verb{
#' increment = min(0.01, 1 / sqrt(batch*batch_length + offset))
#' } 
#' 
#' @return Function returns a list (inheriting from class "\code{adaptControl}") 
#' containing parameter specifications for the adaptive random walk sampling procedure.
#' 
#' @export
adaptControl = function(batch_length = 100, offset = 0){
  
  if(batch_length < 50){
    stop("batch_length must be at least 10")
  }
  if(floor(batch_length)!=batch_length){
    stop("batch_length must be integer")
  }
  
  structure(list(batch_length = batch_length, offset = offset), 
            class = "adaptControl")
}

#' @title Control of Penalized Generalized Linear Mixed Model Fitting
#' 
#' Constructs the control structure for the optimization of the penalized mixed model fit algorithm.
#' 
#' @param conv_EM a non-negative numeric convergence criteria for the convergence of the 
#' EM algorithm. Default is 0.001. EM algorithm is considered to have converge if the average Euclidean 
#' distance between the current coefficient estimates and the coefficient estimates from 
#' \code{t} EM iterations back is less than \code{conv_EM} \code{mcc} times in a row.
#' See \code{t} and \code{mcc} for more details.
#' @param conv_CD a non-negative numeric convergence criteria for the convergence of the 
#' grouped coordinate descent loop within the M step of the EM algorithm. Default 0.0001.
#' @param nMC_burnin positive integer specifying the number of posterior draws to use as
#' burnin for each E step in the EM algorithm. Default 250 when the number of random effects 
#' predictors is less than or equal to 10; default 100 otherwise. Function will not allow \code{nMC_burnin}
#' to be less than 100. 
#' @param nMC_start a positive integer for the initial number of Monte Carlo draws. Default 250 when 
#' the number of random effects predictors is less than 25; default 100 otherwise.
#' @param nMC_max a positive integer for the maximum number of allowed Monte Carlo draws used
#' in each step of the EM algorithm. When the number of random effect predictors is 10 or less, 
#' Default is set to 5000 when no selection is performed and 2500 when selection is performed.
#' Default is set to 1000 when the number of random effect predictors is between 11 and 24;
#' Default is set to 500 when the number of random effect predictors is 25 or more.
#' @param maxitEM a positive integer for the maximum number of allowed EM iterations. 
#' Default equals 50 for the Binomial and Poisson families, 100 for the Gaussian family.
#' @param maxit_CD a positive integer for the maximum number of allowed interations for the
#' coordinate descent algorithms used within the M-step of each EM iteration. Default equals 200.
#' @param M positive integer specifying the number of posterior draws to use within the 
#' Pajor log-likelihood calculation. Default is 10^4; minimum allowed value is 5000.
#' @param t the convergence criteria is based on the average Euclidean distance between 
#' the most recent coefficient estimates and the coefficient estimates from t EM iterations back.
#' Positive integer, default equals 2.
#' @param mcc the number of times the convergence critera must be met before the algorithm is
#' seen as having converged (mcc for 'meet condition counter'). Default set to 2. Value retricted 
#' to be no less than 2.
#' @param sampler character string specifying whether the posterior draws of the random effects
#' should be drawn using Stan (default, from package rstan) or the Metropolis-within-Gibbs procedure 
#' incorporating an adaptive random walk sampler ("random_walk") or an
#' independence sampler ("independence"). If using the random walk sampler, see \code{\link{adaptControl}}
#' for some additional control structure parameters.
#' @param var_start either the character string "recommend" or a positive number specifying the 
#' starting values to initialize the variance of the covariance matrix. Default "recommend" first
#' fits a simple model with a fixed and random intercept only using a Laplace estimate. The 
#' random intercept variance estimate from this model is then multiplied by 2 and used as the 
#' starting variance. 
#' @param max_cores integer describing the number of cores available for computation. If 
#' \code{max_cores} is specified to be greater than 1 and the sampler is specified as "stan", then 
#' parallel computation using multiple cores is used to calculate the Stan MCMC samples within
#' each E step.
#' 
#' @details When the \code{optim_options} arugment in \code{\link{glmm}} and \code{\link{glmmPen}}
#' is set to "recommend", the default settings discussed in the given \code{optimControl} arguments are
#' used. These default settings depend on both the family of the data structure and the number 
#' of random effects predictors specified for use. If the user specifies 
#' \code{optim_options = optimControl()} with any argument specifications, no additional
#' adjustments will be performed on the arguments based on family or random effect predictors. 
#' 
#' @return Function returns a list (inheriting from class "\code{optimControl}") 
#' containing fit and optimization criteria values used in optimization routine.
#' 
#' @export
optimControl = function(conv_EM = 0.001, conv_CD = 0.0001, 
                        nMC_burnin = 250, nMC_start = 250, nMC_max = 5000, nMC_report = 5000,
                        maxitEM = 50, maxit_CD = 100, 
                        M = 10000, t = 2, mcc = 2,
                        sampler = c("stan","random_walk","independence"), 
                        var_start = "recommend", max_cores = 1){
  
  # Acceptable input types and input restrictions - vectors, integers, positive numbers ...
  
  x = c(nMC_burnin, nMC_start, nMC_max, nMC_report, maxitEM, maxit_CD, M, t, mcc)
  if((!all(floor(x)==x)) | (sum(x <= 0) > 0)){ # if any of the above values not positive integers
    stop("the nMC arguments, maxit arguments, M, t, and mcc must be positive integers")
  }
  
  # More restrictive checks
  # if(M < 10^4){
  #   stop("M must be greater than or equal to 10^4")
  # }
  if(M < 5000){
    stop("M must be greater than or equal to 5000")
  }
  if(nMC_burnin < 100){
    warning("nMC_burnin not allowed to be less than 100. Value set to 100", immediate. = T)
    nMC_burnin = 100
  }
  if(nMC_max < nMC_start){
    stop("nMC_max cannot be smaller than nMC_start")
  }
  
  if(is.character(var_start)){
    if(var_start != "recommend"){
      stop("var_start must either be 'recommend' or a positive numeric constant")
    }
  }else{
    if(var_start <= 0){
      stop("var_start must be a positive numeric value")
    }
  }
  
  if((max_cores < 1) | (max_cores) %% 1 != 0){
    stop("max_cores must be a positive integer")
  }
  
  sampler = checkSampler(sampler)
  
  if(mcc < 2){
    stop("mcc must be at least 2")
  }
  
  structure(list(conv_EM = conv_EM, conv_CD = conv_CD, 
                 nMC_burnin = nMC_burnin, nMC = nMC_start, nMC_max = nMC_max, nMC_report = nMC_report,
                 maxitEM = maxitEM, maxit_CD = maxit_CD,  M = M, t = t, mcc = mcc,
                 sampler = sampler, var_start = var_start,
                 max_cores = max_cores),
            class = "optimControl")
}


# q = number of random effects (including random intercept)
# select: TRUE if running the selection algorithm
optim_recommend = function(family, q, select){
  
  optim_options = optimControl()
  # Default options for family = 'binomial' or 'poisson'
  if (family == "gaussian") {
    # increase maxitEM
    optim_options = optimControl()
    optim_options$maxitEM = 100
  }
  
  if (q <= 11) {
    # For selection, decrease nMC_max to improve time
    ## Note: During selection, will initialize with good starting points, so
    ## larger nMC_max not needed for convergence
    ## In order to speed up algorithm, will decrease nMC_max while increasing maxitEM
    # Otherwise, keep everything else the same
    if (select == T) {
      optim_options$nMC_max = 2500
      # optim_options$maxitEM = optim_options$maxitEM + 15
    }
  }else{
    # Decrease nMC
    optim_options$nMC_burnin = 100
    optim_options$nMC = 100
    optim_options$nMC_max = 1500
    optim_options$maxitEM = optim_options$maxitEM + 15
  }
  
  return(optim_options)
  
}