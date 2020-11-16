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
#' small portion of max, range is composed of \code{nlambda} values spread evenly on the log scale.
#' For secondary rounds run using \code{\link{glmmPen_FineSearch}}, the optimal lambda combination 
#' from the initial round (based on the specified \code{BIC_option}) will be used to calculate a 
#' finer grid search. The new max and min lambda values are the lambda values two 
#' positions away from the best lambda in the original lambda sequences for the fixed and random
#' effects.
#' 
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
#' Default "BICh" utilizes the hybrid BIC value described in Delattre, Lavielle, and Poursat (2014).
#' The regular "BIC" option penalty term uses (total non-zero coefficients)*(length(y) = total number
#' observations). The "BICq" option specifies the BIC-ICQ criterion, which requires a fit of 
#' a full model (a small penalty is used for the fixed and random effects). 
#' 
#' @return The *Control functions return a list (inheriting from class "\code{pglmmControl}") 
#' containing penalization parameter values, presented either as an individual set or as a range of
#' possible values.
#' 
#' @export
lambdaControl = function(lambda0 = 0, lambda1 = 0){
  structure(list(lambda0 = lambda0, lambda1 = lambda1), 
            class = c("lambdaControl","pglmmControl"))
}

#' @rdname lambdaControl
#' @export
selectControl = function(lambda0_seq = NULL, lambda1_seq = NULL, nlambda = 10,
                         BIC_option = c("BICh","BIC","BICq")){
  
  # Perform input checks
  
  if(length(BIC_option) > 1){
    BIC_option = BIC_option[1]
  }
  if(!(BIC_option %in% c("BICh","BIC","BICq"))){
    stop("BIC_option must be 'BICh', 'BIC', or 'BICq'")
  }
  if(BIC_option != "BICq"){
    warning("BIC-ICQ will not be calculated since 'BICq' not specified for BIC_option", 
            immediate. = T)
  }
  
  if(!is.null(lambda0_seq)){
    if(!is.numeric(lambda0_seq)){
      stop("lambda0_seq must be numeric")
    }
    if(any(lambda0_seq) < 0){
      stop("lambda0_seq cannot be negative")
    }
  }
  if(!is.null(lambda1_seq)){
    if(!is.numeric(lambda1_seq)){
      stop("lambda1_seq must be numeric")
    }
    if(any(lambda1_seq < 0)){
      stop("lambda1_seq cannot be negative")
    }
  }
  
  # if(!all(floor(x)==x)){ # if any of the above values not integers
  #   stop("batch_length and burnin_batchnum must be integers")
  # }
  if(!(floor(nlambda) == nlambda)){
    stop("nlambda must be an integer")
  }
  if(nlambda < 2){
    stop("nlambda must be at least 2")
  }
  
  structure(list(lambda0_seq = lambda0_seq,
                 lambda1_seq = lambda1_seq,
                 nlambda = nlambda, BIC_option = BIC_option),
            class = c("selectControl", "pglmmControl"))
}

#' @title Control of Metropolis-within-Gibbs Adaptive Random Walk Sampling Procedure
#' 
#' Controls the adaptive random walk Metropolis-within-Gibbs sampling procedure.
#' 
#' @param batch_length positive integer specifying the number of posterior draws to collect
#' before the proposal variance is adjusted based on the acceptance rate of the last 
#' \code{batch_length} accepted posterior draws
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
  
  if(batch_length < 10){
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
#' EM algorithm. Default is 0.001.
#' @param conv_CD a non-negative numeric convergence criteria for the convergence of the 
#' grouped coordinate descent loop within the M step of the EM algorithm. Default 0.0001.
#' @param nMC_burnin positive integer specifying the number of posterior draws to use as
#' burnin for each E step in the EM algorithm. Default 250. Function will not allow \code{nMC_burnin}
#' to be less than 100.
#' @param nMC_start a positive integer for the initial number of Monte Carlo draws. 
#' @param nMC_max a positive integer for the maximum number of allowed Monte Carlo draws
#' @param nMC_report positive integer specifying number of Monte Carlo posterior draws to return
#' when the function ends. Default set to the minimum allowed value: 5,000. 
#' Warning: the returned draws are formatted as a regular
#' matrix (not a big.matrix). Therefore, depending on the number of random effect covariates (q)
#' and the number of groups (d), choose \code{nMC_report} such that a matrix of size 
#' \code{nMC_report} by (q*d) does not cause memory issues on your operating system.
#' @param maxitEM a positive integer for the maximum number of allowed EM iterations. Default equals 100.
#' @param maxit_CD a positive integer for the maximum number of allowed interations for the
#' coordinate descent algorithms used within each EM iteration. Default equals 250.
#' @param M positive integer specifying the number of posterior draws to use within the 
#' Pajor log-likelihood calculation
#' @param t the convergence criteria is based on the average Euclidean distance between 
#' the most recent coefficient estimate and the coefficient estimate from t EM iterations back.
#' Positive integer, default equals 2.
#' @param sampler character string specifying whether the posterior draws of the random effects
#' should be drawn using Stan (default, from package rstan) or the Metropolis-within-Gibbs procedure 
#' incorporating an adaptive random walk sampler ("random_walk") or an
#' independence sampler ("independence"). 
#' @param covar character string specifying whether the covariance matrix should be unstructured
#' ("unstructured") or diagonal with no covariances between variables ("independent").
#' Default is "unstructured", but if the number of random effects (including the intercept) is 
#' greater than or equal to 7 (i.e. high dimensional), then the algorithm automatically assumes an 
#' independent covariance structure (covar switched to "independent").
#' @param var_start positive number specifying the starting values to initialize the variance
#' of the covariance matrix. Default set to 1.0. 
#' 
#' @return Function returns a list (inheriting from class "\code{optimControl}") 
#' containing fit and optimization criteria values used in optimization routine.
#' 
#' @export
optimControl = function(conv_EM = 0.001, conv_CD = 0.0001, 
                        nMC_burnin = 250, nMC_start = 1000, nMC_max = 10^4, nMC_report = 5000,
                        maxitEM = 100, maxit_CD = 200, M = 10000, t = 2,
                        covar = c("unstructured","independent"),
                        sampler = c("stan","random_walk","independence"), 
                        var_start = 1.0, max_cores = 1){
  
  # Acceptable input types and input restrictions - vectors, integers, positive numbers ...
  
  x = c(nMC_burnin, nMC_start, nMC_max, nMC_report, maxitEM, M)
  if(!all(floor(x)==x)){ # if any of the above values not integers
    stop("all nMC arguments must be integers")
  }
  if(sum(x <= 0) > 0){
    stop("nMC_start, nMC_max, and maxitEM must be positive integers")
  }
  # Checks
  if(M < 10^4){
    stop("M must be greater than or equal to 10^4")
  }
  if(nMC_report < 5000){
    stop("nMC_report must be greater than or equal to 5000")
  }
  if(nMC_burnin < 100){
    warning("nMC_burnin not allowed to be less than 100. Value set to 100", immediate. = T)
    nMC_burnin = 100
  }
  if(nMC_max < nMC_start){
    stop("nMC_max should not be smaller than nMC_start \n")
  }
  
  if(var_start <= 0){
    stop("var_start must be a positive number")
  }
  if((max_cores < 1) | (max_cores) %% 1 != 0){
    stop("max_cores must be a positive integer")
  }
  
  if(length(covar) > 1){
    covar = covar[1]
  }
  if(!(covar %in% c("unstructured","independent"))){
    stop("covariance structure, 'covar', must be either 'unstructured' or 'independent'")
  }
  
  if(length(sampler) > 1){
    sampler = sampler[1]
  }
  if(!(sampler %in% c("stan","random_walk","independence"))){
    stop("sampler must be either 'stan', 'random_walk', or 'independence'")
  }
  
  structure(list(conv_EM = conv_EM, conv_CD = conv_CD, 
                 nMC_burnin = nMC_burnin, nMC = nMC_start, nMC_max = nMC_max, nMC_report = nMC_report, 
                 maxitEM = maxitEM, maxit_CD = maxit_CD,  M = M, t = t, covar = covar,
                 sampler = sampler, var_start = var_start,
                 max_cores = max_cores),
            class = "optimControl")
}

