#' @name lambdaControl
#' @aliases selectControl
#' 
#' @title Control of Penalization Parameters and Selection Criteria
#' 
#' Constructs control structures for penalized mixed model fitting.
#' 
#' @inheritParams fit_dat_B
#' @param lambda0_seq a range of non-negative numeric penalty parameter for the fixed 
#' effects parameters. If \code{NULL}, then a range will be calculated as described by (...).
#' @param lambda1_seq a range of non-negative numeric penalty parameter for the grouped random 
#' effects covariance parameters. If \code{NULL}, then a range will be calculated as described 
#' by (...).
#' @param nlambda positive integer specifying number of lambda with which to fit a model.
#' Ignored if \code{lambda0_seq} and \code{lambda1_seq} are specified by the user.
#' @param BIC_option character value specifing the selection criteria used to select the 'best' model.
#' Default "BICh" utilizes the hybrid BIC value described in Delattre, Lavielle, and Poursat (2014).
#' The regular "BIC" option penalty term uses (total non-zero coefficients)*(length(y) = total number
#' observations). The "BICq" option specifies the BIC ICQ criterion, which requires a fit of 
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
selectControl = function(lambda0_seq = NULL, lambda1_seq = NULL, nlambda = 20,
                         BIC_option = c("BICh","BIC","BICq")){
  
  if(length(BIC_option) > 1){
    BIC_option = BIC_option[1]
  }
  if(!(BIC_option %in% c("BICh","BIC","BICq"))){
    stop("BIC_option must be 'BICh', 'BIC', or 'BICq'")
  }
  
  structure(list(lambda0_seq = lambda0_seq,
                 lambda1_seq = lambda1_seq,
                 nlambda = nlambda, BIC_option = BIC_option),
            class = c("selectControl", "pglmmControl"))
}

#' @title Control of Metropolis-within-Gibbs Adaptive Random Walk Sampling Procedure
#' 
#' Constructs control structure for the adaptive random walk Metropolis-within-Gibbs sampling procedure.
#' 
#' @param nMC_burnin positive integer specifying teh number of posterior draws to use as
#' burnin for each E step in the EM algorithm. Default 200. Function will not allow \code{nMC_burnin}
#' to be less than 100.
#' @param batch_length positive integer specifying the number of posterior draws to collect
#' before the proposal variance is adjusted based on the acceptance rate of the last 
#' \code{batch_length} posterior draws
#' @param offset positive integer specifying an offset value for the increment of the proposal
#' variance adjustment. Used to ensure the required diminishing adaptation condition. Default 8500.
#' \verb{
#' increment = min(0.01, 1 / sqrt(batch*batch_length + offset))
#' } 
#' @param burnin_batchnum positive integer specifying the maximum number of batches (of length
#' \code{batch_length}) during which the proposal variance can be adjusted. Default 1,000.
#' 
#' @return The *Control functions return a list (inheriting from class "\code{adaptControl}") 
#' containing parameter specifications for the adaptive random walk sampling procedure.
#' 
#' @export
adaptControl = function(batch_length = 100, offset = 0, burnin_batchnum = 500){
  
  if(batch_length < 10){
    stop("batch_length must be at least 10")
  }
  if(burnin_batchnum < 0){
    stop("burnin_batchnum cannot be negative")
  }
  x = c(batch_length, burnin_batchnum)
  if(!all(floor(x)==x)){ # if any of the above values not integers
    stop("batch_length and burnin_batchnum must be integers")
  }
  
  structure(list(batch_length = batch_length, offset = offset, 
                 burnin_batchnum = burnin_batchnum), 
            class = "adaptControl")
}

#' @title Control of Penalized Generalized Linear Mixed Model Fitting
#' 
#' Constructs the control structure for penalized mixed model fitting.
#' 
#' @inheritParams fit_dat_B
#' @param nMC_start a positive integer for the initial number of Monte Carlo draws. 
#' 
#' @return The *Control functions return a list (inheriting from class "\code{optimControl}") 
#' containing fit and optimization criteria values used in optimization routine.
#' 
#' @export
optimControl = function(conv_EM = 0.001, conv_CD = 0.0001, 
                        nMC_burnin = 500, nMC_start = 5000, nMC_max = 10000, nMC_report = 5000,
                        maxitEM = 100, maxit_CD = 200, M = 10000, t = 2,
                        covar = c("unstructured","independent"),
                        sampler = c("stan","random_walk","independence"), 
                        var_start = 0.5, max_cores = 1){
  
  # Acceptable input types and input restrictions - vectors, integers, positive numbers ...
  if(sum(c(nMC_start, nMC_max, maxitEM) %% 1) > 0 | sum(c(nMC_start, nMC_max, maxitEM) <= 0) > 0){
    stop("nMC_start, nMC_max, and maxitEM must be positive integers")
  }
  # Checks
  if(nMC_burnin < 100){
    warning("nMC_burnin not allowed to be less than 100. Value set to 100", immediate. = T)
    nMC_burnin = 100
  }
  if(nMC_max < nMC_start){
    stop("nMC_max should not be smaller than nMC_start \n")
  }
  x = c(nMC_burnin, nMC_start, nMC_max, nMC_report)
  if(!all(floor(x)==x)){ # if any of the above values not integers
    stop("all nMC arguments must be integers")
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

