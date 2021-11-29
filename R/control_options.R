#' @name lambdaControl
#' @aliases selectControl
#' 
#' @title Control of Penalization Parameters and Selection Criteria
#' 
#' @description Constructs control structures for penalized mixed model fitting.
#' 
#' @details If unspecified, the \code{lambda0_seq} and \code{lambda1_seq} numeric sequences are 
#' automatically calculated. The sequence will be calculated in the same manner as 
#' \code{ncvreg} calculates the range: max penalizes all fixed and random effects to 0, min is a 
#' small portion of max (\code{lambda.min}*(lambda max)), sequence is composed of 
#' \code{nlambda} values spread evenly on the log scale. Unlike \code{ncvreg}, the order of penalty
#' values used in the algorithm must run from the min lambda to the max lambda (as opposed to 
#' running from max lambda to min lambda). The length of the sequence is specified by \code{nlambda}. 
#' By default, these sequences are calculated using \code{\link{LambdaSeq}}.
#' 
#' The \code{lambda0} and \code{lambda1} arguments allow for a user to fit a model with a single 
#' non-zero penalty parameter combination. However, this is generally not recommended. 
#' 
#' Abbreviated grid search: The abbreviated grid search proceeds in two stages. In stage 1, the
#' algorithm fits the following series of models: the fixed effects penalty parameter remains a
#' fixed value evaluated at the minimum of the fixed effects penalty parameters, and a all
#' random effects penalty parameters are examined. The 'best' model from this first stage of models
#' determines the optimum random effect penalty parameter. In stage 2, the algorithm fits the 
#' following series of models: the random effects penalty parameter remains fixed at the value of
#' the optimum random effect penalty parameter (from stage 1) and all fixed effects penalty
#' parameters are considered. The best overall model is the best model from stage 2. This reduces the 
#' number of models considered to length(`lambda0_seq`) + length(`lambda1_seq`). The authors found
#' that this abbreviated grid search worked well in simulations.
#' 
#' @param lambda0 a non-negative numeric penalty parameter for the fixed effects parameters
#' @param lambda1 a non-negative numeric penalty parameter for the (grouped) random effects
#' covariance parameters
#' @param lambda0_seq,lambda1_seq a sequence of non-negative numeric penalty parameters for the fixed 
#' and random effect parameters, respectively. If \code{NULL}, then a sequence will be automatically 
#' calculated. See 'Details' section for more details on these default calculations. 
#' @param nlambda positive integer specifying number of penalty parameters (lambda) 
#' to use for the fixed and random effects penalty parameters. Default set to 10.
#' Ignored if \code{lambda0_seq} and \code{lambda1_seq} are specified by the user.
#' @param BIC_option character string specifing the selection criteria used to select the 'best' model.
#' Default "BICq" option specifies the BIC-ICQ criterion, which requires a fit of 
#' a full model; a small penalty (the minimum of the penalty sequence) 
#' is used for the fixed and random effects. 
#' The "BICh" option utilizes the hybrid BIC value described in Delattre, Lavielle, and Poursat (2014).
#' The regular "BIC" option penalty term uses (total non-zero coefficients)*(length(y) = total number
#' observations). The "BICNgrp" option penalty term uses (total non-zero coefficients)*(nlevels(group) = number
#' groups).
#' @param search character string of "abbrev" (default) or "full_grid" indicating if the search of models over 
#' the penalty parameter space should be the full grid search (total number of models equals
#' `nlambda`^2 or length(`lambda0_seq`)*length(`lambda1_seq`)) or an abbreviated grid search.
#' The abbreviated grid search is described in more detail in the Details section. Te authors
#' highly recommend the abbreviated grid search.
#' @param logLik_calc logical value specifying if the log likelihood (and log-likelihood based 
#' calculations BIC, BICh, and BICNgrp) should be calculated for all of the models in the selection procedure. 
#' If BIC-ICQ is used for selection, the log-likelihood is not needed for each model. 
#' However, if users are interested
#' in comparing the best models from BIC-ICQ and other BIC-type selection criteria, setting
#' \code{logLik_calc} to \code{TRUE} will calculate these other quantities for all of the models.
#' @param lambda.min numeric fraction between 0 and 1. The sequence of the lambda penalty parameters
#' ranges from the maximum lambda where all fixed and random effects are penalized to 0 and 
#' a minimum lambda value, which equals a small fraction of the maximum lambda. The parameter 
#' \code{lambda.min} specifies this fraction. Default value is set to \code{NULL}, which
#' automatically selects \code{lambda.min} to equal 0.01 when p <= 10 and 0.05 when p > 10.
#' @param pre_screen logical value indicating whether pre-screening should be performed before
#' model selection (default \code{TRUE}). If the number of random effects considered less than 5,
#' no pre-screening will be performed. Pre-screening removes random effects from consideration
#' during the model selection process, which can significantly speed up the algorithm.
#' @param lambda.min.presc numeric fraction between 0 and 1. During pre-screening and the full
#' model fit for the BIC-ICQ calculation, the small penalty used on the random effect is
#' the fraction \code{lambda.min.presc} mulitplied by the maximum penalty parameter that penalizes
#' all fixed and random effects to 0. If left as \code{NULL}, the default value is 0.01 when the number
#' of random effects is 10 or less and 0.05 otherwise.
#' 
#' @return The *Control functions return a list (inheriting from class "\code{pglmmControl}") 
#' containing parameter values that determine settings for variable selection.
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
                         search = c("abbrev","full_grid"),
                         BIC_option = c("BICq","BICh","BIC","BICNgrp"), 
                         logLik_calc = switch(BIC_option[1], BICq = F, T), 
                         lambda.min = NULL, pre_screen = T, lambda.min.presc = NULL){
  
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
    if(any(lambda0_seq < 0)){
      stop("lambda0_seq cannot be negative")
    }
    if(length(lambda0_seq) > 1){
      for(i in 2:length(lambda0_seq)){
        if(lambda0_seq[i-1] > lambda0_seq[i]){
          stop("lambda0_seq must be a sequence of ascending order (min value to max value)")
        }
      }
    }
    
  } # End if-else !is.null(lambda0_seq)
  
  if(!is.null(lambda1_seq)){
    if(!is.numeric(lambda1_seq)){
      stop("lambda1_seq must be numeric")
    }
    if(any(lambda1_seq < 0)){
      stop("lambda1_seq cannot be negative")
    }
    if(length(lambda1_seq) > 1){
      for(i in 2:length(lambda1_seq)){
        if(lambda1_seq[i-1] > lambda1_seq[i]){
          stop("lambda1_seq must be a sequence of ascending order (min value to max value)")
        }
      }
    }
  } # End if-else !is.null(lambda1_seq)
  
  if(!(floor(nlambda) == nlambda)){
    stop("nlambda must be an integer")
  }
  if(nlambda < 2){
    stop("nlambda must be at least 2")
  }
  
  if (!is.logical(logLik_calc)) {
    stop("'logLik_calc' must be a logical value (T or F)")
  }
  
  if((BIC_option %in% c("BICh","BIC","BICNgrp")) & (logLik_calc == F)){
    stop("When 'BIC_option' is BICh, BIC, or BICNgroup, 'logLik_calc' must be TRUE")
  }
  
  if (!is.logical(pre_screen)) {
    stop("'pre_screen' must be a logical value (T or F)")
  }
  
  if(!is.null(lambda.min)){
    if(!is.numeric(lambda.min)){
      stop("lambda.min must be numeric")
    }
    if((lambda.min >= 1) | (lambda.min <=0 )){
      stop("lambda.min must be a fraction between 0 and 1")
    }
  }
  
  if(length(search) > 1) search = search[1]
  if(!(search %in% c("full_grid","abbrev"))){
    stop("'search' must be either 'full_grid' or 'abbrev'")
  }
  
  if(!is.null(lambda.min.presc)){
    if(!is.numeric(lambda.min.presc)){
      stop("lambda.min.presc must be numeric")
    }
    if((lambda.min.presc >= 1) | (lambda.min.presc <=0 )){
      stop("lambda.min.presc must be a fraction between 0 and 1")
    }
  }
  
  
  structure(list(lambda0_seq = lambda0_seq, lambda1_seq = lambda1_seq, search = search,
                 nlambda = nlambda, BIC_option = BIC_option, logLik_calc = logLik_calc,
                 lambda.min = lambda.min, pre_screen = pre_screen, 
                 lambda.min.presc = lambda.min.presc),
            class = c("selectControl", "pglmmControl"))
}

#' @title Control of Metropolis-within-Gibbs Adaptive Random Walk Sampling Procedure
#' 
#' Controls the adaptive random walk Metropolis-within-Gibbs sampling procedure.
#' 
#' @param batch_length positive integer specifying the number of posterior samples to collect
#' before the proposal variance is adjusted based on the acceptance rate of the last 
#' \code{batch_length} accepted posterior samples. Default is set to 100. Batch length restricted
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
#' EM algorithm. Default is 0.0015. 
#' EM algorithm is considered to have converge if the average Euclidean 
#' distance between the current coefficient estimates and the coefficient estimates from 
#' \code{t} EM iterations back is less than \code{conv_EM} \code{mcc} times in a row.
#' See \code{t} and \code{mcc} for more details.
#' @param conv_CD a non-negative numeric convergence criteria for the convergence of the 
#' grouped coordinate descent loop within the M step of the EM algorithm. Default 0.0005.
#' @param nMC_burnin positive integer specifying the number of posterior samples to use as
#' burnin for each E step in the EM algorithm. If set to \code{NULL}, the algorithm inputs
#' the following defaults: Default 250 when the number of random effects 
#' predictors is less than or equal to 10; default 100 otherwise. Function will not allow \code{nMC_burnin}
#' to be less than 100. 
#' @param nMC_start a positive integer for the initial number of Monte Carlo draws. If set to
#' \code{NULL}, the algorithm inputs the following defaults: Default 250 when 
#' the number of random effects predictors is less than or equal to 10; default 100 otherwise.
#' @param nMC_max a positive integer for the maximum number of allowed Monte Carlo draws used
#' in each step of the EM algorithm. If set to \code{NULL}, the algorithm inputs the following 
#' defaults: When the number of random effect predictors is 10 or less, 
#' Default is set to 5000 when no selection is performed and 2500 when selection is performed.
#' Default is set to 1000 when the number of random effect predictors is greater than 10.
#' @param nMC_report a positive integer for the number of posterior samples to save from the final
#' model. These posterior samples can be used for diagnostic purposes, see \code{\link{plot_mcmc}}.
#' Default set tp 5000.
#' @param maxitEM a positive integer for the maximum number of allowed EM iterations. 
#' If set to \code{NULL}, then the algorithm inputs the following defaults:
#' Default equals 50 for the Binomial and Poisson families, 100 for the Gaussian family.
#' @param maxit_CD a positive integer for the maximum number of allowed interations for the
#' coordinate descent algorithms used within the M-step of each EM iteration. Default equals 50.
#' @param M positive integer specifying the number of posterior samples to use within the 
#' Pajor log-likelihood calculation. Default is 10^4; minimum allowed value is 5000.
#' @param t the convergence criteria is based on the average Euclidean distance between 
#' the most recent coefficient estimates and the coefficient estimates from t EM iterations back.
#' Positive integer, default equals 2.
#' @param mcc the number of times the convergence critera must be met before the algorithm is
#' seen as having converged (mcc for 'meet condition counter'). Default set to 2. Value retricted 
#' to be no less than 2.
#' @param sampler character string specifying whether the posterior samples of the random effects
#' should be drawn using Stan (default, from package rstan) or the Metropolis-within-Gibbs procedure 
#' incorporating an adaptive random walk sampler ("random_walk") or an
#' independence sampler ("independence"). If using the random walk sampler, see \code{\link{adaptControl}}
#' for some additional control structure parameters.
#' @param var_start either the character string "recommend" or a positive number specifying the 
#' starting values to initialize the variance of the covariance matrix. Default "recommend" first
#' fits a simple model with a fixed and random intercept only using a Laplace approximation. The 
#' random intercept variance estimate from this model is then multiplied by 2 and used as the 
#' starting variance. 
#' 
#' @details Several arguments are set to a default value of \code{NULL}. If these arguments 
#' are left as \code{NULL} by the user, then these values will be filled in with appropriate
#' default values as specified above, which may depend on the number of random effects,
#' the family of the data, and/or whether selection is being performed. If the user
#' specifies particular values for these arguments, no additional modifications to these 
#' arguments will be done within the algorithm.
#' 
#' @return Function returns a list inheriting from class \code{optimControl}
#' containing fit and optimization criteria values used in optimization routine.
#' 
#' @export
optimControl = function(conv_EM = 0.0015, conv_CD = 0.0005, 
                        nMC_burnin = NULL, nMC_start = NULL, nMC_max = NULL, nMC_report = 5000,
                        maxitEM = NULL, maxit_CD = 50, 
                        M = 10000, t = 2, mcc = 2,
                        sampler = c("stan","random_walk","independence"), 
                        var_start = "recommend"){
  
  # Acceptable input types and input restrictions
  ## Arguments with default as NULL
  args_null = list(nMC_burnin = nMC_burnin, nMC_start = nMC_start, 
                   nMC_max = nMC_max, maxitEM = maxitEM)
  ## x = vector of several arguments restricted to be positive integers, will check inputs
  ##    in next steps.
  x = c(nMC_report, maxitEM, maxit_CD, M, t, mcc)
  ## If user set some of the 'args_null' arguments to non-NULL values, check these input values
  ##    as well (should be positive integers).
  for(a in 1:length(args_null)){
    if(!is.null(args_null[[a]])){
      x = c(x, args_null[[a]])
    }
  }
  ## Check arguments that should be positive integers
  if((!all(floor(x)==x)) | (sum(x <= 1) > 0)){ # if any of the above values not positive integers
    stop("M, t, mcc, and all entered nMC and maxit arguments must be positive integers")
  }
  
  # More restrictive checks
  if(M < 5000){
    stop("M must be greater than or equal to 5000")
  }
  if(!is.null(nMC_burnin)){
    if(nMC_burnin < 100){
      warning("nMC_burnin not allowed to be less than 100. Value set to 100", immediate. = T)
      nMC_burnin = 100
    }
  }
  if(!is.null(nMC_max) & !is.null(nMC_start)){
    if(nMC_max < nMC_start){
      stop("nMC_max cannot be smaller than nMC_start")
    }
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
  
  sampler = checkSampler(sampler)
  
  if(mcc < 2){
    stop("mcc must be at least 2")
  }
  
  structure(list(conv_EM = conv_EM, conv_CD = conv_CD, 
                 nMC_burnin = nMC_burnin, nMC = nMC_start, nMC_max = nMC_max, nMC_report = nMC_report,
                 maxitEM = maxitEM, maxit_CD = maxit_CD,  M = M, t = t, mcc = mcc,
                 sampler = sampler, var_start = var_start),
            class = "optimControl")
}

# optim_recommend: For input arguments of optimControl() that are set to NULL by default,
#   recommend appropriate inputs that depend on the family, number of random effects,
#   and whether or not variable selection is being performed.
# q = number of random effects (including random intercept)
# select: TRUE if running the selection algorithm
optim_recommend = function(optim_options, family, q, select){
  
  # If selection and moderate/small q, or if large q (with or without selection):
  
  if (q <= 11) { # q includes random intercept
    # For selection, decrease nMC_max to improve time
    ## Note: During selection, will initialize with good starting points, so
    ## larger nMC_max not needed for convergence
    ## In order to speed up algorithm, will decrease nMC_max while increasing maxitEM
    # Otherwise, keep everything else the same
    if (select == T) {
      if(is.null(optim_options$nMC_max)){
        optim_options$nMC_max = 2500
      }
    }
  }else{
    # Decrease nMC
    if(is.null(optim_options$nMC_burnin)){
      optim_options$nMC_burnin = 100
    }
    if(is.null(optim_options$nMC)){
      optim_options$nMC = 100
    }
    if(is.null(optim_options$nMC_burnin)){
      optim_options$nMC_max = 1000
    }
    
  } # End if-else q <= 11
  
  # If not special case of large q and/or selection, give default values of nMC args and maxitEM
  
  if(is.null(optim_options$nMC_burnin)){
    optim_options$nMC_burnin = 250
  }
  if(is.null(optim_options$nMC)){
    optim_options$nMC = 250
  }
  if(is.null(optim_options$nMC_max)){
    optim_options$nMC_max = 5000
  }
  if(is.null(optim_options$maxitEM)){
    if(family %in% c("binomial","poisson")){
      optim_options$maxitEM = 50
    }else if(family == "gaussian"){
      optim_options$maxitEM = 50
    }
  }
  
  if(optim_options$nMC_max < optim_options$nMC){
    stop("in optimControl, nMC_max, ",optim_options$nMC_max," must be less than nMC_start ", optim_options$nMC)
  }
  
  return(optim_options)
  
}

