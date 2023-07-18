#' @name lambdaControl
#' @aliases selectControl
#' 
#' @title Control of Penalization Parameters and Selection Criteria
#' 
#' @description Constructs control structures for penalized mixed model fitting.
#' 
#' @param lambda0 a non-negative numeric penalty parameter for the fixed effects coefficients
#' @param lambda1 a non-negative numeric penalty parameter for the (grouped) random effects
#' covariance coefficients
#' @param lambda0_seq,lambda1_seq a sequence of non-negative numeric penalty parameters for 
#' the fixed 
#' and random effect coefficients (\code{lambda0_seq} and \code{lambda1_seq}, respectively). 
#' If \code{NULL}, then a sequence will be automatically 
#' calculated. See 'Details' section for more details on these default calculations. 
#' @param nlambda positive integer specifying number of penalty parameters  
#' to use for the fixed and random effects penalty parameters. Default set to 10. Only used
#' if either \code{lambda0_seq} or \code{lambda1_seq} remain unspecified by the user
#' (one or both of these arguments set to \code{NULL}) and, consequently, one or more default sequences need
#' to be calculated.
#' @param BIC_option character string specifying the selection criteria used to select the 'best' model.
#' Default "BICq" option specifies the BIC-ICQ criterion (Ibrahim et al (2011)
#' <doi:10.1111/j.1541-0420.2010.01463.x>),
#' which requires a fit of 
#' a 'minimum penalty' model; a small penalty (the minimum of the penalty sequence) 
#' is used for the fixed and random effects. See "Details" section for what these 
#' small penalties will be.
#' The "BICh" option utilizes the hybrid BIC value described in 
#' Delattre, Lavielle, and Poursat (2014) <doi:10.1214/14-EJS890>.
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
#' automatically selects \code{lambda.min} to equal 0.01 when the number of observations is
#' greater than the number of fixed effects predictors and 0.05 otherwise.
#' Only used
#' if either \code{lambda0_seq} or \code{lambda1_seq} remain unspecified by the user
#' (one or both of these sequence arguments set to \code{NULL}) and, consequently, 
#' one or more default sequences need
#' to be calculated.
#' @param pre_screen logical value indicating whether pre-screening should be performed before
#' model selection (default \code{TRUE}). If the number of random effects covariates
#' considered is 4 or less, then
#' no pre-screening will be performed. Pre-screening removes random effects from consideration
#' during the model selection process, which can significantly speed up the algorithm.
#' See "Details" section for a further discussion.
#' @param lambda.min.presc numeric fraction between 0 and 1. During pre-screening and the 
#' minimal penalty
#' model fit for the BIC-ICQ calculation, the small penalty used on the random effect is
#' the fraction \code{lambda.min.presc} multiplied by the maximum penalty parameter that penalizes
#' all fixed and random effects to 0. If left as \code{NULL}, the default value is 0.01 when the number
#' of random effect covariates is 10 or less and 0.05 otherwise.
#' Only used if \code{lambda1_seq} remains unspecified by the user
#' (this argument set to \code{NULL} so the random effects penalty parameter sequence
#' needs to be automatically calculated) AND either the pre-screening procedure is selected by the argument
#' \code{pre_screen} or the BIC-ICQ is selected as the model selection criteria,
#' i.e., \code{BIC_option} = "BICq". See the "Details" section for a further discussion.
#' 
#' 
#' @details If left as the default \code{NULL} values, 
#' the \code{lambda0_seq} and \code{lambda1_seq} numeric sequences are 
#' automatically calculated. The sequence will be calculated in the same manner as 
#' \code{ncvreg} calculates the range: the max value (let's denote this as \code{lambda_max}) 
#' penalizes all fixed and random effects to 0, the min value is a 
#' small portion of max (\code{lambda.min}*\code{lambda_max}), and the sequence is composed of 
#' \code{nlambda} values ranging from these min and max values spread evenly on the log scale. 
#' Unlike \code{ncvreg}, the order of penalty
#' values used in the algorithm must run from the min lambda to the max lambda (as opposed to 
#' running from max lambda to min lambda). The length of the sequence is specified by \code{nlambda}. 
#' By default, these sequences are calculated using \code{\link{LambdaSeq}}.
#' 
#' The \code{lambda0} and \code{lambda1} arguments used within the \code{\link{glmm}} function 
#' allow for a user to fit a model with a single 
#' non-zero penalty parameter combination. However, this is generally not recommended. 
#' 
#' Abbreviated grid search: The abbreviated grid search proceeds in two stages. In stage 1, the
#' algorithm fits the following series of models: the fixed effects penalty parameter remains a
#' fixed value evaluated at the minimum of the fixed effects penalty parameters, and all
#' random effects penalty parameters are examined. The 'best' model from this first stage of models
#' determines the optimum random effect penalty parameter. In stage 2, the algorithm fits the 
#' following series of models: the random effects penalty parameter remains fixed at the value of
#' the optimum random effect penalty parameter (from stage 1) and all fixed effects penalty
#' parameters are considered. The best overall model is the best model from stage 2. This reduces the 
#' number of models considered to length(`lambda0_seq`) + length(`lambda1_seq`). The authors found
#' that this abbreviated grid search worked well in simulations, and performed considerably
#' faster than the full grid search that examined all possible fixed and random effect penalty
#' parameter combinations.
#' 
#' The arguments \code{nlambda} and \code{lambda.min} are only used
#' if one or both of the \code{lambda0_seq} and \code{lambda1_seq} penalty sequences
#' (corresponding to the fixed and random effects penalty sequences, respectively) 
#' remain unspecified by the user (one or both of these arguments left as \code{NULL}),
#' indicating that the algorithm needs to calculate default penalty sequences.
#' 
#' The argument \code{lambda.min.presc} is only used under the following condition:
#' \code{lambda1_seq} remains unspecified by the user
#' (this argument set to \code{NULL} so the random effects penalty parameter sequence
#' needs to be calculated) AND either the pre-screening procedure is selected by the argument
#' \code{pre_screen} or the BIC-ICQ is selected as the model selection criteria,
#' i.e., \code{BIC_option} = "BICq". 
#' If \code{lambda1_seq} is specified by the user, the minimum
#' value in that sequence will be used as the random effect penalty in the
#' pre-screening procedure and/or the minimal penalty model for the BIC-ICQ calculation.
#' 
#' BIC-ICQ calculation: This model selection criteria requires the fitting of a 'minimal penalty'
#' model, which fits a model with a small penalty on the fixed and random effects.
#' For the fixed effects penalty, the minimal penalty is: (a) 0 if the number of fixed 
#' effects covariates is 4 or less or (b) the minimum fixed effect penalty from the fixed
#' effects penalty sequence (either from the default sequence or from the sequence
#' specified by the user). For the random effects penalty, the minimal penalty
#' is (a) 0 if the number of random 
#' effects covariates is 4 or less; (b) the minimum random effect penalty 
#' from the random effects penalty sequence specified by the user, or 
#' (c) \code{lambda.min.presc} multiplied to the \code{lambda_max} maximum penalty
#' specified above when a default random effects penalty sequence is calculated.
#' 
#' Pre-screening: The minimum fixed effects penalty used in the pre-screening stage 
#' will be the minimum penalty of the fixed effects penalty sequence, \code{lambda0_seq}.
#' The minimum random effects penalty used in the pre-screening stage will be either 
#' (a) the minimum random effects penalty in the sequence \code{lambda1_seq} 
#' if this sequence specified by the user, or (b) \code{lambda.min.pres} x \code{lambda_max},
#' where \code{lambda_max} was described above.
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
                         logLik_calc = switch(BIC_option[1], BICq = FALSE, TRUE), 
                         lambda.min = NULL, pre_screen = TRUE, lambda.min.presc = NULL){
  
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
  
  if((BIC_option %in% c("BICh","BIC","BICNgrp")) & (logLik_calc == FALSE)){
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

#' @title Control of Cox Proportional Hazards Model Fitting
#' 
#' @description Constructs the control structure for additional parameters needed to
#' properly fit survival data using a Piecewise Exponential model
#' 
#' @param cut_num positive integer specifying the number of time intervals to include in
#' the piecewise exponential hazard model assumptions for the sampling step. Default is 8.
#' General recommendation: use between 5 and 10 intervals. See the Details section for
#' additional information.
#' @param interval_type character specifying how the time intervals are calculated.
#' If 'equal' (default),  time intervals
#' are calculated such that there are approximately equal numbers of events per time 
#' interval.
#' If 'manual', the user needs to input
#' their own cut points (see \code{cut_points} for details).
#' If 'group', time intervals are calculated such that there are approximately
#' equal numbers of events per time interval AND there are at least 4 events observed
#' by each group within each time interval. The input number of time intervals \code{cut_num}
#' may be modified to a lower number in order to accomplish this goal. 
#' This option is generally difficult to perform when there are 
#' a large number of groups in the data.
#' @param cut_points numeric vector specifying the value of the cut points to use
#' in the calculation of the time intervals for the piecewise exponential model. 
#' If \code{interval_type} set to 'equal' or 'group', then this argument is not utilized.
#' If \code{interval_type} set to 'manual', then this argument is required.
#' First value must be 0, and all values must be ordered smallest to largest.
#' @param time_scale positive numeric value (greater than 1) used to scale the time variable 
#' in the survival data. In order to calculate the piecewise exponential model,
#' the log of the time a subject survived within a particular interval
#' is used as an offset in the model fit. Sometimes multiplying the time scale
#' by a factor greater than 1 improves the stability of the model fit algorithm.
#' 
#' @return Function returns a list inheriting from class \code{optimControl}
#' containing fit and optimization criteria values used in optimization routine.
#' 
#' @details In the piecewise exponential hazard model assumption, there is an assumption that the 
#' time line of the data can be cut into \code{cut_num}
#' time intervals and the baseline hazard is constant within
#' each of these time intervals. In the fit algorithm, we estimate
#' these baseline hazard values (specifically, we estimate the log of the baseline
#' hazard values). We determine cut points by specifying the total number of cuts
#' to make (\code{cut_num}) and then specifying time values for cut points such
#' that each time interval has an approximately equal number 
#' of events and that each group has at least some (4) events within each time interval. 
#' Having too many cut points could result in not having enough
#' events for each time interval.
#' Additionally, data with few events could result in too few events per time interval
#' even for a small number of cut points. We generally recommend having
#' 8 total time intervals (more broadly, between 5 and 10). Warnings or errors
#' will occur for cases when there are 1 or 0 events for a time interval. 
#' If this happens, either adjust the \code{cut_num} value appropriately,
#' or in the case when the data simply has a very small number of events,
#' consider not using this software for your estimation purposes. 
#' 
#' @export
survivalControl = function(cut_num = 8, interval_type = c("equal","manual","group"), 
                        cut_points = NULL, time_scale = 1){
  
  #########################################################################################
  # Input checks and restrictions
  #########################################################################################
  
  # cut_num
  if((floor(cut_num) != cut_num) | (cut_num < 1)){
    stop("cut_num must be a positive integer")
  }
  
  if(length(interval_type) > 1){
    interval_type = interval_type[1]
  }
  if(!(interval_type %in% c("equal","manual","group"))){
    stop("interval_type must be either 'equal', 'manual, or 'group', see coxphControl() documentation")
  }
  
  if(interval_type == "manual"){
    if(is.null(cut_points)){
      stop("cut_points must be specified if interval_type = 'manual' ")
    }
    if(length(cut_points) != cut_num){
      message("changing 'cut_num' to length of input cut_points vector")
      cut_num = length(cut_points)
    }
    if(cut_points[1] != 0){
      stop("first cut-point value must equal 0")
    }
    if(length(cut_points) > 1){
      for(i in 2:length(cut_points)){
        if(cut_points[i-1] > cut_points[i]){
          stop("cut_points must be a sequence of ascending order (min value to max value)")
        }
      }
    }
  }
  
  if((interval_type != "manual") & (!is.null(cut_points))){
    message("if interval_type set to 'equal' or 'group', then cut_points will be automatically calculated and input cut_points will not be used")
    cut_points = NULL
  }
  
  if((cut_num < 5) | (cut_num > 10)){
    warning("the glmmPen team recommends that you keep cut_num between 5 and 10; 8 is typically a good cut_num value", immediate. = TRUE)
  }
  
  if(time_scale < 1){
    stop("time_scale must be a positive value greater than 1")
  }
  
  # output object
  structure(list(cut_num = cut_num, interval_type = interval_type, 
                 cut_points = cut_points, time_scale = time_scale),
            class = "survivalControl")
  
}

#' @title Control of Penalized Generalized Linear Mixed Model Fitting
#' 
#' @description Constructs the control structure for the optimization of the penalized mixed model fit algorithm.
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
#' burn-in for each E step in the EM algorithm. If set to \code{NULL}, the algorithm inputs
#' the following defaults: Default 250 when the number of random effects 
#' predictors is less than or equal to 10; default 100 otherwise. Function will not allow \code{nMC_burnin}
#' to be less than 100. 
#' @param nMC_start a positive integer for the initial number of Monte Carlo draws. If set to
#' \code{NULL}, the algorithm inputs the following defaults: Default 250 when 
#' the number of random effects predictors is less than or equal to 10; default 100 otherwise.
#' @param nMC_max a positive integer for the maximum number of allowed Monte Carlo draws used
#' in each step of the EM algorithm. If set to \code{NULL}, the algorithm inputs the following 
#' defaults: When the number of random effect covariates is greater than 10,
#' the default is set to 1000; when the number of random effect covariates is
#' 10 or less, the default is set to 2500.
#' @param nMC_report a positive integer for the number of posterior samples to save from the final
#' model. These posterior samples can be used for diagnostic purposes, see \code{\link{plot_mcmc}}.
#' Default set to 5000.
#' @param maxitEM a positive integer for the maximum number of allowed EM iterations. 
#' If set to \code{NULL}, then the algorithm inputs the following defaults:
#' Default equals 50 for the Binomial and Poisson families, 65 for the Gaussian family.
#' @param maxit_CD a positive integer for the maximum number of allowed iterations for the
#' coordinate descent algorithms used within the M-step of each EM iteration. Default equals 50.
#' @param M positive integer specifying the number of posterior samples to use within the 
#' Pajor log-likelihood calculation. Default is 10^4; minimum allowed value is 5000.
#' @param t the convergence criteria is based on the average Euclidean distance between 
#' the most recent coefficient estimates and the coefficient estimates from \code{t} EM iterations back.
#' Positive integer, default equals 2.
#' @param mcc the number of times the convergence criteria must be met before the algorithm is
#' seen as having converged (mcc for 'meet condition counter'). Default set to 2. Value restricted 
#' to be no less than 2.
#' @param sampler character string specifying whether the posterior samples of the random effects
#' should be drawn using Stan (default, from package rstan) or the Metropolis-within-Gibbs procedure 
#' incorporating an adaptive random walk sampler ("random_walk") or an
#' independence sampler ("independence"). If using the random walk sampler, see \code{\link{adaptControl}}
#' for some additional control structure parameters.
#' @param var_start either the character string "recommend" or a positive number specifying the 
#' starting values to initialize the variance of the covariance matrix. For \code{\link{glmmPen}},
#' the default "recommend" first
#' fits a simple model with a fixed and random intercept only using the \link{lme4} package. The 
#' random intercept variance estimate from this model is then multiplied by 2 and used as the 
#' starting variance. For \code{\link{glmmPen_FA}}, the default is set to 0.10 (see \code{B_init_type} 
#' for further information).
#' @param step_size positive numeric value indicating the starting step size to use in the 
#' Majorization-Minimization scheme of the M-step. Only relevant when the distributional assumption
#' used is not Binomial or Gaussian with canonical links (e.g. Poisson with log link)
#' @param standardization logical value indicating whether covariates should
#' standardized (\code{TRUE}, default) or unstandardized (\code{FALSE}) before being
#' used within the algorithm. If \code{standardization = TRUE}, then the standardized covariates
#' will also be used to create the Z matrix used in the estimation of the random effects.
#' @param convEM_type character string indicating the type of convergence criteria to 
#' use within the EM algorithm to determine when a model has converged. The default is "AvgEuclid1",
#' which calculates the average Euclidean distance between the most recent coefficient vector and
#' the coefficient vector \code{t} EM iterations back (Euclidean distance divided by the number
#' of non-zero coefficients \code{t} EM iterations back). Alternative convergence options include
#' "maxdiff", which determines convergence based on the maximum difference between the coefficient vectors; 
#' "AvgEuclid2", which is similar to "AvgEuclid1" except it divides the Euclidean distance by the square-root
#' of the number of non-zero coefficients; and "Qfun", which determines convergence based on
#' the relative difference in the Q-function estimates calculated with the most recent coefficient vector 
#' and the coefficient vector \code{t} EM iterations back.
#' @param B_init_type character string indicating how the B matrix within the \code{\link{glmmPen_FA}}
#' method should be initialized. (This argument is not used within the \code{\link{glmmPen}} function.)
#' The default "deterministic" initializes all non-zero variance and 
#' covariance values of the random effect covariance matrix to the value of \code{var_start},
#' such that each non-zero element of the B matrix is \code{sqrt(var_start / r)} (where \code{r} is
#' the number of latent factors). Option "data" is similar to "deterministic", but the 
#' \code{var_start} value is the default data-driven variance estimate used in \code{\link{glmmPen}} 
#' (see argument \code{var_start} for more details).
#' @param var_restrictions character string indicating how the random effect covariance
#' matrix should be initialized at the beginning of the algorithm
#' when penalties are applied to the coefficients. 
#' If "none" (default), all random effect predictors are initialized to have non-zero variances.
#' If "fixef", the code first examines the initialized fixed effects (initialized using a regular
#' penalized GLM), and only the random effect predictors that are initialized with non-zero fixed effects
#' are initialized with non-zero variances.
#' 
#' @details Several arguments are set to a default value of \code{NULL}. If these arguments 
#' are left as \code{NULL} by the user, then these values will be filled in with appropriate
#' default values as specified above, which may depend on the number of random effects or
#' the family of the data. If the user
#' specifies particular values for these arguments, no additional modifications to these 
#' arguments will be done within the algorithm.
#' 
#' @return Function returns a list inheriting from class \code{optimControl}
#' containing fit and optimization criteria values used in optimization routine.
#' 
#' @export
optimControl = function(var_restrictions = c("none","fixef"),
                        conv_EM = 0.0015, conv_CD = 0.0005, 
                        nMC_burnin = NULL, nMC_start = NULL, nMC_max = NULL, nMC_report = 5000,
                        maxitEM = NULL, maxit_CD = 50, 
                        M = 10000, t = 2, mcc = 2,
                        sampler = c("stan","random_walk","independence"), 
                        var_start = "recommend", step_size = 1.0,
                        standardization = TRUE,
                        convEM_type = c("AvgEuclid1","maxdiff","AvgEuclid2","Qfun"),
                        B_init_type = c("deterministic","data","random")){
  
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
      warning("nMC_burnin not allowed to be less than 100. Value set to 100", immediate. = TRUE)
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
  
  if(step_size <= 0){
    stop("step_size must be positive")
  }
  
  if(!(standardization %in% c(TRUE,FALSE))){
    stop("standardization must be logical, TRUE or FALSE")
  }
  
  convEM_type = convEM_type[1]
  if(!(convEM_type %in% c("AvgEuclid1","AvgEuclid2","Qfun","maxdiff"))){
    stop("convEM_type must be one of 'AvgEuclid1', 'AvgEuclid2', 'Qfun', or 'maxdiff'")
  }
  
  B_init_type = B_init_type[1]
  if(!(B_init_type %in% c("random","deterministic","data"))){
    stop("B_init_type must be either 'random' or 'deterministic' or 'data'")
  }
  
  var_restrictions = var_restrictions[1]
  if(!(var_restrictions %in% c("none","fixef"))){
    stop("var_restrictions must be either 'none' or 'fixef'")
  }
  
  structure(list(conv_EM = conv_EM, conv_CD = conv_CD, 
                 nMC_burnin = nMC_burnin, nMC = nMC_start, nMC_max = nMC_max, nMC_report = nMC_report,
                 maxitEM = maxitEM, maxit_CD = maxit_CD,  M = M, t = t, mcc = mcc,
                 sampler = sampler, var_start = var_start, step_size = step_size,
                 standardization = standardization, convEM_type = convEM_type,
                 B_init_type = B_init_type, var_restrictions = var_restrictions),
            class = "optimControl")
}

# optim_recommend: For input arguments of optimControl() that are set to NULL by default,
#   recommend appropriate inputs that depend on the family, number of random effects,
#   and whether or not variable selection is being performed.
# q = number of random effects (including random intercept) 
#     Alternatively, q = number latent factors for glmm_FA and glmmPen_FA functions
# select: TRUE if running the selection algorithm
optim_recommend = function(optim_options, family, q, select){

  # Default nMC in case with large number of random effects covariates
  if(q <= 11){ # q includes random intercept, "if number random effect covariates <= 10"
    # If not special case of large q, give default values of nMC args
    if(is.null(optim_options$nMC_burnin)){
      optim_options$nMC_burnin = 250
    }
    if(is.null(optim_options$nMC)){
      optim_options$nMC = 250
    }
    if(is.null(optim_options$nMC_max)){
      optim_options$nMC_max = 2500
    }
  }else{ # q includes random intercept, "if number random effect covariates >= 11"
    # Decrease burn-in and starting nMC
    if(is.null(optim_options$nMC_burnin)){
      optim_options$nMC_burnin = 100
    }
    if(is.null(optim_options$nMC)){
      optim_options$nMC = 100
    }
    if(is.null(optim_options$nMC_max)){
      optim_options$nMC_max = 1000
    }
  } # End if-else q <= 11

  # Default maxitEM
  if(is.null(optim_options$maxitEM)){
    if(family %in% c("binomial","poisson")){
      optim_options$maxitEM = 50
    }else if(family == "gaussian"){
      optim_options$maxitEM = 65
    }
  }

  if(optim_options$nMC_max < optim_options$nMC){
    stop("in optimControl, nMC_max, ",optim_options$nMC_max," must be less than nMC_start ", optim_options$nMC)
  }

  return(optim_options)

}



#' @title Control of Latent Factor Model Number Estimation
#' 
#' Constructs the control structure for the estimation of the 
#' number of latent factors (r) for use within the \code{glmmPen_FA} and
#' \code{glmm_FA} estimation procedures.
#' 
#' @param r positive integer specifying number of latent common factors to assume 
#' in the model. If \code{NULL} (default), this value estimated from the data. See 
#' \code{r_est_method} for available estimation procedures, and the Details
#' section for further details on the general estimation procedure.
#' If \code{r} is specified, the no estimation procedure is performed and the algorithm
#' uses the input value or r. All other parameters for this function are
#' relevant for the estimation procedure.
#' @param r_max positive integer specifying maximum number of latent factors to consider.
#' If \code{NULL} (default), this value is automatically calculated.
#' @param r_est_method character string indicating method used to estimate number
#' of latent factors \code{r}. Default "GR" uses the Growth Ratio method of
#' Ahn and Horenstein (2013) (<doi:10.3982/ECTA8968>). 
#' Other available options include "ER" for
#' the Eigenvalue Ratio method of Ahn and Horenstein (2013) (<doi:10.3982/ECTA8968>) 
#' and "BN1" or "BN2",
#' the Bai and Ng (2002) method (<dio:10.1111/1468-0262.00273>) using one of two penalties: 
#' (1) \code{(d + p) / (d p) log(d p/(d+p))} or
#' (2) \code{(d + p) / (d p) log(min(d,p))} where d is the number of groups in
#' the data and p is the number of total random effect covariates (including the intercept)
#' @param size positive integer specifying the total number of pseudo random
#' effect estimates to use in the estimation procedure for the number of latent factors
#' r, which is restricted to be no less than 25. If this \code{size} is greater
#' than the number of groups in the data (i.e.~the number of levels of the grouping
#' variable), then a sampling procedure is used to increase the number of pseudo estimates
#' to the value of \code{size} if the value of \code{sample} is \code{TRUE}.
#' @param sample logical value indicating if the total number of pseudo random effect
#' estimates to use in the estimation procedure for the number of latent common factors r
#' should be larger than the number of unique groups in the data, where the number 
#' of pseudo estimates are increased to the value of \code{size} using a sampling 
#' procedure. Default is \code{FALSE}. If \code{TRUE}, the sampling procedure is only
#' performed if the value of \code{size} is greater than the number of groups in the data.
#' 
#' @details Estimation of \code{r} procedure: For each level of the group variable separately,
#' we identify the observations within that group and 
#' fit a regular penalized generalized linear model where the penalty value is the
#' minimum fixed effect penalty. These group-specific estimates, which we label as 'pseudo random effects',
#' are placed into a matrix \code{G}
#' (rows = number of levels of the grouping variable, columns = number of random effect covariates),
#' and this pseudo random effects matrix is treated as the observed outcome matrix used in
#' the "GR", "ER", and "BN" estimation procedures described above in the description of \code{r_est_method}.
#' 
#' @export
rControl = function(r = NULL, r_max = NULL, r_est_method = "GR",
                    size = 25, sample = FALSE){
  
  # r, r_max, and sample_size must all be positive integers
  int_check = list(r = r, r_max = r_max, size = size)
  var_int_names = names(int_check)
  for(i in 1:length(int_check)){
    x = int_check[[i]]
    if(!is.null(x)){
      if(!(floor(x) == x) | (x <= 0)){
        stop(x, " must be a postive integer")
      }
    }
  }
  
  # Additional size checks: r, r_max, and sample_size must be of sufficient size
  if(!is.null(r)){
    if(r < 2){
      message("number of latent factors r restricted to be no less than 2")
    }
  }
  
  if(size < 25){
    stop("size restricted to be no less than 25")
  }
  
  # Check of r estimation method
  if(!(r_est_method %in% c("GR","ER","BN1","BN2"))){ 
    stop("r_est_method must be one of GR, ER, BN1, or BN2, see rControl() documentation")
  } 
  
  if(!is.logical(sample)){
    stop("sample must be either TRUE or FALSE")
  }
  
  structure(list(r = r, r_max = r_max, r_est_method = r_est_method,
                 size = size, sample = sample),
            class = "rControl")
  
}







## Alternative?
# optim_recommend: For input arguments of optimControl() that are set to NULL by default,
#   recommend appropriate inputs that depend on the family, number of random effects,
#   and whether or not variable selection is being performed.
# q = number of random effects (including random intercept) 
#     Alternatively, q = number latent factors for glmm_FA and glmmPen_FA functions
# select: TRUE if running the selection algorithm
# optim_recommend = function(optim_options, family, q, select){
#   
#   if(is.null(optim_options$nMC_burnin)){
#     optim_options$nMC_burnin = 100
#   }
#   
#   if(select == TRUE){
#     if(is.null(optim_options$nMC)){
#       optim_options$nMC = 100
#     }
#     if(is.null(optim_options$nMC_max)){
#       optim_options$nMC_max = 500
#     }
#     # Default maxitEM
#     if(is.null(optim_options$maxitEM)){
#       if(family %in% c("binomial","poisson")){
#         optim_options$maxitEM = 25
#       }else if(family == "gaussian"){
#         optim_options$maxitEM = 35
#       }
#     }
#   }else if(select == FALSE){
#     if(is.null(optim_options$nMC)){
#       optim_options$nMC = 250
#     }
#     if(is.null(optim_options$nMC_max)){
#       optim_options$nMC_max = 1000
#     }
#     # Default maxitEM
#     if(is.null(optim_options$maxitEM)){
#       if(family %in% c("binomial","poisson")){
#         optim_options$maxitEM = 50
#       }else if(family == "gaussian"){
#         optim_options$maxitEM = 65
#       }
#     }
#   }
#   
#   if(optim_options$nMC_max < optim_options$nMC){
#     stop("in optimControl, nMC_max, ",optim_options$nMC_max," must be less than nMC_start ", optim_options$nMC)
#   }
#   
#   return(optim_options)
#   
# }

