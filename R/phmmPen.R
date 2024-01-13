
#' @title Fit a Proportional Hazards Mixed Model via Monte Carlo Expectation Conditional Minimization (MCECM) using
#' a Piecewise Constant Hazard Mixed Model Approximation
#' 
#' @description \code{phmm_FA} is used to fit a single piecewise constant hazard mixed model as an
#' approximation to a proportional hazards mixed model via Monte Carlo 
#' Expectation Conditional Minimization (MCECM). This piecewise constant hazard mixed model
#' uses a factor model decomposition of
#' the random effects. No model selection is performed.
#' 
#' @inheritParams phmmPen
#' @inheritParams glmmPen 
#' @inheritParams glmmPen_FA
#' @param ... additional arguments that could be passed into \code{phmmPen_FA}. 
#' See \code{\link{phmmPen_FA}}
#' for further details.
#' 
#' @details The \code{phmm_FA} function can be used to approximate a single proportional hazards
#' mixed model using a piecewise constant hazard mixed model.
#' While this approach is meant to be used in the case where the user knows which
#' covariates belong in the fixed and random effects and no penalization is required, one is
#' allowed to specify non-zero fixed and random effects penalties using \code{\link{lambdaControl}}
#' and the (...) arguments. The (...) allow for specification of penalty-related arguments; see
#' \code{\link{phmmPen_FA}} for details. For a high dimensional situation, the user may want to fit a 
#' full model using a small penalty for the fixed and random effects and save the posterior
#' draws from this full model for use in any BIC-ICQ calculations during selection within \code{phmmPen_FA}. 
#' Specifying a file name in the 'BICq_posterior' argument will save the posterior draws from the 
#' \code{phmm_FA} model into a big.matrix with this file name, see the Details section of 
#' \code{\link{phmmPen_FA}} for additional details.
#' 
#' 
#' @return A reference class object of class \code{\link{pglmmObj}} for which many methods are 
#' available (e.g. \code{methods(class = "pglmmObj")})
#' 
#' @export
phmm_FA = function(formula, data = NULL,
                    offset = NULL, r_estimation = rControl(),
                    optim_options = optimControl(), adapt_RW_options = adaptControl(),
                    trace = 0, tuning_options = lambdaControl(),
                    survival_options = survivalControl(),
                    progress = TRUE, ...){
  
  # Specify family
  family = "coxph"
  
  # Check that (...) arguments are subsets of phmmPen arguments
  args_extra = list(...)
  args_avail = c("fixef_noPen","penalty","alpha","gamma_penalty","BICq_posterior")
  if(length(args_extra) >= 1){
    if(!(names(args_extra) %in% args_avail)){
      stop("additional arguments provided in '...' input must match phmmPen arguments \n",
           "allowed extra arguments inclue 'fixef_noPen', 'penalty', 'alpha', 'gamma_penalty', 'BICq_posterior' \n",
           "see phmmPen documentation for details")
    }
  }
  
  # Check that tuning_options has the correct input type
  if(!inherits(tuning_options, "lambdaControl")){
    stop("phmm requires the use of lambdaControl for the tuning_options. \n",
         "  In order to use selectControl for model selection, please use phmmPen function")
  }
  
  call = match.call(expand.dots = FALSE)
  
  # Call glmmPen_FA() (tuning_options = lambdaControl() specifies the fitting of a single model 
  #     and not model selection)
  # If ... arguments empty or only includes a few of the args_avail, 
  # glmmPen_FA will use default arguments 
  output = glmmPen(formula = formula, data = data, family = family,
                   offset = offset, r_estimation = r_estimation,
                   optim_options = optim_options,
                   adapt_RW_options = adapt_RW_options, trace = trace,
                   tuning_options = tuning_options, survival_options = survival_options,
                   progress = progress, ...)
  
  output$call = call
  out_object = pglmmObj$new(output)
  
  return(out_object)
  
  
}

#' @title Fit a Penalized Proportional Hazards Mixed Model via 
#' Monte Carlo Expectation Conditional Minimization (MCECM) using
#' a Piecewise Constant Hazard Mixed Model Approximation
#' 
#' @description \code{phmmPen_FA} is used to approximate penalized proportional hazards mixed models using
#' using a piecewise constant hazard mixed survival model approximation via 
#' Monte Carlo Expectation Conditional Minimization (MCECM) and a factor model decomposition of
#' the random effects. 
#' The purpose of the function is to perform 
#' variable selection on both the fixed and random effects simultaneously for the
#' piecewise constant hazard mixed model. 
#' This function uses a factor model decomposition on the random effects. This assumption
#' reduces the latent space in the E-step (Expectation step) of the algorithm,
#' which reduces the computational complexity of the algorithm. This improves
#' the speed of the algorithm and enables the
#' scaling of the algorithm to larger dimensions.
#' \code{phmmPen_FA} selects the best model using 
#' BIC-type selection criteria (see \code{\link{selectControl}} documentation for 
#' further details). 
#' To improve the speed of the algorithm, consider setting
#' \code{var_restrictions} = "fixef" within the \code{\link{optimControl}} options.
#' 
#' @inheritParams phmmPen
#' @inheritParams glmmPen
#' @inheritParams glmmPen_FA
#' 
#' @details Argument \code{BICq_posterior} details: If the \code{BIC_option} in \code{\link{selectControl}} 
#' (\code{tuning_options}) is specified 
#' to be 'BICq', this requests the calculation of the BIC-ICQ criterion during the selection
#' process. For the BIC-ICQ criterion to be calculated, a full model assuming a small valued 
#' lambda penalty needs to be fit, and the posterior draws from this full model need to be used. 
#' In order to avoid repetitive calculations of
#' this full model (i.e. if the user wants to re-run \code{phmmPen} with a different
#' set of penalty parameters), a \code{big.matrix} of these 
#' posterior draws will be file-backed as two files: a backing file with extention '.bin' and a 
#' descriptor file with extension '.desc'. The \code{BICq_posterior} argument should contain a 
#' path and a filename with no extension of the form "./path/filename" such that the backingfile and
#' the descriptor file would then be saved as "./path/filename.bin" and "./path/filename.desc", respectively.
#' If \code{BICq_posterior} is set to \code{NULL}, then by default, the backingfile and descriptor
#' file are saved in the working directory as "BICq_Posterior_Draws.bin" and "BICq_Posterior_Draws.desc".
#' If the big matrix of posterior draws is already file-backed, \code{BICq_posterior} should
#' specify the path and basename of the appropriate files (again of form "./path/filename"); 
#' the full model 
#' will not be fit again and the big.matrix of 
#' posterior draws will be read using the \code{attach.big.matrix} function of the
#'  \code{bigmemory} package and used in the BIC-ICQ 
#' calcuations. If the appropriate files do not exist or \code{BICq_posterior} 
#' is specified as \code{NULL}, the full model will be fit and the full model posterior
#' draws will be saved as specified above. The algorithm will save 10^4 posterior draws automatically.
#'  
#' Trace details: The value of 0 (default) does not output any extra information. The value of 1
#' additionally outputs the updated coefficients, updated covariance matrix values, and the
#' number of coordinate descent iterations used for the M step for each
#' EM iteration. When pre-screening procedure is used and/or if the BIC-ICQ criterion is
#' requested, trace = 1 gives additional information about the penalties used
#' for the 'full model' fit procedure. If Stan is not used as the E-step sampling mechanism, 
#' the value of 2 outputs all of the above plus gibbs acceptance rate information
#' for the adaptive random walk and independence samplers and the updated proposal standard deviation
#' for the adaptive random walk. 
#' 
#' @return A reference class object of class \code{\link{pglmmObj}} for which many methods are 
#' available (e.g. \code{methods(class = "pglmmObj")}, see ?pglmmObj for additional documentation)
#'  
#' @export
phmmPen_FA = function(formula, data = NULL,
                      offset = NULL, r_estimation = rControl(),
                      fixef_noPen = NULL, penalty = c("MCP","SCAD","lasso"),
                      alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0),
                      optim_options = optimControl(), adapt_RW_options = adaptControl(),
                      trace = 0, tuning_options = selectControl(), 
                      survival_options = survivalControl(),
                      BICq_posterior = NULL,
                      progress = TRUE){
  
  # Specify family
  family = "coxph"
  
  ###########################################################################################
  # Call glmmPen
  ###########################################################################################
  output = glmmPen(formula = formula, data = data, family = family,
                   offset = offset, r_estimation = r_estimation,
                   fixef_noPen = fixef_noPen, penalty = penalty,
                   alpha = alpha, gamma_penalty = gamma_penalty,
                   optim_options = optim_options, adapt_RW_options = adapt_RW_options,
                   trace = trace, tuning_options = tuning_options,
                   survival_options = survival_options,
                   BICq_posterior = BICq_posterior, progress = progress)
  
  
}


#' @title Fit a Proportional Hazards Mixed Model via Monte Carlo Expectation Conditional Minimization (MCECM)
#' using a Piecewise Constant Hazard Mixed Model Approximation
#' 
#' @description \code{phmm} is used to approximate a single proportional hazards mixed model 
#' using a piecewise constant hazard mixed model approximation
#' via Monte Carlo Expectation Conditional Minimization (MCECM). No model selection is performed.
#' 
#' @inheritParams phmmPen
#' @inheritParams glmmPen
#' @param ... additional arguments that could be passed into \code{phmmPen}. 
#' See \code{\link{phmmPen}}
#' for further details.
#' 
#' @details The \code{phmm} function can be used to approximate a single proportional hazards
#' mixed model using a piecewise constant hazard mixed model.
#' While this approach is meant to be used in the case where the user knows which
#' covariates belong in the fixed and random effects and no penalization is required, one is
#' allowed to specify non-zero fixed and random effects penalties using \code{\link{lambdaControl}}
#' and the (...) arguments. The (...) allow for specification of penalty-related arguments; see
#' \code{\link{glmmPen}} and \code{\link{glmmPen_FA}}for details. 
#' For a high dimensional situation, the user may want to fit a 
#' full model using a small penalty for the fixed and random effects and save the posterior
#' draws from this full model for use in any BIC-ICQ calculations during selection within \code{phmmPen}. 
#' Specifying a file name in the 'BICq_posterior' argument will save the posterior draws from the 
#' \code{phmm} model into a big.matrix with this file name, see the Details section of 
#' \code{\link{phmmPen}} for additional details.
#' 
#' 
#' @return A reference class object of class \code{\link{pglmmObj}} for which many methods are 
#' available (e.g. \code{methods(class = "pglmmObj")})
#' 
#' @export
phmm = function(formula, data = NULL, covar = NULL,
                offset = NULL,
                optim_options = optimControl(), adapt_RW_options = adaptControl(),
                trace = 0, tuning_options = lambdaControl(),
                survival_options = survivalControl(),
                progress = TRUE, ...){
  
  # Specify family
  family = "coxph"
  
  # Check that (...) arguments are subsets of phmmPen arguments
  args_extra = list(...)
  args_avail = c("fixef_noPen","penalty","alpha","gamma_penalty","BICq_posterior")
  if(length(args_extra) >= 1){
    if(!(names(args_extra) %in% args_avail)){
      stop("additional arguments provided in '...' input must match phmmPen arguments \n",
           "allowed extra arguments inclue 'fixef_noPen', 'penalty', 'alpha', 'gamma_penalty', 'BICq_posterior' \n",
           "see phmmPen documentation for details")
    }
  }
  
  # Check that tuning_options has the correct input type
  if(!inherits(tuning_options, "lambdaControl")){
    stop("phmm requires the use of lambdaControl for the tuning_options. \n",
         "  In order to use selectControl for model selection, please use phmmPen function")
  }
  
  call = match.call(expand.dots = FALSE)
  
  # Call glmmPen_FA() (tuning_options = lambdaControl() specifies the fitting of a single model 
  #     and not model selection)
  # If ... arguments empty or only includes a few of the args_avail, 
  # glmmPen_FA will use default arguments 
  output = glmmPen(formula = formula, data = data, family = family, covar = covar,
                   offset = offset,
                   optim_options = optim_options,
                   adapt_RW_options = adapt_RW_options, trace = trace,
                   tuning_options = tuning_options, survival_options = survival_options,
                   progress = progress, ...)
  
  output$call = call
  out_object = pglmmObj$new(output)
  
  return(out_object)
  
  
}


#' @title Fit Penalized Proportional Hazards Mixed Models via Monte Carlo Expectation Conditional 
#' Minimization (MCECM) using a Piecewise Constant Hazard Mixed Model Approximation
#' 
#' @description \code{phmmPen_FA} is used to fit penalized proportional hazards mixed models 
#' using a piecewise constant hazard mixed model approximation via 
#' Monte Carlo Expectation Conditional Minimization (MCECM). 
#' The purpose of the function is to perform 
#' variable selection on both the fixed and random effects simultaneously for the
#' piecewise constant hazard mixed model. 
#' \code{phmmPen} selects the best model using 
#' BIC-type selection criteria (see \code{\link{selectControl}} documentation for 
#' further details). 
#' To improve the speed of the algorithm, consider setting
#' \code{var_restrictions} = "fixef" within the \code{\link{optimControl}} options.
#' 
#' @inheritParams glmmPen
#' @param formula a two-sided linear formula object describing both the fixed effects and 
#' random effects part of the model, with the response on the left of a ~ operator and the terms, 
#' separated by + operators, on the right. 
#' The response must be a \code{Surv} object (see
#' \code{\link[survival]{Surv}} from the \code{survival} package).
#' Random-effects terms are distinguished by vertical bars 
#' ("|") separating expression for design matrices from the grouping factor. \code{formula} should be 
#' of the same format needed for \code{\link[lme4]{glmer}} in package \pkg{lme4}. 
#' Only one grouping factor 
#' will be recognized. The random effects covariates need to be a subset of the fixed effects covariates.
#' The offset must be specified outside of the formula in the 'offset' argument.
#' @param survival_options a structure of class "survivalControl" created 
#' from function \code{\link{survivalControl}} that specifies several parameters needed to 
#' properly fit the input survival data using a piecewise constant hazard mixed model. See the 
#' documentation for \code{\link{survivalControl}} for more details on defaults.
#' 
#' @details Argument \code{BICq_posterior} details: If the \code{BIC_option} in \code{\link{selectControl}} 
#' (\code{tuning_options}) is specified 
#' to be 'BICq', this requests the calculation of the BIC-ICQ criterion during the selection
#' process. For the BIC-ICQ criterion to be calculated, a full model assuming a small valued 
#' lambda penalty needs to be fit, and the posterior draws from this full model need to be used. 
#' In order to avoid repetitive calculations of
#' this full model (i.e. if the user wants to re-run \code{phmmPen} with a different
#' set of penalty parameters), a \code{big.matrix} of these 
#' posterior draws will be file-backed as two files: a backing file with extention '.bin' and a 
#' descriptor file with extension '.desc'. The \code{BICq_posterior} argument should contain a 
#' path and a filename with no extension of the form "./path/filename" such that the backingfile and
#' the descriptor file would then be saved as "./path/filename.bin" and "./path/filename.desc", respectively.
#' If \code{BICq_posterior} is set to \code{NULL}, then by default, the backingfile and descriptor
#' file are saved in the working directory as "BICq_Posterior_Draws.bin" and "BICq_Posterior_Draws.desc".
#' If the big matrix of posterior draws is already file-backed, \code{BICq_posterior} should
#' specify the path and basename of the appropriate files (again of form "./path/filename"); 
#' the full model 
#' will not be fit again and the big.matrix of 
#' posterior draws will be read using the \code{attach.big.matrix} function of the
#'  \code{bigmemory} package and used in the BIC-ICQ 
#' calcuations. If the appropriate files do not exist or \code{BICq_posterior} 
#' is specified as \code{NULL}, the full model will be fit and the full model posterior
#' draws will be saved as specified above. The algorithm will save 10^4 posterior draws automatically.
#'  
#' Trace details: The value of 0 (default) does not output any extra information. The value of 1
#' additionally outputs the updated coefficients, updated covariance matrix values, and the
#' number of coordinate descent iterations used for the M step for each
#' EM iteration. When pre-screening procedure is used and/or if the BIC-ICQ criterion is
#' requested, trace = 1 gives additional information about the penalties used
#' for the 'full model' fit procedure. If Stan is not used as the E-step sampling mechanism, 
#' the value of 2 outputs all of the above plus gibbs acceptance rate information
#' for the adaptive random walk and independence samplers and the updated proposal standard deviation
#' for the adaptive random walk. 
#' 
#' @return A reference class object of class \code{\link{pglmmObj}} for which many methods are 
#' available (e.g. \code{methods(class = "pglmmObj")}, see ?pglmmObj for additional documentation)
#'  
#' @export
phmmPen = function(formula, data = NULL, covar = NULL,
                   offset = NULL,
                   fixef_noPen = NULL, penalty = c("MCP","SCAD","lasso"),
                   alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0),
                   optim_options = optimControl(), adapt_RW_options = adaptControl(),
                   trace = 0, tuning_options = selectControl(), 
                   survival_options = survivalControl(),
                   BICq_posterior = NULL,
                   progress = TRUE){
  
  # Specify family
  family = "coxph"
  
  ###########################################################################################
  # Call glmmPen
  ###########################################################################################
  output = glmmPen(formula = formula, data = data, family = family, covar = covar,
                   offset = offset,
                   fixef_noPen = fixef_noPen, penalty = penalty,
                   alpha = alpha, gamma_penalty = gamma_penalty,
                   optim_options = optim_options, adapt_RW_options = adapt_RW_options,
                   trace = trace, tuning_options = tuning_options,
                   survival_options = survival_options,
                   BICq_posterior = BICq_posterior, progress = progress)
  
  
}
