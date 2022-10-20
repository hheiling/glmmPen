
#' @title Fit a Generalized Mixed Model via Monte Carlo Expectation Conditional Minimization (MCECM)
#' 
#' @description \code{glmm} is used to fit a single generalized mixed model via Monte Carlo 
#' Expectation Conditional Minimization (MCECM). Unlike \code{glmmPen}, no model selection 
#' is performed.
#' 
#' @inheritParams glmmPen
#' @param ... additional arguments that could be passed into \code{glmmPen}. See \code{\link{glmmPen}}
#' for further details.
#' 
#' @details The \code{glmm} function can be used to fit a single generalized mixed model.
#' While this approach is meant to be used in the case where the user knows which
#' covariates belong in the fixed and random effects and no penalization is required, one is
#' allowed to specify non-zero fixed and random effects penalties using \code{\link{lambdaControl}}
#' and the (...) arguments. The (...) allow for specification of penalty-related arguments; see
#' \code{\link{glmmPen}} for details. For a high dimensional situation, the user may want to fit a 
#' minimal penalty model using a small penalty for the fixed and random effects and save the posterior
#' samples from this minimal penalty model for use in any BIC-ICQ calculations during selection within \code{glmmPen}. 
#' Specifying a file name in the 'BICq_posterior' argument will save the posterior samples from the 
#' \code{glmm} model into a big.matrix with this file name, see the Details section of 
#' \code{\link{glmmPen}} for additional details.
#' 
#' 
#' @return A reference class object of class \code{\link{pglmmObj}} for which many methods are 
#' available (e.g. \code{methods(class = "pglmmObj")})
#' 
#' @export
glmm = function(formula, data = NULL, family = "binomial", covar = NULL,
                offset = NULL, 
                optim_options = optimControl(), adapt_RW_options = adaptControl(),
                trace = 0, tuning_options = lambdaControl(),
                progress = TRUE, ...){
  
  # Check that (...) arguments are subsets of glmmPen arguments
  args_extra = list(...)
  args_avail = c("fixef_noPen","penalty","alpha","gamma_penalty","BICq_posterior")
  if(length(args_extra) >= 1){
    if(!(names(args_extra) %in% args_avail)){
      stop("additional arguments provided in '...' input must match glmmPen arguments \n",
           "allowed extra arguments inclue 'fixef_noPen', 'penalty', 'alpha', 'gamma_penalty', 'BICq_posterior' \n",
           "see glmmPen documentation for details")
    }
  }
  
  # Check that tuning_options has the correct input type
  if(!inherits(tuning_options, "lambdaControl")){
    stop("glmm requires the use of lambdaControl for the tuning_options. \n",
         "  In order to use selectControl for model selection, please use glmmPen function")
  }
  
  call = match.call(expand.dots = FALSE)
  
  # Call glmmPen() (tuning_options = lambdaControl() specifies the fitting of a single model 
  #     and not model selection)
  # If ... arguments empty or only includes a few of the args_avail, 
  # glmmPen will use default arguments 
  output = glmmPen(formula = formula, data = data, family = family, covar = covar,
                   offset = offset, optim_options = optim_options,
                   adapt_RW_options = adapt_RW_options, trace = trace,
                   tuning_options = tuning_options, progress = progress, ...)
  
  output$call = call
  out_object = pglmmObj$new(output)
  
  return(out_object)
  
  
}


# fD_adj: internal function, used within glmmPen()
# Purpose: (a) run some additional input checks and (b) perform some adjustments to the data
#   to make it easier for the fit algorithm to use the data
# Input: initial output from "glFormula_edit" (an edited version of glFormula() from the lme4 package)
# Output: list with following elements:
## frame: a model frame including all fixed and random covariates, the response, and the 
##    grouping variable
## y, X, Z: response, fixed effects covariates matrix, and random effects covariates model matrix
## group: factor of the grouping variable
## group_num: numeric version of the group factor (e.g. factor with levels "A","B","C" 
##    --> factor with levels 1,2,3
## group_name: character string storing name of variable used for grouping factor
## cnms: vector of names of the random effects covariates
## reTrms: list containing several items relating to the random effects, including: 
##    - Zt: a matrix that will eventually become the Z matrix (after some manipulations)
##    - cnms: see above
##    - flist: a list of grouping factors using inf the random-effects terms
## fixed_vars: vector of variables used for the fixed effects covariates
#' @importFrom stats model.response
fD_adj = function(out){
  frame = out$fr
  y = model.response(frame)
  X = out$X
  reTrms = out$reTrms
  Zt = reTrms$Zt
  cnms = reTrms$cnms
  flist = reTrms$flist
  fixed_vars = out$fixed_vars
  family = out$family
  
  # Check y input
  # If y is not numeric (e.g. a factor), convert to numeric
  if(!is.double(y)) {
    # If y is a character vector, first set y as a factor
    if(is.character(y)){
      y = as.factor(y)
    }
    # Convert y to a numeric vector
    ## If two-level factor, will convert to numeric vector with values {1,2}
    tmp <- try(y <- as.double(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be numeric or able to be coerced to numeric", call.=FALSE)
  }
  
  # Check X input
  if(!is.matrix(X)){
    stop("X must be a matrix \n")
  }else if(typeof(X)=="integer") storage.mode(X) <- "double"
  
  if(nrow(X) != length(y)){
    stop("the dimension of X and y do not match")
  }
  
  # Check that response type matches family
  if(family$family == "binomial"){
    # Check only two response options (binary data)
    if(length(unique(y)) != 2){
      stop("response must be binary for the binomial family; response currently has ", 
           length(unique(y))," different responses")
    }
    # Convert response y to numeric {0,1} if necessary
    if(!all(y %in% c(0,1))){
      # If response input by user is a character or factor, 
      # above "as.double" step should convert response to a
      # factor with levels (1,2). Subtract 1 to get response (0,1)
      if(all(y %in% c(1,2))){
        y = y-1
      }else{ # In case error pops up with glFormula_edit()
        stop("response needs to be coded as a binary variable with levels 0 vs 1")
      }
    }
  }else if(family$family == "poisson"){
    if(!all(floor(y) == y)){
      stop("response must be integer counts for the poisson family")
    }else if(any(y < 0)){
      stop("response must be non-negative integer counts for the poisson family")
    }
  }
  
  # For now, retrict algorithm to only handle single grouping factor
  if(length(flist) > 1){
    stop("procedure can only handle one group")
  }else{
    group = flist[[1]]
    group_name = names(flist)
    cnms = cnms[[1]]
  }
  
  # Convert alpha-numeric group to numeric group (e.g. "A","B","C" to 1,2,3) for ease of 
  # computation within fit_dat() algorithm.
  ## Even if already numeric, convert to levels 1,2,3,... (consecutive integers)
  group_num = factor(group, labels = 1:nlevels(group))
  
  d = nlevels(group)
  numvars = nrow(Zt)/d
  Z = Matrix(0, nrow = ncol(Zt), ncol = nrow(Zt), sparse = TRUE)
  # mkReTrms Zt rows: organized first by level of group, then vars 
  # Want Z columns organized first by vars, then levels of group within vars
  for(lev in 1:numvars){
    Z[,(d*(lev-1)+1):(d*lev)] = Matrix::t(Zt[seq(lev, by = numvars, length.out = d),])
  }
  colnames(Z) = str_c(rep(cnms, each = d), ":", rep(levels(group), length(cnms)))
  
  ## Make sure colnames random effects subset of colnames of fixed effects model matrix X
  if(sum(!(cnms %in% colnames(X))) > 0){
    stop("random effects must be a subset of fixed effects: names of random effect variables must be a subset of names of fixed effect variables; \n",
         "fixef effects: ", paste0(colnames(X), sep=" "), " \n", "random effects: ", paste0(cnms, sep=" "))
  }
  
  if(!any(cnms == "(Intercept)")){
    stop("Model requires a random intercept term. Current random effects: ", cnms)
  }
  
  return(list(frame = frame, y = y, X = X, Z = Z, group = group, 
              group_num = group_num, cnms = cnms,
              group_name = group_name, reTrms = reTrms, fixed_vars = fixed_vars))
  
} # End fD_adj() function

#' @title Fit Penalized Generalized Mixed Models via Monte Carlo Expectation Conditional 
#' Minimization (MCECM)
#' 
#' @description \code{glmmPen} is used to fit penalized generalized mixed models via Monte Carlo Expectation 
#' Conditional Minimization (MCECM). The purpose of the function is to perform 
#' variable selection on both the fixed and random effects simultaneously for the
#' generalized linear mixed model. \code{glmmPen} selects the best model using 
#' BIC-type selection criteria (see \code{\link{selectControl}} documentation for 
#' further details)
#' 
#' @param formula a two-sided linear formula object describing both the fixed effects and 
#' random effects part of the model, with the response on the left of a ~ operator and the terms, 
#' separated by + operators, on the right. Random-effects terms are distinguished by vertical bars 
#' ("|") separating expression for design matrices from the grouping factor. \code{formula} should be 
#' of the same format needed for \code{\link[lme4]{glmer}} in package \pkg{lme4}. Only one grouping factor 
#' will be recognized. The random effects covariates need to be a subset of the fixed effects covariates.
#' The offset must be specified outside of the formula in the 'offset' argument.
#' @param data an optional data frame containing the variables named in \code{formula}. If \code{data} is 
#' omitted, variables will be taken from the environment of \code{formula}.
#' @param family a description of the error distribution and link function to be used in the model. 
#' Currently, the \code{glmmPen} algorithm allows the Binomial, Gaussian, and Poisson families
#' with canonical links only.
#' @param covar character string specifying whether the covariance matrix should be unstructured
#' ("unstructured") or diagonal with no covariances between variables ("independent").
#' Default is set to \code{NULL}. If \code{covar} is set to \code{NULL} and the number of random effects
#' predictors (not including the intercept) is 
#' greater than or equal to 10 (i.e. high dimensional), then the algorithm automatically assumes an 
#' independent covariance structure and \code{covar} is set to "independent". Otherwise if \code{covar}
#' is set to \code{NULL} and the number of random effects predictors is less than 10, then the
#' algorithm automatically assumes an unstructured covariance structure and \code{covar} is set to "unstructured".
#' @param offset This can be used to specify an \emph{a priori} known component to be included in the 
#' linear predictor during fitting. Default set to \code{NULL} (no offset). If the data 
#' argument is not \code{NULL}, this should be a numeric vector of length equal to the 
#' number of cases (the length of the response vector). 
#' If the data argument specifies a data.frame, the offset
#' argument should specify the name of a column in the data.frame.
#' @param fixef_noPen Optional vector of 0's and 1's of the same length as the number of fixed effects covariates
#' used in the model. Value 0 indicates the variable should not have its fixed effect coefficient
#' penalized, 1 indicates that it can be penalized. Order should correspond to the same order of the 
#' fixed effects given in the formula.
#' @param penalty character describing the type of penalty to use in the variable selection procedure.
#' Options include 'MCP', 'SCAD', and 'lasso'. Default is MCP penalty. If the random effect covariance
#' matrix is "unstructured", then a group MCP, group SCAD, or group LASSO penalty is used on the 
#' random effects coefficients. See Breheny and Huang (2011) <doi:10.1214/10-AOAS388> 
#' and Breheny and Huang (2015) <doi:10.1007/s11222-013-9424-2> for details of these penalties.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions 
#' from the MCP/SCAD/LASSO penalty and the ridge, or L2, penalty. \code{alpha=1} is equivalent to 
#' the MCP/SCAD/LASSO penalty, while \code{alpha=0} is equivalent to ridge regression. However,
#' \code{alpha=0} is not supported; \code{alpha} may be arbitrarily small, but not exactly zero
#' @param gamma_penalty The scaling factor of the MCP and SCAD penalties. Not used by LASSO penalty.
#' Default is 4.0 for SCAD and 3.0 for MCP. See Breheny and Huang (2011) <doi:10.1214/10-AOAS388> 
#' and Breheny and Huang (2015) <doi:10.1007/s11222-013-9424-2> for further details.
#' @param optim_options a structure of class "optimControl" created 
#' from function \code{\link{optimControl}} that specifies several optimization parameters. See the 
#' documentation for \code{\link{optimControl}} for more details on defaults.
#' @param adapt_RW_options a list of class "adaptControl" from function \code{\link{adaptControl}} 
#' containing the control parameters for the adaptive random walk Metropolis-within-Gibbs procedure. 
#' Ignored if \code{\link{optimControl}} parameter \code{sampler} is set to "stan" (default) or "independence".
#' @param tuning_options a list of class "selectControl" or "lambdaControl" resulting from 
#' \code{\link{selectControl}} or \code{\link{lambdaControl}} containing additional control parameters.
#' When function \code{glmm} is used,the algorithm may be run using one specific set of
#' penalty parameters \code{lambda0} and \code{lambda1} by specifying such values in \code{lambdaControl()}. 
#' The default for \code{glmm} is to run the model fit with no penalization (\code{lambda0} = \code{lambda1} = 0).
#' When function \code{glmmPen} is run, \code{tuning_options} is specified using \code{selectControl()}. 
#' See the \code{\link{lambdaControl}} and \code{\link{selectControl}} documentation for further details.
#' @param BICq_posterior an optional character string expressing the path and file 
#' basename of a file combination that 
#' will file-back or currently file-backs a \code{big.matrix} of the posterior samples from the 
#' minimal penalty model used for the BIC-ICQ calculation used for model selection. T
#' (BIC-ICQ reference: Ibrahim et al (2011)
#' <doi:10.1111/j.1541-0420.2010.01463.x>). 
#' If this argument is
#' specified as \code{NULL} (default) and BIC-ICQ calculations are requested (see \code{\link{selectControl}})
#' for details), the posterior samples
#' will be saved in the file combination 'BICq_Posterior_Draws.bin' and 'BICq_Posterior_Draws.desc'
#' in the working directory.
#' See 'Details' section for additional details about the required format of \code{BICq_posterior}
#' and the file-backed big matrix. 
#' @param trace an integer specifying print output to include as function runs. Default value is 0. 
#' See Details for more information about output provided when trace = 0, 1, or 2.
#' @param progress a logical value indicating if additional output should be given showing the
#' progress of the fit procedure. If \code{TRUE}, such output includes iteration-level information
#' for the fit procedure (iteration number EM_iter,
#' number of MCMC samples nMC, average Euclidean distance between current coefficients and coefficients
#' from t--defined in \code{\link{optimControl}}--iterations back EM_conv, 
#' and number of non-zero fixed and random effects covariates
#' not including the intercept). Additionally, \code{progress = TRUE}
#' gives some other information regarding the progress of the variable selection 
#' procedure, including the model selection criteria and log-likelihood estimates
#' for each model fit.
#' Default is \code{TRUE}.
#' 
#' @details Argument \code{BICq_posterior} details: If the \code{BIC_option} in \code{\link{selectControl}} 
#' (\code{tuning_options}) is specified 
#' to be 'BICq', this requests the calculation of the BIC-ICQ criterion during the selection
#' process. For the BIC-ICQ criterion to be calculated, a minimal penalty model assuming a small valued 
#' lambda penalty needs to be fit, and the posterior samples from this minimal penalty model need to be used. 
#' In order to avoid repetitive calculations of
#' this minimal penalty model (i.e. if secondary rounds of selection are desired in 
#' \code{\link{glmmPen_FineSearch}} or if the user wants to re-run \code{glmmPen} with a different
#' set of penalty parameters), a \code{big.matrix} of these 
#' posterior samples will be file-backed as two files: a backing file with extension '.bin' and a 
#' descriptor file with extension '.desc'. The \code{BICq_posterior} argument should contain a 
#' path and a filename of the form "./path/filename" such that the backingfile and
#' the descriptor file would then be saved as "./path/filename.bin" and "./path/filename.desc", respectively.
#' If \code{BICq_posterior} is set to \code{NULL}, then by default, the backingfile and descriptor
#' file are saved in the working directory as "BICq_Posterior_Draws.bin" and "BICq_Posterior_Draws.desc".
#' If the big matrix of posterior samples is already file-backed, \code{BICq_posterior} should
#' specify the path and basename of the appropriate files (again of form "./path/filename"); 
#' the minimal penalty model 
#' will not be fit again and the big.matrix of 
#' posterior samples will be read using the \code{attach.big.matrix} function of the
#'  \code{bigmemory} package and used in the BIC-ICQ 
#' calculations. If the appropriate files do not exist or \code{BICq_posterior} 
#' is specified as \code{NULL}, the minimal penalty model will be fit and the minimal penalty model posterior
#' samples will be saved as specified above. The algorithm will save 10^4 posterior samples automatically.
#'  
#' Trace details: The value of 0 (default) does not output any extra information. The value of 1
#' additionally outputs the updated coefficients, updated covariance matrix values, and the
#' number of coordinate descent iterations used for the M step for each
#' EM iteration. When pre-screening procedure is used and/or if the BIC-ICQ criterion is
#' requested, trace = 1 gives additional information about the penalties used
#' for the 'minimal penalty model' fit procedure. If Stan is not used as the E-step sampling mechanism, 
#' the value of 2 outputs all of the above plus gibbs acceptance rate information
#' for the adaptive random walk and independence samplers and the updated proposal standard deviation
#' for the adaptive random walk. 
#' 
#' @return A reference class object of class \code{\link{pglmmObj}} for which many methods are 
#' available (e.g. \code{methods(class = "pglmmObj")}, see ?pglmmObj for additional documentation)
#'  
#' @importFrom stringr str_to_lower str_c str_detect
#' @importFrom Matrix Matrix
#' @importFrom bigmemory write.big.matrix attach.big.matrix
#' @importFrom stats model.offset na.omit
#' @import bigmemory Rcpp
#' @export
glmmPen = function(formula, data = NULL, family = "binomial", covar = NULL,
                   offset = NULL,
                   fixef_noPen = NULL, penalty = c("MCP","SCAD","lasso"),
                   alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0),
                   optim_options = optimControl(), adapt_RW_options = adaptControl(),
                   trace = 0, tuning_options = selectControl(), BICq_posterior = NULL,
                   progress = TRUE){
  
  ###########################################################################################
  # Input argument checks and modifications
  ###########################################################################################
  
  # Input modification and restriction for family
  family_info = family_export(family)
  fam_fun = family_info$family_fun
  
  # Check penalty parameters
  penalty_check = checkPenalty(penalty, gamma_penalty, alpha)
  penalty = penalty_check$penalty
  gamma_penalty = penalty_check$gamma_penalty
  alpha = penalty_check$alpha
  
  # Check covariance matrix specification
  covar = checkCovar(covar)
  
  # Check BICq_posterior path is appropriate
  checkBICqPost(BICq_posterior)
  
  # Check that *_options arguments are the correct types
  if(!inherits(adapt_RW_options, "adaptControl")){ 
    stop("adapt_RW_options must be of class 'adaptControl', see adaptControl documentation")
  }
  
  if(!inherits(tuning_options, "pglmmControl")){
    stop("tuning_option parameter must be a list of type lambdaControl() or selectControl()")
  }
  
  if(!inherits(optim_options, "optimControl")){ 
    stop("optim_options must be of class 'optimControl' (see optimControl documentation)")
  }
  
  out = NULL
  
  # Convert formula and data to useful forms
  # Code glFormula_edit() edited version of glFormula() from lme4 package, 
  #   see code in "formula_data_edits.R"
  # fD_out = 'formula-data output'
  # fD_out0: list with elements formula, fr (model.frame), X (fixed effects matrix), 
  #   reTrms (list with a random effects matrix Zt (needs to be adjusted)), 
  #   cnms (names of random effects), and flist (list of the groups) 
  # Note: restrict na.action to be na.omit
  fD_out0 = glFormula_edit(formula = formula, data = data, family = fam_fun, subset = NULL,
                           weights = NULL, na.action = na.omit, offset = offset)
  
  # Perform additional checks/restrictions/modifications of data input
  # See fD_adj() function earlier in this document
  fD_out = fD_adj(out = fD_out0)
  
  # If offset = NULL, then set offset as arbitrary vector of 0s
  if(is.null(model.offset(fD_out$frame))){
    offset = rep(0, length(fD_out$y))
  }else{
    offset = model.offset(fD_out$frame)
    if((!is.numeric(offset)) | (length(offset) != length(y))){
      stop("offset must be a numeric vector of the same length as y")
    }
  }
  
  # Standardize X and Z - see XZ_std() function later in this document
  std_out = XZ_std(fD_out)
  
  # data_input: list of main data for fit algorithm
  data_input = list(y = fD_out$y, X = std_out$X_std, Z = std_out$Z_std, group = fD_out$group_num)
  
  # Names of covariates for fixed effects, random effects, and group factor
  coef_names = list(fixed = colnames(fD_out$X), random = fD_out$cnms, group = fD_out$group_name)
  
  # Identify fixed effects that should not be penalized in select_tune or fit_dat (if any)
  ## For now, do not allow grouping of fixed effects (i.e. cannot penalize fixed effects as groups of covariates)
  # group_X: 0 indicates intercept or other covariates that should not be penalized
  #   positive integer indicates that the fixed effect can be penalized
  if(is.null(fixef_noPen)){
    group_X = 0:(ncol(data_input$X)-1)
  }else if(is.numeric(fixef_noPen)){
    if(length(fixef_noPen) != (ncol(data_input$X)-1)){
      stop("length of fixef_noPen must match number of fixed effects covariates")
    }
    if(sum(fixef_noPen == 0) == 0){
      group_X = 0:(ncol(data_input$X)-1)
    }else{
      ones = which(fixef_noPen == 1) + 1
      # Sequence: consecutive integers (1 to number of fixed effects to potentially penalize)
      sq = 1:length(ones)
      # Initialize as all 0
      group_X = rep(0, times = ncol(data_input$X))
      # for appropriate fixed effects, set to the appropriate postive integer
      group_X[ones] = sq
    }
  }
  
  # Store call for output object
  call = match.call(expand.dots = F)
  
  # If covar specified as NULL, recommend a covariance structure based on the size of the 
  # random effects.
  if(is.null(covar) & (ncol(data_input$Z)/nlevels(data_input$group) >= 11)){
    message("Setting random effect covariance structure to 'independent'")
    covar = "independent"
  }else if(is.null(covar) & (ncol(data_input$Z)/nlevels(data_input$group) < 11)){
    message("Setting random effect covariance structure to 'unstructured'")
    covar = "unstructured"
  }
  
  if((covar == "unstructured") & (ncol(data_input$Z)/nlevels(data_input$group) >= 11)){
    warning("The random effect covariance matrix is currently specified as 'unstructured'. 
            Due to dimension of sigma covariance matrix, we suggest using covar = 'independent' to simplify computation",
            immediate. = TRUE)
  }
  
  ###########################################################################################
  # Fitting of model
  ###########################################################################################
  
  if(inherits(tuning_options, "lambdaControl")){
    
    ###########################################################################################
    # Fit single model specified by glmm()
    ###########################################################################################
    
    # Default: penalty parameters lambda0 and lambda1 both = 0.
    # Checks of penalty parameters
    lambda0 = tuning_options$lambda0
    lambda1 = tuning_options$lambda1
    if(lambda0 < 0 | lambda1 < 0){
      stop("lambda0 and lambda1 cannot be negative")
    }
    
    # If any optim_option arguments initialized as NULL, give defaults
    ## See bottom of "control_options.R" for function optim_recommend()
    optim_options = optim_recommend(optim_options = optim_options, family = family,
                                    q = ncol(fD_out$Z)/nlevels(data_input$group), select = F)
    
    # Recommend starting variance  for random effects covariance matrix (if necessary)
    if(optim_options$var_start == "recommend"){
      var_start = var_init(data_input, fam_fun)
      optim_options$var_start = var_start
    }
    
    # Extract variables from optimControl
    conv_EM = optim_options$conv_EM
    conv_CD = optim_options$conv_CD
    nMC_burnin = optim_options$nMC_burnin
    nMC = optim_options$nMC
    nMC_max = optim_options$nMC_max
    maxitEM = optim_options$maxitEM
    maxit_CD = optim_options$maxit_CD
    M = optim_options$M
    t = optim_options$t
    mcc = optim_options$mcc
    sampler = optim_options$sampler
    var_start = optim_options$var_start
    
    # Call fit_dat function to fit the model
    # fit_dat function found in "/R/fit_dat.R" file
    
    fit = fit_dat(dat = data_input, lambda0 = lambda0, lambda1 = lambda1, 
                    conv_EM = conv_EM, conv_CD = conv_CD,
                    family = fam_fun, offset_fit = offset, trace = trace, 
                    group_X = group_X, penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                    nMC_burnin = nMC_burnin, nMC = nMC, nMC_max = nMC_max,
                    t = t, mcc = mcc, maxitEM = maxitEM, maxit_CD = maxit_CD,
                    M = M, sampler = sampler, adapt_RW_options = adapt_RW_options,
                    covar = covar, var_start = var_start, logLik_calc = FALSE,
                    checks_complete = TRUE, progress = progress)
    
    selection_results = matrix(NA, nrow = 1, ncol = 1)
    
    ranef_keep = NULL
    
    # If relevant, save posterior samples for later BIC-ICQ calculations
    if(!is.null(BICq_posterior)){
      
      y = data_input$y
      X = data_input$X
      Z = data_input$Z
      group = data_input$group
      d = nlevels(data_input$group)
      
      coef = fit$coef
      J = fit$J
      gamma = fit$Gamma_mat
      # Znew2: For each group k = 1,...,d, calculate Znew2 = Z %*% gamma (calculate for each group individually)
      # Used within E-step
      Znew2 = Z
      for(j in 1:d){
        Znew2[group == j,seq(j, ncol(Z), by = d)] = Z[group == j,seq(j, ncol(Z), by = d)]%*%gamma
      }
      Estep_out = E_step(coef = coef, ranef_idx = sum(diag(fit$sigma) > 0), 
                         y=y, X=X, Znew2=Znew2, group=group, offset_fit = offset,
                         nMC=nMC, nMC_burnin=nMC_burnin, family=fam_fun$family, link=fam_fun$link, 
                         phi=1.0, sig_g=out$sigma_gaus,
                         sampler=sampler, d=d, uold=fit$u_init, proposal_SD=fit$proposal_SD, 
                         batch=fit$updated_batch, batch_length=adapt_RW_options$batch_length, 
                         offset_increment=adapt_RW_options$offset_increment, 
                         trace=trace)
      u0 = attach.big.matrix(Estep_out$u0)
      # File-back the posterior samples of the minimal penalty model
      ufull_big_tmp = as.big.matrix(u0[,],
                                    backingpath = dirname(BICq_posterior),
                                    backingfile = sprintf("%s.bin",basename(BICq_posterior)),
                                    descriptorfile = sprintf("%s.desc",basename(BICq_posterior)))
      rm(ufull_big_tmp)
    }
    
  }else if(inherits(tuning_options, "selectControl")){
    
    ###########################################################################################
    # Perform variable and model selection based on specifications from glmmPen()
    ###########################################################################################
    
    # Extract or calculate penalty parameter values for grid search
    # See LambdaSeq() function later in this document
    ## Fixed effects penalty parameters
    if(is.null(tuning_options$lambda0_seq)){
      lambda0_seq = LambdaSeq(X = data_input$X[,-1,drop=F], y = data_input$y, family = fam_fun$family, 
                                alpha = alpha, nlambda = tuning_options$nlambda, penalty.factor = fixef_noPen,
                                lambda.min = tuning_options$lambda.min)
    }else{
      lambda0_seq = tuning_options$lambda0_seq
    }
    ## Random effects penalty parameters
    if(is.null(tuning_options$lambda1_seq)){
      lambda1_seq = LambdaSeq(X = data_input$X[,-1,drop=F], y = data_input$y, family = fam_fun$family, 
                                alpha = alpha, nlambda = tuning_options$nlambda, penalty.factor = fixef_noPen,
                                lambda.min = tuning_options$lambda.min)
      # Extract or calculate penalty parameters for pre-screening step
      lambda.min.presc = tuning_options$lambda.min.presc
      if(is.null(lambda.min.presc)){ 
        if((ncol(data_input$Z)/nlevels(data_input$group) <= 11)){ # number random effects (including intercept)
          lambda.min.presc = 0.01
        }else{
          lambda.min.presc = 0.05
        }
      }
      
      # Minimum lambda values to use for pre-screening and BIC-ICQ minimal penalty model fit
      # Fixed effects: use minimum fixed effects penalty parameters
      # Random effects: calculate penalty to use using lambda.min.presc 
      ##    (penalty = lambda.min.presc * max random effect penalty parameter)
      full_ranef = LambdaSeq(X = data_input$X[,-1,drop=F], y = data_input$y, family = fam_fun$family, 
                             alpha = alpha, nlambda = 2, penalty.factor = fixef_noPen,
                             lambda.min = lambda.min.presc)
      lambda.min.full = c(min(lambda0_seq), min(full_ranef))
    }else{
      lambda1_seq = tuning_options$lambda1_seq
      # Minimum lambda values to use for pre-screening and BIC-ICQ minimal penalty model fit
      # Fixed and random effects: choose minimum penalty for respective sequence
      lambda.min.full = c(min(lambda0_seq), min(lambda1_seq))
    }
    
    names(lambda.min.full) = c("fixef","ranef")
    
    # Extract additional arguments needed for selection
    BIC_option = tuning_options$BIC_option
    if(lambda.min.full["ranef"] == 0){
      message("minimum random effect penalty set to 0, will consequently skip pre-screening procedure")
      pre_screen = FALSE
    }else{
      pre_screen = tuning_options$pre_screen # TRUE vs FALSE
    }
    logLik_calc = tuning_options$logLik_calc # TRUE vs FALSE
    
    if((is.null(BICq_posterior)) & (BIC_option == "BICq")){
      BICq_post_file = "BICq_Posterior_Draws"
      # Check to see if files already exist in working directory
      # If they exist, remove files. 
      # Since this name is generic for all cases where the input of BICq_posterior was NULL,
      #   safe to assume user forgot that these files were saved previously.
      # Delete to avoid potential issues of using wrong posterior draws for the model
      if(file.exists("BICq_Posterior_Draws.bin") | file.exists("BICq_Posterior_Draws.desc")){
        message("Removing existing BICq_Posterior_Draws files with saved posterior samples")
        if(file.exists("BICq_Posterior_Draws.bin")) file.remove("BICq_Posterior_Draws.bin")
        if(file.exists("BICq_Posterior_Draws.desc")) file.remove("BICq_Posterior_Draws.desc")
      }
    }else{
      BICq_post_file = BICq_posterior
    }
    
    
    if(tuning_options$search == "abbrev"){
      
      ###########################################################################################
      # Two-stage (abbreviated) grid search
      ###########################################################################################
      
      # Fit the following set of models:
      ## fixed effects penalty: fixed at lambda_min
      ## random effects penalty: all lambda1 values (from min to max penalty parameters)
      
      lam_min = min(lambda0_seq)
      
      # Fit first stage of abbreviated grid search
      if(progress == TRUE) message("Start of stage 1 of abbreviated grid search")
      fit_fixfull = select_tune(dat = data_input, offset = offset, family = family,
                                lambda0_seq = lam_min, lambda1_seq = lambda1_seq,
                                penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                                group_X = group_X, trace = trace,
                                adapt_RW_options = adapt_RW_options, 
                                optim_options = optim_options, covar = covar, logLik_calc = logLik_calc,
                                BICq_calc = (BIC_option == "BICq"),
                                BIC_option = BIC_option, BICq_posterior = BICq_post_file, 
                                checks_complete = TRUE, pre_screen = pre_screen, 
                                lambda.min.full = lambda.min.full, stage1 = TRUE,
                                progress = progress)
      if(progress == TRUE) message("End of stage 1 of abbreviated grid search")
      
      # Extract stage 1 results needed for stage 2
      
      res_pre = fit_fixfull$results
      coef_pre = fit_fixfull$coef
      
      # Choose the best model from the above 'fit_fixfull' models
      opt_pre = matrix(res_pre[which.min(res_pre[,BIC_option]),], nrow = 1)
      # optimum penalty parameter on random effects from above first step
      lam_ref = opt_pre[,2]
      # Extract coefficient and posterior results from 'best' model from first step
      coef_old = coef_pre[which(res_pre[,2] == lam_ref),]
      u_init = fit_fixfull$out$u_init
      
      # Extract other relevant info
      ## only consider random effects that are non-zero or above 10^-2 in the optimal random effect model
      # ranef_keep = fit_fixfull$ranef_keep
      vars = diag(fit_fixfull$out$sigma)
      
      # Fit second stage of 'abbreviated' model fit
      if(progress == TRUE) message("Start of stage 2 of abbreviated grid search")
      fit_select = select_tune(dat = data_input, offset = offset, family = family,
                               lambda0_seq = lambda0_seq, lambda1_seq = lam_ref,
                               penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                               group_X = group_X, trace = trace,
                               coef_old = coef_old, u_init = u_init,
                               adapt_RW_options = adapt_RW_options, 
                               optim_options = optim_options, covar = covar, 
                               logLik_calc = logLik_calc, BICq_calc = (BIC_option == "BICq"),
                               BIC_option = BIC_option, BICq_posterior = BICq_post_file, 
                               checks_complete = TRUE, pre_screen = FALSE,
                               ranef_keep = as.numeric((vars > 0)), 
                               lambda.min.full = lambda.min.full, stage1 = FALSE,
                               progress = progress)
      if(progress == TRUE) message("End of stage 2 of abbreviated grid search")
      
      resultsA = rbind(fit_fixfull$results, fit_select$results)
      coef_results = rbind(fit_fixfull$coef, fit_select$coef)
      
    }else if(tuning_options$search == "full_grid"){
      
      ###########################################################################################
      # Full grid search
      ###########################################################################################
      
      fit_select = select_tune(dat = data_input, offset = offset, family = family,
                               lambda0_seq = lambda0_seq, lambda1_seq = lambda1_seq,
                               penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                               group_X = group_X, trace = trace,
                               adapt_RW_options = adapt_RW_options, 
                               optim_options = optim_options, covar = covar, logLik_calc = logLik_calc,
                               BICq_calc = (tuning_options$BIC_option == "BICq"),
                               BIC_option = BIC_option, BICq_posterior = BICq_post_file, 
                               checks_complete = TRUE, pre_screen = pre_screen, 
                               lambda.min.full = lambda.min.full, stage1 = FALSE,
                               progress = progress)
      
      resultsA = fit_select$results
      coef_results = fit_select$coef
      
    } # End if-else tuning_options$search == 'abbrev'
    
    ###########################################################################################
    # Extract output from best model
    ###########################################################################################
    
    # Extract best model
    fit = fit_select$out
    
    # Unstandardize fixed effect coefficient results
    beta_results = matrix(0, nrow = nrow(coef_results), ncol = ncol(data_input$X))
    beta_results[,1] = coef_results[,1] - apply(coef_results[,2:ncol(data_input$X),drop=F], MARGIN = 1, FUN = function(x) sum(x * std_out$X_center / std_out$X_scale))
    for(i in 1:nrow(beta_results)){
      beta_results[,-1] = coef_results[,2:ncol(data_input$X),drop=F] / std_out$X_scale
    }
    
    gamma_results = coef_results[,(ncol(data_input$X)+1):ncol(coef_results),drop=F]
    
    # Combine all relevant selection results
    selection_results = cbind(resultsA,beta_results,gamma_results)
    colnames(selection_results) = c(colnames(resultsA), coef_names$fixed, 
                                    str_c("Gamma",0:(ncol(gamma_results)-1)))
    
    # Find optimum/best model from model selection
    optim_results = matrix(selection_results[which.min(selection_results[,BIC_option]),], nrow = 1)
    
    colnames(optim_results) = colnames(selection_results)
    
    # Record what final optimization arguments were used (could have been modified within
    #   the select_tune() selection procedure if input optimization options were set to NULL)
    optim_options = fit_select$optim_options
    
    if(tuning_options$search == "full_grid"){
      ranef_keep = fit_select$ranef_keep
    }else if(tuning_options$search == "abbrev"){
      ranef_keep = fit_fixfull$ranef_keep
    }
    
    
  } # End if-else inherits(tuning_options)
  
  ###########################################################################################
  # Final output preparation
  ###########################################################################################
  
  sampling = switch(optim_options$sampler, stan = "Stan", 
                    random_walk = "Metropolis-within-Gibbs Adaptive Random Walk Sampler",
                    independence = "Metropolis-within-Gibbs Independence Sampler")
  
  if(penalty == "lasso"){
    gamma_penalty = NULL
  }
  
  # Perform a final E-step
  ## Use results to calculate logLik and posterior samples, and save samples for MCMC diagnostics
  Estep_final_complete = FALSE
  # if problematic model fit, do not perform final calculations for 
  #   logLik and posterior samples and give NA values for appropriate output
  if(!is.null(fit$warnings)){
    if(str_detect(fit$warnings, "coefficient values diverged") | str_detect(fit$warnings, "coefficient estimates contained NA values")){
     
      Estep_out = list(u0 = matrix(NA, nrow = 1, ncol = ncol(data_input$Z)), 
                       post_modes = rep(NA, times = ncol(data_input$Z)), 
                       post_out = matrix(NA, nrow = 1, ncol = ncol(data_input$Z)), 
                       u_init = matrix(NA, nrow = 1, ncol = ncol(data_input$Z)),
                       ll = NA, BICh = NA, BIC = NA, BICNgrp = NA)
      Estep_final_complete = TRUE
    }
  }
  # if model fit completed without serious issues, perform final computations
  if(Estep_final_complete == FALSE){
    Estep_out = E_step_final(dat = data_input, offset_fit = offset, fit = fit, optim_options = optim_options, 
                             fam_fun = fam_fun, extra_calc = TRUE, 
                             adapt_RW_options = adapt_RW_options, trace = trace, 
                             progress = progress)
  }
  
  
  # save optim_results:
  if(inherits(tuning_options, "lambdaControl")){
    # If performed a single model fit, save the final model output as the 'optimal' model
    # in order to standarize the output 
    optim_results = c(fit$lambda0, fit$lambda1, Estep_out$BICh, Estep_out$BIC, fit$BICq,
                      Estep_out$BICNgrp, Estep_out$ll,
                      sum(fit$coef[2:ncol(data_input$X)] != 0),
                      sum(diag(fit$sigma[-1,-1,drop=F]) !=0),
                      sum(fit$coef != 0))
    optim_results = matrix(optim_results, nrow = 1)
    colnames(optim_results) = c("lambda0","lambda1","BICh","BIC","BICq","BICNgrp",
                                "LogLik","Non0 Fixef","Non0 Ranef","Non0 Coef")
    
  }else if(inherits(tuning_options, "selectControl")){
    # Save log-likelihood and several BIC-derived quantities for the best model.
    # Need to save these here because if BIC_option set to "BICq" (BIC-ICQ),
    # these quantities not calculated for each model during the selection process
    # (to improve speed of algorithm)
    optim_results[,c("BICh","BIC","BICNgrp","LogLik")] = c(Estep_out$BICh, Estep_out$BIC, 
                                                           Estep_out$BICNgrp, Estep_out$ll)
    
  }
  
  
  # Create pglmmObj object
  ## output = list object plugged into function used to create pglmmObj object
  output = c(fit, 
             list(Estep_out = Estep_out, formula = formula, y = fD_out$y, fixed_vars = fD_out$fixed_vars,
                  X = fD_out$X, Z_std = std_out$Z_std, group = fD_out$reTrms$flist,
                  coef_names = coef_names, family = fam_fun, covar = covar,
                  offset = offset, frame = fD_out$frame, 
                  sampling = sampling, std_out = std_out, 
                  selection_results = selection_results, optim_results = optim_results,
                  penalty = penalty, gamma_penalty = gamma_penalty, alpha = alpha, 
                  fixef_noPen = fixef_noPen, ranef_keep = ranef_keep,
                  control_options = list(optim_options = optim_options, tuning_options = tuning_options,
                                         adapt_RW_options = adapt_RW_options)))

  # If glmm() function, output the output list to the outer glmm() function
  # If glmmPen() function, create pglmmObj object and output this pglmmObj object directly
  if(inherits(tuning_options, "lambdaControl")){
    return(output)
  }else if(inherits(tuning_options, "selectControl")){
    output$call = call
    out_object = pglmmObj$new(output)
    return(out_object)
  }
  
  
}

###########################################################################################
# Standardization of covariates
###########################################################################################
# Standardize fixed effects covariates, use these standardized values in the 
# Z random effects covariates model matrix, and save standarization information

#' @importFrom ncvreg std
XZ_std = function(fD_out){
  
  group_num = fD_out$group_num
  
  # Standardize X - ncvreg::std method
  X = fD_out$X
  X_noInt_std = std(X[,-1, drop = FALSE])
  X_std = cbind(1, X_noInt_std)
  X_center = attr(X_noInt_std, "center")
  X_scale = attr(X_noInt_std, "scale")
  # Note: X_noInt_std = (X[,-1] - X_center) / X_scale
  
  var_subset = which(colnames(X) %in% fD_out$cnms)
  Z_center = X_center[var_subset]
  Z_scale = X_scale[var_subset]
  
  # Standardize Z using X_std output
  Z_sparse = fD_out$Z
  d = nlevels(group_num)
  num_vars = ncol(Z_sparse) / d
  
  Z_std = Z_sparse
  
  if("(Intercept)" %in% fD_out$cnms){
    v_start = 2
  }else{
    stop("Model requires the assumption of a random intercept")
  }
  
  if(length(fD_out$cnms) > 1){ # If have more than just a random intercept
    if(num_vars > 2){
      idx_seq = v_start:num_vars
    }else{
      idx_seq = v_start
    }
    for(v in idx_seq){ 
      cols = seq(from = (v - 1)*d + 1, to = v*d, by = 1)
      for(k in 1:nlevels(group_num)){
        ids = which(group_num == k)
        Z_std[ids, cols[k]] = (Z_sparse[ids, cols[k]] - Z_center[v-1]) / Z_scale[v-1]
      }
    }
  }
  
  
  colnames(X_std) = colnames(X)
  colnames(Z_std) = colnames(fD_out$Z)
  
  return(list(X_std = X_std, Z_std = Z_std, X_center = X_center, X_scale = X_scale,
              Z_center = Z_center, Z_scale = Z_scale))
}

###########################################################################################
# Penalty parameter calculation
## Note: significant parts of code borrowed from ncvreg
###########################################################################################

#' @title Calculation of Penalty Parameter Sequence (Lambda Sequence)
#' 
#' @description Calculates the sequence of penalty parameters used in the model selection procedure.
#' This function calls functions from package \code{ncvreg}. 
#' 
#' @inheritParams glmmPen
#' @inheritParams lambdaControl
#' @param X matrix of standardized fixed effects (see \code{std} function in \code{ncvreg} 
#' documenation). X should not include intercept.
#' @param y numeric vector of response values
#' @param nlambda positive integer specifying number of penalty parameters (lambda) with 
#' which to fit a model.
#' @param penalty.factor an optional numeric vector equal to the \code{fixef_noPen} argument
#' in \code{\link{glmmPen}}
#' 
#' @return numeric sequence of penalty parameters of length \code{nlambda} ranging from the
#' minimum penalty parameter (first element) equal to fraction \code{lambda.min} multiplied by the 
#' maximum penalty parameter to the maximum penalty parameter (last element)
#' 
#' @importFrom ncvreg setupLambda
#' @export
LambdaSeq = function(X, y, family, alpha = 1, lambda.min = NULL, nlambda = 10,
                       penalty.factor = NULL){
  # Checks
  if(!is.matrix(X)){
    stop("X must be a matrix")
  }
  if(!is.numeric(y)){
    stop("y must be numeric")
  }
  if(nrow(X) != length(y)){
    stop("The number of rows of X must equal the length of y")
  }
  if(!is.null(penalty.factor)){
    if(ncol(X) != length(penalty.factor)){
      stop("The number of columns of X must equal the length of penalty.factor")
    }
    if(any(!(penalty.factor %in% c(0,1)))){
      stop("penalty.factor must be vector of 0 and 1 only")
    }
  }
  
  # Borrowed elements from `ncvreg` function
  n = nrow(X)
  p = ncol(X)
  
  if(family == "gaussian"){
    yy = y - mean(y)
  }else{
    yy = y
  }
  
  if(is.null(lambda.min)){
    lambda.min = ifelse(n>p, 0.01, 0.05)
  }
  
  if(is.null(penalty.factor)){
    penalty.factor = rep(1, p)
  }
  
  # setupLambda from ncvreg package
  ## Order: from max lambda to min lambda
  lambda = setupLambda(X, yy, family, alpha, lambda.min, nlambda, penalty.factor)
  # reverse the order of the lambda - from min lambda to max lambda
  lambda_rev = lambda[order(lambda)]
  
  return(lambda_rev)
  
}





#################################################################################################

## Possible code to only perform final E-step if number of random effects in final model is low
# If final model has relatively low random effect dimensions, perform another E step
## Use results to calculate logLik and posterior samples, and save samples for MCMC diagnostics
## If final model has too many random effects, the calculation of this last E step
## with the necessarily large number of samples may be too computationally burdensome, so
## we will output NA values and allow the user to calculate these values after the model fit.

# Estep_out = list(u0 = NULL, u_init = fit$u_init, 
#                  post_modes = rep(NA, times = ncol(data_input$Z)),
#                  post_out = matrix(NA, nrow = 1, ncol = ncol(data_input$Z)),
#                  ll = NA, BICh = NA, BIC = NA, BICNgrp = NA)
# 
# q_final = sum(diag(fit$sigma) > 0)
# if(q_final <= 51){
#   
#   Estep_out = E_step_final(dat = data_input, offset_fit = offset, fit = fit, optim_options = optim_options, 
#                            fam_fun = fam_fun, extra_calc = T, 
#                            adapt_RW_options = adapt_RW_options, trace = trace)
# }