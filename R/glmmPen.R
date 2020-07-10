
#' Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
#' Minimization (MCECM)
#' 
#' \code{glmmPen} is used to fit penalized generalized mixed models via Monte Carlo Expectation 
#' Conditional Minimization (MCECM)
#' 
#' @inheritParams formulaData
#' @inheritParams fit_dat_B
#' @param offset This can be used to specify an \emph{a priori} known component to be included in the 
#' linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the 
#' number of cases. Currently, the formula does not allow specification of an offset.
#' @param fixef_noPen Optional vector of 0's and 1's of the same length as the number of fixed effect covariates
#' used in the model. Value 0 indicates the variable should not have its fixed effect coefficient
#' penalized, 1 indicates that it can be penalized. Order should correspond to the same order of the 
#' fixed effects given in the formula.
#' @param optim_options a list of class "optimControl" created from function \code{\link{optimControl}}
#' that specifies optimization parameters.
#' @param adapt_RW_options a list of class "adaptControl" from function \code{\link{adaptControl}} 
#' containing the control parameters for the adaptive random walk Metropolis-within-Gibbs procedure. 
#' Ignored if \code{\link{optimControl}} parameter \code{MwG_sampler} is set to "independence"
#' @param tuning_options a list of class selectControl or lambdaControl resulting from 
#' \code{\link{selectControl}} or \code{\link{lambdaControl}} containing additional control parameters.
#' If the user wants to run the algorithm using one specific set of
#' penalty parameters \code{lambda0} and \code{lambda1}, then use \code{lambdaControl()}. 
#' If no penalization is desired, use this setting with \code{lambda0} = \code{lambda1} = 0.
#' If the user wants to run the algorithm over multiple possible \code{lambda0} and \code{lambda1},
#' then use \code{selectControl{}}. See the \code{\link{lambdaControl}} and \code{\link{selectControl}}
#' documentation for details
#' 
## #' @inheritSection fit_dat
#' 
#' @return An reference class object of class \code{\link{pglmmObj}} for which many methods are 
#' available (e.g. \code{methods(class = "pglmmObj")})
#'  
#' @importFrom stringr str_to_lower str_c
#' @importFrom Matrix Matrix
#' @export
glmmPen = function(formula, data = NULL, family = "binomial", na.action = na.omit,
                   offset = NULL, fixef_noPen = NULL, penalty = c("MCP","SCAD","lasso"),
                   alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0),
                   optim_options = optimControl(), adapt_RW_options = adaptControl(),
                   returnMC = T, trace = 0, tuning_options = selectControl()){
  # Things to address / Questions to answer:
  ## Add option for different penalties
  ## Specify what fit_dat output will be
  ## Provide option for offset, weights
  ## gibbs T / F: internally specify T or F depending on X and Z dimensions
  
  # Input modification and restriction for family
  if(is.character(family)){
    family = get(family, mode = "function", envir = parent.frame())
  }
  if(is.function(family)){
    family = family()
  }
  if(class(family) == "family"){
    if(!(family$family %in% c("binomial","poisson","gaussian"))){
      stop("family must be binomial, poisson, or gaussian")
    }
  }
  
  penalty = penalty[1]
  if(!(penalty %in% c("lasso","MCP","SCAD"))){
    stop("penalty ", penalty, " not available, must choose 'lasso', 'MCP', or 'SCAD' \n")
  }

  # Convert formula and data to useful forms to plug into fit_dat
  fD_out = formulaData(formula, data, na.action)
  if(!any(fD_out$cnms == "(Intercept)")){
    print(fD_out$cnms)
    stop("Model requires an intercept term")
  }
  
  ## Convert group to numeric factor - for fit_dat
  ## Even if already numeric, convert to levels 1,2,3,... (consecutive integers)
  group_num = as.factor(as.numeric(fD_out$group))
  # if(any(is.character(fD_out$group))){
  #   group = as.factor(as.numeric(fD_out$group))
  # }else{
  #   group = fD_out$group
  # }
  
  # Standardize X and Z
  std_out = XZ_std(fD_out, group_num)
  
  data_input = list(y = fD_out$y, X = std_out$X_std, Z = std_out$Z_std, group = group_num)
  
  coef_names = list(fixed = colnames(fD_out$X), random = fD_out$cnms, group = fD_out$group_name)
  
  # Identify fixed effects that should not be penalized in select_tune or fit_dat (if any)
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
      zeros = which(fixef_noPen == 0) + 1
      sq = 1:length(ones)
      group_X = rep(0, times = ncol(data_input$X)-1)
      group_X[ones] = sq
    }
  }
  
  # Things that should be included in call:
  ## formula, data, any other items included in glmmPen function call
  call = match.call(expand.dots = F)
  
  # rho = environment()
  if(!is.list(tuning_options) | !inherits(tuning_options, "pglmmControl")){
    stop("tuning_option parameter must be a list of type lambdaControl() or selectControl()")
  }
  
  if(inherits(tuning_options, "lambdaControl")){
    lambda0 = tuning_options$lambda0
    lambda1 = tuning_options$lambda1
    if(lambda0 < 0 | lambda1 < 0){
      stop("lambda0 and lambda1 cannot be negative")
    }
    
    # Extract variables from optimControl
    conv_EM = optim_options$conv_EM
    conv_CD = optim_options$conv_CD
    nMC_burnin = optim_options$nMC_burnin
    nMC = optim_options$nMC
    nMC_max = optim_options$nMC_max
    nMC_report = optim_options$nMC_report
    maxitEM = optim_options$maxitEM
    maxit_CD = optim_options$maxit_CD
    M = optim_options$M
    t = optim_options$t
    covar = optim_options$covar
    sampler = optim_options$sampler
    gibbs = optim_options$gibbs
    fit_type = optim_options$fit_type
    max_cores = optim_options$max_cores
    
    # Call fit_dat function
    # fit_dat function found in "/R/fit_dat.R" file
    fit = fit_dat_B(dat = data_input, lambda0 = lambda0, lambda1 = lambda1, 
                    conv_EM = conv_EM, conv_CD = conv_CD,
                    family = family, offset_fit = offset, trace = trace, 
                    group_X = group_X, penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                    nMC_burnin = nMC_burnin, nMC = nMC, nMC_max = nMC_max, nMC_report = nMC_report,
                    t = t, maxitEM = maxitEM, maxit_CD = maxit_CD,
                    M = M, gibbs = gibbs, sampler = sampler, adapt_RW_options = adapt_RW_options,
                    covar = covar, fit_type = fit_type, returnMC = returnMC, max_cores = max_cores)
    
    selection_results = matrix(NA, nrow = 1, ncol = 1)
    optim_results = matrix(NA, nrow = 1, ncol = 1)
    
  }else if(inherits(tuning_options, "selectControl")){
    if(is.null(tuning_options$lambda0_seq)){
      lambda0_range = LambdaRange(X = data_input$X[,-1], y = data_input$y, family = family$family, 
                                  alpha = alpha, nlambda = tuning_options$nlambda, penalty.factor = fixef_noPen)
    }else{
      lambda0_range = tuning_options$lambda0_seq
      if(!is.numeric(lambda0_range) | any(lambda0_range < 0)){
        stop("lambda0_seq must be a positive numeric sequence")
      }
    }
    if(is.null(tuning_options$lambda1_seq)){
      lambda1_range = lambda0_range
    }else{
      lambda1_range = tuning_options$lambda1_seq
      if(!is.numeric(lambda1_range) | any(lambda1_range < 0)){
        stop("lambda1_seq must be a positive numeric sequence")
      }
    }
    
    fit_select = select_tune(dat = data_input, offset = offset, family = family,
                             lambda0_range = lambda0_range, lambda1_range = lambda1_range,
                             penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                             returnMC = returnMC, trace = trace,
                             adapt_RW_options = adapt_RW_options, 
                             optim_options = optim_options)
    
    if(tuning_options$BIC_option == "BICh"){
      fit = fit_select$out[["BICh"]]
    }else{
      fit = fit_select$out[["BIC"]]
    }
    
    resultsA = fit_select$results
    coef_results = fit_select$coef
    # Unstandardize coefficient results
    beta_results = matrix(0, nrow = nrow(coef_results), ncol = ncol(data_input$X))
    beta_results[,1] = coef_results[,1] - apply(coef_results[,2:ncol(data_input$X)], MARGIN = 1, FUN = function(x) sum(x * std_out$X_center / std_out$X_scale))
    for(i in 1:nrow(beta_results)){
      beta_results[,-1] = coef_results[,2:ncol(data_input$X)] / std_out$X_scale
    }
    colnames(beta_results) = c(coef_names$fixed)
    
    selection_results = cbind(resultsA,beta_results,coef_results[,(ncol(data_input$X)+1):ncol(coef_results)])
    
    if(tuning_options$BIC_option == "BICh"){
      optim_results = matrix(selection_results[which.min(selection_results[,"BICh"]),], nrow = 1)
    }else{
      optim_results = matrix(selection_results[which.min(selection_results[,"BIC"]),], nrow = 1)
    }
    colnames(optim_results) = colnames(selection_results)
  }
  
  # Things that should be included in fit_dat:
  ## (fill in later)
  
  if(optim_options$gibbs){
    sampling = "Metropolis-within-Gibbs Sampling"
  }else{
    if(fit$rej_to_gibbs < 3){
      sampling = "Rejection Sampling"
    }else{
      sampling = "Metropolis-within-Gibbs Sampling"
    }
  }
  
  # Format Output - create pglmmObj object
  output = c(fit, list(call = call, formula = formula, data = data, Y = fD_out$y,
                       X = fD_out$X, Z = fD_out$Z, group = fD_out$flist,
                       coef_names = coef_names, family = family,
                       offset = offset, frame = fD_out$frame, 
                       sampling = sampling, std_out = std_out, 
                       selection_results = selection_results, optim_results = optim_results))

  out_object = pglmmObj$new(output)
  return(out_object)
  
  # return(fit)

}


#' @importFrom ncvreg std
#' @export
XZ_std = function(fD_out, group_num){
  # Standardize X - ncvreg::std method
  X = fD_out$X
  X_noInt_std = std(X[,-1])
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
    stop("Model requires the assumption of an intercept")
    v_start = 1
  }
  
  for(v in v_start:num_vars){ 
    cols = seq(from = (v - 1)*d + 1, to = v*d, by = 1)
    for(k in 1:nlevels(group_num)){
      ids = which(group_num == k)
      Z_std[ids, cols[k]] = (Z_sparse[ids, cols[k]] - Z_center[v-1]) / Z_scale[v-1]
    }
  }
  
  return(list(X_std = X_std, Z_std = Z_std, X_center = X_center, X_scale = X_scale,
              Z_center = Z_center, Z_scale = Z_scale))
}


#' @importFrom ncvreg setupLambda
#' @export
LambdaRange = function(X, y, family, alpha = 1, lambda.min = NULL, nlambda = 20,
                       penalty.factor = NULL){
  # Borrowed elements from `ncvreg` function
  n = nrow(X)
  p = ncol(X)
  
  if(family == "gaussian"){
    yy = y - mean(y)
  }else{
    yy = y
  }
  
  if(is.null(lambda.min)){
    lambda.min = ifelse(n>p, 0.001, 0.05)
  }
  
  if(is.null(penalty.factor)){
    penalty.factor = rep(1, p)
  }
  
  # lambda = calcLambda(X, yy, family, alpha, lambda.min, nlambda, penalty.factor)
  lambda = setupLambda(X, yy, family, alpha, lambda.min, nlambda, penalty.factor)
  
  return(lambda)
  
}