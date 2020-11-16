
#' Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
#' Minimization (MCECM)
#' 
#' \code{glmmPen} is used to fit penalized generalized mixed models via Monte Carlo Expectation 
#' Conditional Minimization (MCECM)
#' 
#' @inheritParams formulaData
#' @param family a description of the error distribution and link function to be used in the model. 
#' Currently, the \code{glmmPen} algorithm allows the binomial, gaussian, and poisson families
#' with canonical links only.
#' @param offset This can be used to specify an \emph{a priori} known component to be included in the 
#' linear predictor during fitting. Default set to \code{NULL} (no offset). If the data 
#' argument is \code{NULL}, this should be a numeric vector of length equal to the 
#' number of cases (the response). If the data argument specifies a data.frame, the offset
#' argument should specify the name of a column in the data.frame. 
#' Currently, the formula does not allow specification of an offset.
#' @param fixef_noPen Optional vector of 0's and 1's of the same length as the number of fixed effects covariates
#' used in the model. Value 0 indicates the variable should not have its fixed effect coefficient
#' penalized, 1 indicates that it can be penalized. Order should correspond to the same order of the 
#' fixed effects given in the formula.
#' @param penalty character describing the type of penalty to use in the variable selection procedure.
#' Options include 'MCP', 'SCAD', and 'lasso'. Default is MCP penalty.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions 
#' from the MCP/SCAD/lasso penalty and the ridge, or L2, penalty. \code{alpha=1} is equivalent to 
#' the MCP/SCAD/lasso penalty, while \code{alpha=0} is equivalent to ridge regression. However,
#' \code{alpha=0} is not supported; \code{alpha} may be arbibrarily small, but not exactly zero
#' @param gamma_penalty The tuning parameter of the MCP and SCAD penalties. Not used by Lasso penalty.
#' Default is 4.0 for SCAD and 3.0 for MCP.
#' @param optim_options a list of class "optimControl" created from function \code{\link{optimControl}}
#' that specifies optimization parameters.
#' @param adapt_RW_options a list of class "adaptControl" from function \code{\link{adaptControl}} 
#' containing the control parameters for the adaptive random walk Metropolis-within-Gibbs procedure. 
#' Ignored if \code{\link{optimControl}} parameter \code{sampler} is set to "stan" or "independence".
#' @param tuning_options a list of class selectControl or lambdaControl resulting from 
#' \code{\link{selectControl}} or \code{\link{lambdaControl}} containing additional control parameters.
#' If the user wants to run the algorithm using one specific set of
#' penalty parameters \code{lambda0} and \code{lambda1}, then use \code{lambdaControl()}. 
#' If no penalization is desired, use this setting with \code{lambda0} = \code{lambda1} = 0.
#' If the user wants to run the algorithm over multiple possible \code{lambda0} and \code{lambda1},
#' then use \code{selectControl{}}. See the \code{\link{lambdaControl}} and \code{\link{selectControl}}
#' documentation for details
#' @param BICq_posterior an optional character string expressing the path and file name of a text file that 
#' the will store or currently stores a \code{big.matrix} of the posterior draws from the full model.
#' These full model posterior draws will be used in BIC-ICQ calculations if these calculations
#' are requested. The name of the file should include a .txt extension. If this argument is
#' specified as \code{NULL} (default) and BIC-ICQ calculations are requested, the posterior draws
#' will be saved in the file 'BICq_Posterior_Draws.txt' in the working directory.
#' See 'Details' section for additional details.
#' @param trace an integer specifying print output to include as function runs. Default value is 0. 
#' See Details for more information about output provided when trace = 0, 1, or 2.
#' 
#' @details 
#' If the \code{BIC_option} in \code{\link{selectControl}} (\code{tuning_options}) is specified 
#' to be 'BICq', this requests the calculation of the BIC-ICQ criterion during the selection
#' process. For the BIC-ICQ criterion to be calculated, a full model assuming a small valued 
#' lambda penalty needs to be fit, and the posterior draws from this full model need to be used. 
#' In order to avoid repetitive calculations of
#' this full model if secondary rounds of selection are desired, a \code{big.matrix} of these 
#' posterior draws will be saved in a .txt file. If the file already exits, the full model 
#' will not be fit again and the big.matrix of 
#' posterior draws will be read using \code{bigmemory::read.big.matrix}
#' (\code{bigmemory::read.big.matrix(filename, sep = ' ', type = 'double')}) and used in the BIC-ICQ 
#' calcuations. If the file does not exist or \code{BICq_posterior} 
#' is specified as \code{NULL}, the full model will be fit and the full model posterior
#' draws will be saved at the specified file location (using \code{bigmemory::write.big.matrix})
#' or in the working directory with generic filename "BICq_Posterior_Draws.txt" (if \code{NULL} 
#' is specified). The algorithm will save 10^4 posterior draws automatically.
#' 
#' Trace details: The value of 0 outputs some general updates for each EM iteration (iteration number EM_iter,
#' number of MCMC draws nMC, average Euclidean distance between current coefficients and coefficients
#' from t iterations back EM_diff, and number of non-zero coefficients Non0 Coef). The value of 1
#' additionally outputs the updated coefficients, updated covariance matrix values, and the
#' number of coordinate descent iterations used for the M step for each
#' EM iteration. The value of 2 outputs all of the above plus gibbs acceptance rate information
#' for the adaptive random walk and independence samplers and the updated proposal standard deviation
#' for the adaptive random walk. 
#' 
#' @return A reference class object of class \code{\link{pglmmObj}} for which many methods are 
#' available (e.g. \code{methods(class = "pglmmObj")})
#'  
#' 
#' @importFrom stringr str_to_lower str_c
#' @importFrom Matrix Matrix
#' @export
glmmPen = function(formula, data = NULL, family = "binomial", na.action = na.omit,
                   offset = NULL, fixef_noPen = NULL, penalty = c("MCP","SCAD","lasso"),
                   alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0),
                   optim_options = optimControl(), adapt_RW_options = adaptControl(),
                   trace = 0, tuning_options = selectControl(), BICq_posterior = NULL){
  # Things to address / Questions to answer:
  ## Specify what fit_dat output will be
  ## Provide option for offset, weights
  
  # Input modification and restriction for family
  family_info = family_export(family)
  fam_fun = family_info$family_fun
  # family = family_info$family
  
  if(length(penalty) > 1){
    penalty = penalty[1]
  }
  if(!(penalty %in% c("lasso","MCP","SCAD"))){
    stop("penalty ", penalty, " not available, must choose 'lasso', 'MCP', or 'SCAD' \n")
  }
  
  if(!is.numeric(gamma_penalty)){
    stop("gamma_penalty must be numeric")
  }
  if(penalty == "MCP" & gamma_penalty <= 1){
    stop("When using MCP penalty, gamma_penalty must be > 1")
  }else if(penalty == "SCAD" & gamma_penalty <= 2){
    stop("When using SCAD penalty, gamma_penalty must be > 2")
  }
  
  if(class(optim_options) != "optimControl"){
    stop("optim_options must be of class 'optimControl', see optimControl documentation")
  }
  if(class(adapt_RW_options) != "adaptControl"){
    stop("adapt_RW_options must be of class 'adaptControl', see adaptControl documentation")
  }
  
  
  # Convert formula and data to useful forms to plug into fit_dat
  fD_out = formulaData(formula, data, na.action)
  if(!any(fD_out$cnms == "(Intercept)")){
    print(fD_out$cnms)
    stop("Model requires an intercept term")
  }
  
  if(is.null(offset)){
    offset = rep(0, length(fD_out$y))
  }else if(!is.null(data) & !is.null(offset)){
    offset_name = offset
    offset = data$offset_name
    if(is.null(offset)){
      stop("offset variable not found in data")
    }
  }else if(is.null(data) & !is.null(offset)){
    # Check offset has appropriate attributes
    if((!is.numeric(offset)) | (length(offset) != length(fD_out$y))){
      stop("offset must be a numeric vector of the same length as y")
    }
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
    max_cores = optim_options$max_cores
    
    # Call fit_dat function
    # fit_dat_B function found in "/R/fit_dat_MstepB.R" file
    
    # Default arguments:
    # dat, lambda0 = 0, lambda1 = 0, conv_EM = 0.001, conv_CD = 0.0001,
    # family = "binomial", offset_fit = NULL,
    # trace = 0, penalty = c("MCP","SCAD","lasso"),
    # alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0), 
    # group_X = 0:(ncol(dat$X)-1),
    # nMC_burnin = 500, nMC = 5000, nMC_max = 20000, t = 2,
    # nMC_report = 5000, u_init = NULL, coef_old = NULL, 
    # ufull_describe = NULL, maxitEM = 100, maxit_CD = 250,
    # M = 10^4, sampler = c("stan","random_walk","independence"),
    # adapt_RW_options = adaptControl(), covar = c("unstructured","independent"),
    # var_start = 0.5, max_cores = 1
    
    fit = fit_dat_B(dat = data_input, lambda0 = lambda0, lambda1 = lambda1, 
                    conv_EM = conv_EM, conv_CD = conv_CD,
                    family = fam_fun, offset_fit = offset, trace = trace, 
                    group_X = group_X, penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                    nMC_burnin = nMC_burnin, nMC = nMC, nMC_max = nMC_max, nMC_report = nMC_report,
                    t = t, maxitEM = maxitEM, maxit_CD = maxit_CD,
                    M = M, sampler = sampler, adapt_RW_options = adapt_RW_options,
                    covar = covar, max_cores = max_cores)
    
    selection_results = matrix(NA, nrow = 1, ncol = 1)
    # c("lambda0","lambda1","BICh","BIC","BICq","LogLik","Non_0_fef","Non_0_ref","Non_0_coef")
    optim_results = c(fit$lambda0, fit$lambda1, fit$BICh, fit$BIC, fit$BICq, fit$ll,
                      sum(fit$coef[2:ncol(data_input$X)] != 0),
                      sum(diag(fit$sigma) != 0),
                      sum(fit$coef != 0))
    optim_results = matrix(optim_results, nrow = 1)
    colnames(optim_results) = c("lambda0","lambda1","BICh","BIC","BICq","LogLik","Non_0_fef","Non_0_ref","Non_0_coef")
    
  }else if(inherits(tuning_options, "selectControl")){
    if(is.null(tuning_options$lambda0_seq)){
      lambda0_range = LambdaRange(X = data_input$X[,-1], y = data_input$y, family = fam_fun$family, 
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
    
    # dat, offset = NULL, family, group_X = 0:(ncol(dat$X)-1),
    # penalty, lambda0_range, lambda1_range,
    # alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0),
    # trace = 0, u_init = NULL, coef_old = NULL, 
    # BICq_calc = T,
    # adapt_RW_options = adaptControl(),
    # optim_options = optimControl(), BICq_posterior
    
    fit_select = select_tune(dat = data_input, offset = offset, family = family,
                             lambda0_range = lambda0_range, lambda1_range = lambda1_range,
                             penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                             group_X = group_X, trace = trace,
                             BICq_calc = (tuning_options$BIC_option == "BICq"),
                             adapt_RW_options = adapt_RW_options, 
                             optim_options = optim_options, BICq_posterior = BICq_posterior)
    
    if(tuning_options$BIC_option == "BICh"){
      fit = fit_select$out[["BICh"]]
    }else if(tuning_options$BIC_option == "BIC"){
      fit = fit_select$out[["BIC"]]
    }else if(tuning_options$BIC_option == "BICq"){
      fit = fit_select$out[["BICq"]]
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
    
    gamma_results = coef_results[,(ncol(data_input$X)+1):ncol(coef_results)]
    colnames(gamma_results) = str_c("Gamma",0:(ncol(gamma_results)-1))
    
    selection_results = cbind(resultsA,beta_results,gamma_results)
    
    if(tuning_options$BIC_option == "BICh"){
      optim_results = matrix(selection_results[which.min(selection_results[,"BICh"]),], nrow = 1)
    }else if(tuning_options$BIC_option == "BIC"){
      optim_results = matrix(selection_results[which.min(selection_results[,"BIC"]),], nrow = 1)
    }else if(tuning_options$BIC_option == "BICq"){
      optim_results = matrix(selection_results[which.min(selection_results[,"BIC"]),], nrow = 1)
    }
    colnames(optim_results) = colnames(selection_results)
  }
  
  # Things that should be included in fit_dat:
  ## (fill in later)
  
  sampling = switch(optim_options$sampler, stan = "Stan", 
                    random_walk = "Metropolis-within-Gibbs Adaptive Random Walk Sampler",
                    independence = "Metropolis-within-Gibbs Independence Sampler")
  
  if(penalty == "lasso"){
    gamma_penalty = NULL
  }
  
  # Format Output - create pglmmObj object
  output = c(fit, list(call = call, formula = formula, y = fD_out$y,
                       X = fD_out$X, Z_std = std_out$Z_std, group = fD_out$flist,
                       coef_names = coef_names, family = fam_fun,
                       offset = offset, frame = fD_out$frame, 
                       sampling = sampling, std_out = std_out, 
                       selection_results = selection_results, optim_results = optim_results,
                       penalty = penalty, gamma_penalty = gamma_penalty, alpha = alpha, fixef_noPen = fixef_noPen,
                       control_options = list(optim_options = optim_options, tuning_options = tuning_options,
                                              adapt_RW_options = adapt_RW_options)))
  
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
  
  colnames(X_std) = colnames(X)
  colnames(Z_std) = colnames(fD_out$Z)
  
  return(list(X_std = X_std, Z_std = Z_std, X_center = X_center, X_scale = X_scale,
              Z_center = Z_center, Z_scale = Z_scale))
}


#' @importFrom ncvreg setupLambda
#' @export
LambdaRange = function(X, y, family, alpha = 1, lambda.min = NULL, nlambda = 10,
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