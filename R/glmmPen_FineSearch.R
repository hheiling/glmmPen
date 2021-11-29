

#' @title Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
#' Minimization (MCECM) using a finer penalty grid search
#' 
#' \code{glmmPen_FineSearch} finds the best model from the selection results of a \code{pglmmObj} object 
#' created by \code{glmmPen}, identifies a more targeted grid search around the optimum lambda penalty
#' values, and performs model selection on this finer grid search.
#' 
#' @param object an object of class \code{pglmmObj} created by \code{glmmPen}. This object must 
#' contain model selection results.
#' @param tuning_options a list of class selectControl resulting from \code{\link{selectControl}} 
#' containing model selection control parameters. See the \code{\link{selectControl}}
#' documentation for details. The user can specify their own fine grid search, or if the 
#' lambda0_seq and lambda1_seq arguments are left as \code{NULL}, the algorithm will automatically
#' select a fine grid search based on the best model from the previous selection. See Details for 
#' more information. Default value set to 1.
#' @param idx_range a positive integer that determines what positions within the sequence of the 
#' fixed and random effect lambda penalty parameters used in the previous coarse grid search
#' will be used as the new fixed and random effect lambda penalty parameter ranges. See Details 
#' for more information.
#' @param optim_options an optional list of class "optimControl" created from function \code{\link{optimControl}}
#' that specifies optimization parameters. If set to the default \code{NULL}, will use the 
#' optimization parameters used for the previous round of selection stored within the 
#' \code{pglmmObj} object.
#' @param adapt_RW_options an optional list of class "adaptControl" from function \code{\link{adaptControl}} 
#' containing the control parameters for the adaptive random walk Metropolis-within-Gibbs procedure. 
#' Ignored if \code{\link{optimControl}} parameter \code{sampler} is set to "stan" or "independence".
#' If set to the default \code{NULL}, will use the adaptive random walk paraters used for the 
#' previous round of selection stored within the \code{pglmmObj} object.
#' @param trace an integer specifying print output to include as function runs. Default value is 0. 
#' See Details for more information about output provided when trace = 0, 1, or 2.
#' @param BICq_posterior an optional character string specifying the file-backed \code{big.matrix} 
#' containing the posterior draws used to calculate the BIC-ICQ selection criterion if such a 
#' \code{big.matrix} was created in the previous round of selection. See \code{\link{glmmPen}} 
#' documentation for further details.
#'   
#' @details 
#' The \code{glmmPen_FineSearch} function extracts the data, the penalty information (penalty type,
#' gamma_penalty, and alpha), the pre-screening results from the initial variable selection
#' procedure, and some other argument specifications from the \code{pglmmObj} object
#' created during a previous round of variable/model selection. In this finer grid search, the user has
#' the ability to make the following adjustments: the user can change the BIC option used for selection,
#' any optimization control parameters, or any adaptive random walk parameters (if the sampler
#' specified in the optimization parameters is "random_walk"). The user could manually specify the 
#' lambda penalty grid to search over within the \code{\link{selectControl}} control parameters,
#' or the user could let the \code{glmmPen_FineSearch} algorithm calculate a finer grid search 
#' automatically (see next paragraph for details).
#' 
#' If the sequences of lambda penalty values are left unspecified in the \code{\link{selectControl}} 
#' tuning options, the \code{glmmPen_FineSearch} algorithm performs the following steps to find
#' the finer lambda grid search: (i) The lambda combination from the best model is identified 
#' from the earlier selection results saved in the \code{pglmmObj} object. (ii) For the fixed and
#' random effects separately, the new max and min lambda values are the lambda values \code{idx_range} 
#' positions away from the best lambda in the original lambda sequences for the fixed and random
#' effects. For instance, suppose we consider a hypothetical lambda sequence of 
#' \{0.1,0.2,0.3,0.4,0.5,0.6,0.7\} for both fixed and random effects, 
#' and the best model was given by the (0.4,0.5) combination. If the \code{idx_lambda} = 2, then
#' the fine search would use the fixed effects sequence would have (min,max) = (0.2,0.6) and
#' the fixed effects sequence would have (min,max) = (0.3,0.7).
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
#' @importFrom stringr str_c
#' @export
glmmPen_FineSearch = function(object, tuning_options = selectControl(), idx_range = 2,
                              optim_options = NULL, adapt_RW_options = NULL, 
                              trace = 0, BICq_posterior = NULL){
  
  ###########################################################################################
  # Input checks and/or extraction of control settings from pglmmObj object
  ###########################################################################################
  
  if(!inherits(object, "pglmmObj")){
    stop("object must be of class 'pglmmObj', output from the glmmPen function")
  }
  
  if(!inherits(tuning_options, "selectControl")){
    stop("tuning_options must be of class 'selectControl', see selectControl documentation")
  }
  
  if(is.null(optim_options)){
    cat("Using optimization parameters stored in pglmmObj object from past selection \n")
    optim_options = object$control_info$optim_options
  }else{
    if(is.character(optim_options)){ # "recommend"
      if(optim_options != "recommend"){
        stop("optim_options must be of class 'optimControl' from optimControl() or the character string 'recommend'")
      }
    }else if(!inherits(optim_options, "optimControl")){
      stop("optim_options must be of class 'optimControl' from optimControl() or the character string 'recommend'")
    }
  }
  
  if(is.null(adapt_RW_options)){
    if(optim_options$sampler == "random_walk"){
      cat("Using adaptive random walk parameters stored in pglmmObj object from past selection \n")
    }
    adapt_RW_options = object$control_info$adapt_RW_options
  }else{
    if(!inherits(adapt_RW_options, "adaptControl")){
      stop("adapt_RW_options must be of class 'adaptControl', see adaptControl documentation")
    }
  }
  
  ###########################################################################################
  # Data and argument extraction from pglmmObj object, checks of input
  ###########################################################################################
  
  y = object$data$y 
  X = object$data$X 
  Z_std = object$data$Z_std
  group = object$data$group[[1]] # Assume only grouping by a single variable
  offset = object$data$offset
  family = object$family
  fixef_names = names(object$fixef)
  
  # Standardize the X matrix
  std_info = object$data$std_info
  X_std = matrix(0, nrow = nrow(X), ncol = ncol(X))
  X_std[,1] = X[,1]
  X_std[,-1] = (X[,-1] - std_info$X_center) / std_info$X_scale
  
  data_input = list(y = y, X = X_std, Z = Z_std, group = as.factor(group))
  
  penalty = object$penalty_info$penalty
  cat("Using penalty from pglmmObj object: ",penalty,"\n")
  if(!(penalty %in% c("MCP","SCAD","lasso"))){
    stop("penalty must be 'MCP','SCAD', or 'lasso'")
  }
  
  gamma_penalty = object$penalty_info$gamma_penalty
  if(!is.numeric(gamma_penalty)){
    stop("gamma_penalty must be numeric")
  }
  if(penalty == "MCP" & gamma_penalty <= 1){
    stop("When using MCP penalty, gamma_penalty must be > 1")
  }else if(penalty == "SCAD" & gamma_penalty <= 2){
    stop("When using SCAD penalty, gamma_penalty must be > 2")
  }
  
  alpha = object$penalty_info$alpha
  if(!is.double(alpha)) {
    tmp <- try(alpha <- as.double(alpha), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("alpha must be numeric or able to be coerced to numeric", call.=FALSE)
  }else if(alpha == 0.0){
    stop("alpha cannot equal 0. Pick a small value > 0 instead (e.g. 0.001) \n");
  }
  
  if((floor(idx_range) != idx_range) | (idx_range <= 0)){
    stop("idx_range must be a positive integer")
  }
  
  select_coarse = object$results_all
  if(!inherits(object$control_info$tuning_options, "selectControl") | (ncol(object$results_all)==1)){
    stop("Input pglmmObj object must output selection results. \n",
         "  Given input only gives results for a single lambda penalty combination")
  }
  
  ###########################################################################################
  # Specify lambda penalty parameters for fine grid search
  ###########################################################################################
  
  # Determine 'best' method by best BIC_option criteria
  # BIC_option from selectControl() argument, not pglmmObj
  crit = tuning_options$BIC_option
  if(length(crit) > 1){
    crit = crit[1]
  }
  
  if((crit == "BICq") & (any(is.na(select_coarse[,crit])))){
    stop("BIC-ICQ not calculated in previous coarse grid search \n",
         "  Selection option BIC_option = 'BICq' not appropriate")
  }else if((crit != "BICq") & (any(is.na(select_coarse[,"LogLik"])))){
    stop("BIC, BICh, and BICNgrps not calculated for all models in previous coarse grid search \n",
         "  These selection options not appropriate")
  }
  
  # 'best' model result
  opt_res = select_coarse[which.min(select_coarse[,crit]),]
  
  # Specify fine lambda grid search - function
  lam_choice = function(select_coarse, opt_res, tuning_options, type = c(0,1)){
    type = type[1]
    if(type == 0){
      coarse = unique(select_coarse[,"lambda0"])
      lam_best = opt_res["lambda0"]
    }else if(type == 1){
      coarse = unique(select_coarse[,"lambda1"])
      lam_best = opt_res["lambda1"]
    }
    idx_best = which(coarse == lam_best)
    lam_seq = exp(seq(from = log(coarse[max(idx_best-idx_range,1)]), 
                      to = log(coarse[min(idx_best+idx_range,length(coarse))]),
                      length.out = tuning_options$nlambda))
    
    return(lam_seq)
    
  }
  
  # If user inputs their own penalty parameter sequence, extract and use
  # If not, determine penalty parameters for fine search using function above
  if(is.null(tuning_options$lambda0_seq)){
    lambda0_seq = lam_choice(select_coarse, opt_res, tuning_options, type = 0)
  }else{
    lambda0_seq = tuning_options$lambda0_seq
  }
  
  if(is.null(tuning_options$lambda1_seq)){
    lambda1_seq = lam_choice(select_coarse, opt_res, tuning_options, type = 1)
  }else{
    lambda1_seq = tuning_options$lambda1_seq
  }
  
  BIC_option = tuning_options$BIC_option
  covar = object$call$covar
  logLik_calc = tuning_options$logLik_calc
  
  # Specify starting coefficient - unstandardized result from final model in coarse search
  coef_old = object$Estep_init$coef
  
  fixef_noPen = object$penalty_info$fixef_noPen
  
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
  
  ranef_keep = object$penalty_info$prescreen_ranef
  names(ranef_keep) = NULL
  
  # Note: no pre-screening should be performed in this fine grid search
  ## pre_screen = F always
  # Use lambda.min.presc for full model penalty parameter specification (if needed)
  lambda.min.presc = tuning_options$lambda.min.presc
  if(is.null(lambda.min.presc)){
    if((ncol(data_input$Z)/nlevels(data_input$group) <= 11)){ # number random effects (including intercept)
      lambda.min.presc = 0.01
    }else{
      lambda.min.presc = 0.05
    }
  }
  # Minimum lambda values to use for pre-screening and BIC-ICQ full model fit
  # Fixed effects: use minimum fixed effects penalty parameters
  # Random effects: calculate penalty to use using lambda.min.presc 
  ##    (penalty = lambda.min.presc * max random effect penalty parameter)
  full_ranef = LambdaSeq(X = data_input$X[,-1,drop=F], y = data_input$y, family = family$family, 
                         alpha = alpha, nlambda = 2, penalty.factor = fixef_noPen,
                         lambda.min = lambda.min.presc)
  lambda.min.full = c(lambda0_seq[1], full_ranef[1])
  names(lambda.min.full) = c("fixef","ranef")
  
  if(tuning_options$search == "abbrev"){
    
    ###########################################################################################
    # Two-stage (abbreviated) grid search
    ###########################################################################################
    
    print("Start of stage 1 of abbreviated grid search")
    fit_fixfull = select_tune(dat = data_input, offset = offset, family = family,
                              covar = covar, group_X = group_X,
                              lambda0_seq = min(lambda0_seq), lambda1_seq = lambda1_seq,
                              penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                              trace = trace, coef_old = coef_old, 
                              adapt_RW_options = adapt_RW_options, 
                              optim_options = optim_options, logLik_calc = logLik_calc,
                              BICq_calc = (BIC_option == "BICq"),
                              BIC_option = BIC_option, BICq_posterior = BICq_posterior, 
                              checks_complete = T, pre_screen = F, ranef_keep = ranef_keep,
                              lambda.min.full = lambda.min.full, stage1 = T)
    print("End of stage 1 of abbreviated grid search")
    
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
    vars = diag(fit_fixfull$out$sigma)
    ## BIC-ICQ full model results to avoid repeat calculation of full model
    if((is.null(BICq_posterior)) & (BIC_option == "BICq")){
      BICq_post_file = "BICq_Posterior_Draws"
    }else{
      BICq_post_file = BICq_posterior
    }
    
    # Fit second stage of 'abbreviated' model fit
    print("Start of stage 2 of abbreviated grid search")
    fit_select = select_tune(dat = data_input, offset = offset, family = family,
                             covar = object$call$covar, group_X = group_X, 
                             lambda0_seq = lambda0_seq, lambda1_seq = lam_ref,
                             penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                             trace = trace,
                             coef_old = coef_old, u_init = u_init,
                             adapt_RW_options = adapt_RW_options, 
                             optim_options = optim_options,
                             logLik_calc = logLik_calc, BICq_calc = (BIC_option == "BICq"),
                             BIC_option = BIC_option, BICq_posterior = BICq_post_file, 
                             checks_complete = T, pre_screen = F,
                             ranef_keep = as.numeric((vars > 0)), 
                             lambda.min.full = lambda.min.full, stage1 = F)
    print("End of stage 2 of abbreviated grid search")
    
    resultsA = rbind(fit_fixfull$results, fit_select$results)
    coef_results = rbind(fit_fixfull$coef, fit_select$coef)
    
  }else if(tuning_options$search == "full_grid"){
    
    ###########################################################################################
    # Full grid search
    ###########################################################################################
    
    fit_select = select_tune(dat = data_input, offset = offset, family = family, 
                             covar = object$call$covar, group_X = group_X, 
                             lambda0_seq = lambda0_seq, lambda1_seq = lambda1_seq,
                             penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                             trace = trace, 
                             coef_old = coef_old, u_init = NULL, 
                             adapt_RW_options = adapt_RW_options, 
                             optim_options = optim_options, 
                             logLik_calc = logLik_calc, BICq_calc = (BIC_option == "BICq"),
                             BIC_option = BIC_option, BICq_posterior = BICq_posterior,
                             checks_complete = T, pre_screen = F, ranef_keep = ranef_keep,
                             lambda.min.full = lambda.min.full, stage1 = F)
    
    resultsA = fit_select$results
    coef_results = fit_select$coef
    
  } # End if-else tuning_options$search == 'abbrev'
  
  ###########################################################################################
  # Extract output from best model, final output preparation
  ###########################################################################################
  
  fit = fit_select$out
  
  # Unstandardize coefficient results
  beta_results = matrix(0, nrow = nrow(coef_results), ncol = ncol(data_input$X))
  beta_results[,1] = coef_results[,1] - apply(coef_results[,2:ncol(data_input$X),drop=F], MARGIN = 1, FUN = function(x) sum(x * std_info$X_center / std_info$X_scale))
  for(i in 1:nrow(beta_results)){
    beta_results[,-1] = coef_results[,2:ncol(data_input$X),drop=F] / std_info$X_scale
  }
  
  gamma_results = coef_results[,(ncol(data_input$X)+1):ncol(coef_results), drop=F]
  
  # Combine all relevant selection results
  selection_results = cbind(resultsA,beta_results,gamma_results)
  colnames(selection_results) = c(colnames(resultsA), fixef_names, 
                                  str_c("Gamma",0:(ncol(gamma_results)-1)))
  
  optim_results = matrix(selection_results[which.min(selection_results[,BIC_option]),], nrow = 1)
  colnames(optim_results) = colnames(selection_results)
  
  # Perform a final E-step
  ## Use results to calculate logLik and posterior draws, and save draws for MCMC diagnostics
  Estep_out = E_step_final(dat = data_input, fit = fit, offset_fit = offset, optim_options = optim_options, 
                           fam_fun = object$family, extra_calc = T, 
                           adapt_RW_options = adapt_RW_options, trace = trace)
  
  # update optim_results:
  optim_results[,c("BICh","BIC","BICNgrp","LogLik")] = c(Estep_out$BICh, Estep_out$BIC, 
                                                         Estep_out$BICNgrp, Estep_out$ll)
  
  optim_options = fit_select$optim_options
  
  sampling = switch(optim_options$sampler, stan = "Stan", 
                    random_walk = "Metropolis-within-Gibbs Adaptive Random Walk Sampler",
                    independence = "Metropolis-within-Gibbs Independence Sampler")
  
  call = match.call(expand.dots = F)
  
  output = c(fit, 
             list(Estep_out = Estep_out, call = call, formula = object$formula, 
                  y = y, fixed_vars = object$fixed_vars,
                   X = X, Z_std = Z_std, group = group,
                   coef_names = list(fixed = names(object$fixef), random = colnames(object$sigma)), 
                   family = object$family, offset = offset, frame = object$data$frame, 
                   sampling = sampling, std_out = object$data$std_info, 
                   selection_results = selection_results, optim_results = optim_results,
                   penalty = penalty, gamma_penalty = gamma_penalty, alpha = alpha, 
                   fixef_noPen = fixef_noPen, ranef_keep = ranef_keep,
                   control_options = list(optim_options = optim_options, tuning_options = tuning_options)))
  
  out_object = pglmmObj$new(output)
  return(out_object)
  
}