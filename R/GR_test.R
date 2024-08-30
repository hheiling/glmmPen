
#' @importFrom stringr str_to_lower str_c str_detect
#' @importFrom Matrix Matrix
#' @importFrom bigmemory write.big.matrix attach.big.matrix
#' @importFrom stats model.offset na.omit
#' @import bigmemory Rcpp 
# @export
GR_test = function(formula, data = NULL, family = "binomial", covar = NULL,
                   offset = NULL,
                   fixef_noPen = NULL, penalty = c("MCP","SCAD","lasso"),
                   alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0),
                   optim_options = optimControl(), adapt_RW_options = adaptControl(),
                   trace = 0, tuning_options = selectControl(), 
                   BICq_posterior = NULL,
                   progress = TRUE, ...){
  
  ###########################################################################################
  # Input argument checks and modifications
  ###########################################################################################
  
  # Check that (...) arguments are subsets of glmmPen_FA arguments (factor assumption version)
  args_extra = list(...)
  args_avail = c("r_estimation","survival_options")
  r_estimation = NULL
  survival_options = NULL
  if(length(args_extra) >= 1){
    if(!all(names(args_extra) %in% args_avail)){
      stop("additional arguments provided in '...' input must match the following glmmPen_FA arguments: \n",
           "'r_estimation' or 'survival_options', see glmmPen_FA and phmmPen_FA documentation for details")
    }
    if(any(names(args_extra) %in% "r_estimation")){
      r_estimation = args_extra$r_estimation
    }
    if(any(names(args_extra) %in% "survival_options")){
      survival_options = args_extra$survival_options
    }
  }
  
  if(is.null(r_estimation)){
    fit_type = "glmmPen"
  }else if(!is.null(r_estimation)){
    fit_type = "glmmPen_FA"
    if(!inherits(r_estimation, "rControl")){
      stop("r_estimation must be of class 'rControl' (see rControl() documentation)")
    }
    # For the record, fix 'covar' to 'unstructured'
    covar = "unstructured"
  }
  
  # Input modification and restriction for family
  family_info = family_export(family)
  fam_fun = family_info$family_fun
  data_type = family_info$data_type
  
  if(data_type == "survival"){
    if(is.null(survival_options)){
      stop("survival_options must be specified for the 'coxph' family, see phmmPen and phmmPen_FA documentation for details")
    }
    if(!inherits(survival_options,"survivalControl")){
      stop("survival_options must be of class 'survivalControl' (see survivalControl() documentation)")
    }
  }
  
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
  
  if(data_type == "survival"){
    if(!inherits(survival_options, "survivalControl")){ 
      stop("survival_options must be of class 'survivalControl' (see survivalControl documentation)")
    }
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
  fD_out = fD_adj(out = fD_out0, data_type = data_type)
  
  # If offset = NULL, then set offset as arbitrary vector of 0s
  if(is.null(model.offset(fD_out$frame))){
    offset = rep(0, length(fD_out$y))
  }else{
    offset = model.offset(fD_out$frame)
    if((!is.numeric(offset)) | (length(offset) != length(fD_out$y))){
      stop("offset must be a numeric vector of the same length as y")
    }
  }
  
  # If survival data, create long-form dataset needed for piecewise exponential
  #   model fit and and update fD_out object
  if(data_type == "survival"){
    # Save original data (input format, one row per observation/subject)
    fD_out_original = fD_out
    offset_original = offset
    # Calculate long-form dataset
    data_surv = survival_data(y = fD_out$y, X = fD_out$X,
                              Z = fD_out$Z, group = fD_out$group_num,
                              offset_fit = offset,
                              survival_options = survival_options)
    # Re-specify offset to include appropriate offset for piecewise exponential model
    offset = data_surv$offset_total
    # Update fD_out object with long-form data elements of X, Z, y, and group
    fD_out$X = data_surv$X
    fD_out$Z = data_surv$Z
    fD_out$y = data_surv$y_status
    fD_out$group_num = data_surv$group
    fD_out$group = fD_out$group[data_surv$IDs]
  }else{
    fD_out_original = NULL
    offset_original = NULL
  }

  # Names of covariates for fixed effects, random effects, and group factor
  coef_names = list(fixed = colnames(fD_out$X), random = fD_out$cnms, group = fD_out$group_name)
  
  # If needed, standardize input covariates
  # data_input: list of main data for fit algorithm
  standardization = optim_options$standardization
  if(standardization == TRUE){
    if(data_type != "survival"){
      # Standardize X and Z - see XZ_std() function later in this document
      std_out = XZ_std(fD_out)
      # Specify data to use in fit procedure
      data_input = list(y = fD_out$y, X = std_out$X_std, Z = std_out$Z_std, 
                         group = fD_out$group_num, coef_names = coef_names)
    }else if(data_type == "survival"){
      # Standardize original X and Z (not long-form dataset)
      std_out0 = XZ_std(fD_out_original)
      # Using standardized values of X and Z, create long-form dataset
      data_surv_std = survival_data(y = fD_out_original$y,
                                    X = std_out0$X_std,
                                    Z = std_out0$Z_std,
                                    group = fD_out_original$group_num,
                                    offset_fit = offset_original,
                                    survival_options = survival_options)
      # Specify data to use in fit procedure
      data_input = list(y = fD_out$y,
                        X = data_surv_std$X,
                        Z = data_surv_std$Z,
                        group = fD_out$group_num,
                        coef_names = coef_names)
      # Update std_out object with long-form data elements
      ## Note: No standardization is applied to the time interval indicator variables
      ##    in the long-form X dataset
      std_out = list(X_std = data_surv_std$X,
                     Z_std = data_surv_std$Z,
                     X_center = c(rep(0,length(data_surv$cut_points)-1), std_out0$X_center),
                     X_scale = c(rep(1,length(data_surv$cut_points)-1), std_out0$X_scale),
                     Z_center = std_out0$Z_center,
                     Z_scale = std_out0$Z_scale) 
    }
  }else if(standardization == FALSE){
    # Put dummy values into std_out list
    std_out = list(X_std = NULL, Z_std = NULL, 
                   X_center = rep(0,ncol(fD_out$X[,-1,drop=FALSE])), 
                   X_scale = rep(1,ncol(fD_out$X[,-1,drop=FALSE])),
                   Z_center = NULL, Z_scale = NULL)
    # Specify data to use in the fit procedure
    data_input = list(y = fD_out$y, X = fD_out$X, Z = fD_out$Z, 
                      group = fD_out$group_num, coef_names = coef_names)
  }
  
  
  # Add cut-points information to data_input list if survival data (piecewise exponential model)
  if(data_type == "survival"){
    data_input$cut_points = data_surv$cut_points
  }
  
  # Identify fixed effects that should not be penalized in select_tune or fit_dat (if any)
  ## For now, do not allow grouping of fixed effects (i.e. cannot penalize fixed effects as groups of covariates)
  # group_X: 0 indicates intercept or other covariates that should not be penalized
  #   positive integer indicates that the fixed effect can be penalized
  if(data_type != "survival"){
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
  }else if(data_type == "survival"){
    p_tmp = ncol(data_input$X)
    cut_num = length(data_input$cut_points)
    if(is.null(fixef_noPen)){
      group_X = c(rep(0,times = cut_num),1:(p_tmp-cut_num))
      fixef_noPen = c(rep(0,times = cut_num-1), rep(1,times=p_tmp-cut_num))
    }else if(is.numeric(fixef_noPen)){
      if(length(fixef_noPen) != (p_tmp-cut_num-1)){
        stop("length of fixef_noPen must match number of fixed effects covariates")
      }
      if(sum(fixef_noPen == 0) == 0){
        group_X = c(rep(0,times = cut_num),1:(p_tmp-cut_num))
      }else{
        ones = which(fixef_noPen == 1) + 1
        # Sequence: consecutive integers (1 to number of fixed effects to potentially penalize)
        sq = 1:length(ones) + cut_num
        # Initialize as all 0
        group_X = rep(0, times = ncol(data_input$X))
        # for appropriate fixed effects, set to the appropriate positive integer
        group_X[ones] = sq
      }
      # Update fixef_noPen to include 0's for the time indicator columns
      fixef_noPen = c(rep(0,times = cut_num-1), fixef_noPen)
    }
  }
  
  
  # Store call for output object
  call = match.call(expand.dots = FALSE)
  
  # If covar specified as NULL, recommend a covariance structure based on the size of the 
  # random effects.
  if(is.null(covar) & (ncol(data_input$Z)/nlevels(data_input$group) >= 11)){
    message("Setting random effect covariance structure to 'independent'")
    covar = "independent"
  }else if(is.null(covar) & (ncol(data_input$Z)/nlevels(data_input$group) < 11)){
    message("Setting random effect covariance structure to 'unstructured'")
    covar = "unstructured"
  }
  
  if((covar == "unstructured") & (ncol(data_input$Z)/nlevels(data_input$group) >= 11) & (fit_type == "glmmPen")){
    warning("The random effect covariance matrix is currently specified as 'unstructured'. 
            Due to dimension of sigma covariance matrix, we suggest using covar = 'independent' to simplify computation",
            immediate. = TRUE)
  }
  
  if(fit_type == "glmmPen_FA"){
    
    ###########################################################################################
    # Estimate number of latent factors, r
    ###########################################################################################
    
    # If input r = NULL, use input data to estimate this value
    if(is.null(r_estimation$r)){
      
      # penalty to use for fixed effects:
      if(inherits(tuning_options, "lambdaControl")){
        lambda0 = tuning_options$lambda0
      }else if(inherits(tuning_options, "selectControl")){
        lambda0_seq = tuning_options$lambda0_seq
        if(is.null(lambda0_seq)){
          lambda.min = tuning_options$lambda.min
          if(is.null(lambda.min)){
            if(ncol(data_input$X) <= 51){ # 50 predictors plus intercept
              lambda.min = 0.01
            }else if(ncol(data_input$X) <= 201){
              lambda.min = 0.05
            }else if(ncol(data_input$X) > 201){
              lambda.min = 0.10
            }
          }
          lambda0 = min(LambdaSeq(X = data_input$X[,-1,drop=FALSE], y = data_input$y, family = fam_fun$family, offset = offset,
                                  alpha = alpha, nlambda = tuning_options$nlambda, penalty.factor = fixef_noPen,
                                  lambda.min = lambda.min))
        }else{
          lambda0 = min(lambda0_seq)
        }
      }
      
      # Estimate number of common factors r
      r_out = estimate_r(dat = data_input, optim_options = optim_options,
                         coef_names = coef_names, 
                         r_est_method = r_estimation$r_est_method,
                         r_max = r_estimation$r_max, 
                         family_info = family_info, offset_fit = offset, 
                         penalty = penalty, lambda0 = lambda0, 
                         gamma_penalty = gamma_penalty, alpha = alpha, 
                         group_X = group_X, 
                         sample = ((r_estimation$size > nlevels(data_input$group)) & (r_estimation$sample)), 
                         size = r_estimation$size, data_type = data_type,
                         trace = trace)
      # Extract estimate of r
      r = r_out$r_est
      # Record maximum considered value of r
      r_max = r_out$r_max
      
      message("estimated r: ", r)
      
    }else{ # Checks on input r
      r = r_estimation$r
      if(r <= 1){
        message("glmmPen_FA restricts the number of latent factors r to be at least 2, ",
                "setting r to 2")
        r = 2
      }
      if(!(floor(r)==r)){
        stop("number of latent factors r must be an integer")
      }
      r_max = NULL
    }
    
  }
  
  return(r)
  
  
}

