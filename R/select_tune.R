
#' Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
#' Minimization (MCECM) 
#' 
#' \code{select_tune} is used to fit penalized generalized mixed models via Monte Carlo Expectation 
#' Conditional Minimization (MCECM) over multiple tuning parameters
#' 
#' @inheritParams glmmPen
#' @inheritParams fit_dat_B
#' @inheritParams selectControl
#' @inheritParams optimControl
#' @param BICq_calc boolean value indicating if the BIC-ICQ criterion should be used to select the
#' best model.
#' 
#' @return A list with the following elements:
#' \item{results}{matrix of summary results for each lambda tuning parameter combination, used
#' to select the 'best' model}
#' \item{out}{list of \code{\link{fit_dat_B}} results for the models with the minimum BICh and 
#' minimum BIC values}
#' \item{coef}{matrix of coefficient results for each lambda tuning parameter combination. 
#' Rows correspond with the rows of the results matrix.}
#'  
#' @importFrom stringr str_c
#' @importFrom bigmemory attach.big.matrix describe write.big.matrix read.big.matrix
select_tune = function(dat, offset = NULL, family, covar = c("unstructured","independent"), 
                       group_X = 0:(ncol(dat$X)-1),
                       penalty, lambda0_range, lambda1_range,
                       alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0),
                       trace = 0, u_init = NULL, coef_old = NULL, 
                       adapt_RW_options = adaptControl(),optim_options = optimControl(), 
                       BIC_option = c("BICh","BIC","BICq","BICNgrp"), BICq_calc = T, 
                       BICq_posterior = NULL, checks_complete = F, pre_screen = T, ranef_keep = NULL){
  
  # Input modification and restriction for family
  family_info = family_export(family)
  fam_fun = family_info$family_fun
  link = family_info$link
  link_int = family_info$link_int # Recoded link as integer
  family = family_info$family
  
  if(length(BIC_option) > 1){
    BIC_option = BIC_option[1]
  }
  if(!(BIC_option %in% c("BICh","BIC","BICq","BICNgrp"))){
    stop("BIC_option must be 'BICh', 'BIC', 'BICq', or 'BICNgrp'")
  }
  
  # # Calculate loglik for saturated model (mu set to y; n parameters, one per observation)
  # if((family == "binomial") & (link == "logit")){
  #   sat_ll = 0
  # }else if((family == "poisson") & (link == "log")){
  #   sat_ll = 0
  #   y = dat$y
  #   for(i in 1:length(y)){
  #     if(y[i] > 0){
  #       sat_ll = sat_ll + dpois(x = y[i], lambda = y[i], log = T)
  #     }
  #   }
  # }else if((family == "gaussian") & (link == "identity")){
  #   # Will not look at model saturation for gaussian family
  #   sat_ll = 0
  # }else{
  #   print(fam_fun)
  #   stop("specified family and link combination not available \n")
  # }
  # # Set null deviance arbitrarily to 1 for the moment
  # nullDev = 1.0
  
  if(is.null(offset)){
    offset = rep(0, length(dat$y))
  }
  
  if((covar == "unstructured") & (ncol(dat$Z)/nlevels(dat$group) >= 11)){
    warning("Due to dimension of sigma covariance matrix, will use covar = 'independent' to simplify computation",
            immediate. = T)
    covar = "independent"
  }
  
  if(is.character(optim_options)){ # "recommend"
    var_start = var_init(data = dat, fam_fun = fam_fun)
    sampler = "stan"
  }else if(inherits(optim_options, "optimControl")){
    var_start = optim_options$var_start
    if(var_start == "recommend"){
      var_start = var_init(data = dat, fam_fun = fam_fun)
    }
    sampler = optim_options$sampler
  }
  
  # Pre-screening step
  ## If number of random effects in model (including intercept) is 11 or more
  ## Note: since random effects are subset of fixed effects, fixed effects necessarily 11 or more
  if((ncol(dat$Z)/nlevels(dat$group) >= 11) & (pre_screen == T)){
    cat("begin prescreening procedure \n")
    out_pre = prescreen(dat = dat, family = fam_fun$family, offset_fit = offset, trace = trace, 
                        penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty, 
                        group_X = group_X, sampler = sampler, 
                        adapt_RW_options = adapt_RW_options, covar = covar,
                        var_start = var_start, max_cores = max_cores, 
                        checks_complete = checks_complete)
    cat("end prescreening procedure \n")
    ranef_keep = out_pre$ranef_keep
    coef_pre = out_pre$coef_pre
    u_pre = out_pre$u_pre
  }else{
    if(is.null(ranef_keep)){
      ranef_keep = rep(1, times = ncol(dat$Z)/nlevels(dat$group))
    }
    coef_pre = NULL 
    u_pre = NULL
  }
  
  if(is.character(optim_options)){
    if(optim_options == "recommend"){
      optim_options = optim_recommend(family, sum(ranef_keep), select = T)
    }
  }
  # Save var_start
  optim_options$var_start = var_start
  
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
  mcc = optim_options$mcc
  sampler = optim_options$sampler
  # var_start = optim_options$var_start
  max_cores = optim_options$max_cores
  
  # If need to calculate full model for BICq calculation
  ## Use minimum lambda in range of lambda
  ufull_describe = NULL
  
  if(BICq_calc == T){
    # If posterior draws from full model not yet created, fit full model and save posterior draws
    # Otherwise, read in the prevoiusly fit full model posterior draw results
    if(!is.null(BICq_posterior)){
      if(file.exists(BICq_posterior)){
        cat("Using saved posterior draws from full model for BIC-ICQ calculation: ",
            BICq_posterior,"\n")
        fitfull_needed = F
      }else{
        cat("Fitting full model and saving posterior draws to ",BICq_posterior, "\n")
        # warning("BICq_posterior file ",BICq_posterior," not found", 
        #         "  Creating new full model posterior draws", immediate. = T)
        fitfull_needed = T
      }
    }else{ # if is.null(BICq_posterior)
      BICq_posterior = "BICq_Posterior_Draws.txt"
      cat("Fitting full model and saving posterior draws to ",BICq_posterior, "\n",
          "  within working directory \n")
      fitfull_needed = T
    }
    
    if(fitfull_needed){
      # For high dimensions (number of fixed or random effects >= 10 not including intercept), 
      # find a small penalty to use for full model: the minimum of the lambda range used by ncvreg
      lam_MaxMin = LambdaRange(dat$X[,-1], dat$y, family = family, nlambda = 2, lambda.min = 0.001)
      lam_min = lam_MaxMin[1]
      if(ncol(dat$X) >= 11) lam0 = lam_min else lam0 = 0
      if(ncol(dat$Z)/nlevels(dat$group) >= 11) lam1 = lam_min else lam1 = 0
      # Fit 'full' model
      ## Note: M restricted to be >= 10^4 in optimControl(). Will report M posterior draws from full model
      out = try(fit_dat_B(dat, lambda0 = lam0, lambda1 = lam1, 
                          nMC_burnin = nMC_burnin, nMC = nMC, nMC_max = nMC_max, nMC_report = nMC_report,
                          family = fam_fun, offset_fit = offset, group_X = group_X,
                          penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                          trace = trace, conv_EM = conv_EM, conv_CD = conv_CD,  
                          coef_old = coef_pre, u_init = u_pre, ufull_describe = NULL,
                          maxitEM = maxitEM + 50, maxit_CD = maxit_CD, t = t, mcc = mcc,
                          M = M, sampler = sampler, adapt_RW_options = adapt_RW_options,
                          covar = covar, var_start = var_start,
                          max_cores = max_cores, checks_complete = checks_complete,
                          ranef_keep = ranef_keep))
      
      if(is.character(out)){
        print("Issue with full model fit, no BICq calculation will be completed")
        ufull_describe = NULL
      }else{
        ufull_big = attach.big.matrix(out$u_big)
        write.big.matrix(ufull_big, filename = BICq_posterior, sep = " ")
        ufull_describe = out$u_big
      }
      
    }else{ # if fitfull_needed = F
      cat("Reading in ",BICq_posterior," for posterior draws for BICq calculation \n")
      ufull_big = read.big.matrix(filename = BICq_posterior, sep = " ", type = 'double')
      ufull_describe = describe(ufull_big)
      # Checks
      if(ncol(ufull_big) != ncol(dat$Z)){
        stop("The number of columns in the saved full model posterior draws do not match \n",
             "  the number of columns in the random effects matrix Z: ",ncol(dat$Z))
      }
      if(nrow(ufull_big) < 10^4){
        warning("The number of posterior draws saved in ",BICq_posterior, "\n",
                "  is less than the recommended 10^4",immediate. = T)
      }
    } # End if-else fitfull_needed
    
  } # End if-else BICq_calc = T
  
  n1 = length(lambda0_range)
  n2 = length(lambda1_range)
  # Keep track of BIC-type selection criteria
  ## BIC_critold: record of largest BIC-type critera so far during selection
  BIC_crit = BIC_critold = Inf
  
  res = matrix(0, nrow = n1*n2, ncol = 11)
  colnames(res) = c("lambda0","lambda1","BICh","BIC","BICq","BICNgrp","LogLik","Non_0_fef","Non_0_ref","Non_0_coef","EM_iter")
  
  # If available, initialize first model with 'full' model from BIC-ICQ calculation
  if(BICq_calc){
    if(fitfull_needed){
      if(!is.character(out)){
        coef_old = out$coef
        coef_oldj0 = coef_old
        uj0 = out$u_init
      }
    }
  }else{
    coef_oldj0 = NULL
    uj0 = NULL
  }
  
  coef = NULL
  
  # saturated = FALSE
  fout = list()
  
  for(j in 1:n2){ # random effects
    for(i in 1:n1){ # fixed effects
      if((i == 1) & (j == 1)){
        # start from the initial values input (or NULL values)
        coef_old0 = coef_old
        uold = u_init # Should be null
      }else if((i == 1) & (j != 1)){ # changed from if(i != 1 & j == 1)
        # if re-starting fixed effects penalty loop but moving on to next random effects
        # penalty parameter, take the previous j result for i = 1
        coef_old0 = coef_oldj0
        uold = uj0
        rm(out)
      }else{
        # take the previous values
        coef_old0 = out$coef
        uold = out$u_init
        rm(out)
      }
      
      # If had divergence issues or NA results in past (lambda0,lambda1) combination,
      # then initialized next combination with NULL results (initialize from scratch
      # instead of from last model)
      if(!is.null(coef_old0)){
        if(any(is.na(coef_old0)) | any(abs(coef_old0) > 10^5)){
          coef_old0 = NULL
          uold = NULL
        }
      }
      
      gc()
      print("------------------------------------------------------------------")
      print(sprintf("lambda0 %f lambda1 %f", lambda0_range[i], lambda1_range[j]))
      print(sprintf("lambda0 i %i lambda1 j %i", i, j))
      print("------------------------------------------------------------------")
      out = try(fit_dat_B(dat, lambda0 = lambda0_range[i], lambda1 = lambda1_range[j], 
                          nMC_burnin = nMC_burnin, nMC = nMC, nMC_max = nMC_max, nMC_report = nMC_report,
                          family = fam_fun, offset_fit = offset, group_X = group_X,
                          penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                          trace = trace, conv_EM = conv_EM, conv_CD = conv_CD,  
                          coef_old = coef_old0, u_init = uold, ufull_describe = ufull_describe,
                          maxitEM = maxitEM, maxit_CD = maxit_CD, t = t, mcc = mcc,
                          M = M, sampler = sampler, adapt_RW_options = adapt_RW_options,
                          covar = covar, var_start = var_start,
                          max_cores = max_cores, checks_complete = checks_complete,
                          ranef_keep = ranef_keep))
      

      if(is.character(out)) next
      
      # for restarting the next j for i = 1
      if(i == 1){
        coef_oldj0 = out$coef
        uj0 = out$u_init
      }
      
      if(BIC_option == "BICh"){
        BIC_crit = out$BICh
      }else if(BIC_option == "BIC"){
        BIC_crit = out$BIC
      }else if(BIC_option == "BICq"){
        BIC_crit = out$BICq
        if(is.na(BIC_crit)){
          warning("BIC-ICQ not calculated, using BICh for selection criteria instead",immediate. = T)
        }
      }else if(BIC_option == "BICNgrp"){
        BIC_crit = out$BICNgrp
      }
      
      if(!is.numeric(BIC_crit)) BIC_crit = Inf
      if(length(BIC_crit) != 1) BIC_crit = Inf
      
      if(BIC_crit < BIC_critold){
        fout = out
        BIC_critold = BIC_crit
      }
      
      BICh = out$BICh
      BIC = out$BIC
      BICq = out$BICq
      BICNgrp = out$BICNgrp
      BIC_report = c(BICh,BIC,BICq,BICNgrp)
      names(BIC_report) = c("BICh","BIC","BICq","BICNgrp")
      print(BIC_report)
      
      res[(j-1)*n1+i,1] = lambda0_range[i]
      res[(j-1)*n1+i,2] = lambda1_range[j]
      res[(j-1)*n1+i,3] = BICh
      res[(j-1)*n1+i,4] = BIC
      res[(j-1)*n1+i,5] = BICq
      res[(j-1)*n1+i,6] = BICNgrp
      res[(j-1)*n1+i,7] = out$ll
      res[(j-1)*n1+i,8] = sum(out$coef[2:ncol(dat$X)] !=0)
      if(length(diag(out$sigma)) > 1){
        res[(j-1)*n1+i,9] = sum(diag(out$sigma[-1]) !=0)
      }else{
        res[(j-1)*n1+i,9] = 0
      }
      res[(j-1)*n1+i,10] = sum(out$coef !=0)
      res[(j-1)*n1+i,11] = out$EM_iter
      coef = rbind(coef, out$coef)
      print(res[(j-1)*n1+i,])
      
      # # Check for model saturation in binomial and poisson case
      # if(family %in% c("binomial","poisson")){
      #   
      #   ## set null deviance as deviance for model with fixed and random intercepts only
      #   ## This assumes lambda.max is the first element of the lambda0 and lambda1 sequences
      #   ## This assumption not valid if users specify their own lambda0/1 sequences
      #   ## Idea: fit random intercept model using lme4 and use the given log-lik for the null model
      #   if((i==1) & (j==1)){
      #     nullDev = 2*(sat_ll - out$ll)
      #   }
      #   
      #   ## Calculate deviance
      #   Dev = 2*(sat_ll - out$ll)
      #   if(is.finite(Dev)){
      #     if(Dev / nullDev < 0.01){
      #       print("Reached model saturation")
      #       saturated = TRUE
      #       break
      #     }
      #   }
      #   
      # } End if-else for model saturation check
      
      
    } # End i loop for lambda0_range
    
    # if(saturated){
    #   break
    # }
    
  } # end j loop for lambda1_range
  
  colnames(coef) = c("(Intercept)",str_c("B",1:(ncol(dat$X)-1)), str_c("Gamma",1:(ncol(coef)-ncol(dat$X))))
  
  return(list(results = res, out = fout, coef = coef, optim_options = optim_options, 
              ranef_keep = ranef_keep))
  
}

