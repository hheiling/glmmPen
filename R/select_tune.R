
#' Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
#' Minimization (MCECM) 
#' 
#' \code{select_tune} is used to fit penalized generalized mixed models via Monte Carlo Expectation 
#' Conditional Minimization (MCECM) over multiple tuning parameters and is called within
#' \code{glmmPen}
#' 
#' @inheritParams glmmPen
#' @inheritParams fit_dat_B
#' @inheritParams selectControl
#' @inheritParams optimControl
#' @param BICq_calc boolean value indicating if the BIC-ICQ criterion should be used to select the
#' best model.
#' @param checks_complete boolean value indicating if several data checks have been completed.
#' @param lambda.min.full a vector of two numeric values that gives the fixed and random effect penalty 
#' values to use in pre-screening and/or the full model fit for the BIC-ICQ calculation 
#' (if applicable)
#' @param stage1 boolean value indicating if the first stage of the abbreviated two-stage grid
#' search in the model selection procedure is being performed.
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
                       penalty, lambda0_seq, lambda1_seq,
                       alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0),
                       trace = 0, u_init = NULL, coef_old = NULL, 
                       adapt_RW_options = adaptControl(),optim_options = optimControl(), 
                       BIC_option = c("BICq","BICh","BIC","BICNgrp"), BICq_calc = T, 
                       logLik_calc = switch(BIC_option[1], BICq = F, T),
                       BICq_posterior = NULL, checks_complete = F, pre_screen = T, 
                       ranef_keep = NULL, lambda.min.full, stage1 = F){
  
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
  
  if(is.null(offset)){
    offset = rep(0, length(dat$y))
  }
  
  if((covar == "unstructured") & (ncol(dat$Z)/nlevels(dat$group) >= 11)){
    warning("Due to dimension of sigma covariance matrix, will use covar = 'independent' to simplify computation",
            immediate. = T)
    covar = "independent"
  }

  if(optim_options$var_start == "recommend"){
    var_start = var_init(dat, fam_fun)
    cat("recommended starting variance: ", var_start, "\n")
    # record starting variance 
    optim_options$var_start = var_start
  }else{
    var_start = optim_options$var_start
  }
  
  sampler = optim_options$sampler
  

  # Pre-screening step
  ## If number of random effects in model (including intercept) is 6 or more
  ## Note: since random effects are subset of fixed effects, fixed effects necessarily 6 or more
  if((ncol(dat$Z)/nlevels(dat$group) >= 6) & (pre_screen == T)){
    cat("begin prescreening procedure \n")
    out_pre = prescreen(dat = dat, family = fam_fun$family, offset_fit = offset, trace = trace, 
                        penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty, 
                        lambda0_min = lambda.min.full[1], lambda1_min = lambda.min.full[2], 
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
  
  # If nMC arguments or maxitEM left as NULL, fill in defaults
  optim_options = optim_recommend(optim_options = optim_options, family = family, 
                                  q = sum(ranef_keep), select = T)
  
  # Extract variables from optimControl
  conv_EM = optim_options$conv_EM
  conv_CD = optim_options$conv_CD
  nMC_burnin = optim_options$nMC_burnin
  nMC = optim_options$nMC
  nMC_max = optim_options$nMC_max
  # nMC_report = optim_options$nMC_report
  maxitEM = optim_options$maxitEM
  maxit_CD = optim_options$maxit_CD
  M = optim_options$M
  t = optim_options$t
  mcc = optim_options$mcc
  sampler = optim_options$sampler
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
      # fixed effects penalty
      if(ncol(dat$X) >= 6) lam0 = lambda.min.full[1] else lam0 = 0
      # random effect penalty - may be larger than fixed effects penalty
      ## There is a general sparsity assumption for random effects, so we may
      ## assume a higher penalty on the random effects for a 'full' model
      if(ncol(dat$Z)/nlevels(dat$group) >= 6) lam1 = lambda.min.full[2] else lam1 = 0
      print("Penalty parameters used to calculate full model for BIC-ICQ calculation:")
      print(sprintf("fixed effects %f, random effects %f", lam0, lam1))
      # Fit 'full' model
      ## Note: M restricted to be >= 10^4 in optimControl(). Will report M posterior draws from full model
      out = try(fit_dat_B(dat, lambda0 = lam0, lambda1 = lam1, 
                          nMC_burnin = nMC_burnin, nMC = nMC, nMC_max = nMC_max,
                          family = fam_fun, offset_fit = offset, group_X = group_X,
                          penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                          trace = trace, conv_EM = conv_EM, conv_CD = conv_CD,  
                          coef_old = coef_pre, u_init = u_pre, ufull_describe = NULL,
                          maxitEM = maxitEM + 15, maxit_CD = maxit_CD, t = t, mcc = mcc,
                          M = M, sampler = sampler, adapt_RW_options = adapt_RW_options,
                          covar = covar, var_start = var_start, logLik_calc = F,
                          max_cores = max_cores, checks_complete = checks_complete,
                          ranef_keep = ranef_keep))
      
      if(is.character(out)){
        print("Issue with full model fit, no BICq calculation will be completed")
        ufull_describe = NULL
      }else{
        
        Estep_out = E_step_final(dat = dat, fit = out, optim_options = optim_options, 
                                 fam_fun = fam_fun, extra_calc = F, 
                                 adapt_RW_options = adapt_RW_options, trace = trace)
        
        ufull_big = attach.big.matrix(Estep_out$u0)
        write.big.matrix(ufull_big, filename = BICq_posterior, sep = " ")
        ufull_describe = Estep_out$u0
        # If pre_screen = T, restrict the next model in the search to only consider random 
        # effects that were not penalized out in this full model
        if(pre_screen){
          ranef_keep = as.numeric(diag(out$sigma) > 0)
        }
        # Restrict lambda1_seq to not include any points smaller than the minimum lambda used
        # in the full model fit (restriction also applied if pre-screening performed)
        lambda1_seq = lambda1_seq[which(lambda1_seq >= lam1)]
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
  
  # Restrict lambda1_seq to not include any points smaller than the minimum lambda used
  # in the 'full' model fit used in pre-screening or BIC-ICQ calculation
  ## Note: this restriction is also done when the full model is fit for the BIC-ICQ calculation
  ## (see code above)
  if(pre_screen){
    lambda1_seq = lambda1_seq[which(lambda1_seq >= lambda.min.full[2])]
  }
  
  n1 = length(lambda0_seq)
  n2 = length(lambda1_seq)
  # Keep track of BIC-type selection criteria
  ## BIC_critold: record of largest BIC-type critera so far during selection
  BIC_crit = Inf
  BIC_critold = Inf
  
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
  
  # If stage 1 of abbreviated grid search, will adjust ranef_keep for subsequent models:
  ## if random effect penalized out in model j, that random effect will be assumed 0 in
  ## subseqent model
  # Therefore, store ranef_keep after pre-screening/full model (if performed)
  ranef_keep_presc = ranef_keep
  
  for(j in 1:n2){ # random effects
    for(i in 1:n1){ # fixed effects
      if((i == 1) & (j == 1)){
        # start from the initial values input (or NULL values)
        coef_old0 = coef_old
        uold = u_init # Should be null
      }else if((i == 1) & (j != 1)){ 
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
      print(sprintf("lambda0 %f lambda1 %f", lambda0_seq[i], lambda1_seq[j]))
      print(sprintf("lambda0 i %i lambda1 j %i", i, j))
      print("------------------------------------------------------------------")
      out = try(fit_dat_B(dat, lambda0 = lambda0_seq[i], lambda1 = lambda1_seq[j], 
                          nMC_burnin = nMC_burnin, nMC = nMC, nMC_max = nMC_max,
                          family = fam_fun, offset_fit = offset, group_X = group_X,
                          penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                          trace = trace, conv_EM = conv_EM, conv_CD = conv_CD,  
                          coef_old = coef_old0, u_init = uold, ufull_describe = ufull_describe,
                          maxitEM = maxitEM, maxit_CD = maxit_CD, t = t, mcc = mcc,
                          M = M, sampler = sampler, adapt_RW_options = adapt_RW_options,
                          covar = covar, var_start = var_start, logLik_calc = logLik_calc,
                          max_cores = max_cores, checks_complete = checks_complete,
                          ranef_keep = ranef_keep))
      

      if(is.character(out)) next
      
      # Abbreviated grid search result: in stage 1, will restrict random effects in subsequent
      ## model to only consider random effects not penalized to 0 in this past model
      if(stage1){
        ranef_keep = as.numeric((diag(out$sigma) > 0))
      }
      
      # for restarting the next j for i = 1
      ## keep coefficients and posterior draws for initialization of next i = 1, j combination
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
      
      res[(j-1)*n1+i,1] = lambda0_seq[i]
      res[(j-1)*n1+i,2] = lambda1_seq[j]
      res[(j-1)*n1+i,3] = BICh
      res[(j-1)*n1+i,4] = BIC
      res[(j-1)*n1+i,5] = BICq
      res[(j-1)*n1+i,6] = BICNgrp
      res[(j-1)*n1+i,7] = out$ll
      res[(j-1)*n1+i,8] = sum(out$coef[2:ncol(dat$X)] !=0) # ignore intercept
      if(length(diag(out$sigma)) > 1){
        res[(j-1)*n1+i,9] = sum(diag(out$sigma[-1,-1,drop=F]) !=0) # ignore random intercept
      }else{
        res[(j-1)*n1+i,9] = 0
      }
      res[(j-1)*n1+i,10] = sum(out$coef != 0)
      res[(j-1)*n1+i,11] = out$EM_iter
      coef = rbind(coef, out$coef)
      print(res[(j-1)*n1+i,])
      
      
    } # End i loop for lambda0_seq
    
  } # end j loop for lambda1_seq
  
  colnames(coef) = c("(Intercept)",str_c("B",1:(ncol(dat$X)-1)), str_c("Gamma",1:(ncol(coef)-ncol(dat$X))))
  
  return(list(results = res, out = fout, coef = coef, optim_options = optim_options, 
              ranef_keep = ranef_keep_presc))
  
}

