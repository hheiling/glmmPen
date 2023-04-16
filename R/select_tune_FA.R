
# @title Fit a Sequence of Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
# Minimization (MCECM) 
# 
# \code{select_tune_FA} is used to fit a sequence of penalized generalized mixed models 
# via Monte Carlo Expectation Conditional Minimization (MCECM) for 
# multiple tuning parameter combinations and is called within
# \code{glmmPen_FA} (cannot be called directly by user)
# 
# @inheritParams glmmPen
# @inheritParams glmmPen_FA
# @inheritParams fit_dat
# @inheritParams fit_dat_FA
# @inheritParams select_tune
# @inheritParams selectControl
# @param beta_init numeric vector used to initialize fixed effects
# @param B_init numeric matrix used to initialize factor loadings matrix for the random effects. 
# Dimensions: total number of columns equal the total number of random effects (including intercept),
# total number of rows equal the total number of assumed latent common factors (\code{r}).
# 
# @return A list with the following elements:
# \item{results}{matrix of summary results for each lambda tuning parameter combination, used
# to select the 'best' model}
# \item{out}{list of \code{\link{fit_dat}} results for the best model}
# \item{coef}{matrix of coefficient results for each lambda tuning parameter combination. 
# Rows correspond with the rows of the results matrix.}
# 
#' @importFrom stringr str_c
#' @importFrom bigmemory attach.big.matrix describe write.big.matrix read.big.matrix
select_tune_FA = function(dat, offset = NULL, family, 
                       group_X = 0:(ncol(dat$X)-1),
                       penalty, lambda0_seq, lambda1_seq,
                       alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0),
                       trace = 0, u_init = NULL, beta_init = NULL, B_init = NULL,
                       adapt_RW_options = adaptControl(),optim_options = optimControl(), 
                       BIC_option = c("BICq","BICh","BIC","BICNgrp"), BICq_calc = TRUE, 
                       logLik_calc = switch(BIC_option[1], BICq = FALSE, TRUE),
                       BICq_posterior = NULL, checks_complete = FALSE, pre_screen = TRUE, 
                       ranef_keep = NULL, lambda.min.full, stage1 = FALSE, 
                       r, progress = TRUE){
  
  #############################################################################################
  # Input Checks
  #############################################################################################
  
  # Input modification and restriction for family
  family_info = family_export(family)
  fam_fun = family_info$family_fun
  link = family_info$link
  link_int = family_info$link_int # Recoded link as integer, see "family_export.R" for details
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
  
  
  #############################################################################################
  # Pre-screening (optional)
  #############################################################################################
  
  # Pre-screening step
  ## If number of random effects in model (including intercept) is 6 or more
  ## Note: since random effects are subset of fixed effects, fixed effects necessarily 6 or more
  if((ncol(dat$Z)/nlevels(dat$group) >= 6) & (pre_screen == TRUE)){
    message("Running prescreening procedure")
    out_pre = prescreen_FA(dat = dat, family = fam_fun$family, offset_fit = offset, trace = trace, 
                        penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty, 
                        lambda0_min = lambda.min.full[1], lambda1_min = lambda.min.full[2], 
                        group_X = group_X, optim_options = optim_options, 
                        adapt_RW_options = adapt_RW_options,
                        checks_complete = checks_complete, r = r, progress = progress)
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
  
  #############################################################################################
  # Optimization parameters
  #############################################################################################
  
  # If nMC arguments or maxitEM left as NULL, fill in defaults
  optim_options = optim_recommend(optim_options = optim_options, family = family, 
                                  q = r, select = TRUE) # q = sum(ranef_keep)
  
  # If need to calculate minimal penalty model for BICq calculation
  ## Use minimum lambda in range of lambda
  ufull_describe = NULL
  
  #############################################################################################
  # Minimal penalty model fit for BIC-ICQ calculation (optional)
  #############################################################################################
  
  if(BICq_calc == T){
    # If posterior draws from minimal penalty model not yet created, fit minimal penalty model and save posterior draws
    # Otherwise, read in the prevoiusly fit minimal penalty model posterior draw results
    if(!is.null(BICq_posterior)){
      if(file.exists(sprintf("%s.desc",BICq_posterior)) & file.exists(sprintf("%s.bin",BICq_posterior))){
        message("Using saved posterior draws from minimal penalty model for BIC-ICQ calculation: ")
        message(sprintf("file-backed big.matrix stored in %s.bin and %s.desc", BICq_posterior, BICq_posterior))
        fitfull_needed = F
      }else{
        message(sprintf("The files %s.bin and %s.desc do not currently exist.", BICq_posterior, BICq_posterior))
        message(sprintf("Fitting minimal penalty model and saving posterior draws to %s.bin and %s.desc",
                        BICq_posterior, BICq_posterior))
        fitfull_needed = T
      }
    }else{ # if is.null(BICq_posterior)
      BICq_posterior = "BICq_Posterior_Draws"
      message(sprintf("Fitting minimal penalty model and saving posterior draws to %s.bin and %s.desc within working directory",
                      BICq_posterior, BICq_posterior))
      fitfull_needed = T
    }
    
    if(fitfull_needed){
      # fixed effects penalty
      if(ncol(dat$X) >= 6) lam0 = lambda.min.full[1] else lam0 = 0
      # random effect penalty - may be larger than fixed effects penalty
      ## There is a general sparsity assumption for random effects, so we may
      ## assume a higher penalty on the random effects for a 'full' model
      if(ncol(dat$Z)/nlevels(dat$group) >= 6) lam1 = lambda.min.full[2] else lam1 = 0
      if(trace >= 1){
        cat("Penalty parameters used to calculate minimal penalty model for BIC-ICQ calculation: \n")
        cat(sprintf("fixed effects %f, random effects %f", lam0, lam1), "\n")
      }
      # Fit 'full' model
      ## Note: M default set to 10^4 in optimControl(). Will report M posterior draws from minimal penalty model
      if(!is.null(coef_pre)){
        beta_old = coef_pre[1:ncol(dat$X)]
        b_old = coef_pre[-c(1:ncol(dat$X))]
      }else{
        beta_old = beta_init
        if(!is.null(B_init)){
          b_old = c(t(B_init))
        }else{
          b_old = NULL
        }
      }
      out = try(fit_dat_FA(dat, lambda0 = lam0, lambda1 = lam1, 
                          family = fam_fun, offset_fit = offset, 
                          optim_options = optim_options, group_X = group_X,
                          penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                          trace = trace,  
                          beta_old = beta_old, b_old = b_old,
                          u_init = u_pre, ufull_describe = NULL,
                          adapt_RW_options = adapt_RW_options,
                          r=r, logLik_calc = FALSE, checks_complete = checks_complete,
                          ranef_keep = ranef_keep, progress = progress))
      
      if(is.character(out)){
        warning("Issue with minimal penalty model fit, no BICq calculation will be completed \n",
                immediate. = TRUE)
        message("Setting model selection criteria to BICh")
        BIC_option = "BICh"
        ufull_describe = NULL
      }else{
        # Using minimal penalty model fit parameters, sample from the posterior (posterior samples used for
        # BIC-ICQ calculation)
        Estep_out = E_step_final_FA(dat = dat, fit = out, offset_fit = offset, 
                                    optim_options = optim_options, 
                                    fam_fun = fam_fun, extra_calc = FALSE, 
                                    adapt_RW_options = adapt_RW_options, trace = trace, 
                                    progress = progress, r = r)
        
        ufull_big = attach.big.matrix(Estep_out$u0)
        ufull_big_tmp = as.big.matrix(ufull_big[,],
                                      backingpath = dirname(BICq_posterior),
                                      backingfile = sprintf("%s.bin",basename(BICq_posterior)),
                                      descriptorfile = sprintf("%s.desc",basename(BICq_posterior)))
        rm(ufull_big_tmp)
        # write.big.matrix(ufull_big, filename = BICq_posterior, sep = " ")
        ufull_describe = Estep_out$u0
        # If pre_screen = T, restrict the next model in the search to only consider random 
        # effects that were not penalized out in this minimal penalty model
        if(pre_screen){
          ranef_keep = as.numeric(diag(out$sigma) > 0)
        }
        # Restrict lambda1_seq to not include any points smaller than the minimum lambda used
        # in the minimal penalty model fit (restriction also applied if pre-screening performed)
        lambda1_seq = lambda1_seq[which(lambda1_seq >= lam1)]
      }
      
    }else{ # if fitfull_needed = F
      if(progress == TRUE){
        cat("Reading in ",BICq_posterior," for posterior draws for BICq calculation \n")
      }
      ufull_big = attach.big.matrix(sprintf("%s.desc",BICq_posterior))
      # ufull_big = read.big.matrix(filename = BICq_posterior, sep = " ", type = 'double')
      ufull_describe = describe(ufull_big)
      # Checks
      if(ncol(ufull_big) != nlevels(dat$group)*r){
        stop("The number of columns in the saved minimal penalty model posterior samples does not equal \n",
             "  the number of groups ", nlevels(dat$group)," times number of common factors ", r, 
             ", ", nlevels(dat$group)*r)
      }
      if(nrow(ufull_big) < 10^4){
        warning("The number of posterior draws saved in ",BICq_posterior, "\n",
                "  is less than the recommended 10^4",immediate. = T)
      }
    } # End if-else fitfull_needed
    
  } # End if-else BICq_calc = T
  
  #############################################################################################
  # Additional selection set-up
  #############################################################################################
  
  # Restrict lambda1_seq to not include any points smaller than the minimum lambda used
  # in the 'full' model fit used in pre-screening or BIC-ICQ calculation
  ## Note: this restriction is also done when the minimal penalty model is fit for the BIC-ICQ calculation
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
  
  res = matrix(0, nrow = n1*n2, ncol = 12)
  colnames(res) = c("lambda0","lambda1","BICh","BIC","BICq","BICNgrp","LogLik",
                    "Non_0_fef","Non_0_ref","Non_0_coef","EM_iter","Converged")
  
  # General initialization - assuming no prescreening or minimal penalty model fit
  beta_old = beta_init
  if(!is.null(B_init)){
    b_old = c(t(B_init))
  }else{
    b_old = NULL
  }
  
  coef_oldj0 = NULL
  uj0 = NULL
  
  # If available, initialize first model with 'full' model from BIC-ICQ calculation or
  # results from pre-screening
  if(BICq_calc){
    if(fitfull_needed){
      if(!is.character(out)){
        coef_old = out$coef
        beta_old = coef_old[1:ncol(dat$X)]
        b_old = coef_old[-c(1:ncol(dat$X))]
        u_init = out$u_init
        coef_oldj0 = coef_old
        uj0 = out$u_init
      }
    }
  }else if(!is.null(coef_pre)){
    coef_old = coef_pre
    beta_old = coef_pre[1:ncol(dat$X)]
    b_old = coef_pre[-c(1:ncol(dat$X))]
    u_init = u_pre
    coef_oldj0 = coef_pre
    uj0 = u_pre
  }
  
  # During selection, will store coefficients. Initialize as NULL.
  coef = NULL
  
  # saturated = FALSE
  fout = list()
  
  # If stage 1 of abbreviated grid search, will adjust ranef_keep for subsequent models:
  ## if random effect penalized out in model j, that random effect will be assumed 0 in
  ## subseqent model
  # Therefore, store ranef_keep after pre-screening/minimal penalty model (if performed)
  ranef_keep_presc = ranef_keep
  
  #############################################################################################
  # Model Selection
  #############################################################################################
  
  # beta_old and beta_old0: fixed effects
  # b_old and b_old0: coefficients for the B matrix, associated with random effects
  # uold: posterior sample (from past model fit) to initialize E-step draws in next model
  for(j in 1:n2){ # random effects
    for(i in 1:n1){ # fixed effects
      if((i == 1) & (j == 1)){
        # start from the initial values input (or NULL values)
        beta_old0 = beta_old
        b_old0 = b_old # output = row 1, row 2, ..., row q
        uold = u_init # Should be null
      }else if((i == 1) & (j != 1)){ 
        # if re-starting fixed effects penalty loop but moving on to next random effects
        # penalty parameter, take the previous j result for i = 1
        beta_old0 = coef_oldj0[1:ncol(dat$X)]
        b_old0 = coef_oldj0[-c(1:ncol(dat$X))]
        uold = uj0
        rm(out)
      }else{
        # take the previous model values (from model (i-1,j))
        coef_old0 = out$coef
        beta_old0 = coef_old0[1:ncol(dat$X)]
        b_old0 = coef_old0[-c(1:ncol(dat$X))]
        uold = out$u_init
        rm(out)
      }
      
      # If had divergence issues or NA results in past (lambda0,lambda1) combination,
      # then initialized next combination with NULL results (initialize from scratch
      # instead of from last model)
      if(!is.null(beta_old)){
        if(any(is.na(beta_old)) | any(abs(beta_old) > 10^5)){
          beta_old0 = NULL
          b_old0 = NULL
          uold = NULL
        }
      }
      
      gc()
      cat("------------------------------------------------------------------ \n")
      cat(sprintf("lambda0 %f lambda1 %f", lambda0_seq[i], lambda1_seq[j]), "\n")
      cat(sprintf("lambda0 i %i lambda1 j %i", i, j), "\n")
      cat("------------------------------------------------------------------ \n")
      out = try(fit_dat_FA(dat, lambda0 = lambda0_seq[i], lambda1 = lambda1_seq[j], 
                          family = fam_fun, offset_fit = offset, 
                          optim_options = optim_options, group_X = group_X,
                          penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                          trace = trace,
                          u_init = uold, beta_old = beta_old0, b_old = b_old0,
                          ufull_describe = ufull_describe,
                          adapt_RW_options = adapt_RW_options,
                          r=r, logLik_calc = logLik_calc,
                          checks_complete = checks_complete,
                          ranef_keep = ranef_keep, progress = progress))
      

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
      
      if(is.na(BIC_crit)){
        stop("Model selection criteria ", BIC_option, " calculated as NA due to issues with model fit,
             stopping variable selection proceedure")
      }
      
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
      res[(j-1)*n1+i,12] = out$converged
      if(progress == TRUE) print(round(res[(j-1)*n1+i,], digits = 5))
      coef = rbind(coef, out$coef)
      
      
    } # End i loop for lambda0_seq
    
  } # end j loop for lambda1_seq
  
  colnames(coef) = c("(Intercept)",str_c("B",1:(ncol(dat$X)-1)), str_c("Gamma",1:(ncol(coef)-ncol(dat$X))))
  
  return(list(results = res, out = fout, coef = coef, optim_options = optim_options, 
              ranef_keep = ranef_keep_presc))
  
}

