
#' Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
#' Minimization (MCECM) 
#' 
#' \code{select_tune} is used to fit penalized generalized mixed models via Monte Carlo Expectation 
#' Conditional Minimization (MCECM) over multiple tuning parameters
#' 
#' @inheritParams fit_dat_B
#' @inheritParams glmmPen
#' 
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
#' @importFrom bigmemory as.big.matrix describe
#' @export
select_tune = function(dat, offset = NULL, family, group_X = 0:(ncol(dat$X)-1),
                       penalty, lambda0_range, lambda1_range,
                       alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0),
                       trace = 0, u_init = NULL, coef_old = NULL, full_model = T,
                       adapt_RW_options = adaptControl(),
                       optim_options = optimControl()){
  
  # Input modification and restriction for family
  if(is.character(family)){
    family = get(family, mode = "function", envir = parent.frame())
  }
  if(is.function(family)){
    family = family()
  }
  if(class(family) == "family"){
    fam = family$family
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
  var_start = optim_options$var_start
  max_cores = optim_options$max_cores
  
  covar = covar[1]
  if(!(covar %in% c("unstructured","independent"))){
    stop("algorithm currently only handles 'unstructured' or 'independent' covariance structure \n")
  }
  
  # Calculate loglik for saturated model (mu set to y; n parameters, one per observation)
  if(fam == "binomial"){
    sat_ll = 0
  }else if(fam == "poisson"){
    stop("poisson family not yet operational \n")
    sat_ll = 0
    y = dat$y
    for(i in 1:length(y)){
      if(y[i] > 0){
        sat_ll = sat_ll + dpois(x = y[i], lambda = y[i], log = T)
      }
    }
  }else if(fam == "gaussian"){
    stop("gaussian family not yet operational \n")
  }else{
    print(family)
    stop("specified family not recognized \n")
  }
  # Set null deviance arbitrarily to 1 for the moment
  nullDev = 1.0
  
  if(is.null(offset)){
    offset = rep(0, length(dat$y))
  }
  
  # If need to calculate full model for BICq calculation
  ## Use minimum lambda in range of lambda
  ufull_describe = NULL
  if(full_model == T){
    # Find a small penalty to use for full model: the minimum of the lambda range used by ncvreg
    lam_MaxMin = LambdaRange(dat$X[,-1], dat$y, family = fam, nlambda = 2)
    lam_min = lam_MaxMin[2]
    # Fit 'full' model
    out = try(fit_dat_B(dat, lambda0 = lam_min, lambda1 = lam_min, 
                        nMC_burnin = nMC_burnin, nMC = nMC, nMC_max = nMC_max, nMC_report = 10^4,
                        family = family, offset_fit = offset, group_X = group_X,
                        penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                        trace = trace, conv_EM = conv_EM, conv_CD = conv_CD,  
                        coef_old = NULL, u_init = NULL, ufull_describe = NULL,
                        maxitEM = maxitEM, maxit_CD = maxit_CD, t = t,
                        M = M, sampler = sampler, adapt_RW_options = adapt_RW_options,
                        covar = covar, var_start = var_start,
                        max_cores = max_cores))
    
    
    if(is.character(out)){
      print("Issue with full model fit, no BICq calculation will be completed")
      ufull_describe = NULL
    }else{
      ufull = out$u
      if(any(is.na(ufull))){
        print("Issue with full model fit, no BICq calculation will be completed")
        ufull_describe = NULL 
      }else{
        ufull_big = bigmemory::as.big.matrix(ufull)
        ufull_describe = bigmemory::describe(ufull_big)
      }
      ufull = NULL
    }
    
    
  }
  
  n1 = length(lambda0_range)
  n2 = length(lambda1_range)
  BICold = BIC = Inf
  BIChold = BICh = Inf
  BICqold = BICq = Inf
  
  res = matrix(0, n1*n2, 9)
  coef = NULL
  coef_oldj0 = NULL
  uj0 = NULL
  outl = list()
  
  saturated = FALSE
  fout = list()
  
  for(j in 1:n2){ # random effects
    for(i in 1:n1){ # fixed effects
      if((i == 1) & (j == 1)){
        # start from the initial values input (or NULL values)
        coef_old0 = coef_old
        uold = u_init
      }else if((i == 1) & (j != 1)){ # changed from if(i != 1 & j == 1)
        # if re-starting fixed effects penalty loop but moving on to next random effects
        # penalty parameter, take the previous j result for i = 1
        coef_old0 = coef_oldj0
        uold = uj0
        rm(out)
      }else{
        # take the previous values
        coef_old0 = out$coef
        uold = out$u
        rm(out)
      }
      
      gc()
      print("------------------------------------------------------------------")
      print(sprintf("lambda0 %f lambda1 %f", lambda0_range[i], lambda1_range[j]))
      print(sprintf("lambda0 i %i lambda1 j %i", i, j))
      print("------------------------------------------------------------------")
      out = try(fit_dat_B(dat, lambda0 = lambda0_range[i], lambda1 = lambda1_range[j], 
                          nMC_burnin = nMC_burnin, nMC = nMC, nMC_max = nMC_max, nMC_report = nMC_report,
                          family = family, offset_fit = offset, group_X = group_X,
                          penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                          trace = trace, conv_EM = conv_EM, conv_CD = conv_CD,  
                          coef_old = coef_old0, u_init = uold, ufull_describe = ufull_describe,
                          maxitEM = maxitEM, maxit_CD = maxit_CD, t = t,
                          M = M, sampler = sampler, adapt_RW_options = adapt_RW_options,
                          covar = covar, var_start = var_start,
                          max_cores = max_cores))
      

      if(is.character(out)) next
      
      # for restarting the next j for i = 1
      if(i == 1){
        coef_oldj0 = out$coef
        uj0 = out$u 
      }
      
      BICh = out$BICh
      print(BICh)
      if(!is.numeric(BICh)) BICh = Inf
      if(length(BICh) != 1) BICh = Inf
      
      if(BICh < BIChold){
        fout[["BICh"]] = out
        BIChold = BICh
      }
      
      BIC = out$BIC
      print(BIC)
      if(!is.numeric(BIC)) BIC = Inf
      if(length(BIC) != 1) BIC = Inf
      
      if(BIC < BICold){
        fout[["BIC"]] = out
        BICold = BIC
      }
      
      BICq = out$BICq
      if(!is.na(BICq)){
        if(!is.numeric(BICq)) BICq = Inf
        if(length(BICq) != 1) BICq = Inf
        
        if(BICq < BICqold){
          fout[["BICq"]] = out
          BICqold = BICq
        }
      }
      
      res[(j-1)*n1+i,1] = lambda0_range[i]
      res[(j-1)*n1+i,2] = lambda1_range[j]
      res[(j-1)*n1+i,3] = BICh
      res[(j-1)*n1+i,4] = BIC
      res[(j-1)*n1+i,5] = BICq
      # res[(j-1)*n1+i,4] = out$llb
      res[(j-1)*n1+i,6] = out$ll
      res[(j-1)*n1+i,7] = sum(out$coef[2:ncol(dat$X)] !=0)
      res[(j-1)*n1+i,8] = sum(diag(out$sigma) !=0)
      res[(j-1)*n1+i,9] = sum(out$coef !=0)
      # if(maxitEM > 1) res[(i-1)*(n2)+j,8] = loglik(dat = dat, coef = out$coef, u0 = out$u, nMC = 100000, J = out$J)
      # if(maxitEM > 1) res[(j-1)*n1+i,8] = out$ll
      # res[(i-1)*(n2)+j,9] = -2*res[(i-1)*(n2)+j,8] + log(length(dat$y))*sum(out$coef!=0)
      outl[[(j-1)*n1+i]] = 1
      coef = rbind(coef, out$coef)
      print(res[(j-1)*n1+i,])
      
      # Check for model saturation
      ## set null deviance as deviance for model with fixed and random intercepts only
      if((i==1) & (j==1)){
        nullDev = 2*(sat_ll - out$ll)
      }
      ## Calculate deviance
      Dev = 2*(sat_ll - out$ll)
      if(Dev / nullDev < 0.01){
        print("Reached model saturation")
        saturated = TRUE
        break
      }
      
    } # End i loop for lambda0_range
    
    if(saturated){
      break
    }
    
  } # end j loop for lambda1_range
  
  colnames(res) = c("lambda0","lambda1","BICh","BIC","BICq","LogLik","Non_0_fef","Non_0_ref","Non_0_coef")
  colnames(coef) = c("(Intercept)",str_c("B",1:(ncol(dat$X)-1)), str_c("Gamma",1:(ncol(coef)-ncol(dat$X))))
  
  return(list(results = res, out = fout, coef = coef))
  
}

