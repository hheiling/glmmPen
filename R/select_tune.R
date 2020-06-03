
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
#' @export
select_tune = function(dat, offset = NULL, family, group_X = 0:(ncol(dat$X)-1),
                       penalty, lambda0_range, lambda1_range,
                       alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0),
                       returnMC = T, trace = 0, ufull = NULL, coeffull = NULL, 
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
  nMC = optim_options$nMC
  nMC_max = optim_options$nMC_max
  maxitEM = optim_options$maxitEM
  maxit_CD = optim_options$maxit_CD
  M = optim_options$M
  t = optim_options$t
  covar = optim_options$covar
  MwG_sampler = optim_options$MwG_sampler
  gibbs = optim_options$gibbs
  fit_type = optim_options$fit_type
  
  covar = covar[1]
  if(!(covar %in% c("unstructured","independent"))){
    stop("algorithm currently only handles 'unstructured' or 'independent' covariance structure \n")
  }
  
  # Calculate loglik for null model (mu set to y; n parameters, one per observation)
  if(fam == "binomial"){
    null_ll = 0
  }else if(fam == "poisson"){
    stop("poisson family not yet operational \n")
    null_ll = 0
    y = dat$y
    for(i in 1:length(y)){
      if(y[i] > 0){
        null_ll = null_ll + dpois(x = y[i], lambda = y[i], log = T)
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
  
  n1 = length(lambda0_range)
  n2 = length(lambda1_range)
  BICold = BIC = BIC2old2 = BIC2= Inf
  BIChold = BICh = Inf
  res = matrix(0, n1*n2, 8)
  coef = NULL
  outl = list()
  
  saturated = FALSE
  fout = list()
  
  for(j in 1:n2){
    for(i in 1:n1){
      if(i == 1 & j == 1){
        # start from the full model values
        coeffull0 = coeffull
        ufullinit = ufull
      }else if(i != 1 & j == 1){
        # take the previous i value for j = 1
        coeffull0 = coeffulli0
        ufullinit = ui0
        rm(out)
      }else{
        # take the previous j values
        coeffull0 = out$coef
        ufullinit = out$u
        rm(out)
      }
      
      gc()
      print("------------------------------------------------------------------")
      print(sprintf("lambda0 %f lambda1 %f", lambda0_range[i], lambda1_range[j]))
      print("------------------------------------------------------------------")
      out = try(fit_dat_B(dat, lambda0 = lambda0_range[i], lambda1 = lambda1_range[j], 
                          nMC = nMC, family = family, offset_fit = offset, group_X = group_X,
                          penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                          trace = trace, conv_EM = conv_EM, conv_CD = conv_CD, nMC_max = nMC_max, 
                          ufull = ufull, coeffull = coeffull0, 
                          gibbs = gibbs, maxitEM = maxitEM, maxit_CD = maxit_CD,
                          returnMC = returnMC, ufullinit = ufullinit, t = t,
                          M = M, MwG_sampler = MwG_sampler, adapt_RW_options = adapt_RW_options,
                          covar = covar, fit_type = fit_type))
      
      
      if(is.character(out)) next
      
      # for restarting the next i for j = 1
      if(j == 1){
        coeffulli0 = out$coef
        ui0 = out$u 
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
      
      res[(j-1)*n1+i,1] = lambda0_range[i]
      res[(j-1)*n1+i,2] = lambda1_range[j]
      res[(j-1)*n1+i,3] = BICh
      res[(j-1)*n1+i,4] = BIC
      # res[(j-1)*n1+i,4] = out$llb
      res[(j-1)*n1+i,5] = out$ll
      res[(j-1)*n1+i,6] = sum(out$coef[2:ncol(dat$X)] !=0)
      res[(j-1)*n1+i,7] = sum(diag(out$sigma) !=0)
      res[(j-1)*n1+i,8] = sum(out$coef !=0)
      # if(maxitEM > 1) res[(i-1)*(n2)+j,8] = loglik(dat = dat, coef = out$coef, u0 = out$u, nMC = 100000, J = out$J)
      # if(maxitEM > 1) res[(j-1)*n1+i,8] = out$ll
      # res[(i-1)*(n2)+j,9] = -2*res[(i-1)*(n2)+j,8] + log(length(dat$y))*sum(out$coef!=0)
      outl[[(j-1)*n1+i]] = 1
      coef = rbind(coef, out$coef)
      print(res[(j-1)*n1+i,])
      
      # Check for model saturation
      ## set null deviance as deviance for model with fixed and random intercepts only
      if((i==1) & (j==1)){
        nullDev = 2*(null_ll - out$ll)
      }
      ## Calculate deviance
      Dev = 2*(null_ll - out$ll)
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
  
  colnames(res) = c("lambda0","lambda1","BICh","BIC","LogLik","Non_0_fef","Non_0_ref", "Non_0_coef")
  colnames(coef) = c("(Intercept)",str_c("B",1:(ncol(dat$X)-1)), str_c("Gamma",1:(ncol(coef)-ncol(dat$X))))
  
  return(list(results = res, out = fout, coef = coef))
  
}

