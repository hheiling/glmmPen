#' @export
sim.data2 = function(n, ptot, pnonzero, nstudies, sd_raneff = 1, family = "binomial", corr = NULL, slopes = F, seed, imbalance = 0, beta = NULL,pnonzerovar = 0, trace = 0){
  set.seed(seed = seed)
  library(mvtnorm)
  # set variables
  p = ptot
  p1 = pnonzero
  d = nstudies
  ok = NULL
  
  if(pnonzero + pnonzerovar > ptot) stop("pnonzero + pnonzerovar > ptot")
  # create fixed effects covariate matrix
  if(is.null(corr)){
    mat = matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n) # 11/15 switching to var = 1, then scaling below
    #mat = matrix(rbinom(n*p, p = 0.5, size = 1), nrow = n) # now switching back to normal to have more resolution to show prediction performance
  }else{
    library(mvtnorm)
    cor = matrix(corr, p, p)
    diag(cor) = 1
    sigma = 0.5*cor
    mat = rmvnorm(n =  n , mean = rep(0,p), sigma = sigma)
  }
  
  # add intercept
  mat = scale(mat)
  X = cbind(rep(1, n), mat)
  
  # create raneff matrix (assuming only 1 random effect with nstudies levels for now)
  drep = factor(rep(1:d, each = n/d))
  if(imbalance == 1){
    first = rep(1, floor(n/3)) ## change to 1/2?
    second = rep(2:d, each = ceiling((2*n/3)/(d-1)))
    if(length(first) + length(second) < n){
      drep = factor(c(first, second, rep(d, length(drep) - length(first) - length(second))))
    }else if(length(first) + length(second) > n){
      drep = factor(c(first, second))
      drep = drep[1:n]
    }else{
      drep = factor(c(first, second))
    }
  }
  Z = model.matrix(~drep-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
  if(slopes == T) Z = model.matrix(~drep:X-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
  
  # simulate random effect and then y
  if(family == "poisson"){
    
    # create coefficient matrix for fixed effects, leaving fixed for now
    beta = matrix(c(0, rep(1, p1), rep(0, p - p1)), ncol = 1)
    mu = exp(X%*%beta + Z%*%matrix(z1, ncol = 1))
    y  = rpois(n, lambda = mu)
  }else if(family == "binomial"){
    
    if(is.null(beta)){
      if(pnonzerovar > 0){
        beta <- c(0, rep(2, p1), rep(0, pnonzerovar))
        X0 = X[,1:(p1+pnonzerovar+1)]
      }else{
        beta <- c(0, rep(2, p1))
        X0 = X[,1:(p1+1)]
      }
    }else{
      if(length(beta) < p1+pnonzerovar+1) beta = c(beta, rep(0, pnonzerovar))
      X0 = X[,1:(p1+pnonzerovar+1)]
    }
    if(slopes == T) Z0 = model.matrix(~drep:X0-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
    z1 = as.numeric(rmvnorm(d, mean = rep(0,ncol(Z0)/d), sigma = diag(rep(sd_raneff^2, ncol(Z0)/d), nrow = ncol(Z0)/d, ncol = ncol(Z0)/d)))
    
    
    # create coefficient matrix for fixed effects, leaving fixed for now
    mu0001 = rowSums(X0 * matrix(beta, ncol = length(beta), nrow = nrow(X0), byrow = T))  #X0%*%matrix(beta, ncol = 1)
    mu0002 = Z0%*%matrix(z1, ncol = 1)
    mu000 = mu0001 + mu0002
    mu00 = as.numeric(mu000) #- apply((X0%*%beta + Z0%*%matrix(z1, ncol = 1)), 1, FUN = function(x){logsumexp(c(0, x))})
    mu0 = exp(mu00)
    mu = mu0/(1+mu0)
    
    #y  = rbinom(n = n, size = 1, prob = mu)
    y = rep(NA, n)
    for(ii in 1:n){
      y[ii] = rbinom(n = 1, size = 1, prob = mu[ii])
    }
    if(any(is.na(y))){
      ok = which(is.na(y))
      print("")
      print(mu[ok])
      print(mu0[ok])
      print(mu00[ok])
      print(mu000[ok])
      print(dim(mu0001))
      print(mu0001)
      print(dim(mu0002))
      print(mu0002[ok])
      print(X0)
      print("")
      #print((Z0%*%matrix(z1, ncol = 1))[ok])
      #print((X0%*%beta + Z0%*%matrix(z1, ncol = 1))[ok])
      #print(exp(X0%*%beta + Z0%*%matrix(z1, ncol = 1))[ok])
    }
    
  }else{
    print(family)
    stop("Family not specifed properly")
  }
  
  ## 10/26/2016 scaling predictors and remaking Z matrix on scaled predictors
  #X = cbind(X[,1], scale(X[,-1]))
  #Z = model.matrix(~drep-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
  if(slopes == T) Z = model.matrix(~drep:X-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
  ##
  
  if(!is.null(ok)){
    dat = list(y = y[-ok], X = X[-ok,], Z = Z[-ok,],  pnonzero =  pnonzero, z1 = matrix(z1, nrow = d), group = drep[-ok], X0 = X0)
  } else{
    dat = list(y = y, X = X, Z = Z,  pnonzero =  pnonzero, z1 = matrix(z1, nrow = d), group = drep, X0 = X0)
  }
  return(dat)
}
