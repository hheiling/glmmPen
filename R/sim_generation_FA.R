
#' @importFrom stringr str_c str_detect
#' @importFrom ncvreg std
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm rpois rbinom rnbinom
#' @rdname sim.data
#' @export
sim.data.FA = function(n, ptot, pnonzero, nstudies, sd_raneff = 0, 
                       family = "binomial", B = NULL, r = 2,
                       corr = NULL, seed, imbalance = 0, beta = NULL, 
                       pnonzerovar = 0, sd_x = 1){
  
  set.seed(seed = seed)
  
  # set variables
  p = ptot
  p1 = pnonzero
  d = nstudies
  ok = NULL
  
  family_info = family_export(family)
  fam_fun = family_info$family_fun
  link = family_info$link
  link_int = family_info$link_int # Recoded link as integer
  family = family_info$family
  
  slopes = TRUE
  
  if(pnonzero + pnonzerovar > ptot) stop("pnonzero + pnonzerovar > ptot")
  # create fixed effects covariate matrix
  if(is.null(corr)){
    mat = matrix(rnorm(n*p, mean = 0, sd = sd_x), nrow = n) # 11/15 switching to var = 1, then scaling below
    #mat = matrix(rbinom(n*p, p = 0.5, size = 1), nrow = n) # now switching back to normal to have more resolution to show prediction performance
  }else if(is.matrix(corr)){
    if((nrow(corr) != p) | (ncol(corr) != p)){
      stop("corr must be either a single numeric value or a matrix of dimension ptot x ptot")
    }
    sigma = corr
    mat = rmvnorm(n =  n , mean = rep(0,p), sigma = sigma)
  }else{
    cor = matrix(corr, p, p)
    diag(cor) = (sd_x)^2
    sigma = cor # 0.5*cor
    mat = rmvnorm(n =  n , mean = rep(0,p), sigma = sigma)
  }
  
  # add intercept
  if(sd_x == 1){
    mat = std(mat)
  }
  X = cbind(rep(1, n), mat)
  colnames(X) = c("(Intercept)",str_c("X", 1:(ncol(X)-1)))
  
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
  # Random effect covariance matrix
  if(is.null(B)){
    B = matrix(rnorm((ncol(Z0)/d)*r), nrow = ncol(Z0)/d, ncol = r)
  }else{
    if((ncol(B) != r) | (nrow(B) != ncol(Z0)/d)){
      stop("dimensions of B not approrpiate, should have ", ncol(Z0)/d," rows and ", r, " columns")
    }
  }
  Sigma = B %*% t(B) + diag(rep(sd_raneff^2, ncol(Z0)/d))
  z1 = as.numeric(rmvnorm(d, mean = rep(0,ncol(Z0)/d),
                          sigma = Sigma))
  
  eta = X0 %*% matrix(beta, ncol = 1) + Z0 %*% matrix(z1, ncol = 1)
  mu = invlink(link_int, eta)
 
  # simulate random effect and then y
  if(family == "poisson"){
    
    y  = rpois(n, lambda = mu)
    
  }else if(family == "binomial"){
    
    y = rep(NA, n)
    for(ii in 1:n){
      y[ii] = rbinom(n = 1, size = 1, prob = mu[ii])
    }
    if(any(is.na(y))){
      ok = which(is.na(y))
      stop("y resulted in NA values")
    }
    
  }else if(family == "gaussian"){
    
    y = rep(NA, n)
    for(ii in 1:n){
      y[ii] = rnorm(n = 1, mean = mu[ii], sd = 0.5)
    }
    if(any(is.na(y))){
      ok = which(is.na(y))
      stop("y resulted in NA values")
    }
    
  }else if(family == "negbin"){
    # Default variance: theta = 2.0, phi = 1/2.0 = 0.5
    # mu + mu^2 / theta = mu + mu^2 * phi
    y = rep(NA, n)
    for(ii in 1:n){
      y[ii] = rnbinom(n = 1, size = 2.0, mu = mu[ii])
    }
    if(any(is.na(y))){
      ok = which(is.na(y))
      stop("y resulted in NA values")
    }
    
  }else{
    print(family)
    stop("Family not specifed properly")
  }
  
  if(slopes == TRUE) Z = model.matrix(~drep:X-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
  
  colnames(Z) = str_c(rep(colnames(X), each = d), ":", rep(1:d, times = length(colnames(X))))
  
  if(!is.null(ok)){
    dat = list(y = y[-ok], X = X[-ok,], Z = Z[-ok,],  
               pnonzero =  pnonzero, z1 = matrix(z1, nrow = d), group = drep[-ok], 
               X0 = X0, B = B)
  }else{
    dat = list(y = y, X = X, Z = Z,  pnonzero =  pnonzero, z1 = matrix(z1, nrow = d), 
               group = drep, X0 = X0, B = B)
  }
  
  return(dat)
  
}
