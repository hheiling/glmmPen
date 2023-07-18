

#' @importFrom stringr str_c str_detect
#' @importFrom ncvreg std
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm rexp runif
#' @rdname sim.data
#' @export
sim.data.piecewise.exp = function(n, ptot, pnonzero, nstudies, sd_raneff = 0, 
                                  B = NULL, r = 2,
                                  cut_points = c(0, 0.5, 1.0, 1.5, 2.0), 
                                  lhaz_vals = c(-1.5,1.0,2.7,3.7,6.8),
                                  cens_type = c("exp","unif"), cens_max = 5, exp_rate = 0.15,
                                  seed, imbalance = 0, beta = NULL, 
                                  pnonzerovar = 0, sd_x = 1){
  
  # Input checks
  if(length(cut_points) != length(lhaz_vals)) stop("length of cut_points and lhaz_vals must be equal")
  if(cut_points[1] != 0) stop("the first cut-point must be 0")
  
  set.seed(seed = seed)
  
  # set variables
  p = ptot
  p1 = pnonzero
  d = nstudies
  
  slopes = TRUE
  
  if(length(cens_type) > 1){
    cens_type = cens_type[1]
  }
  
  if(pnonzero + pnonzerovar > ptot) stop("pnonzero + pnonzerovar > ptot")
  # create fixed effects covariate matrix
  mat = matrix(rnorm(n*p, mean = 0, sd = sd_x), nrow = n)
  
  # standardize covariates if sd_x = 1
  if(sd_x == 1){
    mat = std(mat)
  }
  X = mat
  colnames(X) = str_c("X", 1:(ncol(X)))
  
  # create raneff matrix 
  drep = factor(rep(1:d, each = n/d))
  if(imbalance == 1){
    first = rep(1, floor(n/3))
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
  if(length(drep) != n){ # d not a factor of n, imbalance = 0
    # randomly assign remaining subjects to a group
    drep = c(drep, sample(1:d, size = (n - length(drep)), replace = FALSE))
    drep = factor(drep)
  }
  
  if(is.null(beta)){
    if(pnonzerovar > 0){
      beta <- c(0, rep(2, p1), rep(0, pnonzerovar))
      X0 = X[,1:(p1+pnonzerovar)]
    }else{
      beta <- c(0, rep(2, p1))
      X0 = X[,1:(p1)]
    }
  }else{
    if(length(beta) < p1+pnonzerovar) beta = c(beta, rep(0, pnonzerovar))
    X0 = X[,1:(p1+pnonzerovar)]
  }
  
  # Include random intercept (aka frailty)
  X0_tmp = cbind(rep(1,n), X0)
  if(slopes == TRUE) Z0 = model.matrix(~drep:X0_tmp-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
  
  # Random effect covariance matrix
  if(is.null(B) & (sd_raneff == 0)){
    stop("Either B must be specified, or sd_raneff must be > 0")
  }
  if(!is.null(B)){
    if((ncol(B) != r) | (nrow(B) != ncol(Z0)/d)){
      stop("dimensions of B not approrpiate, should have ", ncol(Z0)/d," rows and ", r, " columns")
    }
    Sigma = B %*% t(B) + diag(rep(sd_raneff^2, ncol(Z0)/d))
  }else{
    Sigma = diag(rep(sd_raneff^2), ncol(Z0)/d)
  }
  
  # Simulate random effects, calculate linear predictor
  z1 = as.numeric(rmvnorm(d, mean = rep(0,ncol(Z0)/d),
                          sigma = Sigma))
  
  eta = X0 %*% matrix(beta, ncol = 1) + Z0 %*% matrix(z1, ncol = 1)
  
  # Simulate survival times using piecewise exponential modeling assumptions
  T_i <- rep(NULL, n)
  
  for(i in 1:n){
    for(j in 1:length(cut_points)){
      int = rexp(1, exp(lhaz_vals[j] + eta[i]))
      if(j < length(cut_points)){
        if(int + cut_points[j] < cut_points[j+1]){
          T_i[i] = int + cut_points[j]
          break
        }
      }else if(j == length(cut_points)){
        T_i[i] = cut_points[j] + int
      }
    }
  }
  
  # Simulate observed times (event vs censoring)
  ## Censoring times generated independently from exponential distribution 
  ## with parameters selected to induce approx 20% censoring
  if(cens_type == "exp"){
    C_i = rexp(n, exp_rate)
  }else if(cens_type == "unif"){
    C_i = runif(n, min = 0, max = cens_max)
  }
  observed_ti = pmin(T_i, C_i)
  status = delta_i = ifelse((T_i <= C_i),1,0) # 1 = event; 0 = censor
  death_rate = sum(delta_i)/n # check if this is close to a desired percentage, e.g. 80%
  censor_i = 1 - delta_i
  censor_rate = sum(censor_i)/n
  
  # Include random intercept
  X_tmp = cbind(rep(1,n), X)
  if(slopes == TRUE) Z = model.matrix(~drep:X_tmp-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
  
  colnames(Z) = str_c(rep(c("Intercept",colnames(X)), each = d), ":", rep(1:d, times = length(colnames(X))+1))
  
  dat = list(y_status = status, y_times = observed_ti, X = X, Z = Z,  
             pnonzero = pnonzero, z1 = matrix(z1, nrow = d), 
             group = drep, X0 = X0, B = B,
             death_rate = death_rate, censor_rate = censor_rate)
  
  return(dat)
  
}



# @importFrom stringr str_c str_detect
# @importFrom ncvreg std
# @importFrom mvtnorm rmvnorm
# @importFrom stats rnorm rexp runif
# @export
# sim.data.coxph = function(n, ptot, pnonzero, nstudies, sd_raneff = 0, 
#                           B = NULL, r = 2,
#                           lhaz_base = 0.5, alpha_PH = 2, exp_rate = 0.15,
#                           cens_type = c("exp","unif"), cens_max = 5,
#                           seed, imbalance = 0, beta = NULL, 
#                           pnonzerovar = 0, sd_x = 1){
#   
#   set.seed(seed = seed)
#   
#   # set variables
#   p = ptot
#   p1 = pnonzero
#   d = nstudies
#   
#   slopes = TRUE
#   
#   if(length(cens_type) > 1){
#     cens_type = cens_type[1]
#   }
#   
#   if(pnonzero + pnonzerovar > ptot) stop("pnonzero + pnonzerovar > ptot")
#   # create fixed effects covariate matrix
#   mat = matrix(rnorm(n*p, mean = 0, sd = sd_x), nrow = n)
#   
#   # standardize covariates if sd_x = 1
#   if(sd_x == 1){
#     mat = std(mat)
#   }
#   X = mat
#   colnames(X) = str_c("X", 1:(ncol(X)))
#   
#   # create raneff matrix 
#   drep = factor(rep(1:d, each = n/d))
#   if(imbalance == 1){
#     first = rep(1, floor(n/3))
#     second = rep(2:d, each = ceiling((2*n/3)/(d-1)))
#     if(length(first) + length(second) < n){
#       drep = factor(c(first, second, rep(d, length(drep) - length(first) - length(second))))
#     }else if(length(first) + length(second) > n){
#       drep = factor(c(first, second))
#       drep = drep[1:n]
#     }else{
#       drep = factor(c(first, second))
#     }
#   }
#   if(length(drep) != n){ # d not a factor of n, imbalance = 0
#     # randomly assign remaining subjects to a group
#     drep = c(drep, sample(1:d, size = (n - length(drep)), replace = FALSE))
#     drep = factor(drep)
#   }
#   
#   if(is.null(beta)){
#     if(pnonzerovar > 0){
#       beta <- c(0, rep(2, p1), rep(0, pnonzerovar))
#       X0 = X[,1:(p1+pnonzerovar)]
#     }else{
#       beta <- c(0, rep(2, p1))
#       X0 = X[,1:(p1)]
#     }
#   }else{
#     if(length(beta) < p1+pnonzerovar) beta = c(beta, rep(0, pnonzerovar))
#     X0 = X[,1:(p1+pnonzerovar)]
#   }
#   
#   # Include random intercept (aka frailty)
#   X0_tmp = cbind(rep(1,n), X0)
#   if(slopes == TRUE) Z0 = model.matrix(~drep:X0_tmp-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
#   
#   # Random effect covariance matrix
#   if(is.null(B) & (sd_raneff == 0)){
#     stop("Either B must be specified, or sd_raneff must be > 0")
#   }
#   if(!is.null(B)){
#     if((ncol(B) != r) | (nrow(B) != ncol(Z0)/d)){
#       stop("dimensions of B not approrpiate, should have ", ncol(Z0)/d," rows and ", r, " columns")
#     }
#     Sigma = B %*% t(B) + diag(rep(sd_raneff^2, ncol(Z0)/d))
#   }else{
#     Sigma = diag(rep(sd_raneff^2), ncol(Z0)/d)
#   }
#   
#   # Simulate random effects, calculate linear predictor
#   z1 = as.numeric(rmvnorm(d, mean = rep(0,ncol(Z0)/d),
#                           sigma = Sigma))
#   
#   eta = X0 %*% matrix(beta, ncol = 1) + Z0 %*% matrix(z1, ncol = 1)
#   mu = exp(lhaz_base + eta)
#   
#   # Simulate survival times - see simsurv vignette for details
#   ## Covariates are related to the Weibull event time t using survival function
#   ##    S(t) = exp(-(t^alpha_PH) * mu)
#   ##    t = (-log(S_t) / mu)^{1/alpha_PH}
#   S_t = runif(n, 0, 1)
#   T_i = (-log(S_t) / mu)^(1/alpha_PH)
#   St_check = exp(-(T_i^alpha_PH) * mu)
#   # Check of simulated data
#   if(!(sum(round(St_check,3) == round(S_t,3)) == n)){
#     print("ERROR: S_t and T_i generated wrong.") # should equal n
#   }
#   
#   # Simulate observed times (event vs censoring)
#   ## Censoring times generated independently from exponential distribution 
#   ## with parameters selected to induce approx 20% censoring
#   if(cens_type == "exp"){
#     C_i = rexp(n, exp_rate)
#   }else if(cens_type == "unif"){
#     C_i = runif(n, min = 0, max = cens_max)
#   }
#   observed_ti = pmin(T_i, C_i)
#   status = delta_i = ifelse((T_i <= C_i),1,0) # 1 = event; 0 = censor
#   death_rate = sum(delta_i)/n # check if this is close to a desired percentage, e.g. 80%
#   censor_i = 1 - delta_i
#   censor_rate = sum(censor_i)/n
#   
#   # Include random intercept
#   X_tmp = cbind(rep(1,n), X)
#   if(slopes == TRUE) Z = model.matrix(~drep:X_tmp-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
#   
#   colnames(Z) = str_c(rep(c("Intercept",colnames(X)), each = d), ":", rep(1:d, times = length(colnames(X))+1))
#   
#   dat = list(y_status = status, y_times = observed_ti, X = X, Z = Z,  
#              pnonzero = pnonzero, z1 = matrix(z1, nrow = d), 
#              group = drep, X0 = X0, B = B,
#              death_rate = death_rate, censor_rate = censor_rate)
#   
#   return(dat)
#   
# }

