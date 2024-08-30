
#' @name sim.data
#' @aliases sim.data.FA sim.data.piecewise.exp sim.data.weibull
#' 
#' @title Simulates data to use for the \code{\link{glmmPen}} package
#' 
#' @description Simulates data to use for testing the \code{\link{glmmPen}} package.
#' \code{sim.data} simulates data for \code{\link{glmmPen}},
#' \code{sim.data.FA} simulates data for \code{\link{glmmPen_FA}},
#' \code{sim.data.piecewise.exp} simulates survival data for \code{\link{phmmPen}} and \code{\link{phmmPen_FA}},
#' and \code{sim.data.weibull} simulates alternative survival data for \code{\link{phmmPen}} and \code{\link{phmmPen_FA}}.
#' Possible parameters to specify includes number of total covariates,
#' number of non-zero fixed and random effects, and the magnitude
#' of the random effect covariance values.
#' 
#' @param n integer specifying total number of samples to generate
#' @param ptot integer specifying total number of covariates to generate 
#' (values randomly generated from the standard normal distribution)
#' @param pnonzero integer specifying how may of the covariates should have
#' non-zero fixed and random effects
#' @param nstudies number of studies/groups to have in the data
#' @param sd_raneff non-negative value specifying the standard deviation of the
#' random effects covariance matrix (applied to the non-zero random effects)
#' @param family character string specifying which family to generate data from.
#' Family options include "binomial" (default), "poisson", and "gaussian".
#' @param corr optional value to specify correlation between covariates
#' in the model matrix. Default \code{NULL}, only available within \code{sim.data}.
#' @param seed integer to use for the setting of a random seed
#' @param imbalance integer of 0 or 1 indicating whether the observations should
#' be equally distributed among the groups (0) or unequally distributed (1).
#' @param beta numeric vector of the fixed effects (including intercept)
#' @param pnonzerovar non-negative integer specifying the number of 
#' covariates with a zero-valued fixed effect but a non-zero random effect.
#' @param sd_x non-negative value specifying the standard deviation of the
#' simulated covariates (drawn from a normal distribution with mean 0,
#' standard deviation \code{sd_x})
#' @param B matrix specifying the factor loadings matrix for the random effects,
#' only used within \code{sim.data.FA}.
#' Dimensions: number of columns is the number of random effects (including the intercept),
#' and number of rows is the number of latent common factors (\code{r}).
#' The random effect covariance matrix is specified as Sigma = B x t(B) + diag(sd_raneff)
#' @param r positive integer specifying number of latent common factors that describe the random effects,
#' only used within \code{sim.data.FA}
#' @param cut_points vector of cut points to use for the time intervals when simulating piecewise 
#' exponential data. Length of cut points must equal length of \code{lhaz_vals}, and the 
#' value of the first cut point must be 0.
#' @param lhaz_vals vector of the log baseline hazard values for each time interval (which
#' correspond to the time intervals defined by the \code{cut_points} argument) within
#' piecewise exponential data. Hazards are assumed to be constant within each time interval.
#' @param cens_type character specifying type of censoring to implement. Default "unif" specifies
#' uniform censoring from 0 to \code{cens_max}. The option "exp" specifies exponential
#' censoring with rate \code{exp_rate}
#' @param cens_max numeric value used to characterize the range of the uniform censoring procedure 
#' (from 0 to \code{cens_max})
#' @param exp_rate numeric value used to characterize the exponential censoring rate (where rate
#' is defined as the rate used in \code{\link[stats]{rexp}})
#' @param lhaz_base numeric value that gives the log of the scale parameter for the Weibull distribution 
#' (for description of Weibull scale parameter without log transformation, 
#' see "lambdas" argument in \code{\link[simsurv]{simsurv}} and the lambda notation in section "Weibull 
#' distribution" in the simsurv vignette 
#' https://cran.r-project.org/web/packages/simsurv/vignettes/simsurv_technical.html)
#' @param alpha_PH numeric value > 0 that gives the shape parameter for the Weibull distribution 
#' (for description of Weibull shape paraemeter,
#' see "gammas" argument in \code{\link[simsurv]{simsurv}} and the gamma notation in section "Weibull
#' distribution" in the simsurv vignette 
#' https://cran.r-project.org/web/packages/simsurv/vignettes/simsurv_technical.html)
#' 
#' @return list containing the following elements:
#' \item{y}{vector of the response}
#' \item{X}{model matrix for the fixed effects}
#' \item{Z}{model matrix for the random effects, organized first by variable
#' and then by group}
#' \item{pnonzero}{number of non-zero fixed effects}
#' \item{z1}{values of the random effects for each variable for each level of 
#' the grouping factor}
#' \item{group}{grouping factor}
#' \item{X0}{model matrix for just the non-zero fixed effects}
#' 
#' @importFrom stringr str_c str_detect
#' @importFrom ncvreg std
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm rpois rbinom rnbinom
#' @export
sim.data = function(n, ptot, pnonzero, nstudies, sd_raneff = 1, family = "binomial", 
                    corr = NULL, seed, imbalance = 0, 
                    beta = NULL, pnonzerovar = 0, sd_x = 1){
  
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
  if(slopes == TRUE) Z = model.matrix(~drep:X-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
  
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
  if(slopes == TRUE) Z0 = model.matrix(~drep:X0-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
  z1 = as.numeric(rmvnorm(d, mean = rep(0,ncol(Z0)/d), 
                          sigma = diag(rep(sd_raneff^2, ncol(Z0)/d), nrow = ncol(Z0)/d, ncol = ncol(Z0)/d)))
  
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
  
  ## 10/26/2016 scaling predictors and remaking Z matrix on scaled predictors
  #X = cbind(X[,1], scale(X[,-1]))
  #Z = model.matrix(~drep-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
  if(slopes == TRUE) Z = model.matrix(~drep:X-1,  contrasts.arg=list(drep=diag(nlevels(drep))))
  ##
  
  colnames(Z) = str_c(rep(colnames(X), each = d), ":", rep(1:d, times = length(colnames(X))))
  
  if(!is.null(ok)){
    dat = list(y = y[-ok], X = X[-ok,], Z = Z[-ok,],  pnonzero =  pnonzero, z1 = matrix(z1, nrow = d), group = drep[-ok], X0 = X0)
  } else{
    dat = list(y = y, X = X, Z = Z,  pnonzero =  pnonzero, z1 = matrix(z1, nrow = d), group = drep, X0 = X0)
  }
  
  return(dat)
  
}
