
# Notes: 
# Eventually add grouped versions of MCP, SCAD, and lasso
# Eventually add lambda0 and lambda1 (for fixed and random -- grouped case)
# Eventually add alpha option (default = 1), for elastic net?

#' @export
M_step = function(y, X, family, link, coef, offset = NULL, 
                  maxit = 250, maxit_CD = 250, conv = 0.0001, fit_type,
                  penalty = c("MCP","SCAD","lasso"), lambda,
                  gamma = switch(penalty, SCAD = 3.7, 3.0), alpha = 1){
  
  if (!is.double(y)) {
    tmp <- try(y <- as.double(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be numeric or able to be coerced to numeric", call.=FALSE)
  }
  
  if(!is.matrix(X)){
    stop("X must be a matrix \n")
  }else if(typeof(X)=="integer") storage.mode(X) <- "double"
  
  # if(is.null(prior)) prior=rep(1,length(y))
  # 
  # if(is.null(prop)) prop=1
  
  p = ncol(X)
  N = length(y)
  
  if(nrow(X) != N){
    stop("the dimension of X and y do not match\n")
  }
  
  if(!is.null(offset)){
    useOffset = 1
    if((! is.numeric(offset)) || length(offset) != N){
      stop("offset must be a numeric vector of the same length as y\n")
    }
  }else{
    useOffset = 0
    offset    = rep(0.0, times = N)
  }
  
  # if(family == "binomial"){
  #   if((! is.numeric(nTotal)) || length(nTotal) != N){
  #     strWarn = "For binomial family, nTotal must be a numeric vector"
  #     stop(strWarn, "of the same length as y\n")
  #   }
  # }
  
  # if(is.null(fitted)){
  #   init = 0
  #   fitted = rep(0, N)
  # }else{
  #   init = 1
  #   if((! is.numeric(fitted)) || length(fitted) != N){
  #     stop("fitted must be a numeric vector of the same length as y\n")
  #   }
  # }
  
  familyR = family$family
  linkR = family$link
  
  if(!(familyR %in% c("binomial","poisson","gaussian","Gamma"))){
    stop("invalid family \n")
  }
  if(!(linkR %in% c("logit","log","identity","inverse"))){
    stop("invalid link \n")
  }
  
  # Re-code family as integer
  # Binomial 1
  # Poisson  2
  # Gaussian 3
  # Gamma    4
  
  # if(familyR=="binomial"){
  #   family_int = 1
  # }else if(familyR=="poisson"){
  #   family_int = 2
  # }else if(familyR=="gaussian"){
  #   family_int = 3
  # }else if(familyR=="Gamma"){
  #   family_int = 4
  # }else{
  #   stop("invalid family\n")
  # }
  
  # Re-code link as integer
  ## All link_int will have two digits
  ## First digit corresponds to family that link is canonical for
  ## 1 = binomial, 2 = poisson, 3 = gaussian, 4 = gamma
  ## Second digit is arbitrary enumeration of links
  if(linkR == "logit"){
    link_int = 10
  }else if(linkR == "probit"){
    link_int = 11
  }else if(linkR == "cloglog"){
    link_int = 12
  }else if(linkR == "log"){
    link_int = 20
  }else if(linkR == "identity"){
    link_int = 30
  }else if(linkR == "inverse"){
    link_int = 40
  }
  
  penalty = penalty[1]
  
  if(!(penalty %in% c("lasso","MCP","SCAD"))){
    stop("penalty ", penalty, " not available, must choose 'lasso', 'MCP', or 'SCAD' \n")
  }
  
  if(penalty == "MCP" & gamma <= 1){
    stop("gamma must be > 1 when using MCP penalty")
  }else if(penalty == "SCAD" & gamma <= 2){
    stop("gama must be > 2 when using SCAD penalty")
  }else if(!is.double(gamma)){
    tmp <- try(gamma <- as.double(gamma), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("gamma must be numeric or able to be coerced to numeric", call.=FALSE)
    
  }
  
  if (!is.double(alpha)) {
    tmp <- try(alpha <- as.double(alpha), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("alpha must be numeric or able to be coerced to numeric", call.=FALSE)
  }
  
  penalty_params = c(gamma, alpha)
  
  dims = c(p, N, conv, maxit, maxit_CD)
  
  coef_new = glm_fit(y, X, dims, coef, family, link_int, fit_type, penalty, lambda, penalty_params)
  
  # Calculate BIC for model (when no random effects)
  if(familyR == "binomial" & linkR == "logit"){
    eta = X %*% coef_new
    p = exp(eta) / (1+exp(eta))
    ll = sum(dbinom(y, size = 1, prob = p, log = T))
  }
  
  BIC = -2*ll + sum(coef!=0)*N
  
  return(list(coef = coef_new, ll = ll, BIC = BIC))
}

#################################################################################################