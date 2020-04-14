
# Notes: 
# Eventually add grouped versions of MCP, SCAD, and lasso
# Eventually add lambda0 and lambda1 (for fixed and random -- grouped case)

# group_X and group_Z: no penalization for covariates with 0 as group value (e.g. intercept)

#' @export
M_step = function(y, X, family, coef, offset = NULL,
                  group_X = 0:(ncol(X)-1),
                  maxit = 250, maxit_CD = 250, conv = 0.0001, fit_type,
                  penalty = c("MCP","SCAD","lasso"), lambda,
                  gamma = switch(penalty, SCAD = 3.7, 3.0), alpha = 1.0){
  
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
    # useOffset = 1
    if((! is.numeric(offset)) || length(offset) != N){
      stop("offset must be a numeric vector of the same length as y\n")
    }
  }else{
    # useOffset = 0
    offset    = rep(0.0, times = N)
  }
  
  # if(family == "binomial"){
  #   if((! is.numeric(nTotal)) || length(nTotal) != N){
  #     strWarn = "For binomial family, nTotal must be a numeric vector"
  #     stop(strWarn, "of the same length as y\n")
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
  
  # Optional:
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
  }else if(alpha == 0.0){
    stop("alpha cannot equal 0. Pick a small value > 0 instead (e.g. 0.001) \n");
  }
  
  penalty_params = c(gamma, alpha)
  
  if(length(group_X) != ncol(X)){
    stop("length of group_X must equal number columns of X \n")
  }else if(!is.integer(group_X)){
    tmp <- try(group_X <- as.integer(group_X), silent = T)
    if(inherits(tmp, "try_error")) stop("group_X must be a vector of integers \n")
  }
  
  # J_X = number groups in X matrix
  J_X = length(unique(group_X))
  if((max(group_X) - min(group_X) + 1) != J_X){ # Assume min(group_X) == 0 (for intercept)
    stop("group_X must be comprised of consecutive integers \n")
  }
  # Size of groups in X
  K_X = numeric(J_X)
  for(j in unique(group_X)){
    K_X[j+1] = sum(group_X == j) # Adjust by one since 0 included
  }
  
  dims = c(p, N, conv, maxit, maxit_CD, J_X)
  
  coef_new = glm_fit(y, X, dims, coef, offset, familyR, link_int, fit_type, group_X, K_X, penalty, lambda, penalty_params)
  
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

#' @export
M_step2 = function(y, X, Z, u, J, group, family, coef, offset = NULL,
                  maxit = 250, maxit_CD = 250, conv = 0.0001,
                  penalty = c("MCP","SCAD","lasso"), lambda,
                  gamma = switch(penalty, SCAD = 4.0, 3.0), alpha = 1.0){
  
  if (!is.double(y)) {
    tmp <- try(y <- as.double(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be numeric or able to be coerced to numeric", call.=FALSE)
  }
  
  if(!is.matrix(X)){
    stop("X must be a matrix \n")
  }else if(typeof(X)=="integer") storage.mode(X) <- "double"
  
  if(!is.matrix(Z)){
    stop("X must be a matrix \n")
  }else if(typeof(X)=="integer") storage.mode(X) <- "double"
  
  p = ncol(X)
  N = length(y)
  M = nrow(u)
  
  if(nrow(X) != N){
    stop("the dimension of X and y do not match\n")
  }
  
  if(!is.null(offset)){
    if((! is.numeric(offset)) || length(offset) != N){
      stop("offset must be a numeric vector of the same length as y\n")
    }
  }else{
    offset    = rep(0.0, times = N)
  }
  
  if(!is.factor(group)){
    group <- factor(group)
  }
  
  familyR = family$family
  linkR = family$link
  
  if(!(familyR %in% c("binomial","poisson","gaussian","Gamma"))){
    stop("invalid family \n")
  }
  if(!(linkR %in% c("logit","log","identity","inverse"))){
    stop("invalid link \n")
  }
  
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
  }else if(alpha == 0.0){
    stop("alpha cannot equal 0. Pick a small value > 0 instead (e.g. 0.001) \n");
  }
  
  penalty_params = c(gamma, alpha)
  
  # Number of groups wrt observations
  d = nlevels(group)
  # Number of random effect variables
  q = ncol(Z) / d
  
  # For now, assume that covariates in X are not grouped
  if(ncol(J) == q){
    XZ_group = 0:(ncol(X)-1+q) # 0 for intercept
    K = rep(1, times = length(XZ_group))
  }else{ # ncol(J) = q(q+1)/2
    XZ_group = 0:(ncol(X)-1)
    K = rep(1, times = ncol(X))
    for(j in 1:q){
      XZ_group = c(XZ_group, rep(ncol(X)-1+j, times = j))
      K = c(K, rep(j, times = j))
    }
  }
  
  # Number of groups wrt covariates (cols of X and Z)
  J_XZ = max(XZ_group)+1
  
  dims = c(p, N, d, q, M, J_XZ, conv, maxit)
  
  coef_new = grp_CD_XZ(y, X, Z, group, u, J, dims, coef, offset, familyR, link_int, XZ_group, K, penalty, lambda, penalty_params)
  
  return(coef_new)
}


#################################################################################################
