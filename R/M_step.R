
# group_X and group_Z: no penalization for covariates with 0 as group value (e.g. intercept)

#' @export
M_step = function(y, X, family, coef, offset = NULL,
                  group_X = 0:(ncol(X)-1), 
                  maxit = 250, maxit_CD = 250, conv = 0.0001, fit_type,
                  penalty = c("MCP","SCAD","lasso"), lambda,
                  gamma = switch(penalty[1], SCAD = 3.7, 3.0), alpha = 1.0){
  
  if (!is.double(y)) {
    tmp <- try(y <- as.double(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be numeric or able to be coerced to numeric", call.=FALSE)
  }
  
  if(!is.matrix(X)){
    stop("X must be a matrix \n")
  }else if(typeof(X)=="integer") storage.mode(X) <- "double"
  
  p = ncol(X)
  N = length(y)
  
  if(nrow(X) != N){
    stop("the dimension of X and y do not match\n")
  }
  
  if(!is.null(offset)){
    if((!is.numeric(offset)) | (length(offset) != N)){
      stop("offset must be a numeric vector of the same length as y\n")
    }
  }else{
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
  
  return(coef_new)
}

#################################################################################################

# Eventually implement the following:
## Add lambda0 and lambda1 (for fixed and random -- grouped case)
## In fit_dat or glmmPen, check that group_X made of subsequent integers
## Change Z to be sparse (fit_dat and glmmPen)

# init: if this is the first attempt at a fit (using initial coef)
#' @export
M_stepB = function(y, X, Z, u_address, M, J, group, family, coef, offset = NULL,
                  maxit_CD = 250, conv_CD = 0.0001,
                  init, group_X = 0:(ncol(X)-1), covgroup,
                  penalty = c("MCP","SCAD","lasso"), lambda0, lambda1,
                  gamma = switch(penalty[1], SCAD = 4.0, 3.0), alpha = 1.0,
                  fit_type = 1){
  
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
  
  if(nrow(X) != N){
    stop("the dimension of X and y do not match\n")
  }
  
  if(!is.null(offset)){
    if((!is.numeric(offset)) | (length(offset) != N)){
      stop("offset must be a numeric vector of the same length as y\n")
    }
  }else{
    offset = rep(0.0, times = N)
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
  
  penalty_params = c(lambda0, lambda1, gamma, alpha)
  
  # index of covariate groups
  XZ_group = c(group_X, (covgroup+max(group_X)))
  # Number of groups wrt covariates (cols of X and Z)
  J_XZ = length(unique(XZ_group))
  if((max(XZ_group) - min(XZ_group) + 1) != J_XZ){ # Assume min(group_X) == 0 (for intercept)
    stop("group_X must be comprised of consecutive integers \n")
  }
  # Size of covariate groups
  K = numeric(J_XZ)
  for(j in unique(XZ_group)){
    idx = which(XZ_group == j)
    K[j+1] = length(idx)
  }
  
  # Number of groups wrt observations
  d = nlevels(group)
  # Number of random effect variables
  q = ncol(Z) / d
  
  dims = c(p, N, d, q, M, J_XZ, conv_CD, maxit_CD)
  
  if(fit_type == 1){
    coef_new = grp_CD_XZ_B1(y, X, Z, group, u_address, J, dims, coef, offset, familyR, link_int, init, XZ_group, K, penalty, penalty_params)
  }else if(fit_type == 2){
    coef_new = grp_CD_XZ_B2(y, X, Z, group, u_address, J, dims, coef, offset, familyR, link_int, init, XZ_group, K, penalty, penalty_params)
  }else if(fit_type == 3){
    coef_new = grp_CD_XZ_B1_std(y, X, Z, group, u_address, J, dims, coef, offset, familyR, link_int, init, XZ_group, K, penalty, penalty_params)
  }else if(fit_type == 4){
    coef_new = grp_CD_XZ_B2_std(y, X, Z, group, u_address, J, dims, coef, offset, familyR, link_int, init, XZ_group, K, penalty, penalty_params)
  }else{
    stop("fit_type", fit_type, "not available \n")
  }
  
  return(as.numeric(coef_new))
}


#################################################################################################
