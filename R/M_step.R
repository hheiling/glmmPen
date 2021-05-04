
#################################################################################################
# Coordinate Descent (ungrouped) for fixed effects only
#################################################################################################

CD = function(y, X, family, link_int, offset,
              coef_init, maxit_CD = 250, conv = 0.0001,
              penalty, lambda, gamma, alpha = 1.0, penalty_factor, trace = 0){
  
  p = ncol(X)
  N = length(y)
  
  # Coding of link_ink:
  # logit = 10, probit = 11, cloglog = 12
  # log = 20, sqrt = 21 identity = 30, inverse = 31
  
  penalty_params = c(lambda, gamma, alpha)
  
  dims = c(p, N, conv, maxit_CD)
  
  coef_new = glm_fit(y, X, dims, coef_init, offset, family, link_int, penalty, penalty_params, penalty_factor, trace)
  
  return(as.numeric(coef_new))
}


#################################################################################################
# Grouped Coordinate Descent for fixed and random effects
#################################################################################################

# Eventually implement the following:
## In fit_dat or glmmPen, check that group_X made of subsequent integers
## Change Z to be sparse (fit_dat and glmmPen)

# init: if this is the first attempt at a fit (using initial coef)
# family: character describing which family
# link_int: integer summarizing with link to use (see coding in fit_dat_B())
M_step = function(y, X, Z, u_address, M, J, group, family, link_int, coef, offset, phi,
                  maxit_CD = 250, conv_CD = 0.0001,
                  init, group_X = 0:(ncol(X)-1), covgroup,
                  penalty, lambda0, lambda1, gamma, alpha = 1.0, trace){
  
  if(!is.matrix(Z)){
    stop("Z must be a matrix \n")
  }else if(typeof(Z)=="integer") storage.mode(Z) <- "double"
  
  p = ncol(X)
  N = length(y)
  
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
  # Number of non-zero random effect variables in previous EM iteration (including random intercept)
  # q_non0 = sum(diag(sigma) > 0)
  
  dims = c(p, N, d, q, M, J_XZ, conv_CD, maxit_CD)
  
  # const arma::vec& y, const arma::mat& X, const arma::mat& Z,
  # const arma::vec& group, 
  # SEXP pBigMat, const arma::sp_mat& J_q, arma::vec dims,
  # arma::vec beta, const arma::vec& offset,
  # const char* family, int link, int init, double phi,
  # const arma::uvec& XZ_group, arma::uvec K, // covariate group index and size of covariate groups
  # const char* penalty, arma::vec params
  coef_new = grp_CD_XZ_fast(y, X, Z, group, u_address, J, dims, coef, offset, family, link_int, init, phi, XZ_group, K, penalty, penalty_params, trace)
  # coef_new = grp_CD_XZ(y, X, Z, group, u_address, J, dims, coef, offset, family, link_int, init, phi, XZ_group, K, penalty, penalty_params, trace)
  
  return(as.numeric(coef_new))
}


#################################################################################################
