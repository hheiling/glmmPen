

#################################################################################################
# Grouped Coordinate Descent for fixed and random effects
# Version: random effects decomposed using factor analysis
#################################################################################################

# init: Indicator if this is the first iteration of the EM algorithm
#   If this is the first M-step of the EM algorithm, allow fixed effect coefficients
#   set to 0 in initialized fixed effects to be properly estimated
# family: character describing which family
# link_int: integer summarizing with link to use (see coding in "family_export.R")
M_step_FA = function(y, X, Z, u_address, M, J, group, family, link_int, coef, offset,
                     sig_g, phi,
                     maxit_CD = 250, conv_CD = 0.0001,
                     init, group_X = 0:(ncol(X)-1),
                     penalty, lambda0, lambda1, gamma, alpha = 1.0, 
                     step_size, trace){
  
  if(!is.matrix(Z)){
    stop("Z must be a matrix \n")
  }else if(typeof(Z)=="integer") storage.mode(Z) <- "double"
  
  p = ncol(X)
  N = length(y)
  
  penalty_params = c(lambda0, lambda1, gamma, alpha)
  
  # Number of groups in the fixed effects covariates
  J_X = length(unique(group_X))
  if((max(group_X) - min(group_X) + 1) != J_X){ # Assume min(group_X) == 0 (for intercept)
    stop("group_X must be comprised of consecutive integers \n")
  }
  # Size of covariate groups
  K = numeric(J_X)
  for(j in unique(group_X)){
    idx = which(group_X == j)
    K[j+1] = length(idx) # Add 1 because smallest group_X value is 0
  }
  
  # Number of groups wrt observations
  d = nlevels(group)
  # Number of random effect variables
  q = ncol(Z) / d
  
  dims = c(p, N, d, q, M, J_X, conv_CD, maxit_CD)
  
  # const arma::vec& y, const arma::mat& X, const arma::mat& Z,
  # const arma::vec& group, 
  # SEXP pBigMat, const arma::sp_mat& J_f, arma::vec dims,
  # arma::vec beta, const arma::vec& offset,
  # double step_size, double sig_g,
  # const char* family, int link, int init, double phi,
  # const arma::uvec& X_group, arma::uvec K, // covariate group index and size of covariate groups
  # const char* penalty, arma::vec params, int trace
  
  output = grp_CD_XZ_FA_step(y, X, Z, group, u_address, J, dims, coef, offset, step_size, sig_g, family, link_int, init, phi, group_X, K, penalty, penalty_params, trace)
  output = as.numeric(output)
  
  coef_new = output[1:length(coef)]
  step_size = output[length(coef)+1]
  phi = output[length(coef)+2]
  
  return(list(coef_new = coef_new, step_size = step_size, phi = phi))
}

# coef_new = grp_CD_XZ_FA(y, X, Z, group, u_address, J, dims, coef, offset, family, link_int, init, phi, group_X, K, penalty, penalty_params, trace)



#################################################################################################
