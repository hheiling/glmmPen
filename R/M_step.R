
#################################################################################################
# Coordinate Descent (ungrouped) for 'naive' model fit
#   Estimates fixed effects at start of EM algorithm in "fit_dat.R" assuming no random effects
#   Naive model fit used when not using a previous model to initialize coefficients (such as
#   in model selection procedure)
#################################################################################################

CD = function(y, X, family, link_int, offset,
              coef_init, maxit_CD = 250, conv = 0.0001,
              penalty, lambda, gamma, alpha = 1.0, penalty_factor, trace = 0){
  
  p = ncol(X)
  N = length(y)
  
  # Coding of link_int: see "family_export.R" for reasoning
  # logit = 10, log = 20, identity = 30
  
  penalty_params = c(lambda, gamma, alpha)
  
  dims = c(p, N, conv, maxit_CD)
  
  coef_new = pglm_fit(y, X, dims, coef_init, offset, family, link_int, penalty, penalty_params, penalty_factor, trace)
  
  return(as.numeric(coef_new))
}


#################################################################################################
# Grouped Coordinate Descent for M-step of fit_dat() function in "fit_dat.R"
#   Estimates both fixed and random effects
#################################################################################################

# init: indicator if this is the first M-step of the EM algorithm
# family: character describing which family
# link_int: integer summarizing with link to use (see coding in "family_export.R")
M_step = function(y, X, Z, u_address, M, J, group, family, link_int, coef, offset, 
                  sig_g, phi,
                  maxit_CD = 250, conv_CD = 0.0001,
                  init, group_X = 0:(ncol(X)-1), covgroup,
                  penalty, lambda0, lambda1, gamma, alpha = 1.0, 
                  step_size, trace){
  
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
  
  dims = c(p, N, d, q, M, J_XZ, conv_CD, maxit_CD)
  
  output = grp_CD_XZ_step(y, X, Z, group, u_address, J, dims, coef, offset, step_size, sig_g, family, link_int, init, phi, XZ_group, K, penalty, penalty_params, trace)
  output = as.numeric(output)
  
  coef_new = output[1:length(coef)]
  step_size = output[length(coef)+1]
  phi = output[length(coef)+2]
  
  return(list(coef_new = coef_new, step_size = step_size, phi = phi))
}

# coef_new = grp_CD_XZ_fast(y, X, Z, group, u_address, J, dims, coef, offset, family, link_int, init, phi, XZ_group, K, penalty, penalty_params, trace)


#################################################################################################
