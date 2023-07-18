

# Functions to estimate number of common factors, r
## estimate_r: Overall function that performs all steps
## pseudo_ranef: Calculate pseudo random effects estimates to use within 
##    r estimation procedure
## rest_ratio: Estimation methods for Growth Ratio (GR, default) and 
##    Eigenvalue Ratio (ER)
## rest_bai_ng: Estimation method of and Bai and Ng (2002) 
##    (BN1, BN2 depending on penalty used)

############################################################################################
# Put all steps together to estimate r
############################################################################################
#' @importFrom stringr str_detect
estimate_r = function(dat, optim_options = optimControl(), coef_names,
                      r_est_method = "GR", r_max, family_info, offset_fit, 
                      penalty, lambda0, gamma_penalty, alpha, group_X, 
                      sample, size, data_type, trace){
  
  # Extract relevant data
  y = dat$y
  X = dat$X 
  group = dat$group
  d = nlevels(group)
  ## Input covariates: Only those specified for random effects
  col_idx = which(coef_names$fixed %in% coef_names$random)
  if(data_type != "survival"){
    X_use = X[,col_idx, drop = FALSE]
  }else if(data_type == "survival"){
    col_idx_surv = c(2:length(dat$cut_points),col_idx)
    col_idx = col_idx_surv[order(col_idx_surv)]
    X_use = X[,col_idx, drop=FALSE]
  }
  
  p = ncol(X_use)
  
  # Extract optimization hyper parameters
  maxit_CD = optim_options$maxit_CD
  conv_CD = optim_options$conv_CD
  
  # penalty_factor
  penalty_factor = numeric(ncol(X))
  penalty_factor[which(group_X == 0)] = 0
  penalty_factor[which(group_X != 0)] = 1
  penalty_factor = penalty_factor[col_idx]
  
  if(is.null(r_max)){
    stop("algorithm not yet set up to recommend an r_max")
    # Bai and Ng (2002) recommendation for r_max:
    # r_max = 8 * round((min(d,p) / 100)^(1/4))
  }else{
    if(r_max > p){
      r_max = p # change to p/2 ?
    }
  }
  
  # Calculate pseudo random effects estimates
  G = pseudo_ranef(y = y, X = X_use, group = group, d = d,
                   family_info = family_info, offset_fit = offset_fit,
                   maxit_CD = maxit_CD, conv_CD = conv_CD,
                   penalty = penalty, lambda0 = lambda0,
                   gamma_penalty = gamma_penalty, alpha = alpha,
                   penalty_factor = penalty_factor, 
                   sample = sample, size = size, trace = trace)
  
  if((r_est_method == "GR") | (r_est_method == "ER")){
    r_est = rest_ratio(G = G, r_max = r_max, type = r_est_method)
  }else if(r_est_method == "BN1"){
    r_est = rest_bai_ng(G = G, r_max = r_max, penalty = 1)
  }else if(r_est_method == "BN2"){
    r_est = rest_bai_ng(G = G, r_max = r_max, penalty = 2)
  }
  
  return(list(r_est = r_est, r_max = r_max))
  
}

######################################################################################
# Calculate pseudo random effects estimates by
# calculating group-specific fixed effects estimates
######################################################################################

# y = numeric vector of outcomes
# X = matrix of fixed effects predictors, including intercept
# group = vector of groups (numeric values, labeled 1:d)
# d = integer, number of groups in model
#' @importFrom mvtnorm rmvnorm
pseudo_ranef = function(y, X, group, d,
                        family_info, offset_fit,
                        maxit_CD, conv_CD, penalty, lambda0, gamma_penalty,
                        alpha, penalty_factor,
                        sample, size, trace){
  
  # Matrix of group-specific fixed effects estimates
  G = matrix(0, nrow = d, ncol = ncol(X))
  # Initialization of fixed effects: intercept-only model
  IntOnly = glm(y ~ 1, family = family_info$family_fun, offset = offset_fit)
  coef_init = c(IntOnly$coefficients, rep(0, length = (ncol(X)-1)))
  
  for(k in 1:d){
    idx_k = which(group == k)
    y_k = y[idx_k]
    X_k = X[idx_k,]
    offset_k = offset_fit[idx_k]
    beta_k = CD(y_k, X_k, family = family_info$family, link_int = family_info$link_int, 
               offset = offset_k, coef_init = coef_init,
               maxit_CD = maxit_CD, conv = conv_CD, penalty = penalty, lambda = lambda0,
               gamma = gamma_penalty, alpha = alpha, penalty_factor = penalty_factor, 
               trace = trace)
    # Checks of NA output results. 
    # If NA, set to intercept only model 
    # Note: NA values sometimes occur when the fit for a single group is bad/divergent
    ## (Occurred in some simulation settings, possible for this to occur in real data
    ## when number of observations per group is small and number of covariates is large)
    if(!any(is.na(beta_k))){
      G[k,] = beta_k
    }else{
      G[k,] = rep(0, ncol(G))
    }
  }
  
  # Subtract mean of the group-specific fixed effects estimates in order to 
  # acquire pseudo estimates of the random effects to use for r estimation
  beta_mean = colMeans(G)
  for(k in 1:d){
    G[k,] = G[k,] - beta_mean
  }
  
  # If number of groups is small, sample additional pseudo random effect estimates
  # from the distribution of available pseudo random effects estimates
  # Total pseudo random effects to use: 25
  d_tot = size
  if(sample == TRUE){
    G_old = G
    G = matrix(0, nrow = d_tot, ncol = ncol(G_old))
    G[1:d,] = G_old
    
    cov_G = var(G_old)
    mu = colMeans(G_old)
    G[(d+1):d_tot,] = rmvnorm(n = (d_tot - d), 
                              mean = mu,
                              sigma = cov_G)
  }
  
  return(G)
  
}

############################################################################################
# Eigenvalue Ratio (ER) and Growth Ratio estimation procedure
# Paper: Ahn and Horenstein (2013)
############################################################################################
# G = matrix of group-specific psuedo random effects estimates
# r_max = maximum r to consider
# type: either eigenvalue ratio (ER) or growth ratio (GR)
rest_ratio = function(G, r_max = NULL, type = c("ER","GR")){
  
  # d = number of groups = nrow(G)
  d = nrow(G)
  # p = number of total predictors plus intercept
  p = ncol(G)
  
  # Calculate equivalent to X X^T / (NT) in Ahn and Horenstein paper
  R = G %*% t(G) / (d * p)
  
  # Perform SVD on above R matrix (svd from base R)
  R_svd = svd(R)
  # Extract eigenvalues, order from largest to smallest, remove zero-valued eigenvalues (?)
  eigen_vals = R_svd$d
  eigen_vals = eigen_vals[order(eigen_vals, decreasing = TRUE)]
  # eigen_vals = eigen_vals[which(eigen_vals > 10^-4)]
  # Maximum considered r (number of common factors)
  if(is.null(r_max)){
    V0 = mean(eigen_vals)
    r_max = sum(eigen_vals > V0)
  }
  # r_max must be less than total number of eigenvalues 
  if(r_max >= length(eigen_vals)){ 
    r_max = length(eigen_vals) - 1
  }
  # Perform Eigenvalue Ratio or Growth Ratio estimation procedure
  if(type == "ER"){
    # Calculate eigenvalue ratios
    ratios = numeric(r_max)
    for(v in 1:r_max){
      ratios[v] = eigen_vals[v] / eigen_vals[v+1]
    }
    r_est = which.max(ratios)
  }else if(type == "GR"){
    gr = numeric(r_max)
    for(v in 1:r_max){
      gr[v] = log(1 + eigen_vals[v]/sum(eigen_vals[-c(1:v)])) / log(1 + eigen_vals[v+1]/sum(eigen_vals[-c(1:(v+1))]))
    }
    r_est = which.max(gr)
  }
  # Minimum: r = 2
  if(r_est == 1){
    r_est = 2
  }
  
  return(r_est)
}

############################################################################################
# Bai and Ng (2002) estimation procedure
############################################################################################
# G = matrix of psuedo group-specific random effects estimates
# r_max = maximum r to consider
# penalty = type of penalty to use in procedure
rest_bai_ng = function(G, r_max, penalty = 1){
  
  # d = number of groups = nrow(G)
  d = nrow(G)
  # p = number of total predictors plus intercept
  p = ncol(G)
  
  # R = equivalent of Y^T * Y matrix in Bai and Ng, where Y = p x d matrix
  R = G %*% t(G)
  # Perform SVD on R matrix
  R_svd = svd(R)
  
  # if(is.null(r_max)){
  #   r_max = r_true * 2
  # }
  
  BN_vec = numeric(r_max)
  if(penalty == 1){
    g_pen = (d + p) / (d*p) * log(d*p/(d+p))
  }else if(penalty == 2){
    g_pen = (d + p) / (d*p) * log(min(d,p))
  }
  for(i in 1:r_max){
    F_r = R_svd$u[,1:i] * sqrt(d)
    X = t(G) - 1/d * t(G) %*% F_r %*% t(F_r)
    XtX = t(X) %*% X
    BN_vec[i] = log(1/(d*p) * sum(diag(XtX))) + i*g_pen
  }
  
  r_est = which.min(BN_vec)
  # Minimum: r = 2
  if(r_est == 1){
    r_est = 2
  }
  
  return(r_est)
  
}





############################################################################################