
# Functions for Harmonic Mean and Modified Harmonic Mean calculations - R only code

logLik_MCMC_R = function(y, X, Z, group, U_old, sigma, coef, family, M, gibbs, 
                         fit, ok_index, Znew2){
  # Useful variables
  d = nlevels(group)
  
  # Get M MCMC samples
  U_list = sample.mc2(fit = fit, cov = sigma, y = y, X = X,
                      Z = Znew2, nMC = M, family = family, 
                      group = group, d = d, okindex = ok_index,
                      nZ = Znew2, gibbs = gibbs, uold = U_old)
  
  U_full = U_list$u0
  
  # Ignore components of U, Z, and sigma corresponding to random effects penalized to 0 variance
  
  # Non-zero variance random effects variables:
  non_zero = (diag(sigma) != 0)
  
  if(sum(non_zero) > 0){
    
    non_zero_ext = rep(non_zero, each = d)
    
    # Only keep non-zero columns of U and Z
    U = U_full[,non_zero_ext]
    U_means = colMeans(U)
    Z_red = Z[,non_zero_ext]
    
    # Reduced sigma: remove rows and columns with diag = 0
    sigma_red = sigma[non_zero,non_zero]
    
    # Gamma = cholesky decomposition of sigma (lower-triangular)
    Gamma = t(chol(sigma_red))
    
    # Calculated fixed effects contribution to eta (linear predictor)
    eta_fef = X %*% coef[1:ncol(X)]
    
    # Chi-square cutoff for truncated normal
    Chi_cutoff = qchisq(0.98, df = sum(non_zero), ncp = 0, lower.tail = T, log.p = F)
    
    ll = list()
    
    ll[[1]] = logLik_mHM_R(U, U_means, sigma_red, M, d, Chi_cutoff, group, y, eta_fef, Z_red,
                           Gamma, family)
    
    ll[[2]] = logLik_HM_R(U, M, d, group, y, eta_fef, Z_red, Gamma, family)
    
  }else{
    
    eta_fef = X %*% coef[1:ncol(X)]
    
    if(family == "binomial"){
      ll = sum(dbinom(y, size = 1, prob = exp(eta_fef) / (1+exp(eta_fef)), log = T))
    }else if(family == "poisson"){
      ll = sum(dpois(y, lambda = exp(eta_fef), log = T))
    }
  }
  
  return(ll)
}

logLik_mHM_R = function(U, U_means, sigma, M, d, Chi_cutoff,
                        group, y, eta_fef, Z, Gamma, family){
  num_vars = ncol(Z) / d
  
  for(k in 1:d){
    # Identify columns corresponding to group k (all vars)
    cols = seq(from = k, by = d, length.out = num_vars)
    # Select appropriate columns of U and elements of U_means (posterior modes)
    samps = U[,cols]
    mv_means = U_means[cols]
    # Identify individuals belonging to group k
    ids = (group == k)
    # Calculate truncated norm weights for samples
    wt_norm = dmvnorm(samps, mean = mv_means, sigma = sigma)
    region = apply(samps, 1, function(x) t(x - mv_means) %*% solve(sigma) %*% (x - mv_means))
    trunc_norm = ifelse(region > Chi_cutoff, 0, wt_norm)
    print("proportion above cut-off:")
    print(sum(region > Chi_cutoff)/M)
    # Calculate prior density values for samples
    prior = dmvnorm(samps)
    
    # Find elements of y and rows of Z for individuals belonging to group k
    Z_k = Z[ids,cols]
    y_k = y[ids]
    
    # Calculate linear preditor (eta) - sum fixed and random effect components
    eta_ref = samps %*% Gamma %*% t(Z_k)
    eta_fef_k = matrix(eta_fef[ids], nrow = 1)
    eta = eta_fef_k[rep(1, times = M),] + eta_ref
    
    # Calculated log-likelihood
    ll = 0
    l_i = numeric(M)
    for(i in 1:length(y_k)){
      # f*(alpha_m) / (f(y|alpha_m) * prior(alpha_m))
      l_i = trunc_norm / 
        (dbinom(y_k[i], size = 1, prob = exp(eta[,i]) / (1+exp(eta[,i])), log = F) * prior)
      # Sum all f*(alpha_m) / (f(y|alpha_m) * prior(alpha_m)), divide by M
      inv_l_i = sum(l_i) / M
      # Take inverse to find likelihood f(y), then take log
      ll_i = log(1 / inv_l_i)
      # Add individual contribution to overall log-lik
      ll = ll + ll_i
    }
    
  }
  
  return(ll)
  
}

logLik_HM_R = function(U, M, d, group, y, eta_fef, Z, Gamma, family){
  
  num_vars = ncol(Z) / d
  
  for(k in 1:d){
    # Identify columns corresponding to group k (all vars)
    cols = seq(from = k, by = d, length.out = num_vars)
    # Select appropriate columns of U 
    samps = U[,cols]
    # Identify individuals belonging to group k
    ids = (group == k)
    
    # Find elements of y and rows of Z for individuals belonging to group k
    Z_k = Z[ids,cols]
    y_k = y[ids]
    
    # Calculate linear preditor (eta) - sum fixed and random effect components
    eta_ref = samps %*% Gamma %*% t(Z_k)
    eta_fef_k = matrix(eta_fef[ids], nrow = 1)
    eta = eta_fef_k[rep(1, times = M),] + eta_ref
    
    # Calculated log-likelihood
    ll = 0
    l_i = numeric(M)
    for(i in 1:length(y_k)){
      # Calculate 1 / f(y|alpha_m)
      l_i = 1 / dbinom(y_k[i], size = 1, prob = exp(eta[,i]) / (1+exp(eta[,i])), log = F)
      # Sum all (1 / f(y|alpha_m)) and divide by M
      inv_l_i = sum(l_i) / M
      # Take inverse to calculate f(y), then take log for log-lik
      ll_i = log(1 / inv_l_i)
      # Add individual components to overall log-lik
      ll = ll + ll_i
    }
    
  }
  
  return(ll)
  
} 

# Functions for Harmonic Mean and Modified Harmonci Mean calculations - R + Rcpp

logLik_MCMC = function(y, X, Z, group, U_old, sigma, coef, family, M, gibbs, 
                       fit, ok_index, Znew2){
  # Useful variables
  d = nlevels(group)
  
  # Get M MCMC samples
  U_list = sample.mc2(fit = fit, cov = sigma, y = y, X = X,
                      Z = Znew2, nMC = M, family = family, 
                      group = group, d = d, okindex = ok_index,
                      nZ = Znew2, gibbs = gibbs, uold = U_old)
  
  U_full = U_list$u0
  
  # Ignore components of U, Z, and sigma corresponding to random effects penalized to 0 variance
  
  # Non-zero variance random effects variables:
  non_zero = (diag(sigma) != 0)
  
  if(sum(non_zero) > 0){
    
    non_zero_ext = rep(non_zero, each = d)
    
    # Only keep non-zero columns of U and Z
    U = U_full[,non_zero_ext]
    U_means = colMeans(U)
    Z_red = Z[,non_zero_ext]
    
    # Reduced sigma: remove rows and columns with diag = 0
    sigma_red = sigma[non_zero,non_zero]
    
    # Remove columns of Z 
    
    # Gamma = cholesky decomposition of sigma (lower-triangular)
    Gamma = t(chol(sigma_red))
    
    # Calculated fixed effects contribution to eta (linear predictor)
    eta_fef = X %*% coef[1:ncol(X)]
    
    # Chi-square cutoff for truncated normal
    Chi_cutoff = qchisq(0.99, df = sum(non_zero), ncp = 0, lower.tail = T, log.p = F)
    
    ll = list()
    
    # See logLik_mHM in logLik_MCMC.cpp
    ll[[1]] = logLik_mHM(U, U_means, sigma_red, M, d, Chi_cutoff, group, y, eta_fef, Z_red,
                         Gamma, family)
    # See logLik_HM in logLik_MCMC.cpp
    ll[[2]] = logLik_HM(U, M, d, group, y, eta_fef, Z_red, Gamma, family)
    
  }else{
    
    eta_fef = X %*% coef[1:ncol(X)]
    
    if(family == "binomial"){
      ll = sum(dbinom(y, size = 1, prob = exp(eta_fef) / (1+exp(eta_fef)), log = T))
    }else if(family == "poisson"){
      ll = sum(dpois(y, lambda = exp(eta_fef), log = T))
    }
  }
  
  return(ll)
}