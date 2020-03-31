
# Pajor Methods for logLik calcuations

indicator = function(bounds, samples){
  inbounds = function(x){mean(x > bounds[1,] & x < bounds[2,]) == 1}
  apply(samples, 1, inbounds)
}

indicator_1 = function(bounds, samples){
  inbounds = function(x){mean(x > bounds[1,] & x < bounds[2,]) == 1}
  apply(samples, 1, inbounds)
}

indicator_2 = function(bounds, samples){
  inbounds = function(x){mean(x > bounds[1] & x < bounds[2]) == 1}
  apply(samples, 1, inbounds)
}

# Pajor IS CAME Method for logLik calcuations
#' @export
CAME_IS = function(posterior, y, X, Z, group, coef, sigma, family, M){
  
  # Define variables
  d = nlevels(group)
  num_var = ncol(Z) / d
  Gamma = t(chol(sigma[which(diag(sigma)!=0),which(diag(sigma)!=0)]))
  eta_fef = X %*% coef[1:ncol(X)]
  post_mean = colMeans(posterior)
  
  # Calculate bounds of the posterior for each dimension
  bounds = apply(posterior, 2, range)
  
  # Thin posterior draws based on acf results
  grp_names = as.numeric(levels(group))
  d = nlevels(group)
  var_names = 1:ncol(X) - 1  # 0 = intercept, 1 to p for p covariates
  var_num = length(var_names)
  grp_index = rep(grp_names, times = var_num)
  var_index = rep(var_names, each = d)
  
  for(j in 1:ncol(posterior)){
    ACF = acf(posterior[,j], plot=F, lag.max = 30)
    ACF_df = with(ACF, data.frame(lag,acf))
    ACF_df$grp_names = grp_index[j]
    ACF_df$var_names = var_index[j]
    if(j == 1){
      ACF_all = ACF_df
    }else{
      ACF_all = rbind(ACF_all, ACF_df)
    }
  }
  
  # Want lag where acf < 0.1
  df_lag = subset(ACF_all, acf > 0.01)
  max_lag = max(max(df_lag$lag) + 1, 10) # At minimum, thin by 10
  thin = seq(from = 1, to = nrow(posterior), by = max_lag)
  post_thin = posterior[thin,]
  
  # Initiate logLik
  ll = 0
  
  for(k in 1:d){
    # Define group-specific individuals
    ids = (group == k)
    y_k = y[ids]
    n_k = sum(ids)
    
    # Identify columns for variables corresponding to group k
    cols_all = seq(from = k, by = d, length.out = num_var)
    # Restrict above colums to those belonging to variables with non-zero random effect variance
    cols = cols_all[which(diag(sigma) != 0)]
    Z_k = Z[ids,cols]
    if(length(cols) == 1){
      Z_k = matrix(Z_k, ncol = 1)
    }
    
    # Sample M samples from the importance function
    ## Importance function: multivariate normal (multivariate t?)
    if(length(cols) > 1){
      post_cov = var(post_thin[,cols])
    }else{
      post_cov = matrix(var(post_thin[,cols]), nrow = 1, ncol = 1)
    }
    
    imp_samp = rmvnorm(M, mean = post_mean[cols], sigma = post_cov)
    # imp_samp = rmvnorm(M, mean = post_mean[cols], sigma = sigma)
    # imp_samp = rmvt(M, df = 10, delta = post_mean[cols], sigma = sigma)
    
    # Determine Importance Weights - f/g
    prior = dmvnorm(imp_samp) # Mean 0, sigma I
    imp_dens = dmvnorm(imp_samp, post_mean[cols], sigma = post_cov)
    # imp_dens = dmvt(imp_samp, delta = post_mean[cols], sigma = sigma, df = 10, log = F)
    wt = prior / imp_dens
    
    # Calculate density given importance samples
    
    ## Fixed-effects component of eta (linear predictor)
    eta_fef_k = matrix(eta_fef[ids], nrow = 1) # 1 X n_k
    ## Random-effects component of eta (linear predictor)
    eta_ref = imp_samp %*% Gamma %*% t(Z_k) # M x n_k
    ## Full eta - fixed + random effects component
    eta = eta_fef_k[rep(1, M),] + eta_ref
    
    # Calculate indicator function
    if(length(cols) > 1){
      indic = indicator_1(bounds[,cols], imp_samp)
    }else{
      indic = indicator_2(bounds[,cols], imp_samp)
    }
    cat("proportion of samples inside bounds:", mean(indic), "\n")
    
    dens_k = rep(1, times = M)
    
    if(family == "binomial"){
      prob_mat = exp(eta) / (1+exp(eta))
      # Columns of prob_mat correspond to individuals
      
      for(i in 1:n_k){
        dens_i = dbinom(y_k[i], size = 1, prob = prob_mat[,i], log = F)
        # Multiply together all elements belonging to same alpha_m
        dens_k = dens_k * dens_i
      } # End i loop
    }else{
      stop("specified family not currently available")
    }
    
    lik_k = mean(dens_k*indic*wt)
    ll_k = log(lik_k)
    ll = ll + ll_k
    
  } # End k loop
  
  return(ll)
} # End CAME_IS function



# Pajor regular CAME Method for logLik calculations

# Sample from prior, no importance sampling
#' @export
CAME = function(posterior, y, X, Z, group, coef, sigma, family, M){
  
  # Define variables
  d = nlevels(group)
  num_var = ncol(Z) / d
  Gamma = t(chol(sigma[which(diag(sigma)!=0),which(diag(sigma)!=0)]))
  eta_fef = X %*% coef[1:ncol(X)]
  # post_mean = colMeans(posterior)
  # post_cov = sigma
  
  # Calculate bounds of the posterior for each dimension
  bounds = apply(posterior, 2, range)
  
  ll = 0
  
  for(k in 1:d){
    # Define group-specific individuals
    ids = (group == k)
    y_k = y[ids]
    n_k = sum(ids)
    
    # Identify columns corresponding to group k
    cols_all = seq(from = k, by = d, length.out = num_var)
    # Restrict above colums to those belonging to variables with non-zero random effect variance
    cols = cols_all[which(diag(sigma) != 0)]
    Z_k = Z[ids,cols]
    
    # Sample M samples from the prior function
    zero = rep(0, num_var)
    I = diag(x = 1, nrow = num_var)
    prior = rmvnorm(M, zero, I) # Mean 0, sigma I
    
    # Calculate density given importance samples
    
    ## Fixed-effects component of eta (linear predictor)
    eta_fef_k = matrix(eta_fef[ids], nrow = 1)
    ## Random-effects component of eta (linear predictor)
    eta_ref = prior %*% Gamma %*% t(Z_k)
    ## Full eta - fixed + random effects component
    eta = eta_fef_k[rep(1, M),] + eta_ref
    
    # Calculate indicator function
    indic = indicator(bounds[,cols], prior)
    cat("proportion of samples inside bounds:", mean(indic), "\n")
    
    dens_k = rep(1, times = M)
    
    if(family == "binomial"){
      prob_mat = exp(eta) / (1+exp(eta))
      # Columns of prob_mat correspond to individuals
      
      for(i in 1:n_k){
        dens_i = dbinom(y_k[i], size = 1, prob = prob_mat[,i], log = F)
        # Multiply together all elements belonging to same alpha_m
        dens_k = dens_k * dens_i
      } # End i loop
    }else{
      stop("specified family not currently available")
    }
    
    lik_k = mean(dens_k*indic)
    ll_k = log(lik_k)
    ll = ll + ll_k
    
  } # End k loop
  
  return(ll)
} # End CAME function