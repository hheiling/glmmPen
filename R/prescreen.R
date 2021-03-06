
# Use small nMC, small maxitEM, and small maxit_CD
# Only want to get a general gist (not interested in accurate final coefficient estimates,
# just interested in which random effects are left in the model)
# Skimp on nMC_report and M because posterior modes and logLik will not be used from this model
prescreen = function(dat, family, offset_fit, trace = 0, 
                     penalty, alpha, gamma_penalty, 
                     lambda0_min, lambda1_min, group_X,
                     sampler, adapt_RW_options, covar,
                     var_start, max_cores, checks_complete){
  
  
  # fixed effects penalty
  lam0 = lambda0_min
  # random effect penalty 
  lam1 = lambda1_min
  
  print(sprintf("Pre-screening penalty parameters: fixed effects %f, random effects %f", lam0, lam1))
  
  # Determine nMC ranges
  q = ncol(dat$Z) / nlevels(dat$group)
  nMC_burnin = 100
  nMC = 100
  nMC_max = 500
  
  # Other convergence criteria
  conv_EM = 0.0015
  conv_CD = 0.0005
  maxitEM = 30
  maxit_CD = 50
  t = 2
  mcc = 2
  
  # Fit 'full' model (small penalty for fixed and random effects)
  out = try(fit_dat_B(dat, lambda0 = lam0, lambda1 = lam1, 
                      nMC_burnin = nMC_burnin, nMC = nMC, nMC_max = nMC_max,
                      family = family, offset_fit = offset_fit, group_X = group_X,
                      penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                      trace = trace, conv_EM = conv_EM, conv_CD = conv_CD,  
                      coef_old = NULL, u_init = NULL, ufull_describe = NULL,
                      maxitEM = maxitEM, maxit_CD = maxit_CD, t = t, mcc = mcc,
                      sampler = sampler, adapt_RW_options = adapt_RW_options,
                      covar = covar, var_start = var_start, logLik_calc = F,
                      max_cores = max_cores, checks_complete = checks_complete))
  
  if(is.character(out)){
    stop("Issue with pre-screening step in model selection procedure")
  }
  
  # Determine which random effects were kept in or penalized out of the full model
  sigma = out$sigma
  vars = diag(sigma)
  
  ranef_keep = numeric(length(vars))
  for(v in 1:length(vars)){
    if(v == 1){
      # Always keep random intercept
      ranef_keep[v] = 1
    }else{
      # Only keep random effect if variance for the random effect sufficiently large
      if(vars[v] >= 10^-2){
        ranef_keep[v] = 1
      }
    }
  }
  
  max_var = max(vars)
  
  return(list(ranef_keep = ranef_keep, coef_pre = out$coef, u_pre = out$u_init, max_var = max_var))
  
}