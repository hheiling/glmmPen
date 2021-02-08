
# Use small nMC, small maxitEM, and small maxit_CD
# Only want to get a general gist (not interested in accurate final coefficient estimates,
# just interested in which random effects are left in the model)
# Skimp on nMC_report and M because posterior modes and logLik will not be used from this model
prescreen = function(dat, family, offset_fit, trace = 0, 
                     penalty, alpha, gamma_penalty, group_X,
                     sampler, adapt_RW_options, covar,
                     var_start, max_cores, checks_complete){
  
  
  lam_MaxMin = LambdaRange(dat$X[,-1], dat$y, family = family, nlambda = 2)
  lam_min = lam_MaxMin[2]
  lam0 = lam_min
  lam1 = lam_min
  
  # Determine nMC ranges (depends on number random effects)
  q = ncol(dat$Z) / nlevels(dat$group)
  if(q <= 26){ # 25 random effect predictors + random intercept
    nMC_burnin = 50
    nMC = 100
    nMC_max = 500
  }else{
    nMC_burnin = 25
    nMC = 50
    nMC_max = 350
  }
  
  # Other convergence criteria
  conv_EM = 0.001
  conv_CD = 0.0001
  maxitEM = 25
  maxit_CD = 100
  t = 2
  mcc = 2
  # Only interested in what random effects are still present in the model
  ## Not otherwise interested in accurate logLik or posterior modes (or accurate coefficient estimates)
  nMC_report = nMC_max
  M = nMC_max
  
  
  # Fit 'full' model
  out = try(fit_dat_B(dat, lambda0 = lam0, lambda1 = lam1, 
                      nMC_burnin = nMC_burnin, nMC = nMC, nMC_max = nMC_max, nMC_report = nMC_report,
                      family = family, offset_fit = offset_fit, group_X = group_X,
                      penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                      trace = trace, conv_EM = conv_EM, conv_CD = conv_CD,  
                      coef_old = NULL, u_init = NULL, ufull_describe = NULL,
                      maxitEM = maxitEM, maxit_CD = maxit_CD, t = t, mcc = mcc,
                      M = M, sampler = sampler, adapt_RW_options = adapt_RW_options,
                      covar = covar, var_start = var_start,
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
      if(vars[v] >= 10^-4){
        ranef_keep[v] = 1
      }
    }
  }
  # ranef_keep[which(diag(sigma) > 0)] = 1
  
  return(list(ranef_keep = ranef_keep, coef_pre = out$coef, u_pre = out$u_init))
  
}