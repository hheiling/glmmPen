
# Use small nMC, small maxitEM, and small maxit_CD
# Only want to get a general gist (not interested in accurate final coefficient estimates,
# just interested in which random effects are left in the model)
# Skimp on nMC_report and M because posterior modes and logLik will not be used from this model
prescreen_FA = function(dat, family, offset_fit, trace = 0, 
                         penalty, alpha, gamma_penalty, 
                         lambda0_min, lambda1_min, group_X,
                         optim_options, adapt_RW_options,
                         checks_complete, r, progress){
  
  
  # fixed effects penalty
  lam0 = lambda0_min
  # random effect penalty 
  lam1 = lambda1_min
  
  if((trace >= 1)) cat(sprintf("Pre-screening penalty parameters: fixed effects %f, random effects %f", lam0, lam1), "\n")
  
  # Overwrite optim control options with small nMC values
  q = ncol(dat$Z) / nlevels(dat$group)
  optim_options$nMC_burnin = 100
  optim_options$nMC = 100
  optim_options$nMC_max = 500
  
  # Overwrite optim control options with smaller maxitEM
  optim_options$maxitEM = 30
  
  # Fit 'full' model (small penalty for fixed and random effects)
  out = try(fit_dat_FA(dat, lambda0 = lam0, lambda1 = lam1, 
                      family = family, offset_fit = offset_fit, 
                      optim_options = optim_options, group_X = group_X,
                      penalty = penalty, alpha = alpha, gamma_penalty = gamma_penalty,
                      trace = trace,  
                      u_init = NULL, beta_old = NULL, b_old = NULL, ufull_describe = NULL,
                      adapt_RW_options = adapt_RW_options,
                      r=r, logLik_calc = FALSE, checks_complete = checks_complete,
                      progress = progress))
  
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
  
  return(list(ranef_keep = ranef_keep, coef_pre = out$coef, u_pre = out$u_init))
  
}