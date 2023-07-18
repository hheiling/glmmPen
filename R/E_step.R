
# E_step(): Function to perform E-step within the EM algorithm for a single model fit
# coef: coefficients from the M-step. Input all fixed and random effect coefficients,
#   only use fixed effects
# ranef_idx: index of which random effects are non-zero as of latest M-step (will not
#   sample from random effects that have been penalized out)
# y: response. If coxph, y = event indicator
# y_times: If coxph, value of observed times (NULL if not coxph family)
# X: fixed effects covariates
# Znew2: For each group k (k = 1,...,d), Znew2 = Z * Gamma (Gamma = t(chol(sigma)))
# group: group index
# offset_fit: offset for the linear predictor
# nMC: how many posterior samples to acquire
# nMC_burnin: how many posterior samples to discard as burn-in
# family, link: family and link for model
# phi, sig_g: additional parameters needed for negative binomial and gaussian families
#   Note: negative binomial family not available in current version of package
# sampler: type of sampler to use (Stan, adaptive random walk ...)
# d: total number groups
# uold: posterior sample to use for initialization of E-step
# proposal_SD, batch, batch_length, offset_increment: arguments for adaptive random walk
# coxph_options: for Cox Proportional Hazards model, additional parameters of interest


#' @importFrom bigmemory attach.big.matrix describe big.matrix
#' @importFrom rstan sampling extract
E_step = function(coef, ranef_idx, y, X, Znew2, group, offset_fit, 
                  nMC, nMC_burnin, family, link, phi = 0.0, sig_g = 1.0,
                  sampler, d, uold, proposal_SD, batch, batch_length,
                  offset_increment, trace){
  
  gibbs_accept_rate = matrix(NA, nrow = d, ncol = nrow(Znew2)/d)
  
  if(sampler == "independence"){
    # No burnin adjustments used
    samplemc_out = sample_mc2_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, 
                                     Z=Znew2, offset_fit=offset_fit,
                                     nMC=nMC, family = family, link = link, group = group, 
                                     d = d, uold = uold, phi = phi, sig_g = sig_g)
    u0 = attach.big.matrix(samplemc_out$u0)
    gibbs_accept_rate = samplemc_out$gibbs_accept_rate
    
    if(trace >= 2){
      print("Gibbs acceptance rate for E step:")
      print(gibbs_accept_rate)
    }
    
  }else if(sampler == "random_walk"){ 
    # First adapt the proposal variance during an initial burnin period
    samplemc_burnin = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, 
                                             Z=Znew2, offset_fit=offset_fit,
                                             nMC=nMC_burnin, family = family, link = link,
                                             group = group, d = d, uold = uold,
                                             proposal_SD = proposal_SD, batch = batch, batch_length = batch_length,
                                             offset = offset_increment, nMC_burnin = nMC_burnin,
                                             phi = phi, sig_g = sig_g)
    
    # Extract updated info from burnin period
    batch = samplemc_burnin$updated_batch
    proposal_SD = samplemc_burnin$proposal_SD
    
    u0 = attach.big.matrix(samplemc_burnin$u0)
    
    if(trace >= 2){
      print("Updated proposal_SD:")
      print(proposal_SD)
    }
    
    # Run the official E step with no adaptation
    uold = as.numeric(u0[nrow(u0),])
    samplemc_out = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, 
                                          Z=Znew2, offset_fit=offset_fit,
                                          nMC=nMC, family = family, link = link, 
                                          group = group, d = d, uold = uold,
                                          proposal_SD = proposal_SD, batch = batch, batch_length = batch_length,
                                          offset = offset_increment, nMC_burnin = 0,
                                          phi = phi, sig_g = sig_g)
    
    gibbs_accept_rate = samplemc_out$gibbs_accept_rate
    u0 = attach.big.matrix(samplemc_out$u0)
    
    if(trace >= 2){
      print("Gibbs acceptance rate for E step:")
      print(gibbs_accept_rate)
    }
    
  }else if(sampler == "stan"){
    
    u0 = big.matrix(nrow = nMC, ncol = ncol(Znew2), init=0) # use ', init = 0' for sampling within EM algorithm
    
    if(family == "binomial" & link == "logit"){ 
      stan_file = stanmodels$binomial_logit_model
    }else if(family == "poisson" & link == "log"){
      stan_file = stanmodels$poisson_log_model
    }else if(family == "gaussian" & link == "identity"){ 
      stan_file = stanmodels$gaussian_identity_model
    }else{
      stop("family and link combination not available")
    }
    
    # Number draws to extract in each chain (after burn-in)
    nMC_chain = nMC 
    
    # For each group, sample from posterior distribution (sample the alpha values)
    for(k in 1:d){
      idx_k = which(group == k)
      cols_k = seq(from = k, to = ncol(Znew2), by = d)
      y_k = y[idx_k]
      X_k = X[idx_k,]
      cols_use = cols_k[ranef_idx]
      if((length(cols_use) == 1) | (length(idx_k) == 1)){
        Z_k = matrix(Znew2[idx_k, cols_use], nrow = length(idx_k), ncol = length(cols_use))
      }else{ # length(cols_use) > 1
        Z_k = Znew2[idx_k, cols_use]
      }
      
      if(length(idx_k) == 1){
        dat_list = list(N = length(idx_k), # Number individuals in group k
                        q = length(ranef_idx), # number random effects (or common factors)
                        eta_fef = as.array(as.numeric(X_k %*% matrix(coef[1:ncol(X)], ncol = 1)) + offset_fit[idx_k]), # fixed effects componenet of linear predictor
                        y = as.array(y_k), # outcomes for group k
                        Z = Z_k) # portion of Z matrix corresponding to group k
      }else{ # length(idx_k) > 1
        dat_list = list(N = length(idx_k), # Number individuals in group k
                        q = length(ranef_idx), # number random effects
                        eta_fef = as.numeric(X_k %*% matrix(coef[1:ncol(X)], ncol = 1) + offset_fit[idx_k]), # fixed effects componenet of linear predictor
                        y = y_k, # outcomes for group k
                        Z = Z_k) # portion of Z matrix corresponding to group k
      }
      
      if(family == "gaussian"){
        
        dat_list$sigma = sig_g # standard deviation of normal dist of y
        
      }
      # else if(family == "negbin"){
      #   
      #   dat_list$theta = 1.0/phi # dispersion parameter, negbin var = (mu + mu^2 / theta) = (mu + mu^2 * phi)
      #   
      # }
      
      # initialize posterior random draws
      init_lst = list()
      if(length(cols_use) == 1){
        init_lst[[1]] = list(alpha = as.array(uold[cols_use]))
      }else{ # length(cols_use) > 1
        init_lst[[1]] = list(alpha = uold[cols_use])
      }
      
      # Sampling step
      # suppressWarnings(): Avoid excessive warnings when nMC_chain is low in early EM iterations
      stan_fit = suppressWarnings(rstan::sampling(stan_file, data = dat_list, init = init_lst, 
                                         iter = nMC_chain + nMC_burnin,
                                         warmup = nMC_burnin, show_messages = FALSE, refresh = 0,
                                         chains = 1, cores = 1))
      
      stan_out = as.matrix(stan_fit)
      
      # Exclude lp__ column of output (log density up to a constant)
      if(length(ranef_idx) > 1){
        draws_mat = stan_out[,1:length(ranef_idx)]
      }else{ # length(ranef_idx) == 1
        draws_mat = matrix(stan_out[,1], ncol = 1)
      }
      
      if(nrow(draws_mat) == nMC){
        u0[,cols_use] = draws_mat
      }else{ # nrow(draws_mat) > nMC due to ceiling function in 'iter' specification
        start_row = nrow(draws_mat) - nMC + 1
        rows_seq = start_row:nrow(draws_mat)
        u0[,cols_use] = draws_mat[rows_seq,]
      }
      
      
    } # End k for loop
    
  } # End if-else sampler
  
  return(list(u0 = describe(u0), proposal_SD = proposal_SD, gibbs_accept_rate = gibbs_accept_rate,
              updated_batch = batch))
  
}