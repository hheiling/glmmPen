
#' @importFrom bigmemory attach.big.matrix describe big.matrix
#' @importFrom rstan sampling extract
#' @export
E_step = function(coef, ranef_idx, y, X, Znew2, group, nMC, nMC_burnin, family, link, phi, sig_g,
                  sampler, d, uold, proposal_SD, batch, batch_length,
                  offset_increment, trace, max_cores){
  
  # Determine number of cores to use - for Stan sampling only
  if(nMC < 2000){
    num_cores = 1
  }else{
    num_cores = 2 + (nMC - 2000) %/% 1000
    if(num_cores > max_cores) num_cores = max_cores
  }
  
  gibbs_accept_rate = matrix(NA, nrow = d, ncol = nrow(Znew2)/d)
  
  if(sampler == "independence"){
    # No burnin adjustments used
    samplemc_out = sample_mc2_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, Z=Znew2, 
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
                                             Z=Znew2, nMC=nMC_burnin, family = family, link = link,
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
    
    # Run the official first E step with no adaptation
    uold = as.numeric(u0[nrow(u0),])
    samplemc_out = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, 
                                          Z=Znew2, nMC=nMC, family = family, link = link, 
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
    nMC_chain = ceiling(nMC/num_cores)
    
    # Record last elements of each chain for initialization of next E step
    last_draws = matrix(0, nrow = num_cores, ncol = ncol(Znew2)) 
    
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
                        q = length(ranef_idx), # number random effects
                        eta_fef = as.array(as.numeric(X_k %*% matrix(coef[1:ncol(X)], ncol = 1))), # fixed effects componenet of linear predictor
                        y = as.array(y_k), # outcomes for group k
                        Z = Z_k) # portion of Z matrix corresonding to group k
      }else{ # length(idx_k) > 1
        dat_list = list(N = length(idx_k), # Number individuals in group k
                        q = length(ranef_idx), # number random effects
                        eta_fef = as.numeric(X_k %*% matrix(coef[1:ncol(X)], ncol = 1)), # fixed effects componenet of linear predictor
                        y = y_k, # outcomes for group k
                        Z = Z_k) # portion of Z matrix corresonding to group k
      }
      
      
      if(family == "gaussian"){
        
        dat_list$sigma = sig_g # standard deviation of normal dist of y
        
      }else if(family == "negbin"){
        
        dat_list$theta = 1.0/phi # dispersion parameter, negbin var = (mu + mu^2 / theta) = (mu + mu^2 * phi)
        
      }
      
      # initialize posterior random draws
      init_lst = list()
      for(cr in 1:num_cores){
        if(length(cols_use) == 1){
          init_lst[[cr]] = list(alpha = as.array(uold[cols_use]))
        }else{ # length(cols_use) > 1
          init_lst[[cr]] = list(alpha = uold[cols_use])
        }
      }
      
      # Avoid excessive warnings when nMC_chain is low in early EM iterations
      stan_fit = suppressWarnings(rstan::sampling(stan_file, data = dat_list, init = init_lst, 
                                         iter = nMC_chain + nMC_burnin,
                                         warmup = nMC_burnin, show_messages = F, refresh = 0,
                                         chains = num_cores, cores = num_cores))
      
      stan_out = as.matrix(stan_fit)
      
      # Exclude lp__ column of output (log density up to a constant)
      if(length(ranef_idx) > 1){
        draws_mat = stan_out[,1:length(ranef_idx)]
      }else{ # length(ranef_idx) == 1
        draws_mat = matrix(stan_out[,1], ncol = 1)
      }
      
      
      # Find last elements for each chain for initialization of next E step
      for(m in 1:num_cores){
        last_draws[m,cols_use] = draws_mat[nMC_chain*m,]
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