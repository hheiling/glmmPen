
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
E_step = function(coef, ranef_idx, y, y_times = NULL, X, Znew2, group, 
                  offset_fit, offset_interval = NULL, interval_mat = NULL,
                  nMC, nMC_burnin, family, link, phi = 0.0, sig_g = 1.0,
                  sampler, d, uold, proposal_SD, batch, batch_length,
                  offset_increment, trace, coxph_options = NULL){
  
  if((family == "coxph") & !(sampler == "stan")){
    stop("'coxph' family currently only supports the 'stan' sampler")
  }
  if((family == "coxph") & (is.null(coxph_options))){
    stop("coxph_options must be of class 'coxphControl' (see coxphControl() documentation) for the 'coxph' family")
  }
  
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
    
  }else if((sampler == "stan") & (family != "coxph")){
    
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
    
  }else if((sampler == "stan") & (family == "coxph")){
    
    stop("coxph family not yet fully implemented")
    
    # Alternative approach: https://rpubs.com/kaz_yos/surv_stan_piecewise1
    
    # if(is.null(y_times) | is.null(offset_interval) | is.null(interval_mat) | is.null(coxph_options)){
    #   stop("y_times, offset_interval, interval_mat, and coxph_options cannot be input as NULL")
    # }
    # 
    # u0 = big.matrix(nrow = nMC, ncol = ncol(Znew2) + ncol(interval_mat), init=0) # use ', init = 0' for sampling within EM algorithm
    # 
    # stan_file = stanmodels$coxph_piecewise_exp_model
    # 
    # # Number draws to extract in each chain (after burn-in)
    # nMC_chain = nMC 
    # 
    # # If necessary, restrict columns of Znew2 matrix to columns associated with non-zero
    # #  latent variables (random effects / latent common factors)
    # # Also determine relevant rows of u0 matrix to save alpha samples
    # cols_analyze = NULL
    # for(k in 1:d){
    #   cols_k = seq(from = k, to = ncol(Znew2), by = d)
    #   cols_analyze = c(cols_analyze,cols_k[ranef_idx])
    # }
    # cols_analyze = cols_analyze[order(cols_analyze)]
    # 
    # # Sample the random effects / latent common factors 'alpha': group-specific values needed
    # # Also sample log-hazard values 'lhaz' for each time interval
    # ## As opposed to other families, sample all (d*q) random effects / (d*r) latent common factors 
    # ##    together instead of sampling by group. Reasoning: want log-hazard values to be 
    # ##    consistent regardless of group identity
    # dat_list = list(N = length(y), # Number of observations (note: subjects could be measured over several time points)
    #                 NT = ncol(interval_mat), # Number of time intervals
    #                 H = length(ranef_idx)*d, # Number groups times number latent variables (latent random effects or latent common factors)
    #                 eta_fef =  as.numeric(X %*% matrix(coef[1:ncol(X)], ncol = 1) + offset_fit), # Fixed effects portion of linear predictor
    #                 offset_interval = offset_interval, # log(t_ij) where t_ij = length of time subject survived during the interval 
    #                 y = y, # event indicator (1 = event, 0 = censor)
    #                 Z = Znew2[,cols_analyze], # Z * Gamma or Z * B matrix, see calculation in fit_dat_coxph
    #                 interval_mat = interval_mat, # Indicator of observed time interval for observation
    #                 lhaz_prior = coxph_options$lhaz_prior) # Specifies standard deviation of normal prior
    # 
    # # initialize posterior random draws
    # alpha_idx = cols_analyze
    # lhaz_idx = (ncol(Znew2)+1):length(uold)
    # # init: See "rstan::stan" documentation
    # ## Set initial values by providing a list equal in length to the number of chains (1).
    # ## The elements of this list should themselves be named lists, where each of these
    # ## named lists has the name of a parameter and is use to specify the initial values for
    # # that parameter for the corresponding chain
    # init_lst = list()
    # init_lst[[1]] = list(alpha = uold[alpha_idx],
    #                      lhaz = uold[lhaz_idx])
    #   
    # # Sampling step
    # # suppressWarnings(): Avoid excessive warnings when nMC_chain is low in early EM iterations
    # stan_fit = suppressWarnings(rstan::sampling(stan_file, data = dat_list, init = init_lst, 
    #                                             iter = nMC_chain + nMC_burnin,
    #                                             warmup = nMC_burnin, show_messages = FALSE, refresh = 0,
    #                                             chains = 1, cores = 1))
    # 
    # stan_out = as.matrix(stan_fit)
    # # Check organization of samples
    # # print(colnames(stan_out)) # first alpha samples, then lhaz samples, then lp__ value
    # # Exclude lp__ column of output (log density up to a constant)
    # samp_idx = 1:(length(cols_analyze) + ncol(interval_mat))
    # draws_mat = stan_out[,samp_idx]
    # # Specify column locations of u0 matrix to save samples from stan_fit object
    # u0_idx = c(cols_analyze, ((1:ncol(interval_mat))+ncol(Znew2)))
    # 
    # if(nrow(draws_mat) == nMC){
    #   u0[,u0_idx] = draws_mat
    # }else{ # nrow(draws_mat) > nMC due to ceiling function in 'iter' specification
    #   start_row = nrow(draws_mat) - nMC + 1
    #   rows_seq = start_row:nrow(draws_mat)
    #   u0[,u0_idx] = draws_mat[rows_seq,]
    # }
    
  } # End if-else sampler
  
  return(list(u0 = describe(u0), proposal_SD = proposal_SD, gibbs_accept_rate = gibbs_accept_rate,
              updated_batch = batch))
  
}


# # Indicator matrix: 
# ## For subject i, determine which columns of the Znew2 matrix are relevant for analyses
# ## In other words, if subject i in group k, indicate which rows of Znew2 matrix associated with group k
# I_mat = matrix(0, nrow = nrow(Znew2), ncol = ncol(Znew2))
# for(k in 1:d){
#   idx_k = which(group == k)
#   cols_k = seq(from = k, to = ncol(Znew2), by = d)
#   I_mat[idx_k,cols_k] = 1
# }