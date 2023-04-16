
#' @importFrom bigmemory big.matrix attach.big.matrix
E_step_final_FA = function(dat, offset_fit, fit, optim_options, 
                           fam_fun, extra_calc, adapt_RW_options, 
                           trace, progress, r){
  
  # Extract terms from dat list object
  y = dat$y
  X = dat$X # standardized X
  Z = Matrix::as.matrix(dat$Z)
  group = dat$group
  d = nlevels(group)
  q = ncol(Z) / d
  
  if(is.null(offset_fit)){
    offset_fit = rep(0, times = length(y))
  }
  
  # Extract family and link information
  family = fam_fun$family
  link = fam_fun$link
  
  # Extract terms from fitted model object
  sigma = fit$sigma
  J = fit$J
  coef = fit$coef
  sig_g = ifelse(family == "gaussian", fit$sigma_gaus, 1.0)
  
  # Extract control values from optim_options
  M = 10^4
  nMC_burnin = optim_options$nMC_burnin
  nMC_report = optim_options$nMC_report
  if(M < nMC_report){
    nMC_report = M
  }
  sampler = optim_options$sampler
  
  
  # For each group k, Znew2_k = Z_k * B
  Znew2 = matrix(0, nrow = nrow(Z), ncol = d*r)
  B = t(matrix(coef[-c(1:ncol(X))], nrow = r, ncol = q))
  for(j in 1:d){
    Znew2[group == j,seq(j, ncol(Znew2), by = d)] = Z[group == j, seq(j, ncol(Z), by = d)] %*% B
  }
  
  if(progress == TRUE) message("Start of sampling from posterior")
  Estep_out = E_step(coef = coef, ranef_idx = 1:r, y = y,
                     X = X, Znew2 = Znew2, group = group, offset_fit = offset_fit,
                     nMC = M, nMC_burnin = nMC_burnin,
                     family = family, link = link, phi = 1.0, sig_g = sig_g,
                     sampler = sampler, d = d, uold = as.numeric(fit$u_init),
                     proposal_SD = fit$proposal_SD, batch = fit$updated_batch, 
                     batch_length = adapt_RW_options$batch_length, 
                     offset_increment = adapt_RW_options$offset, trace = trace)
  if(progress == TRUE) message("Finished sampling from posterior")
  
  if(extra_calc){
    
    u0 = attach.big.matrix(Estep_out$u0)
    
    # Calculate posterior modes
    ## Final output post_U: posterior modes for all possible q random effects
    ## Input u0: posterior samples for the r common factors
    post_U = big.matrix(nrow = nrow(u0), ncol = ncol(Z))
    for(k in 1:d){
      idx = seq(from = k, to = ncol(Z), by = d)
      idx_u = seq(from = k, to = d*r, by = d)
      post_U[,idx] = u0[,idx_u] %*% t(B)
    }
    
    # Take last nMC_report rows of post_U to save for MCMC diagnostics
    r_start = nrow(post_U) - nMC_report + 1
    r_end = nrow(post_U)
    post_out = bigmemory::as.matrix(post_U[c(r_start:r_end),])
    
    post_modes = numeric(ncol(post_U))
    for(j in 1:ncol(post_U)){
      post_modes[j] = mean(post_U[,j])
    }
    
    # Calculate logLik
    ll = CAME_IS_FA(posterior = u0, y = y, X = X, Z = Z, group = group,
                    coef = coef,  B = t(matrix(coef[-c(1:ncol(X))], ncol = q, nrow = r)),
                    family = fam_fun, M = M, 
                    gaus_sig = sig_g, progress = progress, trace = trace)
    
    # Hybrid BIC (Delattre, Lavielle, and Poursat (2014))
    # d = nlevels(group) = number independent subjects/groups
    BICh = -2*ll + sum(coef[-c(1:ncol(X))] != 0)*log(d) + sum(coef[1:ncol(X)] != 0)*log(nrow(X))
    # Usual BIC
    BIC = -2*ll + sum(coef != 0)*log(nrow(X))
    # BIC using N = nlevels(group)
    BICNgrp = -2*ll + sum(coef != 0)*log(d)
    
    return(list(u0 = Estep_out$u0, post_modes = post_modes, post_out = post_out, 
                u_init = matrix(u0[nrow(u0),], nrow = 1),
                ll = ll, BICh = BICh, BIC = BIC, BICNgrp = BICNgrp))
    
  }else{
    
    u0_mat = attach.big.matrix(Estep_out$u0)
    
    return(list(u0 = Estep_out$u0, u_init = u0_mat[nrow(u0_mat),]))
    
  }
  
  
  
  
}