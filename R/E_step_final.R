
#' @importFrom bigmemory big.matrix attach.big.matrix
E_step_final = function(dat, fit, optim_options, fam_fun, extra_calc, adapt_RW_options, trace){
  
  # Extract terms from dat list object
  y = dat$y
  X = dat$X # standardized X
  Z = Matrix::as.matrix(dat$Z)
  group = dat$group
  d = nlevels(group)
  
  # Extract family and link information
  family = fam_fun$family
  link = fam_fun$link
  
  # Extract terms from fitted model object
  sigma = fit$sigma
  J = fit$J
  coef = fit$coef
  sig_g = ifelse(family == "gaussian", fit$sigma_gaus, 1.0)
  
  # Extract control values from optim_options
  q = sum(diag(sigma) != 0)
  M = 10^4
  # if(q <= 21){
  #   M = 10^4
  # }else if(q <= 31){
  #   M = 7500
  # }else{
  #   M = 5000
  # }
  # M = max(M, optim_options$nMC_report) # default nMC_report: 5000
  nMC_burnin = optim_options$nMC_burnin
  nMC_report = optim_options$nMC_report
  if(M < nMC_report){
    nMC_report = M
  }
  sampler = optim_options$sampler
  max_cores = optim_options$max_cores
  
  
  Znew2 = Z
  Gamma_mat = matrix(J%*%matrix(coef[-c(1:ncol(X))], ncol = 1), ncol = ncol(Z)/d)
  finish = 0
  while(finish == 0){
    for(k in 1:d){
      Znew2[group == k,seq(k, ncol(Z), by = d)] = Z[group == k,seq(k, ncol(Z), by = d)]%*%Gamma_mat
    }
    if(!any(is.na(Znew2))) finish = 1
  }
  
  print("Start of sampling from posterior")
  Estep_out = E_step(coef = coef, ranef_idx = which(diag(sigma > 0)), y = y,
                     X = X, Znew2 = Znew2, group = group, nMC = M, nMC_burnin = nMC_burnin,
                     family = family, link = link, phi = 1.0, sig_g = sig_g,
                     sampler = sampler, d = d, uold = as.numeric(fit$u_init),
                     proposal_SD = fit$proposal_SD, batch = fit$updated_batch, 
                     batch_length = adapt_RW_options$batch_length, 
                     offset_increment = adapt_RW_options$offset, trace = trace, 
                     max_cores = max_cores)
  print("Finished sampling from posterior")
  
  if(extra_calc){
    
    u0 = attach.big.matrix(Estep_out$u0)
    
    # Calculate posterior modes
    post_U = big.matrix(nrow = nrow(u0), ncol = ncol(u0))
    for(k in 1:d){
      idx = seq(from = k, to = ncol(Z), by = d)
      post_U[,idx] = u0[,idx] %*% t(Gamma_mat)
    }
    
    # Take last nMC_report rows of post_U to save for MCMC diagnostics
    r_start = nrow(post_U) - nMC_report + 1
    r_end = nrow(post_U)
    post_out = bigmemory::as.matrix(post_U[c(r_start:r_end),])
    
    post_modes = numeric(ncol(u0))
    for(j in 1:ncol(u0)){
      post_modes[j] = mean(post_U[,j])
    }
    
    # Calculate logLik
    ll = CAME_IS(posterior = u0, y = y, X = X, Z = Z, group = group,
                 coef = coef, sigma = sigma, family = fam_fun, M = M, gaus_sig = sig_g, trace = trace)
    
    # Hybrid BIC (Delattre, Lavielle, and Poursat (2014))
    # d = nlevels(group) = number independent subjects/groups
    BICh = -2*ll + sum(coef[-c(1:ncol(X))] != 0)*log(d) + sum(coef[1:ncol(X)] != 0)*log(nrow(X))
    # Usual BIC
    BIC = -2*ll + sum(coef != 0)*log(nrow(X))
    # BIC using N = nlevels(group)
    BICNgrp = -2*ll + sum(coef != 0)*log(d)
    
    return(list(u0 = Estep_out$u0, post_modes = post_modes, post_out = post_out, 
                u_init = matrix(u0[nrow(u0),]),
                ll = ll, BICh = BICh, BIC = BIC, BICNgrp = BICNgrp))
    
  }else{
    
    u0_mat = attach.big.matrix(Estep_out$u0)
    
    return(list(u0 = Estep_out$u0, u_init = u0_mat[nrow(u0_mat),]))
    
  }
  
  
  
  
}