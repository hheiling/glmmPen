 



############################################################################
# Updated for adaptive random walk metropolis w/in gibbs with random scan
# Adapt proposal SD only during burn-in period
# Output a big.matrix for the posterior draws 
############################################################################


#' @importFrom bigmemory big.matrix describe
sample_mc_adapt_BigMat = function(coef, ranef_idx, y, X, Z, nMC, family, link, group, d,
                                  uold, proposal_SD, batch, batch_length = 100, 
                                  offset = 0, nMC_burnin = 500, phi = 0.0, sig_g = 1.0){
  
  # Re-code link as integer
  ## All link_int will have two digits
  ## First digit corresponds to family that link is canonical for
  ## 1 = binomial, 2 = poisson or negative binomial, 3 = gaussian
  ## Second digit: 0 = canonical link, other = arbitrary enumeration of common non-canonical links
  if(link == "logit"){
    link_int = 10
  }else if(link == "probit"){
    link_int = 11
  }else if(link == "cloglog"){
    link_int = 12
  }else if(link == "log"){
    link_int = 20
  }else if(link == "sqrt"){
    link_int = 21
  }else if(link == "identity"){
    link_int = 30
  }else if(link == "inverse"){
    link_int = 31
  }
  
  eta = X %*% matrix(coef, ncol=1) 
  
  # matrix to hold accepted samples
  u0 = big.matrix(nrow = nMC, ncol = ncol(Z), init = 0)
  # u0 = matrix(rnorm(nMC*ncol(Z)) , nMC, ncol(Z))
  
  # fitted
  fitted_mat = as.matrix(X %*% matrix(coef, ncol=1))
  #generate samples for each i
  
  q = ncol(Z) / d
  gibbs_accept_rate = matrix(NA, nrow = d, ncol = q)
  gibbs_output = NULL

  for(i in 1:d){
    select = group == i
    index = seq(i, ncol(Z), by = d)
    ## new code to limit to non-zero z, skipping elements of q where diag(sigma) are 0
    # index = index[which(diag(cov) != 0)]
    index = index[ranef_idx]
    if(length(index) == 0) next
    var_index = ranef_idx
    
    gibbs_output = sample_mc_gibbs_adapt_rw(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                                            matrix(Z[select,index],ncol = length(index), nrow = sum(select)),  
                                            y[select], nMC, uold[index], 
                                            matrix(proposal_SD[i,var_index], nrow = 1), batch, 
                                            batch_length, offset, nMC_burnin,
                                            family, link_int, phi, sig_g)
    
    u0[,index] = gibbs_output[1:nMC,]
    gibbs_accept_rate[i,var_index] = matrix(gibbs_output[(nMC+1),], nrow = 1)
    proposal_SD[i,var_index] = matrix(gibbs_output[(nMC+2),], nrow = 1)
    
  }
  
  
  return(list(u0 = describe(u0), gibbs_accept_rate = gibbs_accept_rate, 
              proposal_SD = proposal_SD, updated_batch = batch))
  
} # End sample_mc_adapt_BigMat()

#' @importFrom bigmemory big.matrix describe
sample_mc2_BigMat = function(coef, ranef_idx, y, X, Z, nMC, family, link, group, d, 
                             uold, phi = 0.0, sig_g = 1.0){
  
  # Re-code link as integer
  ## All link_int will have two digits
  ## First digit corresponds to family that link is canonical for
  ## 1 = binomial, 2 = poisson or negative binomial, 3 = gaussian
  ## Second digit: 0 = canonical link, other = arbitrary enumeration of common non-canonical links
  if(link == "logit"){
    link_int = 10
  }else if(link == "probit"){
    link_int = 11
  }else if(link == "cloglog"){
    link_int = 12
  }else if(link == "log"){
    link_int = 20
  }else if(link == "sqrt"){
    link_int = 21
  }else if(link == "identity"){
    link_int = 30
  }else if(link == "inverse"){
    link_int = 31
  }
  
  uhat  = rep(0, ncol(Z))
  eta = X %*% matrix(coef, ncol=1) 
  
  # matrix to hold accepted samples
  u0 = big.matrix(nrow = nMC, ncol = ncol(Z), init = 0)
  # u0 = matrix(rnorm(nMC*ncol(Z)) , nMC, ncol(Z))
  
  # fitted
  fitted_mat = as.matrix(X %*% matrix(coef, ncol=1))
  #generate samples for each i
  
  error_out = F
  q = ncol(Z) / d
  gibbs_accept_rate = matrix(NA, nrow = d, ncol = q)
  
  
  for(i in 1:d){
    select = group == i
    index = seq(i, ncol(Z), by = d)
    ## new code to limit to non-zero z, skipping elements of q where diag(sigma) are 0
    # index = index[which(diag(cov) != 0)]
    index = index[ranef_idx]
    if(length(index) == 0) next
    
    gibbs_list = sample_mc_inner_gibbs(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                                       matrix(Z[select,index],ncol = length(index), nrow = sum(select)),  
                                       y[select], uhat[index], nMC, uold[index],
                                       family, link_int, phi, sig_g)
    u0[,index] = gibbs_list$u
    gibbs_accept_rate[i,] = matrix(gibbs_list$acc_rate, nrow = 1)
  }
  
  return(list(u0 = describe(u0), gibbs_accept_rate = gibbs_accept_rate))
  
  
} # End sample_mc2_BigMat

# tau calculate for rejection sampling - discontinued
# if(family == "binomial" & gibbs == F){
#   tau = dbinom(y, size = 1, prob = exp(eta)/(1+exp(eta)), log = T)
# }else if(family == "poisson" & gibbs == F){
#   tau = dpois(y, lambda = exp(eta), log = T)
# }else if(family == "gaussian" & gibbs == F){
#   s2 = sum((y-eta)*(y-eta)) / (length(y) - length(coef))
#   tau = dnorm(y, mean = eta, sd = sqrt(s2), log = T)
# }


###################################################################################################