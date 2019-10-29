 


#' Calculate Monte Carlo draws
#' 
#' \code{sample.mc2} calculates the Monte Carlo draws by either Rejection sampling or 
#' Metropolis-within-Gibbs sampling needed for Monte Carlo Expectation Conditional Minimization (MCECM)
#'
#' @inheritParams fit_dat
#' @param fit a grpreg fit object (set \code{\link[grpreg]{grpreg}})
#' @param cov the random effects covariance matrix estimated from the last M step of the EM algorithm
#' @param y a numeric vector of the response variable
#' @param X a model matrix of the fixed effects covariates
#' @param Z a sparse model matrix of the random effects
#' @param group a factor vector of the grouping variable, converted to a factor of 
#' consecutive numeric integers
#' @param d integer, the number of groups present (number factors of the group vector)
#' @param okindex ?
#' @param uold a matrix of the Monte Carlo draws from the last E step of the EM algorithm
#' 
#' @return a list made of the following components:
#' \item{u0}{a matrix of the Monte Carlo draws (from Rejection sampling if gibbs = F, or 
#' from Metropolis-within-Gibbs sampling if gibbs = T). Number rows = nMC, number columns = 
#' (number random effects)x(number groups). Organization of columns: first by random effect variable,
#' then by group within variable (i.e. Var1:Grp1 Var1:Grp2 ... Var1:GrpK Var2:Grp1 ... Varq:GrpK)}
#' \item{switch}{logical, did the sampling scheme switch from Rejection sampling to 
#' Metropolis-within-Gibbs sampling?}
#'
#' @export
sample.mc2 = function(fit, cov, y, X, Z, nMC, trace = 0, family = family, group, d, nZ, okindex,
                      gibbs = F , uold){
  
  f = get(family, mode = "function", envir = parent.frame())
  
  # find tau for rejection sampler (Booth 1999) by first maximizing u
  if(length(okindex) > 0 & gibbs == F){
    fitu = glm(y ~ Z[,okindex]-1, offset = (X %*% fit$coef[1:ncol(X)]), family = f)
    uhat = fitu$coef
    uhat[is.na(uhat)] = 0
    eta = X %*% fit$coef[1:ncol(X)] + Z[,okindex] %*% uhat
  }else{
    uhat  = rep(0, ncol(Z))
    eta = X %*% fit$coef[1:ncol(X)] 
  }
  
  #if(trace == 1) print(uhat)
  # calculate tau
  
  if(family == "binomial" & gibbs == F){
    tau = dbinom(y, size = 1, prob = exp(eta)/(1+exp(eta)), log = T)
  }else if(family == "poisson" & gibbs == F){
    tau = dpois(y, lambda = exp(eta), log = T)
  }
  
  #if(any(tau == 0)) tau = tau - 10^-20
  # matrix to hold accepted samples
  u0 = matrix(rnorm(nMC*ncol(Z)) , nMC, ncol(Z))
  #u0 = Matrix(0 , nMC, ncol(Z), sparse = T)
  
  # fitted
  fitted_mat = as.matrix(X %*% fit$coef[1:ncol(X)])
  #generate samples for each i
  
  error_out = F
  q = ncol(Z) / d
  gibbs_accept_rate = matrix(NA, nrow = d, ncol = q)
  
  if(gibbs == F){ # Rejection sampling
    
    for(i in 1:d){
      
      if(error_out){
        next
      }
      
      select = group == i
      index = seq(i, ncol(Z), by = d)
      ## new code to limit to non-zero z, skipping elements of q where diag(sigma) are 0
      index = index[which(diag(cov) != 0)]
      if(length(index) == 0) next
      
      draws = sample_mc_inner(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                              matrix(Z[select,index],ncol = length(index), nrow = sum(select)), 
                              y[select], tau[select], nMC, trace)
      
      # If sample_mc_inner did not error out (i.e. have too small an acceptance rate), then continue
      if(nrow(draws) == nMC){
        u0[,index] = draws
      }else{ # If too small an acceptance rate
        cat("switched to gibbs sampling (single round) \n")
        error_out = T
        next
      }
      
    }
    
    # If at any point the sample_mc_inner errored out due to too low of acceptance rate, then 
    # switch to gibbs sampling (redo entire u matrix as gibbs sampling)
    if(error_out){
      for(i in 1:d){
        select = group == i
        index = seq(i, ncol(Z), by = d)
        ## new code to limit to non-zero z, skipping elements of q where diag(sigma) are 0
        index = index[which(diag(cov) != 0)]
        if(length(index) == 0) next
        
        gibbs_list = sample_mc_inner_gibbs(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                                           matrix(Z[select,index],ncol = length(index), nrow = sum(select)),  
                                           y[select], uhat[index], nMC, as.numeric((uold[nrow(uold),index, drop = FALSE])),
                                           trace)
        u0[,index] = gibbs_list$u
        gibbs_accept_rate[i,] = matrix(gibbs_list$acc_rate, nrow = 1)
      }
    }
    
  }else{ # Gibbs sampling
    for(i in 1:d){
      select = group == i
      index = seq(i, ncol(Z), by = d)
      ## new code to limit to non-zero z, skipping elements of q where diag(sigma) are 0
      index = index[which(diag(cov) != 0)]
      if(length(index) == 0) next
      
      gibbs_list = sample_mc_inner_gibbs(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                                         matrix(Z[select,index],ncol = length(index), nrow = sum(select)),  
                                         y[select], uhat[index], nMC, as.numeric((uold[nrow(uold),index, drop = FALSE])), 
                                         trace)
      u0[,index] = gibbs_list$u
      gibbs_accept_rate[i,] = matrix(gibbs_list$acc_rate, nrow = 1)
    }
  }
  # for each i, rbind the nMC samples together to make n*nMC x d matrix (d = dim(Z))
  
  if(gibbs == F){
    cat("switch: ", error_out, "\n")
  }
  
  # return(u0)
  # switch: variable indicating if switched from rejection sampling to gibbs sampling
  if(trace == 2 && !anyNA(gibbs_accept_rate)){
    return(list(u0 = u0, gibbs_accept_rate = gibbs_accept_rate, switch = error_out))
  }else{
    return(list(u0 = u0, switch = error_out))
  }
  
}

##################################################################
# Updated for adaptive metropolis w/in gibbs
##################################################################

#' @export
sample.mc3 = function(fit, cov, y, X, Z, nMC, trace = 0, family = family, group, d, nZ, okindex,
                      gibbs = F , uold, proposal_var, batch){
  
  f = get(family, mode = "function", envir = parent.frame())
  
  # find tau for rejection sampler (Booth 1999) by first maximizing u
  if(length(okindex) > 0 & gibbs == F){
    fitu = glm(y ~ Z[,okindex]-1, offset = (X %*% fit$coef[1:ncol(X)]), family = f)
    uhat = fitu$coef
    uhat[is.na(uhat)] = 0
    eta = X %*% fit$coef[1:ncol(X)] + Z[,okindex] %*% uhat
  }else{
    uhat  = rep(0, ncol(Z))
    eta = X %*% fit$coef[1:ncol(X)] 
  }
  
  #if(trace == 1) print(uhat)
  # calculate tau
  
  if(family == "binomial" & gibbs == F){
    tau = dbinom(y, size = 1, prob = exp(eta)/(1+exp(eta)), log = T)
  }else if(family == "poisson" & gibbs == F){
    tau = dpois(y, lambda = exp(eta), log = T)
  }
  
  #if(any(tau == 0)) tau = tau - 10^-20
  # matrix to hold accepted samples
  u0 = matrix(rnorm(nMC*ncol(Z)) , nMC, ncol(Z))
  #u0 = Matrix(0 , nMC, ncol(Z), sparse = T)
  
  # fitted
  fitted_mat = as.matrix(X %*% fit$coef[1:ncol(X)])
  #generate samples for each i
  
  error_out = F
  q = ncol(Z) / d
  gibbs_accept_rate = matrix(NA, nrow = d, ncol = q)
  
  if(gibbs == F){ # Rejection sampling
    
    for(i in 1:d){
      
      if(error_out){
        next
      }
      
      select = group == i
      index = seq(i, ncol(Z), by = d)
      ## new code to limit to non-zero z, skipping elements of q where diag(sigma) are 0
      index = index[which(diag(cov) != 0)]
      if(length(index) == 0) next
      
      draws = sample_mc_inner(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                              matrix(Z[select,index],ncol = length(index), nrow = sum(select)), 
                              y[select], tau[select], nMC, trace)
      
      # If sample_mc_inner did not error out (i.e. have too small an acceptance rate), then continue
      if(nrow(draws) == nMC){
        u0[,index] = draws
      }else{ # If too small an acceptance rate
        cat("switched to gibbs sampling (single round) \n")
        error_out = T
        next
      }
      
    }
    
    # If at any point the sample_mc_inner errored out due to too low of acceptance rate, then 
    # switch to gibbs sampling (redo entire u matrix as gibbs sampling)
    if(error_out){
      for(i in 1:d){
        select = group == i
        index = seq(i, ncol(Z), by = d)
        ## new code to limit to non-zero z, skipping elements of q where diag(sigma) are 0
        index = index[which(diag(cov) != 0)]
        if(length(index) == 0) next
        var_index = which(diag(cov) != 0)
        
        gibbs_list = sample_mc_inner_gibbs(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                                           matrix(Z[select,index],ncol = length(index), nrow = sum(select)),  
                                           y[select], uhat[index], nMC, as.numeric((uold[nrow(uold),index, drop = FALSE])), 
                                           proposal_var[var_index], batch, trace)
        u0[,index] = gibbs_list$u
        gibbs_accept_rate[i,] = matrix(gibbs_list$acc_rate, nrow = 1)
        proposal_var[i,var_index] = gibbs_list$proposal_var
        batch = gibbs_list$batch
      }
    }
    
  }else{ # Gibbs sampling
    for(i in 1:d){
      select = group == i
      index = seq(i, ncol(Z), by = d)
      ## new code to limit to non-zero z, skipping elements of q where diag(sigma) are 0
      index = index[which(diag(cov) != 0)]
      if(length(index) == 0) next
      var_index = which(diag(cov) != 0)
      
      gibbs_list = sample_mc_inner_gibbs(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                                         matrix(Z[select,index],ncol = length(index), nrow = sum(select)),  
                                         y[select], uhat[index], nMC, as.numeric((uold[nrow(uold),index, drop = FALSE])), 
                                         proposal_var[var_index], batch, trace)
      u0[,index] = gibbs_list$u
      gibbs_accept_rate[i,] = matrix(gibbs_list$acc_rate, nrow = 1)
      proposal_var[i,var_index] = gibbs_list$proposal_var
      batch = gibbs_list$batch
    }
  }
  # for each i, rbind the nMC samples together to make n*nMC x d matrix (d = dim(Z))
  
  if(gibbs == F){
    cat("switch from rejection to gibbs: ", error_out, "\n")
  }
  
  # switch: variable indicating if switched from rejection sampling to Metropolis-within-Gibbs sampling
  if(!anyNA(gibbs_accept_rate)){
    return(list(u0 = u0, gibbs_accept_rate = gibbs_accept_rate, 
                switch = error_out, proposal_var = proposal_var))
  }else{
    return(list(u0 = u0, switch = error_out, proposal_var = proposal_var))
  }
  
}

