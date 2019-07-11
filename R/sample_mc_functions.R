#' @export
sample.mc = function(fit, cov, y, X, Z, nMC, trace = 0, family = "binomial", group, d, nZ, okindex ){
  # Things to address:
  ## Will need to make sure family "if else" statements are consistent 
  
  f = get(family, mode = "function", envir = parent.frame())
  
  # find tau for rejection sampler (Booth 1999) by first maximizing u
  if(length(okindex) > 0){
    fitu = glm(y ~ Z[,okindex]-1, offset = (X %*% fit$coef[1:ncol(X)]), family = f)
    uhat = fitu$coef
    uhat[is.na(uhat)] = 0
    eta = X %*% fit$coef[1:ncol(X)] + Z[,okindex] %*% uhat
  }else{
    uhat  = rep(0, ncol(Z))
    eta = X %*% fit$coef[1:ncol(X)] 
  }
  
  if(trace == 1) print(uhat)
  # calculate tau
  
  if(family == "binomial"){
    tau = dbinom(y, size = 1, prob = exp(eta)/(1+exp(eta)), log = T)
  }else if(family == "poisson"){
    tau = dpois(y, lambda = exp(eta), log = T)
  }else{
    print(family)
    stop("family not specified properly")
  }
  
  if(any(tau == 0)) tau = tau - 10^-20
  # matrix to hold accepted samples
  u0 = matrix(0 , nMC, ncol(Z))
  
  # fitted
  fitted = X %*% as.matrix(fit$coef[1:ncol(X)])
  #generate samples for each i
  for(i in 1:d){
    
    # find tau for rejection sampler (Booth 1999) by first maximizing u
    #Z2  = Z[i,] %*% matrix(1, nrow = ncol(Z), ncol = 1) 
    #fitu = glm(y[i] ~ Z2-1, offset = X[i,] %*% fit$coef, family = "poisson")
    #uhat = fitu$coef
    
    # calculate tau
    #eta = X[i,] %*% fit$coef + Z2 %*% uhat
    #tau = dpois(y[i], lambda = exp(eta))
    
    naccept = index = 0
    alreadyok = rep(F, nMC)
    while(naccept < nMC){
      # generate e
      #e = matrix(rmvnorm(n = nMC, mean  = rep(0, ncol(cov)), sigma = cov), ncol = ncol(cov))
      #e = matrix(rmvnorm(n = 1, mean  = rep(0, ncol(Z)), sigma = diag(rep(cov,ncol(Z)))), ncol = ncol(Z))
      select = group == i
      #e = rnorm(n = ceiling((index+1)/(naccept+1))*10, mean  = 0, sd = sqrt(cov))
      e = matrix(rmvnorm(n = ceiling((index+1)/(naccept+1))*10, mean  = rep(0, ncol(cov)), sigma = cov), ncol = ncol(cov))
      
      # vector of MC candidates, may need to revert to initial approach when have more complicated Z
      #etae = matrix(fitted[select], nrow = sum(select), ncol = length(e) ) +  ( matrix(e, nrow = sum(select), ncol = length(e), byrow = T))#sapply(e, FUN = function(ei) return(X[select,] %*% fit$coef + rowSums(Z[select,] * ei)))
      etae = matrix(fitted[select], nrow = sum(select), ncol = nrow(e)) +  Z[select,seq(i, ncol(Z), by = d)]%*%t(e)#sapply(e, FUN = function(ei) return(X[select,] %*% fit$coef + rowSums(Z[select,] * ei)))
      
      
      # generate w
      w = runif(nrow(e))
      
      
      #if(w < dpois(y[i], lambda =  exp(etae))/tau){
      # pick the first e that satisfies the condition
      
      if(family == "binomial"){
        ok = log(w) < apply(etae, 2, FUN = function(etaei) sum(dbinom(y[select], size = 1, prob = exp(etaei)/(1+exp(etaei)), log = T)) - sum(tau[select]))
      }else if(family == "poisson"){
        ok = log(w) < apply(etae, 2, FUN = function(etaei) sum(dpois(y[select], lambda =  exp(etaei), log = T)) - sum(tau[select]))
      }else{
        print(family)
        stop("family not specified properly")
      }
      
      
      if(any(ok)){
        u0[naccept+1,seq(i, ncol(Z), by = d)] = e[which(ok)[1],]
        #alreadyok = rowSums(u0) != 0
        naccept = naccept + 1#sum(alreadyok)
        if(trace == 1) cat(sprintf("mc i:%d naccept: %d index: %d\n",i, naccept, index))
      }
      index = index + 1
    }
    # for each i, rbind the nMC samples together to make n*nMC x d matrix (d = dim(Z))
  }
  return(u0)
}

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
  
  # switch: variable indicating if switched from rejection sampling to gibbs sampling
  switch = F
  k = 0
  
  for(i in 1:d){
    select = group == i
    index = seq(i, ncol(Z), by = d)
    ## new code to limit to non-zero z, skipping elements of q where diag(sigma) are 0
    index = index[which(diag(cov) != 0)]
    if(length(index) == 0) next
    ###
    k = k + 1
    
    if(gibbs == F){
      # If for first group, acceptance rate too low (nrow(output) < nMC), 
      # then switch to gibbs for all other groups
      if(k == 1){
        draws = sample_mc_inner(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                                matrix(Z[select,index],ncol = length(index), nrow = sum(select)), 
                                y[select], tau[select], nMC, trace)
        if(nrow(draws) == nMC){
          u0[,index] = draws
        }else{
          gibbs = T
          u0[,index] = sample_mc_inner_gibbs(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                                             matrix(Z[select,index],ncol = length(index), nrow = sum(select)),  
                                             y[select], uhat[index], nMC, 
                                             as.numeric((uold[nrow(uold),index, drop = FALSE])), trace)
          switch = T
          }
      }else{
        u0[,index] = sample_mc_inner(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                                matrix(Z[select,index],ncol = length(index), nrow = sum(select)), 
                                y[select], tau[select], nMC, trace)
      }
      
    }else{
      u0[,index] = sample_mc_inner_gibbs(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                                         matrix(Z[select,index],ncol = length(index), nrow = sum(select)),  
                                         y[select], uhat[index], nMC, 
                                         as.numeric((uold[nrow(uold),index, drop = FALSE])),
                                         trace) 
    }    
  }
  # for each i, rbind the nMC samples together to make n*nMC x d matrix (d = dim(Z))
  
  # return(u0)
  return(list(u0 = u0, switch = switch))
}

