
#' Fit a penalized generalized mixed model (pGLMM) via Monte Carlo Expectation Conditional 
#' Minimization (MCECM)
#' 
#' @param dat Input list object specifying y (response vector), X (matrix of all covariates), 
#' Z (matrix of covariates with possible random effects), and group (vector whose value indicates 
#' the study, batch, or other group identity to which on observation belongs), with each row an 
#' observation
#' 
#' @export
fit_dat = function(dat,  lambda0 = 0, lambda1 = 0, conv = 0.001, nMC = 1000, 
                   family = "binomial", type = "hmmcov", trace = 0, initcov = NULL, 
                   glmIALconv = 1e-12, d = 2, group, alpha = 1, vartol = 0.1, nMC_max = 5000, 
                   returnMC = F, ufull = NULL, coeffull = NULL, gibbs = F, maxitEM = 100, 
                   pnonzerovar = 0, ufullinit = NULL){
  
  # Things to address:
  ## Already deal with family modification in 'wrapper' function glmmPen
  ## In dat, what should pnonzero be?
  
  ## added nov 4th 2016 for simulation purposes only to set small penalties to 0
  if(lambda0 <=10^-6) lambda0 = 0
  if(lambda1 <=10^-6) lambda1 = 0
  ##
  
  y = dat$y
  X = as.matrix(dat$X)
  Z = as.matrix(dat$Z)
  group = dat$group
  p1 = dat$pnonzero
  lambda1 = lambda1
  lambda0 = lambda0
  
  f = get(family, mode = "function", envir = parent.frame())
  # Deleted commented section below later: already have restriction on family in glmmPen()
  # if(family == "binomial"){
  #   f = binomial()
  # }else if(family == "poisson"){
  #   f = poisson()
  # }else{
  #   print(family)
  #   stop("family not specified properly")
  # }
  
  #group = factor(apply(Z[,1:2], 1, FUN = function(x) which(x == 1))) #as.numeric(Z[,1])
  d = nlevels(factor(group))
  # fit all data using glmer from lme4 package
  ## Note: delete all references to these variables (and references to pnonzervar)
  # if(ncol(X) < 5){
  #   if(ncol(Z) > d){
  #     fit_glmer = try(glmer(y ~ X-1 + (X-1|group), family = f, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))))
  #   }else{
  #     fit_glmer = try(glmer(y ~ X-1 + (1|group), family = f, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))))
  #   }
  #   # fit true (oracle) model
  #   fit_glmer_oracle = fit_glmer #glmer(y ~ X[,1:(p1+1)]-1 + (X[,1:(p1+1)]-1|group), family = f, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
  #   
  # }else{
  #   
  #   # fit true (oracle) model
  #   fit_glmer = fit_glmer_oracle = try(glmer(y ~ X[,1:(p1+1+pnonzerovar)]-1 + (X[,1:(p1+1+pnonzerovar)]-1|group), family = f, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))))
  #   
  # }
  print("end glmer")
  
  if(is.character(fit_glmer_oracle)){
    fit_glmer = fit_glmer_oracle = try(glmer(y ~ X[,1:(p1+1+pnonzerovar)]-1 + (1|group), family = f, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))))
  }
  vc <- VarCorr(fit_glmer_oracle)
  print(summary(fit_glmer)$coefficients[-1,1])
  print(vc,comp=c("Std.Dev."),digits=3)
  
  #initial fit
  if(family == "binomial"){
    nTotal = rep(1, length(y))
  }else{
    nTotal = NULL
  }
  
  
  if(ncol(Z)/d <= 15){
    # create J, q2 x q*(q+1)/2
    J = Matrix(0, (ncol(Z)/d)^2, (ncol(Z)/d)*((ncol(Z)/d)+1)/2, sparse = T) #matrix(0, (ncol(Z)/d)^2, (ncol(Z)/d)*((ncol(Z)/d)+1)/2)
    index = 0
    indexc = 0
    sumy = 0
    sumx = 0
    zeros = 0
    covgroup = NULL
    for(i in 1:(ncol(Z)/d)){
      J[ sumx + zeros + 1:(ncol(Z)/d - (i-1)), sumy + 1:(ncol(Z)/d - (i-1))] = diag((ncol(Z)/d - (i-1)))
      sumy = sumy + (ncol(Z)/d - (i-1))
      sumx = sumy 
      zeros = zeros + i
      covgroup = rbind(covgroup, rep(i, (ncol(Z)/d)))
      #print(c((ncol(Z)/d - (i-1)),sumy, sumx))
    }
    covgroup = covgroup[lower.tri(covgroup, diag = T)]
    #if(ncol(Z)/d > 20){
  }else{
    J = Matrix(0,(ncol(Z)/d)^2, (ncol(Z)/d), sparse = T) #matrix(0, (ncol(Z)/d)^2, (ncol(Z)/d))
    index = 0
    indexc = 0
    sumy = 0
    sumx = 0
    zeros = 0
    covgroup = NULL
    for(i in 1:(ncol(Z)/d)){
      J[ sumx + zeros + 1, i] = 1
      sumy = sumy + (ncol(Z)/d - (i-1))
      sumx = sumy 
      zeros = zeros + i
      #print(c((ncol(Z)/d - (i-1)),sumy, sumx))
    }
    covgroup = rep(1:(ncol(Z)/d))
  }
  #J = Matrix(J, sparse = TRUE)
  #J = as.big.matrix(J)
  # initial value for beta
  #fit = grpreg(X[,-1], y, group=1:(ncol(X)-1), penalty="grMCP", family="binomial",lambda = lambda1, alpha = alpha)
  #fit00 = fit # naive fit  
  
  if(!is.null(ufull) & !is.null(coeffull)){
    fit = list()
    print("using coef from full model to intialize")
    coef = coeffull
    gamma = matrix(J%*%coef[-c(1:ncol(X))], ncol = ncol(Z)/d)
    cov = var = gamma %*% t(gamma)
    fit$coef = coef[c(1:ncol(X))]#as.numeric(fit$beta)
    ok = which(diag(var) > vartol)# & coef[1:ncol(X)] != 0)
    if(length(ok) == 0) ok = 1 # at least include the random intercept
    okindex = NULL
    for(i in 1:(ncol(Z)/d)){
      if(i %in% ok){
        okindex = c(okindex, (i-1)*d + 1:d)
      }
    }
    fit00 = fit
  }else{
    fit = grpreg(X[,-1], y, group=1:(ncol(X)-1), penalty="grMCP", family="binomial",lambda = lambda1, alpha = alpha)###
    fit00 = fit #### naive fit
    
    coef = as.numeric(fit$beta)
    fit$coef = as.numeric(fit$beta)
    if(is.null(initcov)){
      X00 = matrix(0, nrow(Z), ncol(Z)/d)
      for(j in 1:d){
        X00 = X00 + Z[,seq(j, ncol(Z), by = d)] 
      }
      vars = rep(0.0, ncol(X00))
      #vars[1] = as.numeric(unlist((VarCorr(glmer(y ~ (1|group), family = f)))))
      #for(i in 2:ncol(X00)) vars[i] = as.numeric(unlist((VarCorr(glmer(y ~ (X00[,i]-1|group), family = f)))))
    }else{
      vars = rep(initcov, ncol(X))
    }
    
    if(trace == 1) print(coef);print(vars)
    
    cov = var = diag(vars+10^-10)
    gamma = t(chol(var)) # chol outputs upper triangular, so transposing here
    
    ok = which(vars > vartol)# & coef[1:ncol(X)] != 0)
    if(length(ok) == 0) ok = 1 # at least include the random intercept
    okindex = NULL
    for(i in 1:(ncol(Z)/d)){
      if(i %in% ok){
        okindex = c(okindex, (i-1)*d + 1:d) 
      }
    }
  }
  
  
  Znew2 = Z
  finish = 0
  while(finish == 0){
    for(i in 1:d){
      #print(Z[group == i,seq(i, ncol(Z), by = d)]%*% gamma)
      Znew2[group == i,seq(i, ncol(Z), by = d)] = Z[group == i,seq(i, ncol(Z), by = d)]%*% gamma
    }
    if(!any(is.na(Znew2))) finish = 1
  }
  
  if((!is.null(ufull) | !is.null(ufullinit)) & !is.null(coeffull)){
    if(!is.null(ufullinit)){
      print("using u from previous model to intialize")
    }else{
      print("using u from full model to intialize")
      ufullinit = ufull
    }
    u = u0 = sample.mc2(fit=fit, cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group, d = d, okindex = okindex, trace = trace, gibbs = gibbs, uold = ufullinit)
  }else{
    u = u0 = sample.mc2(fit=fit, cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group, d = d, okindex = okindex, trace = trace, gibbs = gibbs, uold = matrix(rnorm(nMC*ncol(Z)), nrow = nMC, ncol = ncol(Z)))
  }
  #u = bmmat(u)
  nMC2 = nrow(u)  
  etae = matrix(X %*% coef[1:ncol(X)], nrow = nrow(X), ncol = nrow(u) ) + Znew2%*%t(u)
  
  
  diff = rep(NA, 100)
  stopcount = 0
  
  #etae = matrix(X %*% coef, nrow = nrow(X), ncol = nrow(u) ) + Z%*%t(u)
  
  ## this needs to get updated to reflect whatever is chosen to evaluate likelihood
  if(family == "poisson"){
    ll = ll0 = ll20 = sum(rowMeans(dpois(matrix(y, nrow = length(y), ncol = ncol(etae)), lambda =  exp(etae), log = T)))
  }else if(family == "binomial"){
    ll = ll0 = ll20  = sum(rowMeans(dbinom(matrix(y, nrow = length(y), ncol = ncol(etae)), size = 1, prob = exp(etae)/(1+exp(etae)), log = T)))
  }  
  Znew = NULL
  
  
  # initialize zero count vectors
  c0 = rep(0, length(coef))
  for(i in 1:maxitEM){
    
    oldll = ll0
    
    if(family == "binomial"){
      nTotal = rep(1, length(y[rep(1:nrow(X), each = nrow(u))]))
    }else{
      nTotal = NULL
    }
    
    print("Znewgen done")
    rm(Znew)
    gc()
    Znew = big.matrix(nrow = nrow(X)*nrow(u), ncol = ncol(J))
    Znew_gen2(u, Z, group, seq(as.numeric(group[1]), ncol(Z), by = d),nrow(Z),ncol(Z)/d,d, Znew@address, J)
    gc()
    
    active0 = rep(1, max(covgroup))
    active1 = rep(1, ncol(X)-1)
    
    oldcoef = coef
    if(type == "hmmcov"){
      # need to find way to suppress intercept here
      start_time <- Sys.time()
      # cat(sprintf("class of Znew: %s \n", class(Znew)))
      # cat(sprintf("class of as.matrix(Znew): %s \n", class(as.matrix(Znew))))
      fit0 = grpreg(Znew, y[rep(1:nrow(X), each = nrow(u))], group=covgroup, 
                    penalty="grMCP", family="binomial",lambda = lambda0, 
                    offset = X[rep(1:nrow(X), each = nrow(u)),] %*% matrix(coef[1:ncol(X)],ncol = 1), alpha = alpha, active = active0, 
                    initbeta = c(0,coef[-c(1:ncol(X))]))
      end_time <- Sys.time()
      print(start_time - end_time)
      gc()
      coef = rep(0,length(covgroup) + ncol(X))
      coef[-c(1:ncol(X))] = fit0$beta[-1]
      c0[-c(1:ncol(X))] = c0[-c(1:ncol(X))] + (fit0$beta[-1] == 0)^2
      start_time <- Sys.time()
      fit1 = grpreg(X[rep(1:nrow(X), each = nrow(u)),-1], y[rep(1:nrow(X), each = nrow(u))], group=1:(ncol(X)-1), penalty="grMCP", family="binomial",lambda = lambda1, offset = bigmemory::as.matrix(Znew %*% matrix(coef[-c(1:ncol(X))],ncol = 1)), alpha = alpha, active = active1, initbeta = coef[c(1:ncol(X))])
      end_time <- Sys.time()
      print(start_time - end_time)
      gc()
      coef[c(1:ncol(X))] = fit1$beta
      c0[c(1:ncol(X))] = c0[c(1:ncol(X))] + (fit1$beta == 0)^2
      fit = fit1
      fit$coef = coef
      
      # need to compile code first before running.  Actives will default to 1 to test, then uncomment the below to skip groups
      # update active set every 5 iterations
      if(floor(i/5) == ceiling(i/5)){
        active1[which(c0[c(2:ncol(X))] == 5)] = 0
        active1[which(c0[c(2:ncol(X))] < 5)] = 1
        
        for(kk in 1:max(covgroup)){
          active0[kk] = (all(c0[-c(1:ncol(X))][covgroup == kk]<5))^2
        }
        print(length(covgroup))
        print(length(c0[-c(1:ncol(X))]))
        print(active1)
        print(active0)
        # reset c0
        c0 = rep(0, length(coef))
      }
      
    }else{
      y2 = y[rep(1:nrow(X), each = nrow(u))]
      X2 =  cbind(as.matrix(X[rep(1:nrow(X), each = nrow(u)),-1]), Znew)
      #off2 = rowSums(Z[rep(1:nrow(X), each = nrow(u)),] *   u[rep(1:nrow(u), nrow(X)),])
      fit = glm(y2 ~ X2, family = binomial())
      #fit = ncvreg(y = y[rep(1:nrow(X), each = nrow(u))], X = X[rep(1:nrow(X), each = nrow(u)),-1], family = "binomial", alpha = 1, lambda = lambda,
      #             offset = rowSums(Z[rep(1:nrow(X), each = nrow(u)),] *   u[rep(1:nrow(u), nrow(X)),]), penalty = "SCAD" )
      coef = as.numeric(fit$coef)#as.numeric(fit$beta)
      fit$coef = coef
      fit$w2use = which(coef[-1] != 0)
    }
    
    problem = F
    if(any(is.na(coef))){
      problem = T
      ll = Inf
    }
    
    u2 = matrix(0, nMC2, ncol(Z))
    for(ii in 1:d){
      u2[,seq(ii, ncol(Z), by = d)] = rmvnorm(n = nMC2,sigma=var)
    }
    etae = as.numeric(X[rep(1:nrow(X), each = nrow(u)),] %*% matrix(coef[1:ncol(X)],ncol = 1) + Znew %*% matrix(coef[-c(1:ncol(X))],ncol = 1))
    etae2 = X %*% matrix(coef[1:ncol(X)],nrow = ncol(X), ncol = nrow(u2)) + Z %*% t(u2)
    if(length(etae) != length(y)*nMC2){
      print(head(etae))
      print(dim(etae))
    }
    
    if(family == "poisson"){
      q = apply(etae, 2, FUN = function(etaei) sum(dpois(y, lambda =  exp(etaei), log = T)))
      ll = sum(rowMeans(dpois(matrix(y, nrow = length(y), ncol = ncol(etae)), lambda =  exp(etae), log = T)))
    }else if(family == "binomial"){
      q = apply(matrix(dbinom( y[rep(1:nrow(X), each = nrow(u))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T), ncol = nrow(u), byrow = T), 2, FUN = function(etaei) sum(dbinom(y, size = 1, prob = exp(etaei)/(1+exp(etaei)), log = T)))     
      #ll = sum(rowMeans(dbinom(matrix(y, nrow = length(y), ncol = ncol(etae)), size = 1, prob = exp(etae)/(1+exp(etae)), log = T)))
      
      ll = (sum((dbinom(y[rep(1:nrow(X), each = nrow(u))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T))) + sum(dmvnorm(u, log = T)))/nrow(u) # calc of norm is fine since cov = I
      ll0 = (sum((dbinom(y[rep(1:nrow(X), each = nrow(u))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T))))/nrow(u)
      ll20 =  sum(log(rowMeans(dbinom(matrix(y, nrow = length(y), ncol = ncol(etae2)), size = 1, prob = exp(etae2)/(1+exp(etae2)), log = F))))
      ## change this 5/23 12:36 AM to check the issue
      #ll = ll0 #commented out to try fixed grpreg function, try again later
      ##
    }  
    
    if(!is.finite(ll)){
      problem = T
      ll = Inf
      print(coef)
    }
    
    if(problem == T){
      stop("Error in M step")
      if(is.null(ufull)){
        BIC = -2*ll+ log(length(y))*sum(d) 
      }else{
        rm(Znew)
        gc()
        Znew = big.matrix(nrow = nrow(X)*nrow(ufull), ncol = ncol(J))
        Znew_gen2(ufull, Z, group, seq(as.numeric(group[1]), ncol(Z), by = d),nrow(Z),ncol(Z)/d,d, Znew@address, J)
        etae = as.numeric(X[rep(1:nrow(X), each = nrow(ufull)),] %*% matrix(coef[1:ncol(X)],ncol = 1) + Znew %*% matrix(coef[-c(1:ncol(X))],ncol = 1))
        if(!is.finite(ll)){
          ll2 = ll ## 11/28 changed to ll0 from 11
        }else{
          BIC = -2*ll + log(d)*sum(coef != 0)  
        }
      }
      out = list(fit = fit, coef = coef, sigma = cov, BIC = BIC, glmer = 
                   fit_glmer, ll = ll, ll0 = ll0,lambda0 = lambda0, lambda1 = lambda1, 
                 fit00 = fit00, covgroup = covgroup, J = J)
      if(returnMC == T) out$u = u0
      return(out)
    }
    if(trace == 1) print(coef)
    
    
    gamma = matrix(J%*%coef[-c(1:ncol(X))], ncol = ncol(Z)/d)
    cov = var = gamma %*% t(gamma)
    
    ok = which(colSums(cov)> 0) #& coef[1:ncol(X)] != 0)
    if(length(ok) == 0) ok = 1 # at least include the random intercept
    okindex = NULL
    for(j in 1:(ncol(Z)/d)){
      if(j %in% ok){
        okindex = c(okindex, (j-1)*d + 1:d) 
      }
    }
    
    Znew2 = Z
    finish = 0
    while(finish == 0){
      for(j in 1:d){
        #print(Z[group == i,seq(i, ncol(Z), by = d)]%*% gamma)
        Znew2[group == j,seq(j, ncol(Z), by = d)] = Z[group == j,seq(j, ncol(Z), by = d)]%*% gamma
      }
      if(!any(is.na(Znew2))) finish = 1
    }
    #print(fit$coef)
    #print(sqrt(cov))
    
    ## increase nMC?
    
    
    #stopping rule
    diff[i] = abs(ll0 - oldll)/abs(ll0) ## if need to change back later update all ll0's for convergence in script to ll
    #diff[i] = max(abs(coef - oldcoef)/(abs(oldcoef)+.1))
    
    if( sum(diff[i:max(i-2, 1)] < conv) >=3 ) break
    
    # if current q is within 95% emp CI of old q, increase nMC
    #if(mean(q) > lim[1] & mean(q) < lim[2]) nMC = min(round(nMC + nMC/3), 10000)
    if(diff[i] > 10^-10){
      nMC = max(round(nMC + min(.15*0.001*nMC/(abs(ll0 - oldll)/abs(ll0)), 250)), 2)+1
    }else{
      nMC = max(round(nMC + .25*0.001*nMC), 5)+1
    }
    
    if(nMC > nMC_max) nMC = nMC_max
    # now update limits  
    print(c(i, nMC , diff[i], ll0, oldll,  ll0 - oldll, sum(coef!=0), coef[2], sqrt(diag(cov)[2])))
    lim = quantile(q, c(0.025, 0.975))
    
    ### E stepf
    #if(any(svd(cov)$d < 10^-15)) cov = cov + diag(rep(10^-10, ncol(cov)))
    u = u0 = sample.mc2(fit=fit, cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group, d = d, okindex = okindex, nZ= ncol(Z), gibbs = gibbs, uold = u0)
    #u = bmmat(u)
    nMC2 = nrow(u)
    
    if(any(is.na(u)) | any(colSums(u) == 0)){
      print("E step: hit iteration limit of 10^10 samples, fit likely inadequate")
      if(is.null(ufull)){
        BIC = -2*ll + log(d)*sum(coef != 0) # switched BIC and BIC0 11/28
        BIC0 = -2*ll0 + log(d)*sum(coef != 0)
        BIC20 = -2*ll20 + log(d)*sum(coef!=0)
      }else{
        rm(Znew)
        gc()
        Znew = big.matrix(nrow = nrow(X)*nrow(ufull), ncol = ncol(J))
        Znew_gen2(ufull, Z, group, seq(as.numeric(group[1]), ncol(Z), by = d),nrow(Z),ncol(Z)/d,d, Znew@address, J)
        etae = as.numeric(X[rep(1:nrow(X), each = nrow(ufull)),] %*% matrix(coef[1:ncol(X)],ncol = 1) + Znew %*% matrix(coef[-c(1:ncol(X))],ncol = 1))
        ll2 = (sum((dbinom(y[rep(1:nrow(X), each = nrow(ufull))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T))) + sum(dmvnorm(ufull,log = T)))/nrow(ufull)
        BIC = -2*ll2 + log(d)*sum(coef != 0)
        ll20b = (sum((dbinom(y[rep(1:nrow(X), each = nrow(ufull))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T))))/nrow(ufull)
        BIC0 = -2*ll20b + log(d)*sum(coef != 0)
        
        BIC20 = -2*ll20 + log(d)*sum(coef != 0)
        #BIC20 is already computed
      }
      out = list(fit = fit, coef = coef, sigma = cov, BIC = BIC, glmer = 
                   fit_glmer, ll = ll, ll0 = ll0, ll2=ll2, ll20b=ll20b,lambda0 = lambda0, 
                 lambda1 = lambda1, fit00 = fit00, BIC0 = BIC0, BIC20 = BIC20, covgroup 
                 = covgroup, J = J)
      if(returnMC == T) out$u = u0
      return(out)
    }
    
    if(trace == 1) print(diag(cov))
    gc()
  }
  
  
  
  ## calculate BIC 
  if(is.null(ufull)){
    BIC = -2*ll0 + log(length(y))*sum(coef != 0) # switched BIC and BIC0 11/28
    BIC0 = -2*ll + log(length(y))*sum(coef != 0)
    BIC20 = -2*ll20 + log(length(y))*sum(coef!=0)
    llb = ll0b = 0
  }else{
    rm(Znew)
    gc()
    BIC20 = -2*ll20  + log(length(y))*sum(coef!=0)
    Znew = big.matrix(nrow = nrow(X)*nrow(ufull), ncol = ncol(J))
    Znew_gen2(ufull, Z, group, seq(as.numeric(group[1]), ncol(Z), by = d),nrow(Z),ncol(Z)/d,d, Znew@address, J)
    etae = as.numeric(X[rep(1:nrow(X), each = nrow(ufull)),] %*% matrix(coef[1:ncol(X)],ncol = 1) + Znew %*% matrix(coef[-c(1:ncol(X))],ncol = 1))
    llb = (sum((dbinom(y[rep(1:nrow(X), each = nrow(ufull))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T))) + sum(dnorm(ufull, 0,1, log = T)))/nrow(ufull)
    BIC = -2*llb + log(length(y))*sum(coef != 0)
    ll0b = (sum((dbinom(y[rep(1:nrow(X), each = nrow(ufull))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T))))/nrow(ufull)
    BIC0 = -2*ll0b + log(length(y))*sum(coef != 0)
    rm(Znew)
    gc()
  }
  
  if(ncol(X) - 1 == p1 | (lambda0 == 0 & lambda1 == 0) ){
    fin =   cbind(summary(fit_glmer)$coefficients[-1,1], coef[2:ncol(X)])
    colnames(fin) = c("glmer_full","MCEM_pen")
  }else{
    fin =   cbind(summary(fit_glmer)$coefficients[-1,1], c(summary(fit_glmer_oracle)$coef[-1,1], rep(0, ncol(X)-p1-1)),coef[2:ncol(X)])
    colnames(fin) = c("glmer_full","glmer_oracle","MCEM_pen")
  }
  rownames(fin) = NULL
  print(fin)
  print(sqrt(diag(cov)[1:3]))
  vc <- VarCorr(fit_glmer_oracle)
  print(summary(fit_glmer)$coefficients[-1,1])
  returnMC
  out = list(fit = fit, coef = coef, sigma = cov, BIC = BIC, glmer = 
               fit_glmer, ll = ll, ll0 = ll0, llb = llb, ll0b = ll0b, ll20 = ll20, 
             lambda0 = lambda0, lambda1 = lambda1, fit00 = fit00, BIC0 = BIC0, BIC = 
               BIC20, covgroup = covgroup, J = J)
  if(returnMC == T) out$u = u
  return(out)
}
