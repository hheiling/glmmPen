
#' Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional
#' Minimization (MCECM)
#' 
#' Description
#' 
#' @inheritParams fit_dat_B
#' 
#' @section Details:
#' Accepted families: binomial 
#' 
#' 
# @export
# fit_dat = function(dat,  lambda0 = 0, lambda1 = 0, conv = 0.001, nMC = 1000, 
#                    family = "binomial", trace = 0, penalty = "grMCP",
#                    alpha = 1, nMC_max = 5000, t = 10,
#                    returnMC = T, ufull = NULL, coeffull = NULL, gibbs = T, maxitEM = 100, 
#                    ufullinit = NULL, M = 10^4, MwG_sampler = c("random_walk","independence"),
#                    adapt_RW_options = adaptControl(), covar = c("unstructured","independent")){
#   
#   # Things to address:
#   ## Eventually, delete this line and following 'ok' references: ok = which(diag(var) > 0)
#   
#   # Set small penalties to zero
#   if(lambda0 <=10^-6) lambda0 = 0
#   if(lambda1 <=10^-6) lambda1 = 0
#   
#   y = dat$y
#   X = as.matrix(dat$X)
#   # Convert sparse Z to dense Z
#   Z = Matrix::as.matrix(dat$Z)
#   group = dat$group
#   
#   if(is.character(family)){
#     family = get(family, mode = "function", envir = parent.frame())
#   }
#   if(is.function(family)){
#     family = family()
#   }
#   if(class(family) == "family"){
#     f = family
#     link = family$link
#     family = family$family
#   }
#   
#   d = nlevels(factor(group))
#   
#   covar = covar[1]
#   if(!(covar %in% c("unstructured","independent"))){
#     stop("algorithm currently only handles 'unstructured' or 'independent' covariance structure \n")
#   }
#   if(covar == "unstructured" & ncol(Z)/d >= 10){
#     warning("Due to dimension of sigma covariance matrix, will use covar = 'independent' to simplify computation \n",
#             immediate. = T)
#     covar == "independent"
#   }
#   
#   #initial fit
#   if(family == "binomial"){
#     nTotal = rep(1, length(y))
#   }else{
#     nTotal = NULL
#   }
#   
#   initial_gibbs = gibbs
#   
#   MwG_sampler = MwG_sampler[1] # Default of random walk
#   if(!(MwG_sampler %in% c("independence", "random_walk"))){
#     stop("MwG_sampler must be specified as either 'independence' or 'random_walk'")
#   }
#   
#   if(covar == "unstructured"){ # Originally: ncol(Z)/d <= 15 
#     # create J, q2 x q*(q+1)/2
#     J = Matrix(0, (ncol(Z)/d)^2, (ncol(Z)/d)*((ncol(Z)/d)+1)/2, sparse = T) #matrix(0, (ncol(Z)/d)^2, (ncol(Z)/d)*((ncol(Z)/d)+1)/2)
#     index = 0
#     indexc = 0
#     sumy = 0
#     sumx = 0
#     zeros = 0
#     covgroup = NULL
#     for(i in 1:(ncol(Z)/d)){
#       J[ sumx + zeros + 1:(ncol(Z)/d - (i-1)), sumy + 1:(ncol(Z)/d - (i-1))] = diag((ncol(Z)/d - (i-1)))
#       sumy = sumy + (ncol(Z)/d - (i-1))
#       sumx = sumy 
#       zeros = zeros + i
#       covgroup = rbind(covgroup, rep(i, (ncol(Z)/d)))
#     }
#     covgroup = covgroup[lower.tri(covgroup, diag = T)]
#   }else{ # covar == "independent". Originally: ncol(Z)/d > 15
#     J = Matrix(0,(ncol(Z)/d)^2, (ncol(Z)/d), sparse = T) #matrix(0, (ncol(Z)/d)^2, (ncol(Z)/d))
#     index = 0
#     indexc = 0
#     sumy = 0
#     sumx = 0
#     zeros = 0
#     covgroup = NULL
#     for(i in 1:(ncol(Z)/d)){
#       J[ sumx + zeros + 1, i] = 1
#       sumy = sumy + (ncol(Z)/d - (i-1))
#       sumx = sumy 
#       zeros = zeros + i
#     }
#     covgroup = rep(1:(ncol(Z)/d))
#   }
#   
#   if(!is.null(ufull) & !is.null(coeffull)){
#     fit = list()
#     print("using coef from full model to intialize")
#     coef = coeffull
#     gamma = matrix(J%*%matrix(coef[-c(1:ncol(X))], ncol = 1), ncol = ncol(Z)/d)
#     cov = var = gamma %*% t(gamma)
#     fit$coef = coef[c(1:ncol(X))]
#     ok = which(diag(var) > 0)# & coef[1:ncol(X)] != 0)
#     if(length(ok) == 0) ok = 1 # at least include the random intercept
#     okindex = NULL
#     for(i in 1:(ncol(Z)/d)){
#       if(i %in% ok){
#         okindex = c(okindex, (i-1)*d + 1:d)
#       }
#     }
#     fit00 = fit
#   }else{
#     if(ncol(X) > 2){
#       fit = grpreg(X[,-1], y, group=1:(ncol(X)-1), penalty = penalty, family=family,lambda = lambda1, alpha = alpha)###
#     }else{
#       fit = grpreg(matrix(X[,-1], nrow = nrow(X)), y, group=1:(ncol(X)-1), penalty = penalty, family=family,lambda = lambda1, alpha = alpha)###
#     }
#     
#     fit00 = fit # naive fit
#     
#     coef = as.numeric(fit$beta)
#     fit$coef = as.numeric(fit$beta)
#     
#     if(trace == 1) print(coef)
#     
#     if(ncol(Z)/d > 1){
#       vars = rep(10^-10, ncol(Z)/d)
#       cov = var = diag(vars)
#       gamma = t(chol(var)) # chol outputs upper triangular, so transposing here
#     }else{
#       vars = 10^-10
#       cov = var = matrix(vars, ncol = 1)
#       gamma = var
#     }
#     
#     
#     ok = which(vars > 0)# & coef[1:ncol(X)] != 0)
#     if(length(ok) == 0) ok = 1 # at least include the random intercept
#     okindex = NULL
#     for(i in 1:(ncol(Z)/d)){
#       if(i %in% ok){
#         okindex = c(okindex, (i-1)*d + 1:d) 
#       }
#     }
#   }
#   
#   
#   Znew2 = Z
#   finish = 0
#   while(finish == 0){
#     for(i in 1:d){
#       Znew2[group == i,seq(i, ncol(Z), by = d)] = Z[group == i, seq(i, ncol(Z), by = d)] %*% gamma
#     }
#     if(!any(is.na(Znew2))) finish = 1
#   }
#   
#   # intitialize switch-from-rejection-sampling-to-gibbs-sampling counter
#   rej_to_gibbs = 0
#   
#   # initialize adaptive Metropolis-within-Gibbs random walk parameters
#   # ignored if use rejection sampling (gibbs = F), but use if gibbs = T or
#   # use if initially rejection sampling but switch to gibbs = T
#   ## initialize proposal standard deviation
#   proposal_SD = matrix(1.0, nrow = d, ncol = ncol(Z)/d)
#   print("Initialized proposal_SD")
#   print(proposal_SD)
#   
#   ## initialize batch number to 0
#   batch = 0.0
#   ## initialize other paramters from adaptControl()
#   batch_length = adapt_RW_options$batch_length
#   offset = adapt_RW_options$offset
#   burnin_batchnum = adapt_RW_options$burnin_batchnum
#   gibbs_accept_rate = matrix(NA, nrow = d, ncol = nrow(Z)/d)
#   
#   if((!is.null(ufull) | !is.null(ufullinit)) & !is.null(coeffull)){
#     if(!is.null(ufullinit)){
#       print("using u from previous model to intialize")
#     }else{
#       print("using u from full model to intialize")
#       ufullinit = ufull
#     }
#     
#     if(MwG_sampler == "independence"){
#       samplemc_out = sample.mc2(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group, 
#                                 d = d, okindex = okindex, trace = trace, gibbs = gibbs, uold = ufullinit)
#     }else{ # MwG_sampler == "random_walk"
#       samplemc_out = sample_mc_adapt(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group,
#                                      d = d, okindex = okindex, trace = trace, gibbs = gibbs, uold = ufullinit,
#                                      proposal_SD = proposal_SD, batch = batch, batch_length = batch_length,
#                                      offset = offset, burnin_batchnum = burnin_batchnum)
#     }
#     
#     u = u0 = samplemc_out$u0
#     
#     # If specified gibbs = T or if specified gibbs = F but switched to gibbs due to low acceptance rates
#     if(gibbs | samplemc_out$switch){
#       # If rejection sampling and switched to gibbs sampling due to low acceptance rate:
#       if(samplemc_out$switch){ 
#         rej_to_gibbs = rej_to_gibbs + 1
#         cat("rej_to_gibbs count: ", rej_to_gibbs, "\n")
#       }
#       
#       if(MwG_sampler == "random_walk"){
#         gibbs_accept_rate = samplemc_out$gibbs_accept_rate
#         batch = samplemc_out$updated_batch
#         proposal_SD = samplemc_out$proposal_SD
#         
#         print("Updated proposal_SD:")
#         print(proposal_SD)
#       }
#       
#     }
#     
#     
#     
#   }else{
#     
#     if(MwG_sampler == "independence"){
#       samplemc_out = sample.mc2(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group, 
#                                 d = d, okindex = okindex, trace = trace, gibbs = gibbs, 
#                                 uold = matrix(rnorm(nMC*ncol(Z)), nrow = nMC, ncol = ncol(Z)))
#     }else{ # MwG_sampler == "random_walk"
#       samplemc_out = sample_mc_adapt(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group,
#                                      d = d, okindex = okindex, trace = trace, gibbs = gibbs,
#                                      uold = matrix(rnorm(nMC*ncol(Z)), nrow = nMC, ncol = ncol(Z)),
#                                      proposal_SD = proposal_SD, batch = batch, batch_length = batch_length,
#                                      offset = offset, burnin_batchnum = burnin_batchnum)
#     }
#     
#     u = u0 = samplemc_out$u0
#     
#     # If specified gibbs = T or if specified gibbs = F but switched to gibbs due to low acceptance rates
#     if(gibbs | samplemc_out$switch){
#       # If rejection sampling and switched to gibbs sampling due to low acceptance rate:
#       if(samplemc_out$switch){ 
#         rej_to_gibbs = rej_to_gibbs + 1
#         cat("rej_to_gibbs count: ", rej_to_gibbs, "\n")
#       }
#       
#       if(MwG_sampler == "random_walk"){
#         gibbs_accept_rate = samplemc_out$gibbs_accept_rate
#         batch = samplemc_out$updated_batch
#         proposal_SD = samplemc_out$proposal_SD
#         
#         print("Updated proposal_SD:")
#         print(proposal_SD)
#       }
#       
#     }
#     
#   }
#   #u = bmmat(u)
#   nMC2 = nrow(u)  
#   etae = matrix(X %*% coef[1:ncol(X)], nrow = nrow(X), ncol = nrow(u) ) + Znew2%*%t(u)
#   
#   if(nrow(cov) == 1){ # Single random intercept
#     cov_record = rep(NA, maxitEM)
#   }else{
#     cov_record = NULL
#   }
#   
#   fit0_record = matrix(NA, nrow = maxitEM, ncol = (length(covgroup) + 1))
#   diff = rep(NA, maxitEM)
#   stopcount = 0
#   
#   ## this needs to get updated to reflect whatever is chosen to evaluate likelihood
#   if(family == "poisson"){
#     ll = ll0 = ll20 = sum(rowMeans(dpois(matrix(y, nrow = length(y), ncol = ncol(etae)), lambda =  exp(etae), log = T)))
#   }else if(family == "binomial"){
#     ll = ll0 = ll20  = sum(rowMeans(dbinom(matrix(y, nrow = length(y), ncol = ncol(etae)), size = 1, prob = exp(etae)/(1+exp(etae)), log = T)))
#   }  
#   Znew = NULL
#   
#   # initialize zero count vectors
#   c0 = rep(0, length(coef))
#   
#   # Record last t coef vectors (each row = coef vector for a past EM iteration)
#   # Initialize with initial coef vector
#   coef = c(coef, rep(0, length(covgroup)))
#   coef_record = matrix(coef, nrow = t, ncol = length(coef), byrow = T)
#   coef_record_all = matrix(NA, nrow = maxitEM, ncol = length(coef), byrow = T)
#   
#   # Start EM Algorithm (M step first)
#   
#   for(i in 1:maxitEM){
#     
#     if(rej_to_gibbs == 3){
#       gibbs = T
#       cat("permanently switched from rejection sampling to gibbs sampling \n")
#       rej_to_gibbs = rej_to_gibbs + 1
#     }
#     
#     oldll = ll0
#     
#     if(family == "binomial"){
#       nTotal = rep(1, length(y[rep(1:nrow(X), each = nrow(u))]))
#     }else{
#       nTotal = NULL
#     }
#     
#     print("Znewgen done")
#     rm(Znew)
#     gc()
#     Znew = big.matrix(nrow = nrow(X)*nrow(u), ncol = ncol(J))
#     Znew_gen2(u, Z, group, seq(as.numeric(group[1]), ncol(Z), by = d),nrow(Z),ncol(Z)/d,d, Znew@address, J)
#     gc()
#     
#     active0 = rep(1, max(covgroup))
#     active1 = rep(1, ncol(X)-1)
#     
#     # oldcoef = coef
#     
#     # M Step
#     
#     fit0 = grpreg(Znew, y[rep(1:nrow(X), each = nrow(u))], group=covgroup, 
#                   penalty="grMCP", family="binomial",lambda = lambda1, 
#                   offset = X[rep(1:nrow(X), each = nrow(u)),] %*% matrix(coef[1:ncol(X)],ncol = 1), 
#                   alpha = alpha, active = active0, 
#                   initbeta = c(0,coef[-c(1:ncol(X))]))
#     gc()
#     coef = rep(0,length(covgroup) + ncol(X))
#     coef[-c(1:ncol(X))] = fit0$beta[-1]
#     c0[-c(1:ncol(X))] = c0[-c(1:ncol(X))] + (fit0$beta[-1] == 0)^2
#     
#     cat("full fit0$beta output: ", fit0$beta, "\n")
#     fit0_record[i,] = fit0$beta
#     
#     if(ncol(X) > 2){
#       fit1 = grpreg(X[rep(1:nrow(X), each = nrow(u)),-1], y[rep(1:nrow(X), each = nrow(u))], 
#                     group=1:(ncol(X)-1), penalty="grMCP", family="binomial",lambda = lambda0, 
#                     offset = bigmemory::as.matrix(Znew %*% matrix(coef[-c(1:ncol(X))],ncol = 1)), 
#                     alpha = alpha, active = active1, initbeta = coef[c(1:ncol(X))])
#     }else{
#       fit1 = grpreg(matrix(X[rep(1:nrow(X), each = nrow(u)),-1], nrow = nrow(X)*nrow(u)), 
#                     y[rep(1:nrow(X), each = nrow(u))], group=1:(ncol(X)-1), penalty="grMCP", 
#                     family="binomial",lambda = lambda0, 
#                     offset = bigmemory::as.matrix(Znew %*% matrix(coef[-c(1:ncol(X))],ncol = 1)), 
#                     alpha = alpha, active = active1, initbeta = coef[c(1:ncol(X))])
#     }
#     gc()
#     coef[c(1:ncol(X))] = fit1$beta
#     c0[c(1:ncol(X))] = c0[c(1:ncol(X))] + (fit1$beta == 0)^2
#     fit = fit1
#     fit$coef = coef
#     
#     # need to compile code first before running.  Actives will default to 1 to test, then uncomment the below to skip groups
#     # update active set every 5 iterations
#     if(floor(i/5) == ceiling(i/5)){
#       active1[which(c0[c(2:ncol(X))] == 5)] = 0
#       active1[which(c0[c(2:ncol(X))] < 5)] = 1
#       
#       for(kk in 1:max(covgroup)){
#         active0[kk] = (all(c0[-c(1:ncol(X))][covgroup == kk]<5))^2
#       }
#       print(length(covgroup))
#       print(length(c0[-c(1:ncol(X))]))
#       print(active1)
#       print(active0)
#       # reset c0
#       c0 = rep(0, length(coef))
#     }
#     
#     problem = F
#     if(any(is.na(coef))){
#       problem = T
#       ll = Inf
#     }
#     
#     u2 = matrix(0, nMC2, ncol(Z))
#     for(ii in 1:d){
#       u2[,seq(ii, ncol(Z), by = d)] = rmvnorm(n = nMC2,sigma=var)
#     }
#     etae = as.numeric(X[rep(1:nrow(X), each = nrow(u)),] %*% matrix(coef[1:ncol(X)],ncol = 1) + Znew %*% matrix(coef[-c(1:ncol(X))],ncol = 1))
#     etae2 = X %*% matrix(coef[1:ncol(X)],nrow = ncol(X), ncol = nrow(u2)) + Z %*% t(u2)
#     if(length(etae) != length(y)*nMC2){
#       print(head(etae))
#       print(dim(etae))
#     }
#     
#     if(family == "poisson"){
#       q = apply(etae, 2, FUN = function(etaei) sum(dpois(y, lambda =  exp(etaei), log = T)))
#       ll = sum(rowMeans(dpois(matrix(y, nrow = length(y), ncol = ncol(etae)), lambda =  exp(etae), log = T)))
#     }else if(family == "binomial"){
#       q = apply(matrix(dbinom( y[rep(1:nrow(X), each = nrow(u))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T), ncol = nrow(u), byrow = T), 2, FUN = function(etaei) sum(dbinom(y, size = 1, prob = exp(etaei)/(1+exp(etaei)), log = T)))     
#       
#       ll = (sum((dbinom(y[rep(1:nrow(X), each = nrow(u))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T))) + sum(dmvnorm(u, log = T)))/nrow(u) # calc of norm is fine since cov = I
#       ll0 = (sum((dbinom(y[rep(1:nrow(X), each = nrow(u))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T))))/nrow(u)
#       ll20 =  sum(log(rowMeans(dbinom(matrix(y, nrow = length(y), ncol = ncol(etae2)), size = 1, prob = exp(etae2)/(1+exp(etae2)), log = F))))
#       
#     }  
#     
#     if(!is.finite(ll)){
#       problem = T
#       ll = Inf
#       print(coef)
#     }
#     
#     if(problem == T){
#       stop("Error in M step")
#       if(is.null(ufull)){
#         BIC = -2*ll+ log(length(y))*sum(d) 
#       }else{
#         rm(Znew)
#         gc()
#         Znew = big.matrix(nrow = nrow(X)*nrow(ufull), ncol = ncol(J))
#         Znew_gen2(ufull, Z, group, seq(as.numeric(group[1]), ncol(Z), by = d),nrow(Z),ncol(Z)/d,d, Znew@address, J)
#         etae = as.numeric(X[rep(1:nrow(X), each = nrow(ufull)),] %*% matrix(coef[1:ncol(X)],ncol = 1) + Znew %*% matrix(coef[-c(1:ncol(X))],ncol = 1))
#         if(!is.finite(ll)){
#           ll2 = ll
#         }else{
#           BIC = -2*ll + log(d)*sum(coef != 0)  
#         }
#       }
#       out = list(fit = fit, coef = coef, sigma = cov, BIC = BIC,  
#                  ll = ll, ll0 = ll0,lambda0 = lambda0, lambda1 = lambda1, 
#                  fit00 = fit00, covgroup = covgroup, J = J)
#       if(returnMC == T) out$u = u0
#       return(out)
#     }
#     if(trace == 1) print(coef)
#     
#     
#     gamma = matrix(J%*%matrix(coef[-c(1:ncol(X))], ncol = 1), ncol = ncol(Z)/d)
#     cov = var = gamma %*% t(gamma)
#     
#     ok = which(colSums(cov)> 0) #& coef[1:ncol(X)] != 0)
#     if(length(ok) == 0) ok = 1 # at least include the random intercept
#     okindex = NULL
#     for(j in 1:(ncol(Z)/d)){
#       if(j %in% ok){
#         okindex = c(okindex, (j-1)*d + 1:d) 
#       }
#     }
#     
#     Znew2 = Z
#     finish = 0
#     while(finish == 0){
#       for(j in 1:d){
#         Znew2[group == j,seq(j, ncol(Z), by = d)] = Z[group == j,seq(j, ncol(Z), by = d)]%*% gamma
#       }
#       if(!any(is.na(Znew2))) finish = 1
#     }
#     
#     # stopping rule
#     # diff[i] = abs(ll0 - oldll)/abs(ll0) ## if need to change back later update all ll0's for convergence in script to ll
#     # stopping rule: based on average Euclidean distance (comparing coef from minus t iterations)
#     if(i <= t){
#       diff[i] = 10^2
#     }else{
#       diff[i] = sqrt(sum((coef - coef_record[1,])^2)) / length(coef)
#     }
#     
#     # Update latest record of coef
#     coef_record = rbind(coef_record[-1,], coef)
#     coef_record_all[i,] = coef
#     
#     if( sum(diff[i:max(i-2, 1)] < conv) >=3 ) break
#     
#     # if current q is within 95% emp CI of old q, increase nMC
#     #if(mean(q) > lim[1] & mean(q) < lim[2]) nMC = min(round(nMC + nMC/3), 10000)
#     if(diff[i] > 10^-10){
#       nMC = max(round(nMC + min(.15*0.001*nMC/(abs(ll0 - oldll)/abs(ll0)), 250)), 2)+1
#     }else{
#       nMC = max(round(nMC + .25*0.001*nMC), 5)+1
#     }
#     
#     if(nMC > nMC_max) nMC = nMC_max
#     # now update limits  
#     print(c(i, nMC , diff[i], ll0, oldll,  ll0 - oldll, sum(coef!=0), coef[2], sqrt(diag(cov)[2])))
#     lim = quantile(q, c(0.025, 0.975))
#     
#     print("cov:")
#     print(cov)
#     if(nrow(cov) == 1){
#       cov_record[i] = cov
#     }
#     
#     # E Step 
#     
#     if(MwG_sampler == "independence"){
#       samplemc_out = sample.mc2(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, trace = trace, family = family, group = group, 
#                                 d = d, okindex = okindex, nZ = ncol(Z), gibbs = gibbs, uold = u0)
#     }else{ # MwG_sampler == "random_walk"
#       samplemc_out = sample_mc_adapt(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, trace = trace, family = family, group = group, 
#                                      d = d, okindex = okindex, gibbs = gibbs, uold = u0,
#                                      proposal_SD = proposal_SD, batch = batch, batch_length = batch_length, 
#                                      offset = offset, burnin_batchnum = burnin_batchnum)
#     }
#     
#     u = u0 = samplemc_out$u0
#     
#     # If specified gibbs = T or if specified gibbs = F but switched to gibbs due to low acceptance rates
#     if(gibbs | samplemc_out$switch){
#       # If rejection sampling and switched to gibbs sampling due to low acceptance rate:
#       if(samplemc_out$switch){ 
#         rej_to_gibbs = rej_to_gibbs + 1
#         cat("rej_to_gibbs count: ", rej_to_gibbs, "\n")
#       }
#       
#       if(MwG_sampler == "random_walk"){
#         gibbs_accept_rate = samplemc_out$gibbs_accept_rate
#         batch = samplemc_out$updated_batch
#         proposal_SD = samplemc_out$proposal_SD
#         
#         print("Updated proposal_SD:")
#         print(proposal_SD)
#         print("Updated batch:")
#         print(batch)
#       }
#       
#     }
#     
#     
#     
#     nMC2 = nrow(u)
#     
#     if(any(is.na(u)) | any(colSums(u) == 0)){
#       print("E step: hit iteration limit of 10^10 samples, fit likely inadequate")
#       if(is.null(ufull)){
#         BIC = -2*ll + log(d)*sum(coef != 0) # switched BIC and BIC0 11/28
#         BIC0 = -2*ll0 + log(d)*sum(coef != 0)
#         BIC20 = -2*ll20 + log(d)*sum(coef!=0)
#       }else{
#         rm(Znew)
#         gc()
#         Znew = big.matrix(nrow = nrow(X)*nrow(ufull), ncol = ncol(J))
#         Znew_gen2(ufull, Z, group, seq(as.numeric(group[1]), ncol(Z), by = d),nrow(Z),ncol(Z)/d,d, Znew@address, J)
#         etae = as.numeric(X[rep(1:nrow(X), each = nrow(ufull)),] %*% matrix(coef[1:ncol(X)],ncol = 1) + Znew %*% matrix(coef[-c(1:ncol(X))],ncol = 1))
#         ll2 = (sum((dbinom(y[rep(1:nrow(X), each = nrow(ufull))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T))) + sum(dmvnorm(ufull,log = T)))/nrow(ufull)
#         BIC = -2*ll2 + log(d)*sum(coef != 0)
#         ll20b = (sum((dbinom(y[rep(1:nrow(X), each = nrow(ufull))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T))))/nrow(ufull)
#         BIC0 = -2*ll20b + log(d)*sum(coef != 0)
#         
#         BIC20 = -2*ll20 + log(d)*sum(coef != 0)
#         #BIC20 is already computed
#       }
#       out = list(fit = fit, coef = coef, sigma = cov, BIC = BIC, 
#                  ll = ll, ll0 = ll0, ll2=ll2, ll20b=ll20b,lambda0 = lambda0, 
#                  lambda1 = lambda1, fit00 = fit00, BIC0 = BIC0, BIC20 = BIC20, covgroup 
#                  = covgroup, J = J)
#       if(returnMC == T) out$u = u0
#       return(out)
#     }
#     
#     if(trace == 1) print(diag(cov))
#     gc()
#   }
#   
#   ## calculate BIC 
#   # if(is.null(ufull)){
#   #   BIC = -2*ll0 + log(length(y))*sum(coef != 0) # switched BIC and BIC0 11/28
#   #   BIC0 = -2*ll + log(length(y))*sum(coef != 0)
#   #   BIC20 = -2*ll20 + log(length(y))*sum(coef!=0)
#   #   llb = ll0b = 0
#   # }else{
#   #   rm(Znew)
#   #   gc()
#   #   BIC20 = -2*ll20  + log(length(y))*sum(coef!=0)
#   #   Znew = big.matrix(nrow = nrow(X)*nrow(ufull), ncol = ncol(J))
#   #   Znew_gen2(ufull, Z, group, seq(as.numeric(group[1]), ncol(Z), by = d),nrow(Z),ncol(Z)/d,d, Znew@address, J)
#   #   etae = as.numeric(X[rep(1:nrow(X), each = nrow(ufull)),] %*% matrix(coef[1:ncol(X)],ncol = 1) + Znew %*% matrix(coef[-c(1:ncol(X))],ncol = 1))
#   #   llb = (sum((dbinom(y[rep(1:nrow(X), each = nrow(ufull))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T))) + sum(dnorm(ufull, 0,1, log = T)))/nrow(ufull)
#   #   BIC = -2*llb + log(length(y))*sum(coef != 0)
#   #   ll0b = (sum((dbinom(y[rep(1:nrow(X), each = nrow(ufull))], size = 1, prob = exp(etae)/(1+exp(etae)), log = T))))/nrow(ufull)
#   #   BIC0 = -2*ll0b + log(length(y))*sum(coef != 0)
#   #   rm(Znew)
#   #   gc()
#   # }
#   
#   
#   print(sqrt(diag(cov)[1:3]))
#   returnMC
#   
#   # Another E step for loglik calculation (number draws = M)
#   samplemc_out = sample_mc_adapt(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=M, trace = trace, family = family, group = group, 
#                                  d = d, okindex = okindex, gibbs = gibbs, uold = u0,
#                                  proposal_SD = proposal_SD, batch = batch, batch_length = batch_length, 
#                                  offset = offset, burnin_batchnum = burnin_batchnum)
#   u = u0 = samplemc_out$u0
#   
#   # If specified gibbs = T or if specified gibbs = F but switched to gibbs due to low acceptance rates
#   if(gibbs | samplemc_out$switch){
#     # If rejection sampling and switched to gibbs sampling due to low acceptance rate:
#     if(samplemc_out$switch){ 
#       rej_to_gibbs = rej_to_gibbs + 1
#       cat("rej_to_gibbs count: ", rej_to_gibbs, "\n")
#     }
#     
#     gibbs_accept_rate = samplemc_out$gibbs_accept_rate
#     batch = samplemc_out$updated_batch
#     proposal_SD = samplemc_out$proposal_SD
#   }
#   
#   print("Updated proposal_SD:")
#   print(proposal_SD)
#   print("Updated batch:")
#   print(batch)
#   
#   # Calculate loglik using Pajor method (see logLik_Pajor.R)
#   if(sum(diag(cov) == 0) == nrow(cov)){
#     if(family == "binomial"){
#       eta = X %*% coef[1:ncol(X)]
#       ll = sum(dbinom(x = y, size = 1, prob = exp(eta)/(1+exp(eta)), log = T))
#     }
#   }else{
#     ll = CAME_IS(posterior = u, y = y, X = X, Z = Z, group = group,
#                  coef = coef, sigma = cov, family = family, M = M)
#   }
#   
#   
#   # Hybrid BIC (Delattre, Lavielle, and Poursat (2014))
#   # d = nlevels(group) = number independent subjects/groups
#   BICh = -2*ll + sum(diag(cov) != 0)*log(d) + sum(coef[1:ncol(X)] != 0)*log(nrow(X))
#   # Usual BIC
#   BIC = -2*ll + sum(coef != 0)*log(nrow(X))
#   
#   if(gibbs){
#     out = list(fit = fit, coef = coef, sigma = cov,  
#                lambda0 = lambda0, lambda1 = lambda1, 
#                covgroup = covgroup, J = J, ll = ll, BICh = BICh, BIC = BIC,
#                extra = list(fit = fit, okindex = okindex, Znew2 = Znew2),
#                gibbs_accept_rate = gibbs_accept_rate, proposal_SD = proposal_SD)
#   }else{
#     out = list(fit = fit, coef = coef, sigma = cov,  
#                lambda0 = lambda0, lambda1 = lambda1, 
#                covgroup = covgroup, J = J, ll = ll, BICh = BICh, BIC = BIC,
#                extra = list(fit = fit, okindex = okindex, Znew2 = Znew2))
#   }
#   
#   if(returnMC == T) out$u = u
#   
#   if((initial_gibbs == F) && rej_to_gibbs > 0){
#     if(rej_to_gibbs <= 3){
#       cat(sprintf("ending rej_to_gibbs count: %i \n", rej_to_gibbs))
#     }else{
#       # To correct for additional rej_to_gibbs + 1 when rej_to_gibbs = 3
#       cat(sprintf("ending rej_to_gibbs count: %i \n", rej_to_gibbs-1))
#     }
#   }
#   
#   if(initial_gibbs == F){
#     out$rej_to_gibbs = rej_to_gibbs
#   }
#   
#   out$fit0_record = fit0_record
#   out$coef_record_all = coef_record_all
#   if(!is.null(cov_record)){
#     out$cov_record = cov_record
#   }
#   
#   return(out)
# }

# Log-likelihood approximation using importance sampling
# @export
# logLik_imp = function(y, X, Z, U, sigma, group, coef, family, df, c = 1, M){
#   
#   # Set-up calculations
#   d = nlevels(group)
#   cols = seq(from = as.numeric(group[1]), to = ncol(Z), by = d)
#   
#   # Ignore columns of U (and Z) corresponding to random effects penalized to 0 variance 
#   U_means_all = colMeans(U)
#   non_zero = (diag(sigma) != 0)
#   
#   if(sum(non_zero) > 0){
#     
#     non_zero_ext = rep(non_zero, each = d)
#     U_means = U_means_all[non_zero_ext]
#     
#     # Reduced sigma: remove rows and columns with diag = 0
#     sigma_red = sigma[non_zero,non_zero]
#     
#     # Gamma = cholesky decomposition of sigma (lower-triangular)
#     Gamma = t(chol(sigma_red))
#     
#     # Calculated fixed effects contribution to eta (linear predictor)
#     eta_fef = X %*% coef[1:ncol(X)]
#     
#     ll = logLik_cpp(U_means, c*sigma_red, M, group, d, df, y, eta_fef, Z[,non_zero_ext], 
#                     Gamma, family)
#   }else{
#     
#     eta_fef = X %*% coef[1:ncol(X)]
#     
#     if(family == "binomial"){
#       ll = sum(dbinom(y, size = 1, prob = exp(eta_fef) / (1+exp(eta_fef)), log = T))
#     }else if(family == "poisson"){
#       ll = sum(dpois(y, lambda = exp(eta_fef), log = T))
#     }
#   }
#   
#   return(ll)
#   
# }
