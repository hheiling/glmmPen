
#' Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
#' Minimization (MCECM)
#' 
#' Description
#' 
#' @param dat a list object specifying y (response vector), X (model matrix of all covariates), 
#' Z (model matrix for the random effects), and group (vector whose value indicates 
#' the study, batch, or other group identity to which on observation belongs), with each row an 
#' observation
#' @param lambda0 a non-negative numeric penalty parameter for the fixed effects parameters
#' @param lambda1 a non-negative numeric penalty parameter for the grouped random effects
#' covariance parameters
#' @param conv_EM a non-negative numeric convergence criteria for the convergence of the 
#' EM algorithm. Default is 0.001.
#' @param conv_IRLS a non-negative numeric convergence criteria for the convergenc of the 
#' outer IRLS loop within the grouped coordinate descent algorithm. Default 0.0001
#' @param conv_CD a non-negative numeric convergence criteria for the convergence of the 
#' inner coordinate descent loop within the grouped coordinate descent algorithm. Default 0.0001.
#' @param family a description of the error distribution and link function to be used in the model. 
#' For \code{fit_dat_B} this can be a character string naming a family function or a family function. 
#' See \code{family} options in the Details section.
#' @param offset_fit This can be used to specify an a priori known component to be included in the 
#' linear predictor during fitting. This should be NULL or a numeric vector of length equal to the 
#' number of cases. 
#' @param trace an integer specifying print output to include as function runs. Default value is 0 
#' (details ...)
#' @param penalty character descripting the type of penalty to use in the variable selection procedure.
#' Options include 'MCP', 'SCAD', and 'lasso'. Default is MCP penalty.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions 
#' from the MCP/SCAD/lasso penalty and the ridge, or L2, penalty. \code{alpha=1} is equivalent to 
#' the MCP/SCAD/lasso penalty, while \code{alpha=0} is equivalent to ridge regression. However,
#' \code{alpha=0} is not supported; \code{alpha} may be arbibrarily small, but not exactly zero
#' @param gamma_penalty The tuning parameter of the MCP and SCAD penalties. Not used by Lasso penalty.
#' Default is 4.0 for SCAD and 3.0 for MCP.
#' @param nMC a positive integer for the initial number of Monte Carlo draws
#' @param nMC_max a positive integer for the maximum number of allowed Monte Carlo draws
#' @param t the convergence criteria is based on the average Euclidean distance between 
#' the most recent coefficient estimate and the coefficient estimate from t EM iterations back.
#' Positive integer, default equals 2.
#' @param returnMC logical, should \code{nMC_report} Monte Carlo draws be returned? Default \code{TRUE}.
#' @param nMC_report positive integer specifying number of Monte Carlo posterior draws to return
#' when the function ends. Default 5,000. Warning: the returned draws are formatted in a regular
#' matrix (not a big.matrix). Therefore, depending on the number of random effect covariates (q)
#' and the number of groups (d), choose \code{nMC_report} such that a matrix of size 
#' \code{nMC_report} by (q*d) does not cause memory issues on your operating system.
#' @param maxitEM a positive integer for the maximum number of allowed EM iterations. Default equals 100.
#' @param maxit_CD a positive integer for the maximum number of allowed interations for the
#' coordinate descent algorithms used within each EM iteration. Default equals 250.
#' @param M positive integer specifying the number of posterior draws to use within the 
#' Pajor log-likelihood calculation
#' @param gibbs logical, should Metropolis-within-Gibbs sampling be used for the Monte Carlo draws 
#' (if \code{TRUE}) or should Rejection sampling be used for the draws (if \code{FALSE})? 
#' Default of gibbs = T. Gibbs = F not recommended if greater than 5 random effects
#' @param MwG_sampler character value specifying whether the Metropolis-within-Gibbs procedure 
#' should incorporate an adaptive random walk sampler (default, "random_walk") or an
#' independence sampler ("independence"). 
#' @param adapt_RW_options a list of class "adaptControl" from function \code{\link{adaptControl}} 
#' containing the control parameters for the adaptive random walk Metropolis-within-Gibbs procedure. 
#' Ignored if \code{\link{optimControl}} parameter \code{MwG_sampler} is set to "independence"
#' @param covar character value specifying whether the covariance matrix should be unstructured
#' ("unstructured") or diagonal with no covariances between variables ("independent").
#' Default is "unstructured", but if the number of random effects (including the intercept) is 
#' greater than 7 (i.e. high dimensional), then the algorithm automatically assumes an 
#' independent covariance structure (covar switched to "independent").
#' @param fit_type integer specifying the verion of the grouped coordinate descent algorithm to use.
#' Default equals 2, the 'naive fit', assuming that the random effects covariates do not need any
#' additional standardization or adjustments to account for absence of orthogonalization.
#' For the time being, this version is recommended. Version 4 adds standardization of the random
#' effects covariates compared to version 2. Version 1 adds an adjustment to account for absence
#' of orthogonalization, and version 3 adds standardization compared to version 1.
#' 
#' @section Details: 
#' Accepted families: binomial, poisson, gaussian. Currently, the algorithm is only equipted to handle
#' the canonical links associated with these families (logit link for binomial, log link for poisson,
#' identity link for gaussian).
#' 
#' @return a list with the following elements:
#' \item{coef}{a numeric vector of coefficients of fixed effects estimates and 
#' non-zero estimates of the lower-triangular cholesky decomposition of the random effects
#' covariance matrix (in vector form)}
#' \item{sigma}{random effects covariance matrix}
#' \item{lambda0, lambda1}{the penalty parameters input into the function}
#' \item{covgroup}{Organization of how random effects coefficients are grouped.}
#' \item{J}{a sparse matrix that transforms the non-zero elements of the lower-triangular cholesky 
#' decomposition of the random effects covariance matrix into a vector. For unstructured
#' covariance matrices, dimension of dimension q^2 x (q(q+1)/2) (where q = number of random effects).
#' For independent covariance matrices, q^2 x q.}
#' \item{ll}{estimate of the log likelihood, calculated using the Pajor method}
#' \item{BICh}{the hybrid BIC estimate described in Delattre, Lavielle, and Poursat (2014)}
#' \item{BIC}{Regular BIC estimate}
#' \item{u}{a matrix of the Monte Carlo draws. Organization of columns: first by random effect variable,
#' then by group within variable (i.e. Var1:Grp1 Var1:Grp2 ... Var1:GrpK Var2:Grp1 ... Varq:GrpK)
#' Output if \code{returnMC} = \code{TRUE}}
#' \item{gibbs_accept_rate}{a matrix of the ending gibbs acceptance rates for each variable (columns)
#' and each group (rows)}
#' \item{proposal_SD}{a matrix of the ending proposal standard deviations (used in the adaptive
#' random walk version of the Metropolis-within-Gibbs sampling) for each variable (columns) and
#' each group (rows)}
#' \item{dev}{deviance value}
#' \item{rej_to_gibbs}{logical, did the Monte Carlo sampling switch from Rejection sampling 
#' to Metropolis-within-Gibbs sampling due to unacceptably small acceptance rates in the Rejection sampling?
#' Output only if started initially with Rejection sampling}
#'  
#' 
#' @importFrom bigmemory attach.big.matrix describe
#' @importFrom ncvreg ncvreg
#' @export
fit_dat_B = function(dat, lambda0 = 0, lambda1 = 0, conv_EM = 0.001, 
                     conv_IRLS = 0.0001, conv_CD = 0.0001,
                     family = "binomial", offset_fit = NULL,
                     trace = 0, penalty = c("MCP","SCAD","lasso"),
                     alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0), 
                     group_X = 0:(ncol(dat$X)-1),
                     nMC = 1000, nMC_max = 10000, t = 2,
                     returnMC = T, nMC_report = 5000, ufull = NULL, coeffull = NULL, ufullinit = NULL, 
                     maxitEM = 100, maxit_CD = 250,
                     M = 10^4, gibbs = T, MwG_sampler = c("random_walk","independence"),
                     adapt_RW_options = adaptControl(), covar = c("unstructured","independent"),
                     fit_type = 1){
  
  
  penalty = penalty[1]
  if(!(penalty %in% c("lasso","MCP","SCAD"))){
    stop("penalty ", penalty, " not available, must choose 'lasso', 'MCP', or 'SCAD' \n")
  }
  
  # Set small penalties to zero
  if(lambda0 <=10^-6) lambda0 = 0
  if(lambda1 <=10^-6) lambda1 = 0
  
  y = dat$y
  X = as.matrix(dat$X)
  # Convert sparse Z to dense Z
  Z = Matrix::as.matrix(dat$Z)
  group = dat$group
  
  if(is.character(family)){
    family = get(family, mode = "function", envir = parent.frame())
  }
  if(is.function(family)){
    family = family()
  }
  if(class(family) == "family"){
    fam_fun = family
    link = family$link
    family = family$family
    
    # Re-code link as integer
    ## All link_int will have two digits
    ## First digit corresponds to family that link is canonical for
    ## 1 = binomial, 2 = poisson, 3 = gaussian, 4 = gamma
    ## Second digit is arbitrary enumeration of links
    if(link == "logit"){
      link_int = 10
    }else if(link == "probit"){
      link_int = 11
    }else if(link == "cloglog"){
      link_int = 12
    }else if(link == "log"){
      link_int = 20
    }else if(link == "identity"){
      link_int = 30
    }else if(link == "inverse"){
      link_int = 40
    }
  }
  
  if(is.null(offset_fit)){
    offset_fit = rep(0.0, length(y))
  }
  
  d = nlevels(factor(group))
  
  covar = covar[1]
  if(!(covar %in% c("unstructured","independent"))){
    stop("algorithm currently only handles 'unstructured' or 'independent' covariance structure \n")
  }
  if(covar == "unstructured" & ncol(Z)/d >= 7){
    warning("Due to dimension of sigma covariance matrix, will use covar = 'independent' to simplify computation \n",
            immediate. = T)
    covar = "independent"
  }
  
  #initial fit
  if(family == "binomial"){
    nTotal = rep(1, length(y))
  }else{
    nTotal = NULL
  }
  
  initial_gibbs = gibbs
  
  MwG_sampler = MwG_sampler[1] # Default of random walk
  if(!(MwG_sampler %in% c("independence", "random_walk"))){
    stop("MwG_sampler must be specified as either 'independence' or 'random_walk'")
  }
  
  # Create J matrix (from t(alpha_k kronecker z_ki) * J from paper) and 
  # covgroup = indexing of random effect covariate group
  if(covar == "unstructured"){ # Originally: ncol(Z)/d <= 15 
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
    }
    covgroup = covgroup[lower.tri(covgroup, diag = T)]
  }else{ # covar == "independent". Originally: ncol(Z)/d > 15
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
    }
    covgroup = rep(1:(ncol(Z)/d))
  }
  
  # Initialize cov and coef for start of EM algorithm
  if(!is.null(ufull) & !is.null(coeffull)){
    
    print("using coef from full model to intialize")
    gamma = matrix(J%*%matrix(coeffull[-c(1:ncol(X))], ncol = 1), ncol = ncol(Z)/d)
    cov = var = gamma %*% t(gamma)
    
    ok = which(diag(var) > 0)# & coef[1:ncol(X)] != 0)
    if(length(ok) == 0) ok = 1 # at least include the random intercept
    okindex = NULL
    for(i in 1:(ncol(Z)/d)){
      if(i %in% ok){
        okindex = c(okindex, (i-1)*d + 1:d)
      }
    }
    
    coef = coeffull[1:ncol(X)]
    
  }else{
    
    # Coordinate descent ignoring random effects: naive fit
    # int_only = glm(y ~ 1, family = fam_fun, offset = offset_fit)
    # coef_init = c(int_only$coefficients, rep(0, times = ncol(X)-1))
    # coef = M_step(y, X, fam_fun, coef_init, offset_fit, group_X = group_X,
    #               maxit = 200, maxit_CD = 200, conv = conv_IRLS, fit_type = 2, # 3 if allow grouped fixed effects
    #               penalty = penalty, lambda = lambda0, gamma = gamma_penalty, alpha = alpha)
    
    penalty_factor = numeric(ncol(X)-1)
    penalty_factor[which(group_X[-1] == 0)] = 0
    penalty_factor[which(group_X[-1] != 0)] = 1
    
    fit_naive = ncvreg(X[,-1], y, family = family, penalty = penalty, gamma = gamma_penalty,
                       alpha = alpha, lambda = lambda0, penalty.factor = penalty_factor)
    
    coef = as.numeric(fit_naive$beta)
    
    if(trace == 1) print(coef)
    
    if(ncol(Z)/d > 1){
      vars = rep(10^-10, ncol(Z)/d)
      cov = var = diag(vars)
      gamma = t(chol(var)) # chol outputs upper triangular, so transposing here
    }else{
      vars = 10^-10
      cov = var = matrix(vars, ncol = 1)
      gamma = var
    }
    
    
    ok = which(vars > 0)# & coef[1:ncol(X)] != 0)
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
      Znew2[group == i,seq(i, ncol(Z), by = d)] = Z[group == i, seq(i, ncol(Z), by = d)] %*% gamma
    }
    if(!any(is.na(Znew2))) finish = 1
  }
  
  # intitialize switch-from-rejection-sampling-to-gibbs-sampling counter
  rej_to_gibbs = 0
  
  # initialize adaptive Metropolis-within-Gibbs random walk parameters
  # ignored if use rejection sampling (gibbs = F), but use if gibbs = T or
  # use if initially rejection sampling but switch to gibbs = T
  ## initialize proposal standard deviation
  proposal_SD = matrix(1.0, nrow = d, ncol = ncol(Z)/d)
  print("Initialized proposal_SD")
  print(proposal_SD)
  
  ## initialize batch number to 0
  batch = 0.0
  ## initialize other paramters from adaptControl()
  batch_length = adapt_RW_options$batch_length
  offset_increment = adapt_RW_options$offset
  burnin_batchnum = adapt_RW_options$burnin_batchnum
  gibbs_accept_rate = matrix(NA, nrow = d, ncol = nrow(Z)/d)
  
  if((!is.null(ufull) | !is.null(ufullinit)) & !is.null(coeffull)){
    if(!is.null(ufullinit)){
      print("using u from previous model to intialize")
    }else{
      print("using u from full model to intialize")
      ufullinit = ufull
    }
    
    if(MwG_sampler == "independence"){
      samplemc_out = sample_mc2_BigMat(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group, 
                                       d = d, okindex = okindex, trace = trace, gibbs = gibbs, uold = as.numeric(ufullinit[nrow(ufullinit),]))
    }else{ # MwG_sampler == "random_walk"
      samplemc_out = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group,
                                            d = d, okindex = okindex, trace = trace, gibbs = gibbs, uold = as.numeric(ufullinit[nrow(ufullinit),]),
                                            proposal_SD = proposal_SD, batch = batch, batch_length = batch_length,
                                            offset = offset_increment, burnin_batchnum = burnin_batchnum)
    }
    
    u0 = attach.big.matrix(samplemc_out$u0)
    
    # If specified gibbs = T or if specified gibbs = F but switched to gibbs due to low acceptance rates
    if(gibbs | samplemc_out$switch){
      # If rejection sampling and switched to gibbs sampling due to low acceptance rate:
      if(samplemc_out$switch){ 
        rej_to_gibbs = rej_to_gibbs + 1
        cat("rej_to_gibbs count: ", rej_to_gibbs, "\n")
      }
      
      if(MwG_sampler == "random_walk"){
        gibbs_accept_rate = samplemc_out$gibbs_accept_rate
        batch = samplemc_out$updated_batch
        proposal_SD = samplemc_out$proposal_SD
        
        print("Updated proposal_SD:")
        print(proposal_SD)
      }
      
    }
    
    
    
  }else{
    
    if(MwG_sampler == "independence"){
      samplemc_out = sample_mc2_BigMat(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group, 
                                       d = d, okindex = okindex, trace = trace, gibbs = gibbs, 
                                       uold = rnorm(n = ncol(Z)))
    }else{ # MwG_sampler == "random_walk"
      samplemc_out = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group,
                                            d = d, okindex = okindex, trace = trace, gibbs = gibbs,
                                            uold = rnorm(n = ncol(Z)),
                                            proposal_SD = proposal_SD, batch = batch, batch_length = batch_length,
                                            offset = offset_increment, burnin_batchnum = burnin_batchnum)
    }
    
    u0 = attach.big.matrix(samplemc_out$u0)
    
    # If specified gibbs = T or if specified gibbs = F but switched to gibbs due to low acceptance rates
    if(gibbs | samplemc_out$switch){
      # If rejection sampling and switched to gibbs sampling due to low acceptance rate:
      if(samplemc_out$switch){ 
        rej_to_gibbs = rej_to_gibbs + 1
        cat("rej_to_gibbs count: ", rej_to_gibbs, "\n")
      }
      
      if(MwG_sampler == "random_walk"){
        gibbs_accept_rate = samplemc_out$gibbs_accept_rate
        batch = samplemc_out$updated_batch
        proposal_SD = samplemc_out$proposal_SD
        
        print("Updated proposal_SD:")
        print(proposal_SD)
      }
      
    }
    
  }
  
  nMC2 = nrow(u0)  
  etae = matrix(X %*% coef[1:ncol(X)], nrow = nrow(X), ncol = nrow(u0)) + Znew2%*%t(u0[,])
  mu = matrix(0, nrow = length(y), ncol = ncol(etae))
  for(i in 1:nrow(etae)){
    mu[i,] = invlink(link_int, matrix(etae[i,], ncol=1))
  }
  
  if(nrow(cov) == 1){ # Single random intercept
    cov_record = rep(NA, maxitEM)
  }else{
    cov_record = NULL
  }
  
  diff = rep(NA, maxitEM)
  stopcount = 0
  
  if(family == "poisson"){
    ll0 = sum(rowMeans(dpois(matrix(y, nrow = length(y), ncol = ncol(etae)), lambda = mu, log = T)))
  }else if(family == "binomial"){
    ll0 = sum(rowMeans(dbinom(matrix(y, nrow = length(y), ncol = ncol(etae)), size = 1, prob = mu, log = T)))
  } 
  
  # Record last t coef vectors (each row = coef vector for a past EM iteration)
  # Initialize with initial coef vector
  coef = c(coef, rep(0, length(covgroup)))
  coef_record = matrix(coef, nrow = t, ncol = length(coef), byrow = T)
  # coef_record_all = matrix(NA, nrow = maxitEM, ncol = length(coef), byrow = T)
  
  # Start EM Algorithm (M step first)
  
  for(i in 1:maxitEM){
    
    if(rej_to_gibbs == 3){
      gibbs = T
      cat("permanently switched from rejection sampling to gibbs sampling \n")
      rej_to_gibbs = rej_to_gibbs + 1
    }
    
    oldll = ll0
    
    if(family == "binomial"){
      nTotal = rep(1, length(y[rep(1:nrow(X), each = nrow(u0))]))
    }else{
      nTotal = NULL
    }
    
    # M Step
    coef = M_stepB(y, X, Z, u0@address, nrow(u0), J, group, family=fam_fun, coef, offset=offset_fit,
                   maxit_outer=maxit_CD, maxit_inner=maxit_CD, conv_outer=conv_IRLS, 
                   conv_inner=conv_CD, init=(i == 1), group_X, covgroup,
                   penalty, lambda0, lambda1, gamma=gamma_penalty, alpha,
                   fit_type=fit_type)
    
    problem = F
    if(any(is.na(coef))){
      problem = T
      ll = Inf
    }
    
    u2 = matrix(0, nMC2, ncol(Z))
    for(ii in 1:d){
      u2[,seq(ii, ncol(Z), by = d)] = rmvnorm(n = nMC2,sigma=var)
    }
    
    # Q-function version of log-likelihoood
    ## Need to fix (incorporate big.matrix into this Rcpp function)
    ll0 = Qfun(y, X, Z, u0@address, group, J, matrix(coef, ncol = 1), offset_fit, c(d,ncol(Z)/d), family, link_int)
    
    if(!is.finite(ll0)){
      problem = T
      ll0 = Inf
      print(coef)
    }
    
    if(problem == T){
      stop("Error in M step")
      out = list(coef = coef, sigma = cov, lambda0 = lambda0, lambda1 = lambda1, 
                 covgroup = covgroup, J = J)
      if(returnMC == T) out$u = describe(u0)
      return(out)
    }
    
    if(trace == 1) print(coef)
    
    gamma = matrix(J%*%matrix(coef[-c(1:ncol(X))], ncol = 1), ncol = ncol(Z)/d)
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
        Znew2[group == j,seq(j, ncol(Z), by = d)] = Z[group == j,seq(j, ncol(Z), by = d)]%*%gamma
      }
      if(!any(is.na(Znew2))) finish = 1
    }
    
    
    
    # stopping rule: based on average Euclidean distance (comparing coef from minus t iterations)
    if(i <= t){
      diff[i] = 10^2
    }else{
      diff[i] = sqrt(sum((coef - coef_record[1,])^2)) / length(coef)
    }
    
    # Update latest record of coef
    coef_record = rbind(coef_record[-1,], t(coef))
    # coef_record_all[i,] = coef
    
    if( sum(diff[i:max(i-2, 1)] < conv_EM) >=3 ) break
    
    # if current q is within 95% emp CI of old q, increase nMC
    if(diff[i] > 10^-10){
      nMC = max(round(nMC + min(.15*0.001*nMC/(abs(ll0 - oldll)/abs(ll0)), 250)), 2)+1
    }else{
      nMC = max(round(nMC + .25*0.001*nMC), 5)+1
    }
    
    if(nMC > nMC_max) nMC = nMC_max
    # now update limits  
    print(c(i, nMC , diff[i], sum(coef!=0)))
    
    print("cov:")
    print(cov)
    if(nrow(cov) == 1){
      cov_record[i] = cov
    }
    
    # E Step 
    
    # if(i %% 5 == 0){
    #   uold = rnorm(n = ncol(u0), mean = 0, sd = 1)
    # }else{
    #   uold = as.numeric(u0[nrow(u0),])
    # }
    
    uold = as.numeric(u0[nrow(u0),])
    
    if(MwG_sampler == "independence"){
      samplemc_out = sample_mc2_BigMat(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, trace = trace, family = family, group = group, 
                                       d = d, okindex = okindex, nZ = ncol(Z), gibbs = gibbs, uold = uold)
    }else{ # MwG_sampler == "random_walk"
      samplemc_out = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, trace = trace, family = family, group = group, 
                                           d = d, okindex = okindex, gibbs = gibbs, uold = uold,
                                           proposal_SD = proposal_SD, batch = batch, batch_length = batch_length, 
                                           offset = offset_increment, burnin_batchnum = burnin_batchnum)
    }
    
    u0 = attach.big.matrix(samplemc_out$u0)
    
    # If specified gibbs = T or if specified gibbs = F but switched to gibbs due to low acceptance rates
    if(gibbs | samplemc_out$switch){
      # If rejection sampling and switched to gibbs sampling due to low acceptance rate:
      if(samplemc_out$switch){ 
        rej_to_gibbs = rej_to_gibbs + 1
        cat("rej_to_gibbs count: ", rej_to_gibbs, "\n")
      }
      
      if(MwG_sampler == "random_walk"){
        gibbs_accept_rate = samplemc_out$gibbs_accept_rate
        batch = samplemc_out$updated_batch
        proposal_SD = samplemc_out$proposal_SD
        
        print("Updated proposal_SD:")
        print(proposal_SD)
        print("Updated batch:")
        print(batch)
      }
      
    }
    
    nMC2 = nrow(u0)
    
    if(trace == 1) print(diag(cov))
    gc()
  }
  
  print(sqrt(diag(cov)))
  returnMC
  
  # Another E step for loglik calculation (number draws = M)
  if(MwG_sampler == "independence"){
    samplemc_out = sample_mc2_BigMat(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, trace = trace, family = family, group = group, 
                                     d = d, okindex = okindex, nZ = ncol(Z), gibbs = gibbs, uold = as.numeric(u0[nrow(u0),]))
  }else{ # MwG_sampler == "random_walk"
    samplemc_out = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC, trace = trace, family = family, group = group, 
                                          d = d, okindex = okindex, gibbs = gibbs, uold = as.numeric(u0[nrow(u0),]),
                                          proposal_SD = proposal_SD, batch = batch, batch_length = batch_length, 
                                          offset = offset_increment, burnin_batchnum = burnin_batchnum)
  }
  
  u0 = attach.big.matrix(samplemc_out$u0)
  
  # If specified gibbs = T or if specified gibbs = F but switched to gibbs due to low acceptance rates
  if(gibbs | samplemc_out$switch){
    # If rejection sampling and switched to gibbs sampling due to low acceptance rate:
    if(samplemc_out$switch){ 
      rej_to_gibbs = rej_to_gibbs + 1
      cat("rej_to_gibbs count: ", rej_to_gibbs, "\n")
    }
    
    gibbs_accept_rate = samplemc_out$gibbs_accept_rate
    batch = samplemc_out$updated_batch
    proposal_SD = samplemc_out$proposal_SD
  }
  
  print("Updated proposal_SD:")
  print(proposal_SD)
  print("Updated batch:")
  print(batch)
  
  # Calculate loglik using Pajor method (see logLik_Pajor.R)
  if(sum(diag(cov) == 0) == nrow(cov)){
    eta = X %*% coef[1:ncol(X)]
    if(family == "binomial"){
      ll = sum(dbinom(x = y, size = 1, prob = exp(eta)/(1+exp(eta)), log = T))
    }else if(family == "poisson"){
      ll = sum(dpois(x = y, lambda = log(eta), log = T))
    }else if(family == "gaussian"){
      s2 = sum((y - eta)^2)/(length(y) - (sum(coef[1:ncol(X)] > 0)))
      ll = sum(dnorm(x = y, mean = eta, sd = sqrt(s2), log = T))
    }
  }else{
    ll = CAME_IS(posterior = u0, y = y, X = X, Z = Z, group = group,
                 coef = coef, sigma = cov, family = family, M = M)
  }
  
  
  # Hybrid BIC (Delattre, Lavielle, and Poursat (2014))
  # d = nlevels(group) = number independent subjects/groups
  BICh = -2*ll + sum(coef[-c(1:ncol(X))] != 0)*log(d) + sum(coef[1:ncol(X)] != 0)*log(nrow(X))
  # Usual BIC
  BIC = -2*ll + sum(coef != 0)*log(nrow(X))
  
  if(gibbs){
    out = list(coef = coef, sigma = cov,  
               lambda0 = lambda0, lambda1 = lambda1, 
               covgroup = covgroup, J = J, ll = ll, BICh = BICh, BIC = BIC,
               extra = list(okindex = okindex, Znew2 = Znew2),
               gibbs_accept_rate = gibbs_accept_rate, proposal_SD = proposal_SD)
  }else{
    out = list(coef = coef, sigma = cov,  
               lambda0 = lambda0, lambda1 = lambda1, 
               covgroup = covgroup, J = J, ll = ll, BICh = BICh, BIC = BIC,
               extra = list(okindex = okindex, Znew2 = Znew2))
  }
  
  if(returnMC == T){
    
    if(nrow(u0) >= nMC_report){
      # Take last nMC_report rows
      r_start = nrow(u0) - nMC_report + 1
      r_end = nrow(u0)
      u0_out = bigmemory::as.matrix(u0[c(r_start:r_end),])
    }else{
      if(MwG_sampler == "independence"){
        samplemc_out = sample_mc2(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC_report, trace = trace, family = family, group = group, 
                                  d = d, okindex = okindex, nZ = ncol(Z), gibbs = gibbs, matrix(u0[nrow(u0),], nrow = 1))
      }else{ # MwG_sampler == "random_walk"
        samplemc_out = sample_mc_adapt(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC_report, trace = trace, family = family, group = group, 
                                       d = d, okindex = okindex, gibbs = gibbs, uold = matrix(u0[nrow(u0),], nrow = 1),
                                       proposal_SD = proposal_SD, batch = batch, batch_length = batch_length, 
                                       offset = offset_increment, burnin_batchnum = burnin_batchnum)
      }
      
      u0_out = samplemc_out$u0
    }
    
    out$u = u0_out
    
  }else{
    out$u = matrix(NA, nrow = 1, ncol = 1)
  } 
  
  if((initial_gibbs == F) && rej_to_gibbs > 0){
    if(rej_to_gibbs <= 3){
      cat(sprintf("ending rej_to_gibbs count: %i \n", rej_to_gibbs))
    }else{
      # To correct for additional rej_to_gibbs + 1 when rej_to_gibbs = 3
      cat(sprintf("ending rej_to_gibbs count: %i \n", rej_to_gibbs-1))
    }
  }
  
  if(initial_gibbs == F){
    out$rej_to_gibbs = rej_to_gibbs
  }
  
 # out$coef_record_all = coef_record_all
  if(!is.null(cov_record)){
    out$cov_record = cov_record
  }
  
  return(out)
}


