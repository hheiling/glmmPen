
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
#' @param conv_CD a non-negative numeric convergence criteria for the convergence of the 
#' grouped coordinate descent loop within the M step of the EM algorithm. Default 0.0001.
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
#' @param group_X vector ...
#' @param nMC a positive integer for the initial number of Monte Carlo draws
#' @param nMC_max a positive integer for the maximum number of allowed Monte Carlo draws
#' @param t the convergence criteria is based on the average Euclidean distance between 
#' the most recent coefficient estimate and the coefficient estimate from t EM iterations back.
#' Positive integer, default equals 2.
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
#' @param sampler character value specifying whether the Metropolis-within-Gibbs procedure 
#' should incorporate an adaptive random walk sampler (default, "random_walk") or an
#' independence sampler ("independence"). 
#' @param adapt_RW_options a list of class "adaptControl" from function \code{\link{adaptControl}} 
#' containing the control parameters for the adaptive random walk Metropolis-within-Gibbs procedure. 
#' Ignored if \code{\link{optimControl}} parameter \code{sampler} is set to "independence"
#' @param covar character value specifying whether the covariance matrix should be unstructured
#' ("unstructured") or diagonal with no covariances between variables ("independent").
#' Default is "unstructured", but if the number of random effects (including the intercept) is 
#' greater than or equal to 7 (i.e. high dimensional), then the algorithm automatically assumes an 
#' independent covariance structure (covar switched to "independent").
#' 
#' @section Details: 
#' Accepted families: binomial, poisson, and gaussian. Currently, only canonical links for each 
#' family are available.
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
#' then by group within variable (i.e. Var1:Grp1 Var1:Grp2 ... Var1:GrpK Var2:Grp1 ... Varq:GrpK)}
#' \item{gibbs_accept_rate}{a matrix of the ending gibbs acceptance rates for each variable (columns)
#' and each group (rows)}
#' \item{proposal_SD}{a matrix of the ending proposal standard deviations (used in the adaptive
#' random walk version of the Metropolis-within-Gibbs sampling) for each variable (columns) and
#' each group (rows)}
#' 
#'  
#' @useDynLib glmmPen
#' @importFrom bigmemory attach.big.matrix describe as.big.matrix
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom Matrix Matrix
#' @export
fit_dat_B = function(dat, lambda0 = 0, lambda1 = 0, conv_EM = 0.001, conv_CD = 0.0001,
                     family = "binomial", offset_fit = NULL,
                     trace = 0, penalty = c("MCP","SCAD","lasso"),
                     alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0), 
                     group_X = 0:(ncol(dat$X)-1),
                     nMC_burnin = 500, nMC = 5000, nMC_max = 20000, t = 2,
                     nMC_report = 5000, u_init = NULL, coef_old = NULL, 
                     ufull_describe = NULL, maxitEM = 100, maxit_CD = 250,
                     M = 10^4, sampler = c("stan","random_walk","independence"),
                     adapt_RW_options = adaptControl(), covar = c("unstructured","independent"),
                     var_start = 0.5, max_cores = 1){
  
  ############################################################################################
  # Data input checks
  ############################################################################################
  
  y = dat$y
  X = base::as.matrix(dat$X)
  # Convert sparse Z to dense ZF
  Z = Matrix::as.matrix(dat$Z)
  group = dat$group
  d = nlevels(factor(group))
  
  if (!is.double(y)) {
    tmp <- try(y <- as.double(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be numeric or able to be coerced to numeric", call.=FALSE)
  }
  
  if(!is.matrix(X)){
    stop("X must be a matrix \n")
  }else if(typeof(X)=="integer") storage.mode(X) <- "double"
  
  if(nrow(X) != length(y)){
    stop("the dimension of X and y do not match")
  }
  
  if(is.null(offset_fit)){
    offset_fit = rep(0.0, length(y))
  }else{
    if((!is.numeric(offset_fit)) | (length(offset_fit) != length(y))){
      stop("offset must be a numeric vector of the same length as y")
    }
  }
  
  if(!is.factor(group)){
    group <- factor(group)
  }
  
  # Potential code to use if want to constrict random effects such that zero-valued fixed effects
  #   implies zero-valued random effects
  # if(is.null(XZ_map)){
  #   # Organization of Z: First d columns correspond to all group levels within Variable 1
  #   ## In similar fashion, sequential batches of d columns correspond to a particular variable
  #   # column names of Z: Variable:Group
  #   Z_names = str_split_fixed(colnames(Z), pattern = ":", n = 2)
  #   Z_vars = Z_names[seq(from = 1, to = ncol(Z), by = d), 1]
  #   X_vars = colnames(X)
  #   XZ_map = numeric(ncol(X))
  #   for(j in 1:ncol(X)){
  #     XZ_map[j] = ifelse(X_vars[j] %in% Z_vars, which(Z_vars == X_vars[j]), 0)
  #   }
  # }
  
  family_info = family_export(family)
  fam_fun = family_info$family_fun
  link = family_info$link
  link_int = family_info$link_int # Recoded link as integer
  family = family_info$family
  
  if(length(penalty) > 1){
    penalty = penalty[1]
  }
  if(!(penalty %in% c("lasso","MCP","SCAD"))){
    stop("penalty ", penalty, " not available, must choose 'lasso', 'MCP', or 'SCAD' \n")
  }
  
  # Set small penalties to zero
  if(lambda0 <=10^-6) lambda0 = 0
  if(lambda1 <=10^-6) lambda1 = 0
  
  if(penalty == "MCP" & gamma_penalty <= 1){
    stop("gamma_penalty must be > 1 when using MCP penalty")
  }else if(penalty == "SCAD" & gamma_penalty <= 2){
    stop("gamma_penalty must be > 2 when using SCAD penalty")
  }else if(!is.double(gamma_penalty)){
    tmp <- try(gamma_penalty <- as.double(gamma_penalty), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("gamma_penalty must be numeric or able to be coerced to numeric", call.=FALSE)
  }
  
  if (!is.double(alpha)) {
    tmp <- try(alpha <- as.double(alpha), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("alpha must be numeric or able to be coerced to numeric", call.=FALSE)
  }else if(alpha == 0.0){
    stop("alpha cannot equal 0. Pick a small value > 0 instead (e.g. 0.001) \n");
  }
  
  covar = covar[1]
  if(!(covar %in% c("unstructured","independent"))){
    stop("algorithm currently only handles 'unstructured' or 'independent' covariance structure \n")
  }
  if(covar == "unstructured" & ncol(Z)/d >= 7){
    warning("Due to dimension of sigma covariance matrix, will use covar = 'independent' to simplify computation \n",
            immediate. = T)
    covar = "independent"
  }
  
  sampler = sampler[1] # Default of stan
  if(!(sampler %in% c("independence","random_walk","stan"))){
    stop("sampler must be specified as either 'independence' or 'random_walk' or 'stan'")
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
  
  ############################################################################################
  # Initialization
  ############################################################################################
  
  # Initialize cov and coef for start of EM algorithm
  
  if(!is.null(coef_old)){
    
    print("using coef from past model to intialize")
    gamma = matrix(J%*%matrix(coef_old[-c(1:ncol(X))], ncol = 1), ncol = ncol(Z)/d)
    cov = var = gamma %*% t(gamma)
    # cov = var = Matrix(data = gamma %*% t(gamma), sparse = T)
    
    # If the random effect for a variable penalized out in a past model, still possible for that
    # variable to not be penalized out this next model. Therefore, initialize the variance of
    # these penalized-out random effects to have a non-zero variance (specified by var_start)
    ok = which(diag(var) > 0)
    if(length(ok) != length(diag(var))){
      ranef0 = which(diag(var) == 0)
      for(j in ranef0){
        var[j,j] = var_start
      }
      
      # Re-define cov and gamma with new variances
      cov = var
      gamma = t(chol(cov))
    }
    
    coef = coef_old[1:ncol(X)]
    
  }else{
    
    # Coordinate descent ignoring random effects: naive fit
    
    penalty_factor = numeric(ncol(X))
    penalty_factor[which(group_X == 0)] = 0
    penalty_factor[which(group_X != 0)] = 1
    
    # if(family == "negbin"){
    #   # Initialize coefficients by assuming poisson distribution (assume no overdispersion to start)
    #   family0 = "poisson"
    #   fam_fun0 = poisson(link = link)
    # }else{
    #   family0 = family
    #   fam_fun0 = fam_fun
    # }
    
    # Initialize coef as intercept-only for input into CD function
    IntOnly = glm(y ~ 1, family = fam_fun, offset = offset_fit)
    coef_init = c(IntOnly$coefficients, rep(0, length = (ncol(X)-1)))
    
    fit_naive = CD(y, X, family = family, link = link_int, offset = offset_fit, coef_init = coef_init,
                   maxit_CD = maxit_CD, conv = conv_CD, penalty = penalty, lambda = lambda0,
                   gamma = gamma_penalty, alpha = alpha, penalty_factor = penalty_factor, trace = trace)
    
    coef = fit_naive
    
    
    if(trace >= 1){
      print("initialized fixed effects:")
      print(coef)
    } 
    
    if(any(is.na(coef))){
      print(coef)
      stop("Error in initial coefficient fit: NAs produced")
    }
    
    # Initialize covariance matrix (cov)
    if(ncol(Z)/d > 1){
      vars = rep(var_start, ncol(Z)/d)
      cov = var = diag(vars)
      gamma = t(chol(var)) # chol outputs upper triangular, so transposing here to get lower triangular
    }else{
      vars = var_start
      cov = var = matrix(vars, ncol = 1)
      gamma = var
    }
    
    if(trace >= 1){
      print("initialized covariance matrix:")
      print(cov)
    }
    
    
  } # End if-else !is.null(coef_old)
  
  # Initialize additional parameters
  
  # if(family == "negbin"){
  #   # Initialize phi
  #   # variance of negbin: mu + phi * mu^2
  #   # arma::vec y, arma::mat eta, int link, int limit, double eps
  #   eta = X %*% coef
  #   phi = phi_ml_init(y, eta, link_int, 200, conv_CD)
  #   cat("phi: ", phi, "\n")
  # }else{
  #   phi = 0.0
  # }
  phi = 0.0
  
  ## Initialize residual standard error for gaussian family
  if(family == "gaussian"){
    
    if(is.null(u_init)){
      # Find initial estimate of the variance of the error term using initial fixed effects only
      eta = X %*% coef
      s2_g = sum((y - invlink(link_int, eta))^2)
      sig_g = sqrt(s2_g / length(y))
    }else{
      # Find initial estimate of variance of the error term using fixed and random effects from 
      # last round of selection
      u_big = as.big.matrix(u_init)
      sig_g = sig_gaus(y, X, Z, u_big@address, group, J, c(coef, gamma), offset_fit, c(d, ncol(Z)/d), link_int)
    }
    
    if(trace >= 1){
      cat("initial residual error SD: ", sig_g, "\n")
    }
    
  }else{
    sig_g = 1.0 # specify an arbitrary placeholder (will not be used in calculations)
  }
  
  
  
  
  Znew2 = Z
  finish = 0
  while(finish == 0){
    for(i in 1:d){
      Znew2[group == i,seq(i, ncol(Z), by = d)] = Z[group == i, seq(i, ncol(Z), by = d)] %*% gamma
    }
    if(!any(is.na(Znew2))) finish = 1
  }
  
  # initialize adaptive Metropolis-within-Gibbs random walk parameters
  ## initialize proposal standard deviation
  proposal_SD = matrix(1.0, nrow = d, ncol = ncol(Z)/d)
  if(sampler == "random_walk" & trace >= 2){
    print("initialized proposal_SD:")
    print(proposal_SD)
  }
  
  ## initialize batch number to 0
  batch = 0.0
  ## initialize other paramters from adaptControl()
  batch_length = adapt_RW_options$batch_length
  offset_increment = adapt_RW_options$offset
  gibbs_accept_rate = matrix(NA, nrow = d, ncol = nrow(Z)/d)
  
  # Determine initial number of cores to use
  if(nMC < 2000){
    num_cores = 1
  }else{
    num_cores = 2 + (nMC - 2000) %/% 1000
    if(num_cores > max_cores) num_cores = max_cores
  }
  
  # At start of EM algorithm, acquire posterior draws from all random effect variables
  ranef_idx = 1:length(diag(cov))
  
  if(!is.null(u_init)){
    print("using u from previous model to initialize")
    
    if(is.matrix(u_init)){
      uold = as.numeric(u_init[nrow(u_init),])
    }else if(is.vector(u_init)){
      uold = u_init
    }else{
      stop("u_init must be either a matrix of posterior draws or a vector of a single set of posterior draws")
    }
    
    Estep_out = E_step(coef = coef, ranef_idx = ranef_idx, y=y, X=X, Znew2=Znew2, group=group, 
                       nMC=nMC, nMC_burnin=nMC_burnin, family=family, link=link, phi=phi, sig_g=sig_g,
                       sampler=sampler, d=d, uold=uold, proposal_SD=proposal_SD, 
                       batch=batch, batch_length=batch_length, offset_increment=offset_increment, 
                       trace=trace, num_cores=num_cores)
    
    u0 = attach.big.matrix(Estep_out$u0)
    proposal_SD = Estep_out$proposal_SD
    gibbs_accept_rate = Estep_out$proposal_SD
    batch = Estep_out$updated_batch
    
  }else{ # u_init is null, initialize posterior with random draw from N(0,1)
    
    Estep_out = E_step(coef = coef, ranef_idx = ranef_idx, y=y, X=X, Znew2=Znew2, group=group, 
                       nMC=nMC, nMC_burnin=nMC_burnin, family=family, link=link, phi=phi, sig_g=sig_g,
                       sampler=sampler, d=d, uold=rnorm(n = ncol(Z)), proposal_SD=proposal_SD, 
                       batch=batch, batch_length=batch_length, offset_increment=offset_increment, 
                       trace=trace, num_cores=num_cores)
    
    u0 = attach.big.matrix(Estep_out$u0)
    proposal_SD = Estep_out$proposal_SD
    gibbs_accept_rate = Estep_out$proposal_SD
    batch = Estep_out$updated_batch
    
  } # end if-else is.null(u_init)
  
  nMC2 = nrow(u0)
  
  if(nrow(cov) == 1){ # Single random intercept
    cov_record = rep(NA, maxitEM)
  }else{
    cov_record = NULL
  }
  
  diff = rep(NA, maxitEM)
  stopcount = 0
  
  # Record last t coef vectors (each row = coef vector for a past EM iteration)
  # Initialize with initial coef vector
  coef = c(coef, rep(0, length(covgroup)))
  coef_record = matrix(coef, nrow = t, ncol = length(coef), byrow = T)
  # coef_record_all = matrix(NA, nrow = maxitEM, ncol = length(coef), byrow = T)
  
  problem = FALSE # Determining if issue with EM algorithm results
  
  ############################################################################################
  # EM Algorithm
  ############################################################################################
  
  # Start EM Algorithm
  EM_converged = 0
  
  for(i in 1:maxitEM){
    
    # M Step
    # y, X, Z, u_address, M, J, group, family, link_int, coef, offset, phi,
    # maxit_CD = 250, conv_CD = 0.0001,
    # init, group_X = 0:(ncol(X)-1), covgroup,
    # penalty, lambda0, lambda1, gamma, alpha = 1.0
    coef = M_step(y=y, X=X, Z=Z, u_address=u0@address, M=nrow(u0), J=J, 
                  group=group, family=family, link_int=link_int, coef=coef, offset=offset_fit,
                  phi=phi, maxit_CD=maxit_CD, conv_CD=conv_CD, init=(i == 1), group_X=group_X, 
                  covgroup=covgroup, penalty=penalty, lambda0=lambda0, lambda1=lambda1, 
                  gamma=gamma_penalty, alpha=alpha, trace=trace)
    
    if(family == "gaussian"){
      # if family = 'gaussian', calculate sigma of error term (standard deviation)
      sig_g = sig_gaus(y, X, Z, u0@address, group, J, coef, offset_fit, c(d, ncol(Z)/d), link_int)
      cat("sig_g: ", sig_g, "\n")
    }
    
    # if(family == "negbin"){
    #   phi = coef[length(coef)]
    #   coef = coef[-length(coef)]
    # }
    
    if(trace >= 1){
      print("Updated coef:")
      print(coef)
    }
    
    # If fixed effects penalized to 0, make sure associated random effects also penalized to 0
    
    u2 = matrix(0, nMC2, ncol(Z))
    for(ii in 1:d){
      u2[,seq(ii, ncol(Z), by = d)] = rmvnorm(n = nMC2,sigma=var)
    }
    
    # Q-function estimate
    ll0 = Qfun(y, X, Z, u0@address, group, J, matrix(coef, ncol = 1), offset_fit, c(d,ncol(Z)/d), family, link_int, phi)
    
    if(!is.finite(ll0) | any(is.na(coef)) | any(abs(coef) > 10^5)){
      # For now, assume divergent in any abs(coefficient) > 10^5
      problem = T
      ll0 = Inf
      print("Updated coef:")
      print(coef)
    }
    
    if(problem == T){
      warning("Error in M step", immediate. = T)
      out = list(coef = coef, sigma = cov,  
                 lambda0 = lambda0, lambda1 = lambda1, 
                 covgroup = covgroup, J = J, ll = -Inf, BIC = Inf, BICh = Inf,
                 extra = list(Znew2 = Znew2),
                 gibbs_accept_rate = gibbs_accept_rate, proposal_SD = proposal_SD,
                 EM_iter = i)
      if(is.null(coef_old)){
        out$coef_naive = fit_naive
      }
      
      if(nrow(u0) >= nMC_report){
        # Take last nMC_report rows
        r_start = nrow(u0) - nMC_report + 1
        r_end = nrow(u0)
        u0_out = bigmemory::as.matrix(u0[c(r_start:r_end),])
      }else{
        # Use all available u0 entries
        u0_out = bigmemory::as.matrix(u0)
      }
      out$u = u0_out
      
      return(out)
    } # End problem == T
    
    gamma = matrix(J%*%matrix(coef[-c(1:ncol(X))], ncol = 1), ncol = ncol(Z)/d)
    cov = var = gamma %*% t(gamma)
    
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
    
    if( sum(diff[i:max(i-2, 1)] < conv_EM) >=3 ){
      EM_converged = 1
      break
    } 
    
    # Increase nMC in nMC below nMC_max
    if(nMC < nMC_max){
      # Set f at constant 1.05
      nMC = nMC * 1.05
    }
    # If after update nMC exceeds nMC_max, reduce to nMC_max value
    if(nMC > nMC_max) nMC = nMC_max
    
    # now update limits  
    update_limits = c(i, nMC , diff[i], sum(coef!=0))
    names(update_limits) = c("Iter","nMC","EM diff","Non0 Coef")
    print(update_limits)
    
    if(trace >= 1){
      if(nrow(cov) <= 5){
        print("cov:")
        print(cov)
      }else{
        print("cov diagonal:")
        print(diag(cov))
      }
    }
    
    if(nrow(cov) == 1){
      cov_record[i] = cov
    }
    
    # E Step 
    
    # Initial points for Metropolis within Gibbs E step algorithms
    uold = as.numeric(u0[nrow(u0),])
    ranef_idx = which(diag(cov) > 0)
    
    Estep_out = E_step(coef = coef, ranef_idx = ranef_idx, y=y, X=X, Znew2=Znew2, group=group, 
                       nMC=nMC, nMC_burnin=nMC_burnin, family=family, link=link, phi=phi, sig_g=sig_g,
                       sampler=sampler, d=d, uold=uold, proposal_SD=proposal_SD, 
                       batch=batch, batch_length=batch_length, offset_increment=offset_increment, 
                       trace=trace, num_cores=num_cores)
    
    u0 = attach.big.matrix(Estep_out$u0)
    proposal_SD = Estep_out$proposal_SD
    gibbs_accept_rate = Estep_out$proposal_SD
    batch = Estep_out$updated_batch
    
    nMC2 = nrow(u0)
    
  }
  
  if(EM_converged == 0){
    warning("glmmPen algorithm did not converge within maxit_EM iterations", immediate. = T)
  }
  
  # Another E step for loglik calculation (number draws = M)
  
  # Note: Assume that at end of EM algorithm, which(diag(cov) > 0) corresponds with the
  # sampling restrictions specified within EM algorithm:
    # Restrict sampling such that variables are NOT sampled if:
    ## random effect for the variable was penalized to zero in a previous M step
    ## fixed effect for the variable was penalized to zero in a prevoius M step
    ## Make sure to always keep random intercept
  
  Estep_out = E_step(coef=coef, ranef_idx=which(diag(cov) > 0), y=y, X=X, Znew2=Znew2, group=group, 
                     nMC=nMC, nMC_burnin=nMC_burnin, family=family, link=link, phi=phi, sig_g = sig_g,
                     sampler=sampler, d=d, uold=as.numeric(u0[nrow(u0),]), proposal_SD=proposal_SD, 
                     batch=batch, batch_length=batch_length, offset_increment=offset_increment, 
                     trace=trace, num_cores=num_cores)
  
  u0 = attach.big.matrix(Estep_out$u0)
  proposal_SD = Estep_out$proposal_SD
  gibbs_accept_rate = Estep_out$proposal_SD
  batch = Estep_out$updated_batch
 
  if(sampler == "random_walk" & trace >= 2){
    print("Updated proposal_SD:")
    print(proposal_SD)
    print("Updated batch:")
    print(batch)
  }
  
  # Calculate loglik using Pajor method (see "logLik_Pajor.R")
  ll = CAME_IS(posterior = u0, y = y, X = X, Z = Z, group = group,
               coef = coef, sigma = cov, family = fam_fun, M = M, gaus_sig = sig_g)
  
  # Hybrid BIC (Delattre, Lavielle, and Poursat (2014))
  # d = nlevels(group) = number independent subjects/groups
  BICh = -2*ll + sum(coef[-c(1:ncol(X))] != 0)*log(d) + sum(coef[1:ncol(X)] != 0)*log(nrow(X))
  # Usual BIC
  BIC = -2*ll + sum(coef != 0)*log(nrow(X))
  
  # Calculated BICq
  ## Using posterior draws from the full model (ufull) and the coefficients from the
  ## current penalized model, calculate the Q function
  if(!is.null(ufull_describe)){
    ufull = attach.big.matrix(ufull_describe)
    q1 = Qfun(y, X, Z, ufull@address, group, J, matrix(coef, ncol = 1), offset_fit, c(d,ncol(Z)/d), family, link_int, phi)
    q2 = 0
    for(k in 1:d){
      cols_idx = seq(from = k, to = ncol(Z), by = d)
      post_k = ufull[,cols_idx]
      q2 = q2 + sum(dmvnorm(post_k, log=T)) / nrow(ufull)
    }
    llq = q1 + q2
    BICq = -2*llq + sum(coef != 0)*log(nrow(X))
  }else{
    BICq = NA
  }
  
  if(sampler %in% c("random_walk","independence")){
    out = list(coef = coef, sigma = cov,  
               lambda0 = lambda0, lambda1 = lambda1, 
               covgroup = covgroup, J = J, ll = ll, BICh = BICh, BIC = BIC, BICq = BICq,
               extra = list(Znew2 = Znew2),
               gibbs_accept_rate = gibbs_accept_rate, proposal_SD = proposal_SD,
               EM_iter = i)
  }else if(sampler == "stan"){
    out = list(coef = coef, sigma = cov,  
               lambda0 = lambda0, lambda1 = lambda1, 
               covgroup = covgroup, J = J, ll = ll, BICh = BICh, BIC = BIC, BICq = BICq,
               extra = list(Znew2 = Znew2), EM_iter = i)
  }
  
  if((is.null(coef_old)) & (trace >= 1)){
    out$coef_naive = fit_naive
  }
  
  if(nrow(u0) >= nMC_report){
    # Take last nMC_report rows
    r_start = nrow(u0) - nMC_report + 1
    r_end = nrow(u0)
    u0_out = bigmemory::as.matrix(u0[c(r_start:r_end),])
  }else{
    
    Estep_out = E_step(coef=coef, ranef_idx=which(diag(cov) > 0), y=y, X=X, Znew2=Znew2, group=group, 
                       nMC=nMC, nMC_burnin=nMC_burnin, family=family, link=link, phi=phi, sig_g=sig_g,
                       sampler=sampler, d=d, uold=as.numeric(u0[nrow(u0),]), proposal_SD=proposal_SD, 
                       batch=batch, batch_length=batch_length, offset_increment=offset_increment, 
                       trace=trace, num_cores=num_cores)
    
    u0 = attach.big.matrix(Estep_out$u0)
    u0_out = bigmemory::as.matrix(u0)
    proposal_SD = Estep_out$proposal_SD
    gibbs_accept_rate = Estep_out$proposal_SD
    batch = Estep_out$updated_batch
    
  }
  
  out$u = u0_out
  
  # if(family == "negbin"){
  #   out$phi = phi
  # }
  
  if(family == "gaussian"){
    out$sigma_gaus = sig_g
  }
  
 # out$coef_record_all = coef_record_all
  if(!is.null(cov_record)){
    out$cov_record = cov_record
  }
  
  return(out)
}


# negative binomial. In order to select the
# negative binomial distribution, either set family = "negbin" (which assumes a canonical log link)
# or family = \code{MASS::negative.binomial} using any arbitrary \code{theta} value and the desired
# link function.