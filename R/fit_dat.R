
#' @title Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
#' Minimization (MCECM)
#' 
#' \code{fit_dat} is used to fit a penalized generalized mixed model
#' via Monte Carlo Expectation Conditional Minimization (MCECM) for 
#' a single tuning parameter combinations and is called within
#' \code{glmmPen} or \code{glmm} (cannot be called directly by user)
#' 
#' @inheritParams optimControl 
#' @inheritParams lambdaControl
#' @inheritParams glmmPen
#' @param dat a list object specifying y (response vector), X (model matrix of all covariates), 
#' Z (model matrix for the random effects), and group (numeric factor vector whose value indicates 
#' the study, batch, or other group identity to which on observation belongs)
#' @param offset_fit This can be used to specify an a priori known component to be included in the 
#' linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the 
#' number of cases. 
#' @param group_X vector describing the grouping of the covariates in the model matrix.
#' @param nMC a positive integer for the initial number of Monte Carlo draws. See the \code{nMC_start}
#' argument in \code{\link{optimControl}} for more details.
#' @param u_init matrix giving values to initialize samples from the posterior. If 
#' Binomial or Poisson families, only need a single row to initialize samples from
#' the posterior; if Gaussian family, multiple rows needed to initialize the estimate
#' of the residual error (needed for the E-step). Columns correspond to the 
#' columns of the Z random effect model matrix.
#' @param coef_old vector giving values to initialized the coefficients (both fixed
#' and random effects)
#' @param ufull_describe output from \code{bigmemory::describe} (which returns a list 
#' of the information needed to attach to a big.matrix object) applied to the
#' big.matrix of posterior samples from the 'full' model. The big.matrix 
#' described by the object is used to calculate the BIC-ICQ value for the model.
#' @param ranef_keep vector of 0s and 1s indicating which random effects should 
#' be considered as non-zero at the start of the algorithm. For each random effect,
#' 1 indicates the random effect should be considered non-zero at start of algorithm,
#' 0 indicates otherwise. The first element for the random intercept should always be 1.
#' @param checks_complete boolean value indicating whether the function has been called within
#' \code{glmm} or \code{glmmPen} or whether the function has been called by itself. 
#' Used for package testing purposes (user cannot directly call \code{fit_dat}). If true,
#' performs additional checks on the input data. If false, assumes data input checks have 
#' already been performed. 
#' @param conv_type integer specifying which type of convergence criteria to use. Default 1 specifies
#' using the average Eucledian distance, and 2 specifies using relative change in the Q-function
#' estimate. For now, all calls to \code{fit_dat} within the \code{glmmPen} framework
#' restrict this convergence type to be the average Euclidean distance. However,
#' we keep this argument in case we decide to allow multiple convergence type options in
#' future versions of the package.
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
#' \item{BICNgrps}{BIC estimate with N = number of groups in penalty term instead of N = number
#' of total observations.}
#' \item{BICq}{BIC-ICQ estimate}
#' \item{u}{a matrix of the Monte Carlo draws. Organization of columns: first by random effect variable,
#' then by group within variable (i.e. Var1:Grp1 Var1:Grp2 ... Var1:GrpK Var2:Grp1 ... Varq:GrpK)}
#' \item{gibbs_accept_rate}{a matrix of the ending gibbs acceptance rates for each variable (columns)
#' and each group (rows) when the sampler is either "random_walk" or "independence"}
#' \item{proposal_SD}{a matrix of the ending proposal standard deviations (used in the adaptive
#' random walk version of the Metropolis-within-Gibbs sampling) for each variable (columns) and
#' each group (rows)}
#' 
#' @useDynLib glmmPen
#' @importFrom bigmemory attach.big.matrix describe as.big.matrix
#' @importFrom mvtnorm dmvnorm
#' @importFrom Matrix Matrix
fit_dat = function(dat, lambda0 = 0, lambda1 = 0, conv_EM = 0.001, conv_CD = 0.0001,
                     family = "binomial", offset_fit = NULL,
                     trace = 0, penalty = c("MCP","SCAD","lasso"),
                     alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0), 
                     group_X = 0:(ncol(dat$X)-1),
                     nMC_burnin = 250, nMC = 250, nMC_max = 5000, t = 2, mcc = 2,
                     u_init = NULL, coef_old = NULL, 
                     ufull_describe = NULL, maxitEM = 50, maxit_CD = 250,
                     M = 10^4, sampler = c("stan","random_walk","independence"),
                     adapt_RW_options = adaptControl(), covar = c("unstructured","independent"),
                     var_start = 1.0, logLik_calc = F, checks_complete = F,
                     ranef_keep = rep(1, times = (ncol(dat$Z)/nlevels(dat$group))),
                     conv_type = 1){
  
  ############################################################################################
  # Data input checks
  ############################################################################################
  
  # Extract data from dat list object
  
  y = dat$y
  X = base::as.matrix(dat$X)
  # Convert sparse Z to dense ZF
  Z = Matrix::as.matrix(dat$Z)
  group = dat$group
  d = nlevels(factor(group))
  
  # Testing purposes: if calling fit_dat directly instead of within glmm or glmmPen, perform data checks
  
  if(!checks_complete){
    
    if(!is.double(y)) {
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
      group <- as.factor(as.numeric(group))
    }
    
    # Check penalty parameters
    penalty_check = checkPenalty(penalty, gamma_penalty, alpha)
    penalty = penalty_check$penalty
    gamma_penalty = penalty_check$gamma_penalty
    alpha = penalty_check$alpha
    
    # Check covariance matrix specification
    covar = checkCovar(covar)
    
    # Check sampler specification
    sampler = checkSampler(sampler)
    
  } # End checks of input
  
  
  # Extract relevant family information
  # See family_export() in "family_export.R"
  family_info = family_export(family)
  fam_fun = family_info$family_fun
  link = family_info$link
  link_int = family_info$link_int # Recoded link as integer, see "family_export.R" for details
  family = family_info$family
  
  # Set small penalties to zero
  if(lambda0 <=10^-6) lambda0 = 0
  if(lambda1 <=10^-6) lambda1 = 0
  
  
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
  ## coef: Initializes fixed effect coefficients
  ## gamma: Initializes random effects coefficients
  ## cov: Initializes covariance matrix. If variance for a predictor (the diagonal) is
  ##    greater than 0, E step will sample from posterior for this predictor. Otherwise,
  ##    no posterior samples will be drawn for that predictor
  
  if(!is.null(coef_old)){
    
    cat("using coef from past model to intialize \n")
    gamma = matrix(J%*%matrix(coef_old[-c(1:ncol(X))], ncol = 1), ncol = ncol(Z)/d)
    cov = var = gamma %*% t(gamma)
    
    # If the random effect for a variable penalized out in a past model, still possible for that
    # variable to not be penalized out this next model. 
    # Whether or not the random effect should or should not be considered for the model fit
    # is based on the 'ranef_keep' variable (see select_tune() and glmmPen() code
    # for logic on restrictions for random effects)
    # If the j-th element of ranef_keep = 1 but the random effect was penalized
    # out in a past model, initialize the variance of
    # these penalized-out random effects to have a non-zero variance (specified by var_start)
    # However, if ranef_keep = 0, then keep this
    # effect with a 0 variance
    # In E step, only random effects with non-zero variance in covariance matrix have MCMC
    # samples from the posterior distribution
    # Also need to update gamma (gamma used in E step)
    ok = which(diag(var) > 0)
    if(length(ok) != length(diag(var))){
      ranef0 = which(diag(var) == 0)
      for(j in ranef0){
        if(ranef_keep[j] == 1){
          var[j,j] = var_start
        }
      }
      cov = var
      # Update gamma matrix
      ## add small constant to diagonals so that chol() operation works
      gamma = t(chol(var + diag(10^-6,nrow(var))))
    }
    
    if(covar == "unstructured"){
      gamma_vec = c(gamma)[which(lower.tri(matrix(0,nrow=nrow(cov),ncol=ncol(cov)),diag=T))]
    }else if(covar == "independent"){
      gamma_vec = diag(gamma)
    }
    
    # Initialize fixed and random effects 
    coef = c(coef_old[1:ncol(X)], gamma_vec)
    
  }else{ # Will not initialize with past model results
    
    # Coordinate descent ignoring random effects: naive fit
    
    # penalty_factor: indicaor of whether the fixed effect can be penalized in this naive fit
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
    
    # Coordinate descent ignoring random effects: naive fit
    ## See "Mstep.R" for CD function
    fit_naive = CD(y, X, family = family, link = link_int, offset = offset_fit, coef_init = coef_init,
                   maxit_CD = maxit_CD, conv = conv_CD, penalty = penalty, lambda = lambda0,
                   gamma = gamma_penalty, alpha = alpha, penalty_factor = penalty_factor, trace = trace)
    
    coef = fit_naive
    
    if(trace >= 1){
      cat("initialized fixed effects: \n")
      cat(coef, "\n")
    } 
    
    if(any(is.na(coef))){
      cat(coef, "\n")
      stop("Error in initial coefficient fit: NAs produced")
    }
    
    # Initialize covariance matrix (cov)
    if(ncol(Z)/d > 1){
      vars = rep(var_start, ncol(Z)/d) * ranef_keep
      cov = var = diag(vars)
      gamma = diag(sqrt(vars)) 
    }else{
      vars = var_start
      cov = var = matrix(vars, ncol = 1)
      gamma = matrix(sqrt(var), ncol = 1)
    }
    
    if(trace >= 1){
      if(covar == "unstructured"){
        cat("initialized covariance matrix: \n")
        print(cov)
      }else if(covar == "independent"){
        cat("initialized covariance matrix diagonal: \n", diag(cov), "\n")
      }
    }
    
    if(covar == "unstructured"){
      gamma_vec = c(gamma)[which(lower.tri(matrix(0,nrow=nrow(cov),ncol=ncol(cov)),diag=T))]
    }else if(covar == "independent"){
      gamma_vec = diag(gamma)
    }
    coef = c(coef[1:ncol(X)], gamma_vec)
    
    
  } # End if-else !is.null(coef_old)
  
  ############################################################################################
  # Initialization of additional parameters
  ############################################################################################
  
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
      eta = X %*% coef[1:ncol(X)]
      s2_g = sum((y - invlink(link_int, eta))^2)
      sig_g = sqrt(s2_g / length(y))
    }else{
      # Find initial estimate of variance of the error term using fixed and random effects from 
      # last round of selection
      u_big = as.big.matrix(u_init)
      sig_g = sig_gaus(y, X, Z, u_big@address, group, J, coef, offset_fit, c(d, ncol(Z)/d), link_int)
      
    }
   
    if(trace >= 1){
      cat("initial residual error SD: ", sig_g, "\n")
    }
    
  }else{
    sig_g = 1.0 # specify an arbitrary placeholder (will not be used in calculations)
  } # End if-else family == "gaussian"
  
  # initialize adaptive Metropolis-within-Gibbs random walk parameters
  if(sampler == "random_walk"){
    ## initialize proposal standard deviation
    proposal_SD = matrix(1.0, nrow = d, ncol = ncol(Z)/d)
    ## initialize batch number to 0
    batch = 0.0
    ## initialize other paramters from adaptControl()
    batch_length = adapt_RW_options$batch_length
    offset_increment = adapt_RW_options$offset
  }else{
    proposal_SD = NULL
    batch = NULL
    batch_length = NULL
    offset_increment = NULL
  }
  
  # Record Gibbs acceptance rate for both sampler = "random_walk" and sampler = "independence"
  if(sampler %in% c("random_walk","independence")){
    gibbs_accept_rate = matrix(NA, nrow = d, ncol = nrow(Z)/d)
  }
  
  if(sampler == "random_walk" & trace >= 2){
    cat("initialized proposal_SD: \n")
    print(proposal_SD)
  }
  
  # Znew2: For each group k = 1,...,d, calculate Znew2 = Z %*% gamma (calculate for each group individually)
  # Used within E-step
  Znew2 = Z
  for(i in 1:d){
    Znew2[group == i,seq(i, ncol(Z), by = d)] = Z[group == i, seq(i, ncol(Z), by = d)] %*% gamma
  }

  # Initialization of posteiror samples
  
  # At start of EM algorithm, acquire posterior draws from all relevant random effect variables
  # Note: restriction of which random effects to use are based on the ranef_keep variable
  # (see gamma and cov initialization in above code for details)
  ranef_idx = which(diag(cov) != 0)
  
  if(!is.null(u_init)){
    cat("using results from previous model to initialize posterior samples \n")
    if(is.matrix(u_init)){
      uold = as.numeric(u_init[nrow(u_init),])
    }else{
      stop("u_init must be a matrix of posterior draws")
    }
    
    Estep_out = E_step(coef = coef, ranef_idx = ranef_idx, y=y, X=X, Znew2=Znew2, group=group, offset_fit = offset_fit,
                       nMC=nMC, nMC_burnin=nMC_burnin, family=family, link=link, phi=phi, sig_g=sig_g,
                       sampler=sampler, d=d, uold=uold, proposal_SD=proposal_SD, 
                       batch=batch, batch_length=batch_length, offset_increment=offset_increment, 
                       trace=trace)
    
    u0 = attach.big.matrix(Estep_out$u0)
    proposal_SD = Estep_out$proposal_SD
    gibbs_accept_rate = Estep_out$gibbs_accept_rate
    batch = Estep_out$updated_batch
    
  }else{ # u_init is null, initialize posterior with random draw from N(0,1)
    
    Estep_out = E_step(coef = coef, ranef_idx = ranef_idx, y=y, X=X, Znew2=Znew2, group=group, offset_fit = offset_fit,
                       nMC=nMC, nMC_burnin=nMC_burnin, family=family, link=link, phi=phi, sig_g=sig_g,
                       sampler=sampler, d=d, uold=rnorm(n = ncol(Z)), proposal_SD=proposal_SD, 
                       batch=batch, batch_length=batch_length, offset_increment=offset_increment, 
                       trace=trace)
    
    u0 = attach.big.matrix(Estep_out$u0)
    proposal_SD = Estep_out$proposal_SD
    gibbs_accept_rate = Estep_out$gibbs_accept_rate
    batch = Estep_out$updated_batch
    
  } # end if-else is.null(u_init)
  
  
  ############################################################################################
  # EM Algorithm
  ############################################################################################
  
  # Start EM Algorithm
  EM_converged = 0
  # Determining if issue with EM algorithm results
  problem = FALSE 
  
  # Record convergence criteria value (average Euclidean distance) for each EM iteration
  diff = rep(NA, maxitEM)
  
  # Record last t coef vectors (each row = coef vector for a past EM iteration)
  # Initialize with initial coef vector
  coef_record = matrix(coef, nrow = t, ncol = length(coef), byrow = T)
  
  # If using Q-function as convergence criteria, save old Q-function value
  Q_est = rep(NA, maxitEM)
  
  for(i in 1:maxitEM){
    
    ############################################################################################
    # M Step
    ############################################################################################
    if(sum(diag(cov) > 0) > 15){
      # Alternative maxit_CD for when q is large (to speed up start of algorithm when
      # few random effects penalized out of model)
      maxit_CD_use = 25
    }else{
      maxit_CD_use = maxit_CD
    }
    
    # cat("Start M-step \n")
    coef = M_step(y=y, X=X, Z=Z, u_address=u0@address, M=nrow(u0), J=J, 
                  group=group, family=family, link_int=link_int, coef=coef, offset=offset_fit,
                  phi=phi, maxit_CD=maxit_CD_use, conv_CD=conv_CD, init=(i == 1), group_X=group_X, 
                  covgroup=covgroup, penalty=penalty, lambda0=lambda0, lambda1=lambda1, 
                  gamma=gamma_penalty, alpha=alpha, trace=trace)
    # cat("End M-step \n")
    
    if(trace >= 1){
      cat("Fixed effects (scaled X): \n")
      cat(coef[1:ncol(X)], "\n")
    }
    
    # Re-calculate random effects covariance matrix from M step coefficients
    gamma = matrix(J%*%matrix(coef[-c(1:ncol(X))], ncol = 1), ncol = ncol(Z)/d)
    cov = var = gamma %*% t(gamma)
    
    if(trace >= 1){
      if(nrow(cov) <= 5){
        cat("random effect covariance matrix: \n")
        print(cov)
      }else{
        cat("random effect covariance matrix diagonal: \n")
        cat(diag(cov), "\n")
      }
    }
    
    if(family == "gaussian"){
      # if family = 'gaussian', calculate sigma of residual error term (standard deviation)
      sig_g = sig_gaus(y, X, Z, u0@address, group, J, coef, offset_fit, c(d, ncol(Z)/d), link_int)
      if(trace >= 1){
        cat("sig_g: ", sig_g, "\n")
      }
     
    }
    
    # if(family == "negbin"){
    #   phi = coef[length(coef)]
    #   coef = coef[-length(coef)]
    # }
    
    if(any(is.na(coef)) | any(abs(coef) > 10^5)){ # !is.finite(ll0) | 
      # For now, assume divergent in any abs(coefficient) > 10^5
      problem = T
      cat("Updated coef: \n")
      cat(coef, "\n")
    }
    
    # Check for errors in M step
    ## If errors, output warnings and enough relevant output to let select_tune() 
    ## (or glmm()) continue 
    if(problem == T){
      warning("Error in M step, see optinfo$warnings in output for details", immediate. = T)
      out = list(coef = coef, sigma = cov, Gamma_mat = gamma,
                 lambda0 = lambda0, lambda1 = lambda1, 
                 covgroup = covgroup, J = J, ll = -Inf, 
                 BICh = Inf, BIC = Inf, BICq = Inf, BICNgrp = Inf,
                 EM_iter = i, EM_conv = NA, converged = 0, 
                 u_init = NULL)
      # !is.finite(ll0) | any(is.na(coef)) | any(abs(coef) > 10^5)
      if(any(is.na(coef))){
        out$warnings = sprintf("coefficient estimates contained NA values at iteration %i",i)
      }else if(any(abs(coef) > 10^5)){
        if(family == "gaussian"){
          warning("Error in M step: coefficient values diverged. Consider increasing 'var_start' value in optimControl", immediate. = T)
          out$warnings = "Error in M step: coefficient values diverged. Consider increasing 'var_start' value in optimControl"
        }else{
          out$warnings = "Error in M step: coefficient values diverged"
        }
      }
      
      if(sampler %in% c("random_walk","independence")){
        out$gibbs_accept_rate = gibbs_accept_rate
        out$proposal_SD = proposal_SD
      }
      
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
      
      if(family == "gaussian"){
        out$sigma_gaus = sig_g
      }
      
      return(out)
    } # End problem == T
    
    
    
    ############################################################################################
    # Convergence Check
    ############################################################################################
    # stopping rule: based on average Euclidean distance (comparing coef from minus t iterations)
    if(i <= t){
      diff[i] = 10^2
      if(conv_type == 2){
        Q_est[i] = Qfun(y, X, Z, u0@address, group, J, matrix(coef, ncol = 1), offset_fit, c(d,ncol(Z)/d), family, link_int, sig_g, phi)
      }
    }else{
      if(conv_type == 1){
        diff[i] = sqrt(sum((coef - coef_record[1,])^2)) / sum(coef_record[1,] != 0)
      }else if(conv_type == 2){
        Q_est[i] = Qfun(y, X, Z, u0@address, group, J, matrix(coef, ncol = 1), offset_fit, c(d,ncol(Z)/d), family, link_int, sig_g, phi)
        diff[i] = abs(Q_est[i] - Q_est[i-t]) / abs(Q_est[i])
      }
    }
    
    # Update latest record of coef
    coef_record = rbind(coef_record[-1,], t(coef))
    
    # now update EM iteration information  
    update_limits = c(i, nMC , round(diff[i],6), (sum(coef[1:ncol(X)] !=0)), (sum(diag(cov) != 0))) # (sum(coef[-c(1:ncol(X))] != 0))
    names(update_limits) = c("Iter","nMC","EM conv","Non0 Fixef", "Non0 Ranef")
    print(update_limits)
    
    if(sum(diff[max(i-mcc+1,1):i] < conv_EM) >= mcc){
      EM_converged = 1
      break
    } 
    
    ############################################################################################
    # Update number of posterior samples, nMC
    ############################################################################################
    # Increase nMC in nMC below nMC_max
    ## Increase nMC by a multiplicative factor. 
    ## The size of this factor depends on the EM iteration (similar to mcemGLM package recommendations)
    if(i <= 15){
      nMC_fact = 1.1
    }else{
      nMC_fact = 1.20
    }
    # nMC_fact = 1.2
    
    if(nMC < nMC_max){
      nMC = round(nMC * nMC_fact)
      # If after update nMC exceeds nMC_max, reduce to nMC_max value
      if(nMC > nMC_max) nMC = nMC_max
    }
    
    ############################################################################################
    # E Step
    ############################################################################################
    
    # Update Znew2: For each group k = 1,...,d, calculate Znew2 = Z %*% gamma (calculate for each group individually)
    # Used within E-step
    Znew2 = Z
    for(j in 1:d){
      Znew2[group == j,seq(j, ncol(Z), by = d)] = Z[group == j,seq(j, ncol(Z), by = d)]%*%gamma
    }
    
    # Initial points for Metropolis within Gibbs E step algorithms
    uold = as.numeric(u0[nrow(u0),])
    # if random effect penalized out in past model / in previous M-step, do not
    # collect posterior samples for this random effect
    ranef_idx = which(diag(cov) > 0)
    # cat("Start E-step \n")
    Estep_out = E_step(coef = coef, ranef_idx = ranef_idx, y=y, X=X, Znew2=Znew2, group=group, offset_fit = offset_fit,
                       nMC=nMC, nMC_burnin=nMC_burnin, family=family, link=link, phi=phi, sig_g=sig_g,
                       sampler=sampler, d=d, uold=uold, proposal_SD=proposal_SD, 
                       batch=batch, batch_length=batch_length, offset_increment=offset_increment, 
                       trace=trace)
    # cat("End E-step \n")
    u0 = attach.big.matrix(Estep_out$u0)
    proposal_SD = Estep_out$proposal_SD
    gibbs_accept_rate = Estep_out$gibbs_accept_rate
    batch = Estep_out$updated_batch
    
  }
  
  ############################################################################################
  # End of EM algorithm
  ############################################################################################
  
  if(EM_converged == 0){
    warning("glmmPen algorithm did not converge within maxitEM iterations, conv = ", round(diff[i],6), 
            "\n Consider increasing maxitEM iterations or nMC_max in optimControl()", immediate. = T)
  }
  
  ############################################################################################
  # Log-likelihood and BIC calculations
  ############################################################################################
  
  
  # Another E step
  ## if logLik_calc = T, use M draws to calculate the log-likelihood
  ## else if logLik_calc = F but family == "gaussian", calculate 10^3 draws if necessary and 
  ##  use these draws in next round of selection to initialize sig_g value
  if(logLik_calc | ((family == "gaussian") & (nrow(u0) < 10^3))){
    Estep_out = E_step(coef=coef, ranef_idx=which(diag(cov) > 0), y=y, X=X, Znew2=Znew2, group=group, offset_fit = offset_fit,
                       nMC=ifelse(logLik_calc, M, 10^3), 
                       nMC_burnin=nMC_burnin, family=family, link=link, phi=phi, sig_g = sig_g,
                       sampler=sampler, d=d, uold=as.numeric(u0[nrow(u0),]), proposal_SD=proposal_SD, 
                       batch=batch, batch_length=batch_length, offset_increment=offset_increment, 
                       trace=trace)
    
    u0 = attach.big.matrix(Estep_out$u0)
    proposal_SD = Estep_out$proposal_SD
    gibbs_accept_rate = Estep_out$gibbs_accept_rate
    batch = Estep_out$updated_batch
    
    if(sampler == "random_walk" & trace >= 2){
      cat("Updated proposal_SD: ", proposal_SD, "\n")
      cat("Updated batch: ", batch, "\n")
    }
  }
  
  
  # Calculate loglik using Pajor corrected arithmetic mean estimator with
  # importance sampling weights (see "logLik_Pajor.R")
  if(logLik_calc){
    ll = CAME_IS(posterior = u0, y = y, X = X, Z = Z, group = group,
                 coef = coef, sigma = cov, family = fam_fun, M = M, gaus_sig = sig_g, trace = trace)
    
    # Hybrid BIC (Delattre, Lavielle, and Poursat (2014))
    # d = nlevels(group) = number independent subjects/groups
    BICh = -2*ll + sum(coef[-c(1:ncol(X))] != 0)*log(d) + sum(coef[1:ncol(X)] != 0)*log(nrow(X))
    # Usual BIC
    BIC = -2*ll + sum(coef != 0)*log(nrow(X))
    # BIC using N = nlevels(group)
    BICNgrp = -2*ll + sum(coef != 0)*log(d)
  }else{
    ll = NA
    BICh = NA
    BIC = NA
    BICNgrp = NA
  }
  
  # Calculate BICq
  if(!is.null(ufull_describe)){
    ## Using posterior draws from the full model (ufull) and the coefficients from the
    ## current penalized model, calculate the Q function
    ufull = attach.big.matrix(ufull_describe)
    q1 = Qfun(y, X, Z, ufull@address, group, J, matrix(coef, ncol = 1), offset_fit, c(d,ncol(Z)/d), family, link_int, sig_g, phi)
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
  
  ############################################################################################
  # Final output object
  ############################################################################################
  
  out = list(coef = coef, sigma = cov, Gamma_mat = gamma,
             lambda0 = lambda0, lambda1 = lambda1, 
             covgroup = covgroup, J = J, ll = ll, BICh = BICh, BIC = BIC, BICq = BICq, 
             BICNgrp = BICNgrp, EM_iter = i, EM_conv = diff[i], converged = EM_converged) 
  
  # If gaussian family, take 1000 draws of last u0 for u_init in next round of selection 
  #  (to use for sig_g calculation)
  # Else, take last row for initializing E step in next round of selection
  if(family == "gaussian"){
    out$u_init = u0[(nrow(u0)-999):nrow(u0),]
  }else{
    out$u_init = matrix(u0[nrow(u0),], nrow = 1)
  }
  
  if(sampler %in% c("random_walk","independence")){
    out$gibbs_accept_rate = gibbs_accept_rate
    out$proposal_SD = proposal_SD
    out$updated_batch = batch
  }else if(sampler == "stan"){
    out$gibbs_accept_rate = NULL
    out$proposal_SD = NULL
    out$updated_batch = NULL
  }
  
  if((is.null(coef_old)) & (trace >= 1)){
    out$coef_naive = fit_naive
  }
  if(EM_converged == 0){
    out$warnings = "glmmPen algorithm did not converge within maxit_EM iterations, consider increasing maxitEM or nMC_max in optimControl()"
  }
  
  # if(family == "negbin"){
  #   out$phi = phi
  # }
  
  if(family == "gaussian"){
    out$sigma_gaus = sig_g
  }
  
  
  # cat(Q_est, "\n")
  return(out)
}





#####################################################################################################


# negative binomial. In order to select the
# negative binomial distribution, either set family = "negbin" (which assumes a canonical log link)
# or family = \code{MASS::negative.binomial} using any arbitrary \code{theta} value and the desired
# link function.

