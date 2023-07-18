

# @title Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
# Minimization (MCECM) assuming a factor model assumption on the random effects
# 
# \code{fit_dat_FA} is used to fit a penalized generalized mixed model
# via Monte Carlo Expectation Conditional Minimization (MCECM) for 
# a single tuning parameter combinations and is called within
# \code{glmmPen_FA} or \code{glmm_FA} (cannot be called directly by user)
# 
# @inheritParams lambdaControl
# @inheritParams glmmPen
# @inheritParams glmmPen_FA
# @inheritParams fit_dat
# @param beta_old vector giving values to initialize the fixed effects coefficients 
# @param b_old vector giving values to initialize the random effects coefficients, which 
# are values of the B matrix
# @param r positive integer specifying number of latent common factors to assume 
# in the model.
# 
# @return a list with the following elements:
# \item{coef}{a numeric vector of coefficients of fixed effects estimates and 
# non-zero estimates of the lower-triangular cholesky decomposition of the random effects
# covariance matrix (in vector form)}
# \item{sigma}{random effects covariance matrix}
# \item{lambda0, lambda1}{the penalty parameters input into the function}
# \item{covgroup}{Organization of how random effects coefficients are grouped.}
# \item{J}{a sparse matrix that transforms the non-zero elements of the lower-triangular cholesky 
# decomposition of the random effects covariance matrix into a vector. For unstructured
# covariance matrices, dimension of dimension q^2 x (q(q+1)/2) (where q = number of random effects).
# For independent covariance matrices, q^2 x q.}
# \item{ll}{estimate of the log likelihood, calculated using the Pajor method}
# \item{BICh}{the hybrid BIC estimate described in Delattre, Lavielle, and Poursat (2014)}
# \item{BIC}{Regular BIC estimate}
# \item{BICNgrps}{BIC estimate with N = number of groups in penalty term instead of N = number
# of total observations.}
# \item{BICq}{BIC-ICQ estimate}
# \item{u}{a matrix of the Monte Carlo draws. Organization of columns: first by random effect variable,
# then by group within variable (i.e. Var1:Grp1 Var1:Grp2 ... Var1:GrpK Var2:Grp1 ... Varq:GrpK)}
# \item{gibbs_accept_rate}{a matrix of the ending gibbs acceptance rates for each variable (columns)
# and each group (rows) when the sampler is either "random_walk" or "independence"}
# \item{proposal_SD}{a matrix of the ending proposal standard deviations (used in the adaptive
# random walk version of the Metropolis-within-Gibbs sampling) for each variable (columns) and
# each group (rows)}
#' @useDynLib glmmPen
#' @importFrom bigmemory attach.big.matrix describe as.big.matrix
#' @importFrom mvtnorm dmvnorm
#' @importFrom Matrix Matrix
fit_dat_FA = function(dat, lambda0 = 0, lambda1 = 0, 
                     family = "binomial", offset_fit = NULL,
                     optim_options = optimControl(nMC_start = 250, nMC_max = 1000, nMC_burnin = 100),
                     trace = 0, penalty = c("MCP","SCAD","lasso"),
                     alpha = 1, gamma_penalty = switch(penalty[1], SCAD = 4.0, 3.0), 
                     group_X = 0:(ncol(dat$X)-1),
                     u_init = NULL, beta_old = NULL, b_old = NULL, # coef_old = NULL, 
                     ufull_describe = NULL, 
                     adapt_RW_options = adaptControl(), 
                     r, logLik_calc = FALSE, checks_complete = FALSE,
                     ranef_keep = rep(1, times = (ncol(dat$Z)/nlevels(dat$group))), 
                     progress = TRUE){
  
  ############################################################################################
  # Extract optimization parameters
  ############################################################################################
  
  # Extract variables from optimControl
  conv_EM = optim_options$conv_EM
  conv_CD = optim_options$conv_CD
  nMC_burnin = optim_options$nMC_burnin
  nMC = optim_options$nMC
  nMC_max = optim_options$nMC_max
  maxitEM = optim_options$maxitEM
  maxit_CD = optim_options$maxit_CD
  M = optim_options$M
  t = optim_options$t
  mcc = optim_options$mcc
  sampler = optim_options$sampler
  step_size = optim_options$step_size
  convEM_type = optim_options$convEM_type
  B_init_type = optim_options$B_init_type
  var_restrictions = optim_options$var_restrictions
  
  ############################################################################################
  # Data input checks
  ############################################################################################
  
  # if calling fit_dat_FA directly instead of within glmm or glmmPen, perform data checks
  
  y = dat$y
  X = base::as.matrix(dat$X)
  # Convert sparse Z to dense ZF
  Z = Matrix::as.matrix(dat$Z)
  group = dat$group
  coef_names = dat$coef_names # list with items 'fixed', 'random', and 'group'
  d = nlevels(factor(group))
  q = ncol(Z) / d
  
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
    
    # Check sampler specification
    sampler = checkSampler(sampler)
    
  } # End checks of input
  
  
  # Extract relevant family information
  family_info = family_export(family)
  fam_fun = family_info$family_fun
  link = family_info$link
  link_int = family_info$link_int # Recoded link as integer, see "family_export.R" for details
  family = family_info$family
  
  if(family == "gaussian"){
    stop("gaussian familiy not yet available for glmmPen_FA function, please use glmmPen function")
  }
  
  # Set small penalties to zero
  if(lambda0 <=10^-6) lambda0 = 0
  if(lambda1 <=10^-6) lambda1 = 0
  
  # Create J_b matrix (to use in t(f kronecker z_ki) * J_b part of linear predictor)
  # vec(B) = J_b b = J_b vec(t(B))
  # Formula from https://math.stackexchange.com/questions/307299/kronecker-product-and-the-commutation-matrix
  ## Note: Taking transpose of matrix described in link. Formula given in link:
  ## Consider K(q,r) vec(B) = vec(t(B)) for matrix B of dim q x r
  ## For row k, column ((k-1) mod r)*q + floor((k-1)/r) + 1 will have entry 1
  # Transpose the above to get desired J_b = K(r,q) such that vec(B) = J_b vec(t(B))
  ## Note: Transpose of matrix = inverse of matrix in this case
  J = Matrix(0,q*r, q*r, sparse = TRUE)
  for(i in 1:(q*r)){
    idx = ((i-1) %% r)*q + floor((i-1)/r) + 1
    J[idx,i] = 1
  }
  
  ############################################################################################
  # Initialization
  ############################################################################################
  
  # Initialize coef for start of EM algorithm
  ## coef: Initializes coefficients for M step: c(beta0, b0)
  ## beta0: Initializes fixed effects coefficients
  ## b0: Initializes coefficients of the B matrix (where cov = B %*% t(B))
  
  # Messages
  if((!is.null(beta_old)) & (!is.null(b_old))){
    if(progress == TRUE) message("using coefficients from past model to intialize")
  }else if(!is.null(beta_old)){
    if(progress == TRUE) message("using fixed effects coefficients from past model to intialize")
  }else if(!is.null(b_old)){
    if(progress == TRUE) message("using random effects coefficients from past model to intialize")
  }
  
  if(!is.null(beta_old)){
    if(length(beta_old) != ncol(X)){
      stop("issues with dimensions of initialized fixed effects: length of beta_old is not equal to number columns of X \n")
    }
    beta0 = beta_old
  }else{
    # Coordinate descent ignoring random effects: naive fit
    
    penalty_factor = numeric(ncol(X))
    penalty_factor[which(group_X == 0)] = 0
    penalty_factor[which(group_X != 0)] = 1
    
    # Initialize coef as intercept-only for input into CD function
    IntOnly = glm(y ~ 1, family = fam_fun, offset = offset_fit)
    coef_init = c(IntOnly$coefficients, rep(0, length = (ncol(X)-1)))
    
    fit_naive = CD(y, X, family = family, link_int = link_int, offset = offset_fit, coef_init = coef_init,
                   maxit_CD = maxit_CD, conv = conv_CD, penalty = penalty, lambda = lambda0,
                   gamma = gamma_penalty, alpha = alpha, penalty_factor = penalty_factor, trace = trace)
    beta0 = fit_naive
    
    
    if(trace >= 1){
      cat("initialized fixed effects: \n")
      cat(beta0, "\n")
    } 
    
    if(any(is.na(beta0))){
      cat(beta0)
      stop("Error in initial coefficient fit: NAs produced")
    }
    
  } # end if-else !is.null(beta_old)
  
  # Initialization of B matrix (ranef covariance = B %*% t(B))
  if(!is.null(b_old)){
    if(length(b_old) != q*r){
      stop("length of b_old equals ", length(b_old), " but should equal q*r: ", q*r)
    }
    b0 = b_old
    B = t(matrix(b0, nrow = r, ncol = q))
    
    # If element j of ranef_keep equals 1, this indicates that the corresponding random effect j
    # has not been eliminated from consideration.
    # Consequently, the jth row of B should not be penalized to 0 (should be initialized as non-zero)
    # If set to zero, initialized with set deterministic values
    # Alternatively, if pre-screening step eliminates a random effect from consideration but
    # the row of B is not zero, set to 0 (would occur if row of B was non-zero after pre-screening
    # but the resulting variance in the covariance matrix was < 10^-2)
    for(j in 1:q){
      if((ranef_keep[j] == 1) & (all(B[j,] == 0))){
        if(B_init_type == "random"){
          b_init_vec = rnorm(n = r, mean = 0, sd = 0.25)
        }else if(B_init_type == "deterministic"){
          # Set values of B matrix so that all values of covariance matrix = 0.10 OR initialized value of var_start
          if(optim_options$var_start == "recommend"){
            var_start = 0.10
          }else{
            var_start = optim_options$var_start
          }
          b0_val = sqrt(var_start / r)
          b_init_vec = rep(b0_val, r)
        }else if(B_init_type == "data"){
          # Set values of B matrix so that all values of covariance matrix = value from var_init()
          var_start = var_init(dat, fam_fun)
          b0_val = sqrt(var_start / r)
          b_init_vec = rep(b0_val, r)
        }
        B[j,] = b_init_vec 
      }else if(ranef_keep[j] == 0){
        B[j,] = rep(0, r)
      }
    }
    b0 = c(t(B))
  }else{
    if(B_init_type == "random"){
      B = matrix(rnorm(n=q*r, mean = 0, sd = 0.25), nrow = q, ncol = r)
    }else if(B_init_type == "deterministic"){
      # Set values of B matrix so that all values of covariance matrix = 0.10 OR initialized value of var_start
      if(optim_options$var_start == "recommend"){
        var_start = 0.10
      }else{
        var_start = optim_options$var_start
      }
      b0_val = sqrt(var_start / r)
      B = matrix(b0_val, nrow = q, ncol = r)
    }else if(B_init_type == "data"){
      # Set values of B matrix so that all values of covariance matrix = value from var_init()
      var_start = var_init(dat, fam_fun)
      b0_val = sqrt(var_start / r)
      B = matrix(b0_val, nrow = q, ncol = r)
    }
    
    # apply variance restrictions if necessary based on initialized fixed effects
    if(var_restrictions == "fixef"){
      fixed_keep = coef_names$fixed[c(1,which(beta0[-1] != 0)+1)]
      restrict_idx = which(!(coef_names$random %in% fixed_keep))
      for(j in restrict_idx){
        B[j,] = rep(0, r)
      }
    }
   
    b0 = c(t(B))
  } # End if-else !is.null(b_old)
  
  # Calculate covariance matrix from initialized B matrix
  cov = B %*% t(B)
  
  # Initialize coefficient vector:
  coef = c(beta0, b0)
  
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
      eta = X %*% coef[1:ncol(X)]
      s2_g = sum((y - invlink(link_int, eta))^2)
      sig_g = sqrt(s2_g / length(y)) # Divide by N or not?
    }else{
      # Find initial estimate of variance of the error term using fixed and random effects from
      # last round of selection
      u_big = as.big.matrix(u_init)
      sig_g = sig_gaus_FA(y, X, Z, u_big@address, group, J, coef, offset_fit, c(d, ncol(Z)/d, r), link_int)

    }

    if(trace >= 1){
      cat("initial residual error SD: ", sig_g, "\n")
    }
    
  }else{
    sig_g = 1.0 # specify an arbitrary placeholder (will not be used in calculations)
  } # End if-else family == "gaussian"
  
  # Znew2: For each group k = 1,...,d, calculate Znew2 = Z %*% B (calculate for each group individually)
  # Used within E-step
  Znew2 = matrix(0, nrow = nrow(Z), ncol = d*r)
  B = t(matrix(coef[-c(1:ncol(X))], nrow = r, ncol = q))
  for(j in 1:d){
    Znew2[group == j,seq(j, ncol(Znew2), by = d)] = Z[group == j, seq(j, ncol(Z), by = d)] %*% B
  }
  
  # initialize adaptive Metropolis-within-Gibbs random walk parameters
  ## initialize proposal standard deviation
  proposal_SD = matrix(1.0, nrow = d, ncol = ncol(Z)/d)
  if(sampler == "random_walk" & trace >= 2){
    cat("initialized proposal_SD: \n")
    print(proposal_SD)
  }
  
  ## initialize batch number to 0
  batch = 0.0
  ## initialize other paramters from adaptControl()
  batch_length = adapt_RW_options$batch_length
  offset_increment = adapt_RW_options$offset
  gibbs_accept_rate = matrix(NA, nrow = d, ncol = nrow(Z)/d)
  
  # At start of EM algorithm, acquire posterior draws from all random effect variables
  # Note: restricted to just which variables were not selected out in previous selection models
  # ranef_idx = which(diag(cov) != 0)
  
  # For now, do not remove common factors of random effects from EM algorithm
  ranef_idx = 1:r
  
  if(!is.null(u_init)){
    if(progress == TRUE) message("using results from previous model to initialize posterior samples \n")
    if(is.matrix(u_init)){
      u_big = as.big.matrix(u_init)
      uold = as.numeric(u_big[nrow(u_big),])
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
    gibbs_accept_rate = Estep_out$proposal_SD
    batch = Estep_out$updated_batch
    
  }else{ # u_init is null, initialize posterior with random draw from N(0,1)
    
    Estep_out = E_step(coef = coef, ranef_idx = ranef_idx, y=y, X=X, Znew2=Znew2, group=group, offset_fit = offset_fit,
                       nMC=nMC, nMC_burnin=nMC_burnin, family=family, link=link, phi=phi, sig_g=sig_g,
                       sampler=sampler, d=d, uold=rnorm(n = ncol(Znew2)), proposal_SD=proposal_SD, 
                       batch=batch, batch_length=batch_length, offset_increment=offset_increment, 
                       trace=trace)
    
    u0 = attach.big.matrix(Estep_out$u0)
    proposal_SD = Estep_out$proposal_SD
    gibbs_accept_rate = Estep_out$proposal_SD
    batch = Estep_out$updated_batch
    
  } # end if-else is.null(u_init)
  
  if(nrow(cov) == 1){ # Single random intercept
    cov_record = rep(NA, maxitEM)
  }else{
    cov_record = NULL
  }
  
  diff = rep(NA, maxitEM)
  stopcount = 0
  
  # Record last t coef vectors (each row = coef vector for a past EM iteration)
  # Initialize with initial coef vector
  coef_record = matrix(coef, nrow = t, ncol = length(coef), byrow = TRUE)
  
  problem = FALSE # Determining if issue with EM algorithm results
  randInt_issue = 0 # Determine if issue where random-intercept only model has too small of variance
  
  # Remove initialized big matrix of posterior draws
  u_big = NULL
  
  # If using Q-function as convergence criteria, save 'old' Q-function value
  Q_est = rep(NA, maxitEM)
  
  ############################################################################################
  # EM Algorithm
  ############################################################################################
  
  # Start EM Algorithm
  EM_converged = 0
  
  for(i in 1:maxitEM){
    
    ############################################################################################
    # M Step
    ############################################################################################
    
    M_out = M_step_FA(y=y, X=X, Z=Z, u_address=u0@address, M=nrow(u0), J=J, 
                  group=group, family=family, link_int=link_int, coef=coef, offset=offset_fit,
                  sig_g = sig_g, phi=phi, maxit_CD=maxit_CD, conv_CD=conv_CD, 
                  init=(i == 1), group_X=group_X, 
                  penalty=penalty, lambda0=lambda0, lambda1=lambda1, 
                  gamma=gamma_penalty, alpha=alpha, step_size=step_size, trace=trace)
    
    coef = M_out$coef_new
    step_size = M_out$step_size
    phi = M_out$phi # relevant when family = "negbin"
    
    if((trace >= 1) & (family == "poisson")){
      cat("step size: ", step_size,"\n")
    }
    
    if(family == "gaussian"){
      # if family = 'gaussian', calculate sigma of error term (standard deviation)
      sig_g = sig_gaus_FA(y, X, Z, u0@address, group, J, coef, offset_fit, c(d, ncol(Z)/d, r), link_int)
      if(trace >= 1){
        cat("sig_g: ", sig_g, "\n")
      }
    }
    
    B = t(matrix(coef[-c(1:ncol(X))], nrow = r, ncol = q))
    cov = B %*% t(B)
    
    # If trace >= 1, output most recent coefficients
    if(trace >= 1){
      cat("Fixed effects (scaled X): \n")
      cat(round(coef[1:ncol(X)],3), "\n")
    }
    
    if(trace >= 1){
      if(nrow(cov) <= 5){
        cat("random effect covariance matrix: \n")
        print(round(cov,3))
      }else{
        cat("random effect covariance matrix diagonal: \n")
        cat(round(diag(cov),3), "\n")
      }
    }
    
    
    # Check of size of random intercept variance in models with only random intercept,
    #   no random slopes. Cases: pre-specified random effects only includes random
    #   intercept, OR all random slopes have been penalized out of model
    if((nrow(cov) == 1) | all(diag(cov[-1,-1,drop=FALSE]) == 0)){
      if(cov[1,1] < 0.001){ # 0.01
        randInt_issue = 1
      }
    }else{
      randInt_issue = 0
    }
    
    if(any(is.na(coef)) | any(abs(coef) > 10^5) | randInt_issue == 1){  
      # For now, assume divergent in any abs(coefficient) > 10^5
      problem = TRUE
      cat("Updated coef: ", coef, "\n")
    }
    
    # Check for errors in M step
    ## If errors, output warnings and enough relevant output to let select_tune() 
    ## (or glmm()) continue 
    if(problem == TRUE){
      out = list(coef = coef, sigma = cov, B = B,
                 lambda0 = lambda0, lambda1 = lambda1, 
                 J = J, ll = -Inf, 
                 BICh = Inf, BIC = Inf, BICq = Inf, BICNgrp = Inf,
                 EM_iter = i, EM_conv = NA, converged = 0, 
                 u_init = NULL)
      if(any(is.na(coef))){
        warning(sprintf("coefficient estimates contained NA values at iteration %i",i), immediate. = TRUE)
        out$warnings = sprintf("coefficient estimates contained NA values at iteration %i",i)
      }else if(any(abs(coef) > 10^5)){
        if(family == "gaussian"){
          warning("Error in M step: coefficient values diverged. Consider increasing 'var_start' value in optimControl", immediate. = TRUE)
          out$warnings = "Error in M step: coefficient values diverged. Consider increasing 'var_start' value in optimControl"
        }else{
          warning("Error in M step: coefficient values diverged", immediate. = TRUE)
          out$warnings = "Error in M step: coefficient values diverged"
        }
      }else if(randInt_issue == 1){
        warning("Error in model fit: Random intercept variance is too small, indicating either that 
                there are high correlations among the covariates (if so, consider reducing these correlations 
                or changing the Elastic Net alpha value) or that this model should be fit 
                using traditional generalized linear model techniques.", immediate. = TRUE)
        out$warnings = "Error in model fit: random intercept variance becomes too small, model should be fit using regular generalized linear model techniques"
      }
      
      if(sampler %in% c("random_walk","independence")){
        out$gibbs_accept_rate = gibbs_accept_rate
        out$proposal_SD = proposal_SD
      }
      
      if(is.null(beta_old)){
        out$coef_naive = fit_naive
      }
      
      if(family == "gaussian"){
        out$sigma_gaus = sig_g
      }
      
      return(out)
      
      
      return(out)
    } # End problem == T
    
    
    
    if(nrow(cov) == 1){
      cov_record[i] = cov
    }
    
    # Znew2: For each group k = 1,...,d, calculate Znew2 = Z %*% B (calculate for each group individually)
    # Used within E-step
    ## B calculated above
    Znew2 = matrix(0, nrow = nrow(Z), ncol = d*r)
    for(j in 1:d){
      Znew2[group == j,seq(j, ncol(Znew2), by = d)] = Z[group == j, seq(j, ncol(Z), by = d)] %*% B
    }
    
    ############################################################################################
    # Convergence Check
    ############################################################################################
    # stopping rule: based on average Euclidean distance (comparing coef from minus t iterations)
    if(i <= t){
      diff[i] = 10^2
      if(convEM_type == "Qfun"){
        Q_est[i] = Qfun_FA(y, X, Z, u0@address, group, J, coef, offset_fit, c(d, q, r), family, link_int, sig_g, phi)
      }
    }else{
      if(convEM_type == "AvgEuclid1"){
        diff[i] = sqrt(sum((coef - coef_record[1,])^2)) / sum(coef_record[1,] != 0)
      }else if(convEM_type == "AvgEuclid2"){
        diff[i] = sqrt(sum((coef - coef_record[1,])^2)) / sqrt(sum(coef_record[1,] != 0))
      }else if(convEM_type == "Qfun"){
        Q_est[i] = Qfun_FA(y, X, Z, u0@address, group, J, coef, offset_fit, c(d, q, r), family, link_int, sig_g, phi)
        diff[i] = abs(Q_est[i] - Q_est[i-t]) / abs(Q_est[i])
      }else if(convEM_type == "maxdiff"){
        diff[i] = max(abs((coef - coef_record[1,])))
      }
    }
    
    # Update latest record of coef
    coef_record = rbind(coef_record[-1,], t(coef))
    
    # now update EM iteration information  
    if(progress == TRUE){
      update_limits = c(i, nMC , round(diff[i],6), 
                        max((sum(coef[1:ncol(X)] !=0))-1,0), 
                        max((sum(diag(cov) != 0))-1,0))
      names(update_limits) = c("Iter","nMC","EM conv","Non0 Fixef", "Non0 Ranef")
      print(update_limits)
    }
    
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
    
    # Initial points for E step sampling algorithms
    uold = as.numeric(u0[nrow(u0),])
    # ranef_idx = which(diag(cov) > 0)
    ranef_idx = 1:r
    # print("Start E-step")
    Estep_out = E_step(coef = coef, ranef_idx = ranef_idx, y=y, X=X, Znew2=Znew2, group=group, offset_fit = offset_fit,
                       nMC=nMC, nMC_burnin=nMC_burnin, family=family, link=link, phi=phi, sig_g=sig_g,
                       sampler=sampler, d=d, uold=uold, proposal_SD=proposal_SD, 
                       batch=batch, batch_length=batch_length, offset_increment=offset_increment, 
                       trace=trace)
    # print("End E-step")
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
            "\n Consider increasing maxitEM iterations or nMC_max in optimControl()", immediate. = TRUE)
  }
  
  ############################################################################################
  # Log-likelihood and BIC calculations
  ############################################################################################
  
  # Another E step
  ## if logLik_calc = T, use M draws to calculate the log-likelihood
  ## else if logLik_calc = F but family == "gaussian", calculate 10^3 draws if necessary and 
  ##  use these draws in next round of selection to initialize sig_g value
  
  
  # Note: Assume that at end of EM algorithm, which(diag(cov) > 0) corresponds with the
  # sampling restrictions specified within EM algorithm:
    # Restrict sampling such that variables are NOT sampled if:
    ## random effect for the variable was penalized to zero in a previous M step
    ## fixed effect for the variable was penalized to zero in a prevoius M step
    ## Make sure to always keep random intercept
  
  if(logLik_calc | ((family == "gaussian") & (nrow(u0) < 10^3))){
    Estep_out = E_step(coef=coef, ranef_idx=1:r, y=y, X=X, Znew2=Znew2, 
                       group=group, offset_fit = offset_fit,
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
      cat("Updated proposal_SD: \n")
      print(proposal_SD)
      cat("Updated batch: ", batch, "\n")
    }
  }
  
  
  # Calculate loglik using Pajor method (see "logLik_Pajor.R")
  if(logLik_calc){
    # stop("log-likelihood calculation not yet finalized for factor analysis version")
    ll = CAME_IS_FA(posterior = u0, y = y, X = X, Z = Z, group = group,
                    coef = coef, B = t(matrix(coef[-c(1:ncol(X))], ncol = q, nrow = r)), 
                    family = fam_fun, M = M, gaus_sig = sig_g, 
                    trace = trace, progress = progress)

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
  
  # Calculated BICq
  ## Using posterior draws from the full model (ufull) and the coefficients from the
  ## current penalized model, calculate the Q function
  if(!is.null(ufull_describe)){
    ufull = attach.big.matrix(ufull_describe)
    q1 = Qfun_FA(y, X, Z, ufull@address, group, J, coef, offset_fit, c(d, q, r), family, link_int, sig_g, phi)
    q2 = 0
    for(k in 1:d){
      cols_idx = seq(from = k, to = ncol(ufull), by = d)
      post_k = ufull[,cols_idx]
      q2 = q2 + sum(dmvnorm(post_k, log=TRUE)) / nrow(ufull)
    }
    llq = q1 + q2
    BICq = -2*llq + sum(coef != 0)*log(nrow(X))
    if(is.na(BICq)){
      warning("BICq value calculated as NA due to divergent coefficient values. 
              Consider checking correlations in the covariates and/or adjusting
              model fit parameters.", immediate. = TRUE)
      cat("Fixed effects (scaled X): \n")
      cat(round(coef[1:ncol(X)],3), "\n")
      
      if(nrow(cov) <= 5){
        cat("random effect covariance matrix: \n")
        print(round(cov,3))
      }else{
        cat("random effect covariance matrix diagonal: \n")
        cat(round(diag(cov),3), "\n")
      }
    }
  }else{
    BICq = NA
  }
  # message("BICq value outside if statement: ", BICq)
  
  ############################################################################################
  # Final output object
  ############################################################################################
  
  out = list(coef = coef, sigma = cov,
             lambda0 = lambda0, lambda1 = lambda1, 
             J = J, ll = ll, BICh = BICh, BIC = BIC, BICq = BICq, 
             BICNgrp = BICNgrp, EM_iter = i, EM_conv = diff[i],
             converged = EM_converged) 
  
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
  
  if((is.null(beta_old)) & (trace >= 1)){
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
  
 # out$coef_record_all = coef_record_all
  if(!is.null(cov_record)){
    out$cov_record = cov_record
  }
  
  out$B = t(matrix(coef[-c(1:ncol(X))], nrow = r, ncol = q))
  
  # cat(Q_est, "\n")
  return(out)
}


