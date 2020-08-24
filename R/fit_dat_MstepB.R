
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
#' greater than 7 (i.e. high dimensional), then the algorithm automatically assumes an 
#' independent covariance structure (covar switched to "independent").
#' @param fit_type integer specifying the verion of the grouped coordinate descent algorithm to use.
#' Default equals 1, the 'naive fit', assuming that the random effects covariates do not need any
#' additional standardization or adjustments to account for absence of orthogonalization.
#' For the time being, this version is recommended. Version 2 adds standardization of the random
#' effects covariates compared to version 1. 
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
#' @useDynLib glmmPen
#' @importFrom bigmemory attach.big.matrix describe big.matrix
#' @importFrom ncvreg ncvreg
#' @importFrom rstan sampling extract
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
  
  penalty = penalty[1]
  if(!(penalty %in% c("lasso","MCP","SCAD"))){
    stop("penalty ", penalty, " not available, must choose 'lasso', 'MCP', or 'SCAD' \n")
  }
  
  # Set small penalties to zero
  if(lambda0 <=10^-6) lambda0 = 0
  if(lambda1 <=10^-6) lambda1 = 0
  
  y = dat$y
  X = base::as.matrix(dat$X)
  # Convert sparse Z to dense Z
  Z = Matrix::as.matrix(dat$Z)
  group = dat$group
  d = nlevels(factor(group))
  
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
  
  sampler = sampler[1] # Default of random walk
  if(!(sampler %in% c("independence", "random_walk","stan"))){
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
    
    ok = which(diag(var) > 0)
    okindex = NULL
    for(i in 1:(ncol(Z)/d)){
      if(i %in% ok){
        okindex = c(okindex, (i-1)*d + 1:d)
      }
    }
    
    coef = coef_old[1:ncol(X)]
    
  }else{
    
    # Coordinate descent ignoring random effects: naive fit
    
    penalty_factor = numeric(ncol(X)-1)
    penalty_factor[which(group_X[-1] == 0)] = 0
    penalty_factor[which(group_X[-1] != 0)] = 1
    
    fit_naive = ncvreg(as.matrix(X[,-1]), y, family = family, penalty = penalty, gamma = gamma_penalty,
                       alpha = alpha, lambda = lambda0, penalty.factor = penalty_factor)
    
    coef = as.numeric(fit_naive$beta)
    
    if(trace >= 1) print(coef)
    
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
    
    ok = which(vars > 0) 
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
  
  # initialize adaptive Metropolis-within-Gibbs random walk parameters
  ## initialize proposal standard deviation
  proposal_SD = matrix(1.0, nrow = d, ncol = ncol(Z)/d)
  if(sampler != "stan" & trace >= 1){
    print("initialized proposal_SD:")
    print(proposal_SD)
  }
  
  ## initialize batch number to 0
  batch = 0.0
  ## initialize other paramters from adaptControl()
  batch_length = adapt_RW_options$batch_length
  offset_increment = adapt_RW_options$offset
  burnin_batchnum = adapt_RW_options$burnin_batchnum
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
    
    if((sampler == "independence")){
      # No burnin adjustments used
      samplemc_out = sample_mc2_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group, 
                                       d = d, okindex = okindex, trace = trace, uold = uold)
      u0 = attach.big.matrix(samplemc_out$u0)
      
    }else if(sampler == "random_walk"){ 
      # First adapt the proposal variance during an initial burnin period
      samplemc_burnin = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, Z=Znew2, nMC=nMC_burnin, family = family, group = group,
                                               d = d, okindex = okindex, trace = trace, uold = uold,
                                               proposal_SD = proposal_SD, batch = batch, batch_length = batch_length,
                                               offset = offset_increment, burnin_batchnum = burnin_batchnum)
      
      # Extract updated info from burnin period
      batch = samplemc_burnin$updated_batch
      proposal_SD = samplemc_burnin$proposal_SD
      
      u0 = attach.big.matrix(samplemc_burnin$u0)
      
      if(trace >= 1){
        print("Updated proposal_SD:")
        print(proposal_SD)
      }
      
      # Run the official first E step with no adaptation
      uold = as.numeric(u0[nrow(u0),])
      samplemc_out = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group,
                                            d = d, okindex = okindex, trace = trace, uold = uold,
                                            proposal_SD = proposal_SD, batch = batch, batch_length = batch_length,
                                            offset = offset_increment, burnin_batchnum = 0)
      
      gibbs_accept_rate = samplemc_out$gibbs_accept_rate
      u0 = attach.big.matrix(samplemc_out$u0)
      
    }else if(sampler == "stan"){
      
      u0 = big.matrix(nrow = nMC, ncol = ncol(Z)) # use ', init = 0' for sampling within EM algorithm
      
      if(family == "binomial" & link == "logit"){
        stan_file = stanmodels$binomial_logit_model
      }
      
      # Number draws to extract in each chain (after burn-in)
      nMC_chain = ceiling(nMC/num_cores)
      
      # Record last elements of each chain for initialization of next E step
      last_draws = matrix(0, nrow = num_cores, ncol = ncol(Z)) 
      
      # For each group, sample from posterior distribution (sample the alpha values)
      for(k in 1:d){
        idx_k = which(group == k)
        cols_k = seq(from = k, to = ncol(Z), by = d)
        y_k = y[idx_k]
        X_k = X[idx_k,]
        if(length(cols_k) == 1){
          Z_k = matrix(Znew2[idx_k, cols_k], nrow = length(idx_k), ncol = length(cols_k))
        }else{ # length(cols_k) > 1
          Z_k = Znew2[idx_k, cols_k]
        }
        
        dat_list = list(N = length(idx_k), # Number individuals in group k
                        q = length(ranef_idx), # number random effects
                        eta_fef = as.numeric(X_k %*% matrix(coef[1:ncol(X)], ncol = 1)), # fixed effects componenet of linear predictor
                        y = y_k, # outcomes for group k
                        Z = Z_k) # portion of Z matrix corresonding to group k
        
        # initialize posterior random draws
        init_lst = list()
        for(cr in 1:num_cores){
          init_lst[[cr]] = list(alpha = uold[cols_k])
        }
        
        stan_fit = rstan::sampling(stan_file, data = dat_list, init = init_lst, 
                                   iter = nMC_chain + nMC_burnin,
                                   warmup = nMC_burnin, show_messages = F, refresh = 0,
                                   chains = num_cores, cores = num_cores)
        
        stan_out = as.matrix(stan_fit)
        # Exclude lp__ column of output (log density up to a constant)
        draws_mat = stan_out[,1:length(ranef_idx)]
        
        # Find last elements for each chain for initialization of next E step
        for(m in 1:num_cores){
          last_draws[m,cols_k] = draws_mat[nMC_chain*m,]
        }
        
        if(nrow(draws_mat) == nMC){
          u0[,cols_k] = draws_mat
        }else{ # nrow(draws_mat) > nMC due to ceiling function in 'iter' specification
          start_row = nrow(draws_mat) - nMC + 1
          rows_seq = start_row:nrow(draws_mat)
          u0[,cols_k] = draws_mat[rows_seq,]
        }
        
        
      } # End k for loop
      
    } # End if-else sampler
   
    
  }else{ # u_init is null, initialize posterior with random draw from N(0,1)
    
    if(sampler == "independence"){
      samplemc_out = sample_mc2_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group, 
                                       d = d, okindex = okindex, trace = trace, 
                                       uold = rnorm(n = ncol(Z)))
      
      u0 = attach.big.matrix(samplemc_out$u0)
      
    }else if(sampler == "random_walk"){ 
      # First adapt the proposal variance during an initial burnin period
      samplemc_burnin = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, Z=Znew2, nMC=nMC_burnin, family = family, group = group,
                                               d = d, okindex = okindex, trace = trace,
                                               uold = rnorm(n = ncol(Z)),
                                               proposal_SD = proposal_SD, batch = batch, batch_length = batch_length,
                                               offset = offset_increment, burnin_batchnum = burnin_batchnum)
      
      # Extract updated proposal_SD info from burnin period
      batch = samplemc_burnin$updated_batch
      proposal_SD = samplemc_burnin$proposal_SD
      
      u0 = attach.big.matrix(samplemc_burnin$u0)
      
      if(trace >= 1){
        print("Updated proposal_SD:")
        print(proposal_SD)
      }
      
      # Run official first E step with no adaptation
      uold = as.numeric(u0[nrow(u0),])
      samplemc_out = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, Z=Znew2, nMC=nMC, family = family, group = group,
                                            d = d, okindex = okindex, trace = trace, uold = uold,
                                            proposal_SD = proposal_SD, batch = batch, batch_length = batch_length,
                                            offset = offset_increment, burnin_batchnum = 0)
      
      gibbs_accept_rate = samplemc_out$gibbs_accept_rate
      
      u0 = attach.big.matrix(samplemc_out$u0)
      
    }else if(sampler == "stan"){
      
      # stop("Stan not yet fully integrated")
      
      u0 = big.matrix(nrow = nMC, ncol = ncol(Z)) # use ', init = 0' for sampling within EM algorithm
      
      if(family == "binomial" & link == "logit"){
        stan_file = stanmodels$binomial_logit_model
      }
      
      # Number draws to extract in each chain (after burn-in)
      nMC_chain = ceiling(nMC/num_cores)
      
      # Record last elements for each chain for initialization of next E step
      last_draws = matrix(0, nrow = num_cores, ncol = ncol(Z))
      
      # For each group, sample from posterior distribution (sample the alpha values)
      for(k in 1:d){
        idx_k = which(group == k)
        cols_k = seq(from = k, to = ncol(Z), by = d)
        y_k = y[idx_k]
        X_k = X[idx_k,]
        if(length(cols_k) == 1){
          Z_k = matrix(Znew2[idx_k, cols_k], nrow = length(idx_k), ncol = length(cols_k))
        }else{ # length(cols_k) > 1
          Z_k = Znew2[idx_k, cols_k]
        }
        
        dat_list = list(N = length(idx_k), # Number individuals in group k
                        q = length(ranef_idx), # number random effects
                        eta_fef = as.numeric(X_k %*% matrix(coef[1:ncol(X)], ncol = 1)), # fixed effects componenet of linear predictor
                        y = y_k, # outcomes for group k
                        Z = Z_k) # portion of Z matrix corresonding to group k
        
        # initialize posterior random draws
        init_lst = list()
        for(cr in 1:num_cores){
          init_lst[[cr]] = list(alpha = rnorm(length(ranef_idx), 0, 1))
        }
        
        stan_fit = rstan::sampling(stan_file, data = dat_list, iter = nMC_chain + nMC_burnin,
                                   warmup = nMC_burnin, show_messages = F, refresh = 0,
                                   init = init_lst, chains = num_cores, cores = num_cores)
        
        # Extract all chains for use in M step
        stan_out = as.matrix(stan_fit)
        # Exclude lp__ column of output
        draws_mat = as.matrix(stan_out[,1:length(ranef_idx)])
        
        # Find last elements for each chain for initialization of next E step
        for(m in 1:num_cores){
          last_draws[m,cols_k] = draws_mat[nMC_chain*m,]
        }
        
        if(nrow(draws_mat) == nMC){
          u0[,cols_k] = draws_mat
        }else{ # nrow(draws_mat) > nMC due to ceiling function in 'iter' specification
          start_row = nrow(draws_mat) - nMC + 1
          rows_seq = start_row:nrow(draws_mat)
          u0[,cols_k] = draws_mat[rows_seq,]
        }
        
      } # End k for loop
      
    } # End if-else sampler
    
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
  
  # Start EM Algorithm
  
  for(i in 1:maxitEM){
    
    if(family == "binomial"){
      nTotal = rep(1, length(y[rep(1:nrow(X), each = nrow(u0))]))
    }else{
      nTotal = NULL
    }
    
    # M Step
    coef = M_stepB(y, X, Z, u0@address, nrow(u0), J, group, family=fam_fun, coef, offset=offset_fit,
                   maxit_CD=maxit_CD, conv_CD=conv_CD, init=(i == 1), group_X, covgroup,
                   penalty, lambda0, lambda1, gamma=gamma_penalty, alpha,
                   fit_type=fit_type)
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
    ll0 = Qfun(y, X, Z, u0@address, group, J, matrix(coef, ncol = 1), offset_fit, c(d,ncol(Z)/d), family, link_int)
    
    if(!is.finite(ll0) | any(is.na(coef))){
      problem = T
      ll0 = Inf
      print(coef)
    }
    
    if(problem == T){
      warning("Error in M step", immediate. = T)
      out = list(coef = coef, sigma = cov,  
                 lambda0 = lambda0, lambda1 = lambda1, 
                 covgroup = covgroup, J = J, ll = -Inf, BIC = Inf, BICh = Inf,
                 extra = list(okindex = okindex, Znew2 = Znew2),
                 gibbs_accept_rate = gibbs_accept_rate, proposal_SD = proposal_SD)
      
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
    
    if(trace >= 1) print(coef)
    
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
      print("cov:")
      print(cov)
    }
    
    if(nrow(cov) == 1){
      cov_record[i] = cov
    }
    
    # E Step 
    
    # Initial points for Metropolis within Gibbs E step algorithms
    uold = as.numeric(u0[nrow(u0),])
    
    # ranef_keep = which(diag(cov) > 0)
    # fixef_keep = XZ_map[which(coef[1:ncol(X)] != 0)]
    # ranef_idx = intersect(ranef_keep, fixef_keep)
    # if(!(any(ranef_idx == 1))) ranef_idx = c(1, ranef_idx)
    ranef_idx = which(diag(cov) > 0)
    
    if(sampler == "independence"){
      samplemc_out = sample_mc2_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, Z=Znew2, nMC=nMC, trace = trace, family = family, group = group, 
                                       d = d, okindex = okindex, uold = uold)
      u0 = attach.big.matrix(samplemc_out$u0)
      
    }else if(sampler == "random_walk"){
      
     # First adapt the proposal variance during an initial burnin period
      samplemc_burnin = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, Z=Znew2, nMC=nMC_burnin, trace = trace, family = family, group = group, 
                                            d = d, okindex = okindex, uold = uold,
                                            proposal_SD = proposal_SD, batch = batch, batch_length = batch_length, 
                                            offset = offset_increment, burnin_batchnum = burnin_batchnum)
      
      # Extract updated proposal_SD info from burnin period
      batch = samplemc_burnin$updated_batch
      proposal_SD = samplemc_burnin$proposal_SD
      
      u0 = attach.big.matrix(samplemc_burnin$u0)
      
      if(trace >= 1){
        print("Updated proposal_SD:")
        print(proposal_SD)
        cat("Updated batch number: ", batch, "\n")
      }
      
      # Run official E step with no adaptation
      uold = as.numeric(u0[nrow(u0),])
      samplemc_out = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], ranef_idx=ranef_idx, y=y, X=X, Z=Znew2, nMC=nMC, trace = trace, family = family, group = group, 
                                           d = d, okindex = okindex, uold = uold,
                                           proposal_SD = proposal_SD, batch = batch, batch_length = batch_length, 
                                           offset = offset_increment, burnin_batchnum = 0)
      
      gibbs_accept_rate = samplemc_out$gibbs_accept_rate
      u0 = attach.big.matrix(samplemc_out$u0)
      
    }else if(sampler == "stan"){
      
      # Determine number of cores to use
      if(num_cores < max_cores){
        if(nMC < 2000){
          num_cores = 1
        }else{
          num_cores = 2 + (nMC - 2000) %/% 1000
          if(num_cores > max_cores) num_cores = max_cores
        }
      }
      
      u0 = big.matrix(nrow = nMC, ncol = ncol(Z), init=0)
      
      if(family == "binomial" & link == "logit"){
        stan_file = stanmodels$binomial_logit_model
      }
      
      # Number draws to extract in each chain (after burn-in)
      nMC_chain = ceiling(nMC/num_cores)
      
      # Re-name last_draws matrix so it can be over-written in upcoming E step
      init_draws = last_draws
      # Record last elements for each chain for initialization of next E step
      last_draws = matrix(0, nrow = num_cores, ncol = ncol(Z))
      
      # For each group, sample from posterior distribution (sample the alpha values)
      for(k in 1:d){
        idx_k = which(group == k)
        cols_k = seq(from = k, to = ncol(Z), by = d)
        cols_k = cols_k[ranef_idx]
        y_k = y[idx_k]
        X_k = X[idx_k,]
        if(length(cols_k) == 1){
          Z_k = matrix(Znew2[idx_k, cols_k], nrow = length(idx_k), ncol = length(cols_k))
        }else{ # length(cols_k) > 1
          Z_k = Znew2[idx_k, cols_k]
        }
        
        dat_list = list(N = length(idx_k), # Number individuals in group k
                        q = length(ranef_idx), # number random effects
                        eta_fef = as.numeric(X_k %*% matrix(coef[1:ncol(X)], ncol = 1)), # fixed effects componenet of linear predictor
                        y = y_k, # outcomes for group k
                        Z = Z_k) # portion of Z matrix corresonding to group k
        
        # initialize posterior random draws from last E step
        # Note: If number of random effects to sample = 1, need to set the initial values to an array
        ## Stan interprets initial value as a scalar otherwise, when we need it to be a vector of length 1
        
        init_lst = list()
        if(nrow(init_draws) == num_cores){
          for(m in 1:num_cores){
            init_lst[[m]] = list(alpha = as.array(init_draws[m,cols_k]))
          }
        }else{ # nrow(init_draws) < num_cores
          for(m in 1:nrow(init_draws)){
            init_lst[[m]] = list(alpha = as.array(init_draws[m,cols_k]))
          }
          # If current number of cores exceeds previously used number of cores,
          # recycle initial points from first chain
          for(m in (nrow(init_draws)+1):num_cores){
            init_lst[[m]] = list(alpha = as.array(init_draws[1,cols_k]))
          }
        }
        
        stan_fit = rstan::sampling(stan_file, data = dat_list, iter = nMC_chain + nMC_burnin,
                                   warmup = nMC_burnin, show_messages = F, refresh = 0,
                                   init = init_lst, chains = num_cores, cores = num_cores)
        
        stan_out = as.matrix(stan_fit)
        # Exclude lp__ column of output
        draws_mat = as.matrix(stan_out[,1:length(ranef_idx)])
        
        # Find last elements for each chain for initialization of next E step
        for(m in 1:num_cores){
          last_draws[m,cols_k] = draws_mat[nMC_chain*m,]
        }
        
        if(nrow(draws_mat) == nMC){
          u0[,cols_k] = draws_mat
        }else{ # nrow(draws_mat) > nMC due to ceiling function in 'iter' specification
          start_row = nrow(draws_mat) - nMC + 1
          rows_seq = start_row:nrow(draws_mat)
          u0[,cols_k] = draws_mat[rows_seq,]
        }
        
      } # End k for loop
      
    } # End if-else sampler
    
    nMC2 = nrow(u0)
    
    if(trace >= 1) print(diag(cov))
    gc()
  }
  
  # Another E step for loglik calculation (number draws = M)
  
  # Note: Assume that at end of EM algorithm, which(diag(cov) > 0) corresponds with the
  # sampling restrictions specified within EM algorithm:
    # Restrict sampling such that variables are NOT sampled if:
    ## random effect for the variable was penalized to zero in a previous M step
    ## fixed effect for the variable was penalized to zero in a prevoius M step
    ## Make sure to always keep random intercept
  
  if(sampler == "independence"){
    samplemc_out = sample_mc2_BigMat(coef=coef[1:ncol(X)], ranef_idx=which(diag(cov) > 0), y=y, X=X, Z=Znew2, nMC=nMC, trace = trace, family = family, group = group, 
                                     d = d, okindex = okindex, uold = as.numeric(u0[nrow(u0),]))
    u0 = attach.big.matrix(samplemc_out$u0)
  }else if(sampler == "random_walk"){ 
    # No adaptation at convergence point
    samplemc_out = sample_mc_adapt_BigMat(coef=coef[1:ncol(X)], ranef_idx=which(diag(cov) > 0), y=y, X=X, Z=Znew2, nMC=nMC, trace = trace, family = family, group = group, 
                                          d = d, okindex = okindex, uold = as.numeric(u0[nrow(u0),]),
                                          proposal_SD = proposal_SD, batch = batch, batch_length = batch_length, 
                                          offset = offset_increment, burnin_batchnum = 0)
    u0 = attach.big.matrix(samplemc_out$u0)
  }else if(sampler == "stan"){
    
    # Determine number of cores to use
    if(num_cores < max_cores){
      if(nMC < 2000){
        num_cores = 1
      }else{
        num_cores = 2 + (nMC - 2000) %/% 1000
        if(num_cores > max_cores) num_cores = max_cores
      }
    }
    
    u0 = big.matrix(nrow = nMC, ncol = ncol(Z), init=0)
    ranef_idx=which(diag(cov) > 0)
    uold = as.numeric(u0[nrow(u0),])
    
    if(family == "binomial" & link == "logit"){
      stan_file = stanmodels$binomial_logit_model
    }
    
    # Number draws to extract in each chain (after burn-in)
    nMC_chain = ceiling(nMC/num_cores)
    
    # Re-name last_draws matrix so it can be over-written in upcoming E step
    init_draws = last_draws
    
    # For each group, sample from posterior distribution (sample the alpha values)
    for(k in 1:d){
      idx_k = which(group == k)
      cols_k = seq(from = k, to = ncol(Z), by = d)
      cols_k = cols_k[ranef_idx]
      y_k = y[idx_k]
      X_k = X[idx_k,]
      if(length(cols_k) == 1){
        Z_k = matrix(Znew2[idx_k, cols_k], nrow = length(idx_k), ncol = length(cols_k))
      }else{ # length(cols_k) > 1
        Z_k = Znew2[idx_k, cols_k]
      }
      
      dat_list = list(N = length(idx_k), # Number individuals in group k
                      q = length(ranef_idx), # number random effects
                      eta_fef = as.numeric(X_k %*% matrix(coef[1:ncol(X)], ncol = 1)), # fixed effects componenet of linear predictor
                      y = y_k, # outcomes for group k
                      Z = Z_k) # portion of Z matrix corresonding to group k
      
      # initialize posterior random draws from last E step
      init_lst = list()
      if(nrow(init_draws) == num_cores){
        for(m in 1:num_cores){
          init_lst[[m]] = list(alpha = as.array(init_draws[m,cols_k]))
        }
      }else{ # nrow(init_draws) < num_cores
        for(m in 1:nrow(init_draws)){
          init_lst[[m]] = list(alpha = as.array(init_draws[m,cols_k]))
        }
        # If current number of cores exceeds previously used number of cores,
        # recycle initial points from first chain
        for(m in (nrow(init_draws)+1):num_cores){
          init_lst[[m]] = list(alpha = as.array(init_draws[1,cols_k]))
        }
      }
      
      stan_fit = rstan::sampling(stan_file, data = dat_list, iter = nMC_chain + nMC_burnin,
                                 warmup = nMC_burnin, show_messages = F, refresh = 0,
                                 init = init_lst, chains = num_cores, cores = num_cores)
      
      stan_out = as.matrix(stan_fit)
      # Exclude lp__ column of output
      draws_mat = as.matrix(stan_out[,1:length(ranef_idx)])
      
      if(nrow(draws_mat) == nMC){
        u0[,cols_k] = draws_mat
      }else{ # nrow(draws_mat) > nMC due to ceiling function in 'iter' specification
        start_row = nrow(draws_mat) - nMC + 1
        rows_seq = start_row:nrow(draws_mat)
        u0[,cols_k] = draws_mat[rows_seq,]
      }
      
    } # End k for loop
    
  } # End if-else sampler
  
 
  if(sampler == "random_walk"){
    print("Updated proposal_SD:")
    print(proposal_SD)
    print("Updated batch:")
    print(batch)
  }
  
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
  
  # Calculated BICq
  ## Using posterior draws from the full model (ufull) and the coefficients from the
  ## current penalized model, calculate the Q function
  if(!is.null(ufull_describe)){
    ufull = attach.big.matrix(ufull_describe)
    q1 = Qfun(y, X, Z, ufull@address, group, J, matrix(coef, ncol = 1), offset_fit, c(d,ncol(Z)/d), family, link_int)
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
               extra = list(okindex = okindex, Znew2 = Znew2),
               gibbs_accept_rate = gibbs_accept_rate, proposal_SD = proposal_SD,
               EM_iter = i)
  }else if(sampler == "stan"){
    out = list(coef = coef, sigma = cov,  
               lambda0 = lambda0, lambda1 = lambda1, 
               covgroup = covgroup, J = J, ll = ll, BICh = BICh, BIC = BIC, BICq = BICq,
               extra = list(okindex = okindex, Znew2 = Znew2), EM_iter = i)
  }
  
  if(nrow(u0) >= nMC_report){
    # Take last nMC_report rows
    r_start = nrow(u0) - nMC_report + 1
    r_end = nrow(u0)
    u0_out = bigmemory::as.matrix(u0[c(r_start:r_end),])
  }else{
    if(sampler == "independence"){
      samplemc_out = sample_mc2(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC_report, trace = trace, family = family, group = group, 
                                d = d, okindex = okindex, matrix(u0[nrow(u0),], nrow = 1))
      u0_out = samplemc_out$u0
      
    }else if(sampler == "random_walk"){ 
      # No adaptation at convergence point
      samplemc_out = sample_mc_adapt(coef=coef[1:ncol(X)], cov=cov, y=y, X=X, Z=Znew2, nMC=nMC_report, trace = trace, family = family, group = group, 
                                     d = d, okindex = okindex, uold = matrix(u0[nrow(u0),], nrow = 1),
                                     proposal_SD = proposal_SD, batch = batch, batch_length = batch_length, 
                                     offset = offset_increment, burnin_batchnum = 0)
      u0_out = samplemc_out$u0
      
    }else if(sampler == "stan"){
      
      # Note: only use one chain here so that trace plots are an accurate representation 
      # of a single chain
      
      u0 = matrix(0, nrow = nMC_report, ncol = ncol(Z))
      ranef_idx=which(diag(cov) > 0)
      uold = as.numeric(u0[nrow(u0),])
      
      if(family == "binomial" & link == "logit"){
        stan_file = stanmodels$binomial_logit_model
      }
      
      # For each group, sample from posterior distribution (sample the alpha values)
      for(k in 1:d){
        idx_k = which(group == k)
        cols_k = seq(from = k, to = ncol(Z), by = d)
        cols_k = cols_k[ranef_idx]
        y_k = y[idx_k]
        X_k = X[idx_k,]
        if(length(cols_k) == 1){
          Z_k = matrix(Znew2[idx_k, cols_k], nrow = length(idx_k), ncol = length(cols_k))
        }else{ # length(cols_k) > 1
          Z_k = Znew2[idx_k, cols_k]
        }
        
        dat_list = list(N = length(idx_k), # Number individuals in group k
                        q = length(ranef_idx), # number random effects
                        eta_fef = as.numeric(X_k %*% matrix(coef[1:ncol(X)], ncol = 1)), # fixed effects componenet of linear predictor
                        y = y_k, # outcomes for group k
                        Z = Z_k) # portion of Z matrix corresonding to group k
        
        # initialize posterior random draws
        init_lst = list()
        init_lst[[1]] = list(alpha = as.array(uold[cols_k]))
        
        stan_fit = rstan::sampling(stan_file, data = dat_list, chains = 1, init = init_lst,
                                   iter = nMC_report + nMC_burnin, warmup = nMC_burnin, 
                                   show_messages = F, refresh = 0)
        list_of_draws = extract(stan_fit)
        u0[,cols_k] = list_of_draws$alpha
        
        u0_out = u0
        
      } # End k for loop
      
    } # End if-else sampler
    
  } # End if-else nrow(u0)
  
  out$u = u0_out
  
  
 # out$coef_record_all = coef_record_all
  if(!is.null(cov_record)){
    out$cov_record = cov_record
  }
  
  return(out)
}


