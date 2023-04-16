# Utility functions for the Cox Proportional Hazards family procedure

###############################################################################################
# Adjustments to data for use in the coxph family
###############################################################################################
# Note: must input non-sparse Z for this function to work as-is
## If want to apply to non-sparse Z, need to make adjustments
# @export
coxph_data = function(y, X, Z, group){ # , na.action = na.omit
  
  if(!inherits(y,"Surv")){
    stop("response in formula must be of class 'Surv', use survival::Surv()")
  }
  # Extract observed times
  y_times = y[,1]
  # Extract status (event = 1, censoring = 0)
  y_status = y[,2]
  
  # Checks to y_times: Cannot have negative values
  if(any(y_times <= 0)){
    stop("survival times must be positive, check or remove times that are negative or equal to 0")
  }
  
  # Checks to y_status: 
  ## Make sure values are 0 or 1 
  ## See survival::Surv() for possible options for events, ideas for checks
  if(!all(unique(y_status) %in% c(0,1))){
    stop("event variable must be 0 or 1 (censoring vs event)")
  }
  
  # Order y, X, Z, and group by observed times
  ## Note: In glFormula_edit, removed NA observations
  order_idx = order(y_times)
  y_times = y_times[order_idx]
  y_status = y_status[order_idx]
  X = X[order_idx,,drop=FALSE]
  Z = Z[order_idx,,drop=FALSE]
  group = group[order_idx]
  
  # If intercept term in X, remove
  if(ncol(X) >= 2){
    if(all(X[,1] == 1)){
      X = X[,-1,drop=FALSE]
    }
  }else if(ncol(X) == 1){
    if(all(X  == 1)){
      stop("Fixed-effect intercept not calculated for the 'coxph' model fit")
    }
  }
  
  
  # Calculate unique times (event or censoring)
  times_unique = unique(y_times)
  if(length(times_unique) == length(y_times)){ # No ties
    nevent = y_status
    nsubject = rep(1, length(y_times))
    t_loc = 1:length(y_times)
    t_group = 1:length(y_times)
  }else if(length(times_unique) < length(y_times)){ # At least one tie
    t_loc = numeric(length(times_unique))
    nsubject = numeric(length(times_unique))
    nevent = numeric(length(times_unique))
    t_group = numeric(length(y_times))
    for(t in 1:length(times_unique)){
      # Determine which subjects have an observed time at unique time t
      idx = which(y_times == times_unique[t])
      # Determine index locations of where the unique times occur in the original dataset
      ## (First observed location of this unique time)
      t_loc[t] = idx[1]
      # Calculate number of subjects observed at each unique time
      nsubject[t] = length(idx)
      # Calculate number of events observed at each unique time
      nevent[t] = sum(y_status[idx])
      # Grouping within y_times
      t_group[idx] = t
    }
  }
  
  
  return(list(y_status = y_status, y_times = y_times, nsubject = nsubject,
              nevent = nevent, t_loc = t_loc, t_group = t_group,
              X = X, Z = Z, group = group))
  
}

###############################################################################################
# Initializations of random effect covariates 
# Used within fit_dat_coxph
###############################################################################################

# Initialization of the B matrix and b vector within the glmmPen_FA formulation
## cov = B %*% t(B)
ranef_init_B = function(dat, fam_fun, beta0, b_old, q, r, ranef_keep, B_init_type, 
                        var_start, var_restrictions, coef_names){
  
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
          if(var_start == "recommend"){
            var_start = 0.10
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
      if(var_start == "recommend"){
        var_start = 0.10
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
  # coef = c(beta0, b0)
  
  return(list(b0 = b0, cov = cov))
  
}

# Initialization of the Gamma matrix and gamma vector within the glmmPen formulation
## cov = Gamma %*% t(Gamma)
ranef_init_Gamma = function(dat, gamma_old, covar, q, ranef_keep, J, 
                            var_start, var_restrictions, coef_names, trace, progress){
  
  if(var_start == "recommend"){
    var_start = var_init(data = dat, fam_fun = list(family = "coxph"))
  }
  
  if(!is.null(gamma_old)) {
    
    if(progress == TRUE) message("using coef from past model to intialize")
    gamma = matrix(J%*%matrix(gamma_old, ncol = 1), ncol = q)
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
      gamma_vec = c(gamma)[which(lower.tri(matrix(0,nrow=nrow(cov),ncol=ncol(cov)),diag=TRUE))]
    }else if(covar == "independent"){
      gamma_vec = diag(gamma)
    }
    
  }else{ # Will not initialize with past model results
    
    # apply variance restrictions if necessary based on initialized fixed effects
    if(var_restrictions == "fixef"){
      fixed_keep = coef_names$fixed[c(1,which(coef[-1] != 0)+1)]
      restrict_idx = which(!(coef_names$random %in% fixed_keep))
      for(j in 1:q){
        if(j %in% restrict_idx){
          ranef_keep[j] = 0
        }
      }
    }
    # Initialize covariance matrix (cov)
    if(q > 1){
      vars = rep(var_start, q) * ranef_keep
      cov = var = diag(vars)
      gamma = diag(sqrt(vars)) 
    }else{
      vars = var_start
      cov = var = matrix(vars, ncol = 1)
      gamma = matrix(sqrt(vars), ncol = 1)
    }
    
    if(trace >= 1){
      cat("initialized covariance matrix diagonal: \n", diag(cov), "\n")
    }
    
    if(covar == "unstructured"){
      gamma_vec = c(gamma)[which(lower.tri(matrix(0,nrow=nrow(cov),ncol=ncol(cov)),diag=T))]
    }else if(covar == "independent"){
      gamma_vec = diag(gamma)
    }
    
    
  } # End if-else !is.null(coef_old)
  
  return(list(gamma_vec = gamma_vec, cov = cov))
  
}

# setupLambdaCox copied and slightly modified from ncvreg code
# maxprod also copied
# y: y_times vector (ordered)
# Delta: y_status vector (ordered)
# @importFrom survival coxph
# setupLambdaCox_copy <- function(X, y, Delta, alpha, lambda.min, nlambda, penalty.factor) {
#   n <- nrow(X)
#   p <- ncol(X)
#   
#   # Fit to unpenalized covariates
#   ind <- which(penalty.factor!=0)
#   if (length(ind)!=p) {
#     nullFit <- survival::coxph(survival::Surv(y, Delta) ~ X[, -ind, drop=FALSE])
#     eta <- nullFit$linear.predictors
#     rsk <- rev(cumsum(rev(exp(eta))))
#     s <- Delta - exp(eta)*cumsum(Delta/rsk)
#   } else {
#     w <- 1/(n-(1:n)+1)
#     s <- Delta - cumsum(Delta*w)
#   }
#   
#   # Determine lambda.max
#   zmax <- .Call("maxprod", X, s, ind, penalty.factor) / n
#   lambda.max <- zmax/alpha
#   
#   if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)), 0)
#   else lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), len=nlambda))
#   lambda
# }

cut_points_calc = function(coxph_options, y, y_times, group, d){
  # Cox Proportional Hazards family: Calculate cut-points to use for time intervals
  ## Divide the timeline into J = cut_num (or less) intervals such that there are
  ##  at least 3 events for each group within each cut-point
  ## Note: must have at least one event in each interval (preferably > 2) to be identifiable
  # Starting suggested number of cut-points; want at least 5 total
  cut_num = coxph_options$cut_num 
  # Total number of events in dataset
  event_total = sum(y)
  # Goal: if random intercept variance small, would expect we could equally split events
  #   into each time interval, and would expect round(event_total/cut_num) events per interval
  goal_events = floor(event_total / cut_num)
  # Create table of events by group for each time point
  ## Note: Need to figure out way to adjust for ties in y_times
  event_mat = matrix(0, nrow = length(y_times), ncol = d)
  for(i in 1:length(y_times)){
    if(y[i] == 1){
      k = group[i]
      event_mat[i,k] = 1
    }
  }
  
  cut_points = rep(NA, times=cut_num)
  event_mat_tmp = event_mat
  event_cumsum = matrix(0, nrow = nrow(event_mat), ncol = d)
  for(j in 1:cut_num){
    for(k in 1:d){
      event_cumsum[,k] = cumsum(event_mat_tmp[,k])
    }
    tot_cumsum = rowSums(event_cumsum)
    if(max(tot_cumsum) < goal_events){
      tmp_idx = nrow(event_mat)
    }else{
      tmp_idx = min(which(tot_cumsum >= goal_events))
    }
    grp_issue = TRUE
    message("initial tmp_idx: ", tmp_idx)
    while(grp_issue & (tmp_idx < nrow(event_mat))){
      print(event_cumsum[tmp_idx,])
      if(any(event_cumsum[tmp_idx,] < 4)){
        tmp_idx = tmp_idx + 1
      }else if(all(event_cumsum[tmp_idx,] >= 4)){
        grp_issue = FALSE
      }
      print("tmp_idx update")
    }
    message("ending tmp_idx: ", tmp_idx)
    if(grp_issue & (tmp_idx == nrow(event_mat))){
      cut_points[j-1] = max(y_times) + 1
      break
    }
    if(tmp_idx < nrow(event_mat)){
      cut_points[j] = mean(y_times[tmp_idx], y_times[tmp_idx + 1])
    }else if(tmp_idx == nrow(event_mat)){
      cut_points[j] = max(y_times) + 1
    }
    
    # Reset
    print(tmp_idx)
    event_mat_tmp[1:tmp_idx,] = 0
    
  }
  cut_points = cut_points[which(!is.na(cut_points))]
  print(cut_points)
  if(max(cut_points) < max(y_times)){
    cut_points[length(cut_points)] = max(y_times)
  }
  print(cut_points)
  message("Number intervals used: ", length(cut_points))
  for(j in 1:length(cut_points)){
    if(j == 1){
      cut_min = 0
    }else{
      cut_min = cut_points[j-1]
    }
    int_idx = which((y_times >= cut_min) & (y_times < cut_points[j]))
    y_tmp = y[int_idx]
    grp_tmp = group[int_idx]
    message("Total events in interval ", j, ": ", sum(y_tmp))
    message("Total events by group for interval:")
    print(table(y_tmp, grp_tmp))
  }
  
  # Adjust cut-points so that first cut-point is 0 and final cut-point is the
  #   current 2nd-to-last cut-point. (In survSplit function, do not want the
  #   final cut-point to be the max value of y)
  cut_points = c(0, cut_points[-length(cut_points)])
  
  return(cut_points)
}


cut_points_calc_simple = function(coxph_options, y, y_times){
  
  # Initial cut-point estimation - equal observations per time period, not taking into
  #   account the number of observations for each group within each time period
  ##############################################################################################
  # Cox Proportional Hazards family: Calculate cut-points to use for time intervals
  ## Divide the timeline into J = cut_num intervals such that there are an equal
  ##  (or approximately equal) number of events in each interval
  ## Note: must have at least one event in each interval (preferably > 2) to be identifiable
  cut_num = coxph_options$cut_num
  event_total = sum(y)
  event_idx = which(y == 1)
  # Determine number of events per time interval, event_j
  if((event_total %% cut_num) == 0){ # event_total is a factor of cut_num
    event_cuts = rep(event_total / cut_num, times = cut_num)
  }else{
    tmp = event_total %/% cut_num
    event_cuts = rep(tmp, times = cut_num)
    for(j in 1:(event_total - tmp*cut_num)){
      event_cuts[j] = event_cuts[j] + 1
    }
  }

  # warning if only 1 event for an interval, stop if 0 events for an interval
  if(any(event_cuts == 1)){
    warning("at least one time interval for the piecewise exponential hazard model has only 1 event, ",
            "please see the coxphControl() documentation for details and tips on how to fix the issue",
            immediate. = TRUE)
  }else if(any(event_cuts == 0)){
    stop("at least one time interval for the piecewise exponential hazard model has 0 events, ",
         "please see the coxphControl() documentation for details and tips on how to fix the issue")
  }

  cut_pts_idx = numeric(cut_num)
  for(j in 1:cut_num){
    cut_pts_idx[j] = event_idx[sum(event_cuts[1:j])]
  }

  cut_points = numeric(cut_num)
  for(j in 1:(cut_num-1)){
    cut_points[j] = mean(y_times[cut_pts_idx[j]], y_times[cut_pts_idx[j]+1])
  }
  cut_points[cut_num] = max(y_times) + 1
  
  # re-do cut-points so that first value = 0 and the max value is ignored
  cut_points = c(0, cut_points[-cut_num])
  
  return(cut_points)
  
}