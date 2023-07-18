# Utility functions for the Cox Proportional Hazards family procedure

###############################################################################################
# Adjustments to data for use in the piecewise exponential model fitting of survival data
###############################################################################################
# Note: must input non-sparse Z for this function to work as-is
## If want to apply to non-sparse Z, need to make adjustments
#' @title Convert Input Survival Data Into Long-Form Data Needed for Fitting a Piecewise
#' Exponential Model
#' 
#' @description Converts the input survival data with one row or element corresponding to a 
#' single observation or subject into a long-form dataset where one observation or subject
#' contributes \code{j} rows, where \code{j} is the number of time intervals that 
#' a subject survived at least part-way through.
#' 
#' @inheritParams phmmPen
#' @param y response, which must be a \code{Surv} object (see
#' \code{\link[survival]{Surv}} from the \code{survival} package)
#' @param X matrix of fixed effects covariates
#' @param Z matrix of random effects covariates
#' @param group vector specifying the group assignment for each subject
#' @param offset_fit vector specifying the offset. 
#' This can be used to specify an \emph{a priori} known component to be included in the 
#' linear predictor during fitting. Default set to \code{NULL} (no offset). If the data 
#' argument is not \code{NULL}, this should be a numeric vector of length equal to the 
#' number of cases (the length of the response vector). 
#'
#' @importFrom stringr str_c
#' @export
survival_data = function(y, X, Z, group, offset_fit = NULL, survival_options){ # , na.action = na.omit
  
  interval_type = survival_options$interval_type
  cut_points = survival_options$cut_points
  cut_num = survival_options$cut_num
  time_scale = survival_options$time_scale
  
  if(!inherits(y,"Surv")){
    stop("response in formula must be of class 'Surv', use survival::Surv()")
  }
  # Extract observed times
  y_times = y[,1]
  # Extract status (event = 1, censoring = 0)
  y_status = y[,2]
  
  # Multiply y_times by specified time_scale (default = 1)
  y_times = y_times * time_scale
  
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
  
  # Check of offset_fit
  if(is.null(offset_fit)){
    offset_fit = rep(0.0, length(y))
  }else{
    if((!is.numeric(offset_fit)) | (length(offset_fit) != length(y))){
      stop("offset must be a numeric vector of the same length as y")
    }
  }
  
  # Order y, X, Z, and group by observed times
  ## Note: In glFormula_edit, removed NA observations
  order_idx = order(y_times)
  y_times = y_times[order_idx]
  y_status = y_status[order_idx]
  X = X[order_idx,,drop=FALSE]
  Z = Z[order_idx,,drop=FALSE]
  group = group[order_idx]
  offset_fit = offset_fit[order_idx]
  
  # If intercept term in X, remove
  if(ncol(X) >= 2){
    if(all(X[,1] == 1)){
      X = X[,-1,drop=FALSE]
    }
  }else if(ncol(X) == 1){
    if(all(X  == 1)){
      stop("Intercept-only model not allowed for the piecewise exponential (survival) model fit")
    }
  }
  
  # save this dataset - used in initialization of random effect variance values
  data_ordered = list(y_status = y_status, y_times = y_times, X = X, Z = Z, group = group, offset_fit = offset_fit)
  
  ############################################################################################
  # Calculate time interval cut-points to use in E-step
  # E-step: Piecewise Exponential approximation of the Cox Proportional Hazards model
  # For now, assume no ties in times (will need to adjust for ties eventually)
  ############################################################################################
  
  # Create long-form dataset
  if(interval_type == "equal"){
    cut_points = cut_points_calc_simple(cut_num = cut_num, 
                                        y = y_status, y_times = y_times)
  }else if(interval_type == "group"){
    cut_points = cut_points_calc(cut_num = cut_num, 
                                 y = y_status, y_times = y_times, 
                                 group = group, d = nlevels(group))
  }else if(interval_type == "manual"){
    if(is.null(cut_points)){
      stop("cut_points must be specified if interval_type = 'manual' ")
    }
  }
  
  
  ############################################################################################
  # Create long-form dataset
  # Note: PE stands for Piecewise-Exponential, referring to our use of the 
  #   Piecewise-Exponential approximation to the Cox Proportional Hazards model
  ############################################################################################
  
  df_init = data.frame(id = 1:length(y_status), y_status=y_status, y_times=y_times)
  df_PE = survival::survSplit(formula = Surv(y_times, y_status) ~ .,
                              data = df_init, cut = cut_points)
  # Re-define data variables in the long format
  IDs = df_PE$id
  y_status = df_PE$y_status
  y_times = df_PE$y_times
  offset_fit = offset_fit[IDs]
  group = group[IDs]
  interval = factor(df_PE$tstart)
  interval_mat = model.matrix(y_status ~ interval) # Reference coding for time interval indicator
  colnames(interval_mat) = c("(Intercept)",str_c("Time_Interval",2:ncol(interval_mat)))
  interval_length = y_times - df_PE$tstart
  offset_interval = log(interval_length)
  X = cbind(interval_mat,X[IDs,])
  Z = Z[IDs,]
  offset_total = offset_fit + offset_interval
  
  return(list(IDs = IDs, y_status = y_status, y_times = y_times,
              X = X, Z = Z, group = group, offset_total = offset_total,
              cut_points = cut_points,
              data_ordered = data_ordered))
  
}


cut_points_calc = function(cut_num, y, y_times, group, d){
  # Cox Proportional Hazards family: Calculate cut-points to use for time intervals
  ## Divide the timeline into J = cut_num (or less) intervals such that there are
  ##  at least 3 events for each group within each cut-point
  ## Note: must have at least one event in each interval (preferably > 2) to be identifiable
  # Starting suggested number of cut-points; want at least 5 total
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
    # message("initial tmp_idx: ", tmp_idx)
    while(grp_issue & (tmp_idx < nrow(event_mat))){
      # print(event_cumsum[tmp_idx,])
      if(any(event_cumsum[tmp_idx,] < 4)){
        tmp_idx = tmp_idx + 1
      }else if(all(event_cumsum[tmp_idx,] >= 4)){
        grp_issue = FALSE
      }
      # print("tmp_idx update")
    }
    # message("ending tmp_idx: ", tmp_idx)
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
    # print(tmp_idx)
    event_mat_tmp[1:tmp_idx,] = 0
    
  }
  cut_points = cut_points[which(!is.na(cut_points))]
  if(max(cut_points) < max(y_times)){
    cut_points[length(cut_points)] = max(y_times)
  }
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


cut_points_calc_simple = function(cut_num, y, y_times){
  
  # Initial cut-point estimation - equal observations per time period, not taking into
  #   account the number of observations for each group within each time period
  ##############################################################################################
  # Cox Proportional Hazards family: Calculate cut-points to use for time intervals
  ## Divide the timeline into J = cut_num intervals such that there are an equal
  ##  (or approximately equal) number of events in each interval
  ## Note: must have at least one event in each interval (preferably > 2) to be identifiable
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
            "please see the survivalControl() documentation for details and tips on how to fix the issue",
            immediate. = TRUE)
  }else if(any(event_cuts == 0)){
    stop("at least one time interval for the piecewise exponential hazard model has 0 events, ",
         "please see the survivalControl() documentation for details and tips on how to fix the issue")
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