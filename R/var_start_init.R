
# data: list object containing response y and group vector 'group'
# fam_fun: family function output from 'family_export' function
#' @importFrom lme4 lmer glmer VarCorr
var_init = function(data, fam_fun){
  y = data$y
  grp = data$group
  if(fam_fun$family == "gaussian"){
    int_only = lmer(formula = y ~ 1 + (1 | grp))
  }else if(fam_fun$family %in% c("binomial","poisson")){
    int_only = glmer(formula = y ~ 1 + (1 | grp), family = fam_fun)
  }
  
  var_start = VarCorr(int_only)$grp[[1]] * 2.0
  
  # Alternative starting variance specification: 
  # if(fam_fun$family == "gaussian"){
  #   var_start = VarCorr(int_only)$grp[[1]] * 2.0
  # }else{
  #   var_start = VarCorr(int_only)$grp[[1]]
  # }
  
  # Protect against case when variance estimate above is 0 or very small
  var_start = max(0.5, var_start) 
  message(sprintf("recommended starting variance: %f", var_start))
  
  return(var_start)
}

#' @importFrom lme4 lmer glmer VarCorr
var_init_survival = function(data){
  y = data$y
  grp = data$group
  cut_points = data$cut_points
  cut_num = length(cut_points)
  # Note: X1 inlcudes the reference coding for the time point indicators
  X1 = data$X[,c(1:length(cut_points)), drop=FALSE]
  # Since X1 includes an intercept term, remove automatic intercept from coding
  int_only = glmer(formula = y ~ -1 + X1 + (1 | grp), family = poisson())
  
  var_start = VarCorr(int_only)$grp[[1]]
  
  # Protect against case when variance estimate above is 0 or very small
  var_start = max(0.5, var_start) 
  message(sprintf("recommended starting variance: %f", var_start))
  
  return(var_start)
}

