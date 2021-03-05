
# Perform checks of input arguments

# check penalization input
checkPenalty = function(penalty, gamma_penalty, alpha){
  
  if(length(penalty) > 1){
    penalty = penalty[1]
  }
  if(!(penalty %in% c("lasso","MCP","SCAD"))){
    stop("penalty ", penalty, " not available, must choose 'lasso', 'MCP', or 'SCAD' \n")
  }
  
  if(!is.double(gamma_penalty)) {
    tmp <- try(gamma_penalty <- as.double(gamma_penalty), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("gamma_penalty must be numeric or able to be coerced to numeric", call.=FALSE)
  }
  
  if(penalty == "MCP" & gamma_penalty <= 1){
    stop("gamma_penalty must be > 1 when using MCP penalty")
  }else if(penalty == "SCAD" & gamma_penalty <= 2){
    stop("gamma_penalty must be > 2 when using SCAD penalty")
  }else if(!is.double(gamma_penalty)){
    tmp <- try(gamma_penalty <- as.double(gamma_penalty), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("gamma_penalty must be numeric or able to be coerced to numeric", call.=FALSE)
  }
  
  if(!is.double(alpha)) {
    tmp <- try(alpha <- as.double(alpha), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("alpha must be numeric or able to be coerced to numeric", call.=FALSE)
  }else if(alpha == 0.0){
    stop("alpha cannot equal 0. Pick a small value > 0 instead (e.g. 0.001) \n");
  }
  
  return(list(penalty = penalty, gamma_penalty = gamma_penalty, alpha = alpha))
}

# check covariance matrix structure 
#' @importFrom stringr str_c
checkCovar = function(covar, acceptable = c("unstructured","independent")){
  
  if(length(covar) > 1){
    covar = covar[1]
  }
  if(!(covar %in% acceptable)){
    stop("covariance structure 'covar' must be ", str_c(acceptable, collapse = " or "))
  }
  return(covar)
}

# check that sampler is specified correctly
#' @importFrom stringr str_c
checkSampler = function(sampler, acceptable = c("stan","random_walk","independence")){
  if(length(sampler) > 1){
    sampler = sampler[1]
  }
  if(!(sampler %in% acceptable)){
    stop("sampler must be specified as ", str_c(acceptable, collapse = " or "))
  }
  return(sampler)
}

# check that BICq_posterior can be saved appropriately
#' @importFrom stringr str_detect
checkBICqPost = function(BICq_posterior){
  if(!is.null(BICq_posterior)){
    
    file_name = basename(BICq_posterior)
    path_name = dirname(BICq_posterior)
    # Check file is specified as a .txt file
    if(!str_detect(file_name,".txt$")){
      stop("BICq_posterior file name must end in a .txt extension")
    }
    # Check that path to file exists
    if(!dir.exists(path_name)){
      stop("The path ", path_name, " specified for the 'BICq_posterior' does not exist")
    }
    
  }
  
}