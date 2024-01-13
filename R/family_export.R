
#' @importFrom stringr str_detect
#' @importFrom MASS negative.binomial
#' @importFrom stats poisson
family_export = function(family){
  
  data_type = NULL
  
  if(is.character(family)){
    if(family == "coxph"){
      data_type = "survival"
    }
  }
  
  if(is.null(data_type)){
    data_type = "glmm"
  }
  
  if(is.character(family)){
    if(family != "coxph"){
      family = get(family, mode = "function", envir = parent.frame())
    }else if(family == "coxph"){
      family = poisson(link="log")
    }
  }
  if(is.function(family)){
    family = family()
  }
  if(inherits(family, "family")){ # class(family) == "family"
    fam_fun = family
    link = family$link
    family = family$family
  }
  
  ## Possible code for future inclusion of negative binomial family
  # if(is.character(family)){
  #   if(family == "negbin"){
  #     fam_fun = negative.binomial(theta = 1) # Give negbin dist arbitrary theta for now
  #     link = "log"
  #     family = "negbin"
  #   }else{
  #     family = get(family, mode = "function", envir = parent.frame())
  #   }
  # }
  # if(is.function(family)){
  #   family = family()
  # }
  # if(class(family) == "family"){
  #   fam_fun = family
  #   link = family$link
  #   family = family$family
  #   if(str_detect(family, "Negative Binomial")){
  #     family = "negbin"
  #   }
  # }
  
  ## Currently, only allow the following family-link combinations
  if(!(family %in% c("binomial","poisson","gaussian"))){
    stop("Invalid family. Available families: 'binomial', 'poisson', 'gaussian', or 'coxph'")
  }
  if(!(link %in% c("logit","log","identity"))){
    stop("Invalid link. Available link functions: 'logit', 'log', 'identity'")
  }
  if((family == "binomial" & link != "logit") | (family == "poisson" & link != "log") | (family == "gaussian" & link != "identity")){
    stop("family and link combination not available")
  }
  # if(!(link %in% c("logit","probit","cloglog","log","sqrt","identity","inverse"))){
  #   stop("Invalid link. Available link functions: 'logit', 'probit', 'cloglog', 'log', 'sqrt', 'identity', 'inverse'")
  # }
  
  # Re-code link as integer
  ## All link_int will have two digits
  ## First digit corresponds to family that link generally associated with
  ## 1 = binomial, 2 = poisson or negative binomial, 3 = gaussian
  ## Second digit: 0 = canonical link, other = arbitrary enumeration of common non-canonical links
  if(link == "logit"){
    link_int = 10
  }else if(link == "log"){
    link_int = 20
  }else if(link == "identity"){
    link_int = 30
  }
  ## For future: add additional family-link combinations
  # if(link == "logit"){
  #   link_int = 10
  # }else if(link == "probit"){
  #   link_int = 11
  # }else if(link == "cloglog"){
  #   link_int = 12
  # }else if(link == "log"){
  #   link_int = 20
  # }else if(link == "sqrt"){
  #   link_int = 21
  # }else if(link == "identity"){
  #   link_int = 30
  # }else if(link == "inverse"){
  #   link_int = 31
  # }
  
  return(list(family_fun = fam_fun, family = family, link = link, link_int = link_int, data_type = data_type))
  
}