# Functions to define object's S3 methods

# Object name: pglmmObj
# Properties: 
## Reference Class object

# fixef = function(object){
#   UseMethod("fixef")
# }

#' @export
fixef.pglmmObj = function(object){
  object$fixef
}

# ranef = function(object){
#   UseMethod("ranef")
# }

#' @export
ranef.pglmmObj = function(object){
  object$ranef
}

coefGlmmPen = function(object){
  # Combined coefficients
  fixef = object$fixef
  ranef = object$ranef
  coefficients = data.frame()
  for(j in length(fixef)){
    var = names(fixef[j])
    if(var %in% colnames(ranef)){
      coefficients$var = coefficients$var + fixef[j]
    }else{
      coefficients$var = fixef[j]
    }
  }
  rownames(coefficients) = levels(object$group)
  
  return(coefficients)
}

# coef = function(object){
#   UseMethod("coef")
# }

#' @importFrom stats coef
#' @export
coef.pglmmObj = coefGlmmPen

# family = function(object){
#   UseMethod{"family"}
# }

#' @importFrom stats family
#' @export
family.pglmmObj = function(object){
  object$family
}

# nobs = function(object){
#   UseMethod("nobs")
# }

##' @export
nobs.pglmmObj = function(object){
  nrows(object$X)
}

# ngrps = function(){
#   UseMethod("ngrps")
# }
# 
# ngprs = function(object){
#   
# }

# fitted = function(object){
#   UseMethod("fitted")
# }

#' @importFrom stats fitted
#' @export
fitted.pglmmObj = function(object){
  X = object$X
  Z = object$Z
  beta = object$fixef
  U = object$gibbs_mcmc
  if(ncol(U) == 1){
    gamma = mean(U)
  }else{
    gamma = colMeans(object$gibbs_mcmc)
  }
  
  mu = X %*% beta + Z %*% gamma
  return(mu)
}

# Functions to add:
## offset() - will not have $offset option
## weights() - will not have $weights option