# Functions to define object's S3 methods

# Object name: pglmmObj
# Properties: 
## Reference Class object

#' @export
fixef.pglmmObj = function(object){
  object$fixef
}

#' @export
ranef.pglmmObj = function(object){
  object$ranef
}

# coefGlmmPen = function(object){
#   # Combined coefficients
#   fixef = object$fixef
#   ranef = object$ranef
#   coefficients = data.frame()
#   for(j in length(fixef)){
#     var = names(fixef[j])
#     if(var %in% colnames(ranef)){
#       coefficients$var = coefficients$var + fixef[j]
#     }else{
#       coefficients$var = fixef[j]
#     }
#   }
#   rownames(coefficients) = levels(object$group)
#   
#   return(coefficients)
# }

coefMer <- function(object, ...){
  # Direct copy from lme4 source code
  if (length(list(...)))
    warning('arguments named "', paste(names(list(...)), collapse = ", "),
            '" ignored')
  fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
  ref <- ranef(object)
  ## check for variables in RE but missing from FE, fill in zeros in FE accordingly
  refnames <- unlist(lapply(ref,colnames))
  nmiss <- length(missnames <- setdiff(refnames,names(fef)))
  if (nmiss > 0) {
    fillvars <- setNames(data.frame(rbind(rep(0,nmiss))),missnames)
    fef <- cbind(fillvars,fef)
  }
  val <- lapply(ref, function(x)
    fef[rep.int(1L, nrow(x)),,drop = FALSE])
  for (i in seq(a = val)) {
    refi <- ref[[i]]
    row.names(val[[i]]) <- row.names(refi)
    nmsi <- colnames(refi)
    if (!all(nmsi %in% names(fef)))
      stop("unable to align random and fixed effects")
    for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[,nm]
  }
  class(val) <- "coef.mer"
  val
} 

#' @importFrom stats coef
#' @export
coef.pglmmObj = coefMer

#' @importFrom stats family
#' @export
family.pglmmObj = function(object){
  object$family
}

#' @export
nobs.pglmmObj = function(object){
  nrow(object$X)
}

# ngprs = function(object){
#   
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

#' @importFrom stats formula
#' @export
formula.pglmmObj = function(object){
  object$formula
}

# Functions to add:
## offset() - will not have $offset option
## weights() - will not have $weights option
### Ideas for offset() and weights() - extract from model.frame ?