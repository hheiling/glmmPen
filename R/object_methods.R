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

coefGlmmPen = function(object){
  # Find combined coefficients
  ## Borrowed code elements from coef.merMod in lme4
  
  fixef = data.frame(rbind(object$fixef), check.names = F)
  ranef = object$ranef
  group = object$group
  
  # check for variables in RE but missing from FE, fill in zeros in FE accordingly
  refnames = unlist(lapply(ranef,colnames))
  nmiss = length(missnames <- setdiff(refnames,names(fixef)))
  if (nmiss > 0) {
    fillvars = setNames(data.frame(rbind(rep(0,nmiss))),missnames)
    fixef = cbind(fillvars,fixef)
  }
  
  output = lapply(ranef, function(j) fixef[rep.int(1L, nrow(j)), , drop = FALSE])
  
  for (i in seq(a = output)){
    refi <- ranef[[i]]
    row.names(output[[i]]) <- row.names(refi)
    nmsi <- colnames(refi)
    if (!all(nmsi %in% names(fixef))){
      stop("unable to align random and fixed effects")
    }
    for (nm in nmsi) {
      output[[i]][[nm]] <- output[[i]][[nm]] + refi[,nm]
    }
  }
  
  return(output)
  
}

#' @importFrom stats coef
#' @export
coef.pglmmObj = coefGlmmPen

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

reOnly <- function(f, response=FALSE) {
  # Copied from lme4 source code (predict.R)
  reformulate(paste0("(", vapply(findbars(f), safeDeparse, ""), ")"),
              response = if(response && length(f)==3L) f[[2]])
}

#' @importFrom stats formula
#' @export
formula.pglmmObj = function(object, fixed.only = F, random.only = F){
  form = object$formula
  frame = object$frame
  if(fixed.only && random.only){
    stop("cannot specify both 'fixed.only' and 'random.only")
  }
  if(fixed.only){
    formula = nobars(form)
  }
  if(random.only){
    formula = reOnly(form, response = T)
  }else{
    formula = form
  }
  
  return(formula)
  
}

#' @importFrom stats model.frame
#' @export
model.frame.pglmmObj = function(object, fixed.only = F){
  frame = object$frame
  # Borrowed some code from lme4
  if(fixed.only){
    form = formula(object, fixed.only = T)
    vars.fixed = rownames(attr(terms.formula(form), "factors"))
    frame = frame[vars.fixed]
  }
  
  return(frame)
}

randomListRaw = function(object){
  bars = findbars(object$formula)
  ref = ranef(object)
  vars.random = lapply(ref, function(j) colnames(j))
  # Assumption: only one group allowed
  frame_rand = lapply(bars, function(j) object$frame[vars.random[[j]]])
  return(frame_rand)
}

#' @importFrom stats model.matrix
#' @export
model.matrix.pglmmObj = function(object, type = c("fixed", "random", "randomListRaw")) {
  # if(length(type) > 1){
  #   warning("only the first element of 'type' will be used")
  # }
  switch(type[1],
         "fixed" = object$X,
         "random" = Matrix(object$Z, sparse = T),
         "randomListRaw" = randomListRaw(object))
  }

# Functions to add:
## offset() - will not have $offset option
## weights() - will not have $weights option
### Ideas for offset() and weights() - extract from model.frame ?