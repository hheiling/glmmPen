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
  }else if(random.only){
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

# randomListRaw = function(object){
#   # Assumption: only one group allowed
#   bars = findbars(object$formula)
#   ref = ranef(object)
#   vars.random = lapply(ref, function(j) colnames(j))
#   frame_rand = lapply(bars, function(j) object$frame[as.factor(vars.random[[j]])])
#   return(frame_rand)
# }

#' @importFrom stats model.matrix
#' @export
model.matrix.pglmmObj = function(object, type = c("fixed", "random", "randomListRaw")) {
  # if(length(type) > 1){
  #   warning("only the first element of 'type' will be used")
  # }
  switch(type[1],
         "fixed" = object$X,
         "random" = Matrix(object$Z, sparse = T),
         "randomListRaw" = stop("'randomListRaw' option not available at this time"))
}

etaCalc = function(X, Z, beta, U){

  if(class(U) == "numeric"){ # U = vector instead of matrix; highly unlikely
    gamma = mean(U)
  }else{
    gamma = colMeans(U)
  }
  
  eta = X %*% beta + Z %*% gamma
  
  return(eta)
}

invLink = function(family, eta){
  # Inverse link functions (canonical links only)
  if(family == "binomial"){
    mu = exp(eta) / (1+exp(eta))
  }else if(family == "poisson"){
    mu = exp(eta)
  }else if(family == "gaussian"){
    mu = eta
  }else{
    stop("fitted means for given family not available")
  }
  
  return(mu)
}

#' @importFrom stats fitted
#' @export
fitted.pglmmObj = function(object){
  ## ToDo: add names/rownames (from X?)
  ## ToDo: Take into account potential offset (in etaCalc)
  
  eta = etaCalc(X = object$X, Z = object$Z, beta = fixef(object), U = object$gibbs_mcmc)
  mu = invLink(family = family(object), eta)
  
  return(mu)
  
}

#' @importFrom stats predict
#' @export
predict.pglmmObj = function(object, newdata = NULL, type = c("link","response"),
                            re.form = NULL, na.action = na.pass){
  ## Other arguments used by lme4: re.form, random.only = F, allow.new.levels = F, newparams = NULL
  
  if(!is.null(newdata) && class(newdata) != data.frame){
    stop("newdata must be a dataframe")
  }
  if(!is.null(re.form)){
    stop("option re.form not available at this time")
  }

  if(is.null(newdata) && is.null(re.form)){
    pred = switch(type[1], # if unspecified, default = link output (linear predictor)
                  response = fitted(object),
                  link = etaCalc(object))
  }else{
    # if(is.null(newdata)){ # Use original data (X matrix and offset)
      # In future, code for is.null(newdata) && !is.null(re.form)
      # X = object$X
      # offset = model.offset(model.frame(object))
      # if(is.null(offset)){offset = 0}
      ## for is.null(newdata), only option not yet covered is !is.null(re.form)
      
    # }else{ # Use newdata
      ## Get fixed effects - fixed vars be same in both original and new data
      form_fixef = formula(object, fixed.only = T)
      
      
      # Deal with NAs
      if(na.action == na.omit){
        data = na.omit(newdata[,colnames(model.frame(formula(object), newdata))])
      }else if(na.action == na.pass){
        data = newdata
      }else{
        print(na.action)
        stop("specified na.action not recognized by function")
      }
      
      ## Get new model frame - full model frame
      mf = model.frame(formula(object), data)
      
      ## Get new fixed effects model.matrix
      X = model.matrix(form_fixef, data)
      ## for !is.null(newdata), both is.null(re.form) and !is.null(re.form) not covered
      ## for now, only have is.null(re.form) option
      if(is.null(re.form)){
        reExprs = findbars(formula(object))
        reTrms = mkReTrms(reExprs, mf)
        # t(Zt) from mkReTrms: columns organized by group level, then vars within group level
        Zt = reTrms$Zt
        # ASsume only one group
        group = reTrms$flist[[1]]
        
        d = nlevels(group)
        Z = Matrix(0, nrow = ncol(Zt), ncol = nrow(Zt), sparse = T)
        # Want Z columns organized by vars, then levels of group within vars
        for(lev in 1:d){
          Z[,(d*(lev-1)+1):(d*lev)] = Matrix::t(Zt[seq(lev, by = d, length.out = nrow(Zt)/d),])
        }
        eta = etaCalc(X, Z, beta = fixef(object), U = object$gibbs_mcmc)
        pred = switch(type, 
                      link = eta,
                      response = invLink(family = family(object), eta = eta))
        if(class(pred) %in% c("dgeMatrix","matrix")){
          if(is.null(rownames(data))){
            rownames(pred) = seq_along(nrow(data))
          }else{
            rownames(pred) = rownames(data)
          }
        }else if(class(pred) == "numeric"){
          if(is.null(rownames(data))){
            names(pred) = seq_along(nrow(data))
          }else{
            names(pred) = rownames(data)
          }
        }
      }
      # In future, create code for option else(!is.null(re.form))
  }
  return(pred)
}

#' @importFrom stats predict
#' @export
residuals.pglmmObj = function(object, type = "response"){
  # Add more type options?
  Y = object$Y
  mu = fitted(object)
  res = Y = mu
  return(res)
}

# Functions to add:
## offset() - will not have $offset option
## weights() - will not have $weights option
### Ideas for offset() and weights() - extract from model.frame ?