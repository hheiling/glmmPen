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

#' @export
etaCalc = function(X, Z, beta, U){

  if(class(U) == "numeric"){ # U = vector instead of matrix; highly unlikely
    gamma = mean(U)
  }else{
    gamma = colMeans(U)
  }
  
  eta = X %*% beta + Z %*% gamma
  
  return(eta)
}

#' @export
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
  
  if(!is.null(newdata) && class(newdata) != "data.frame"){
    stop("newdata must be a dataframe")
  }
  if(!is.null(re.form)){
    stop("option re.form not available at this time")
  }

  if(is.null(newdata) && is.null(re.form)){
    pred = switch(type[1], # if unspecified, default = link output (linear predictor)
                  response = fitted(object),
                  link = etaCalc(X = object$X, Z = object$Z, beta = fixef(object), U = object$gibbs_mcmc))
  }else if(!is.null(newdata) && is.null(re.form)){
      fD_out = formulaData(formula = formula(object), data = newdata, na.action = na.action)
      # Assume only one group
      if(nlevels(fD_out$group) != nlevels(object$group[[1]])){
        levs_orig = levels(object$group[[1]])
        levs_new = levels(fD_out$group)
        # Columns of Z and object$gibbs_mcmc organized first by vars, then by group within vars
        present_levs = as.numeric(levs_new %in% levs_orig)
        keep_cols = rep(present_levs, each = ncol(object$Z)/nlevels(object$group[[1]]))
        U = object$gibbs_mcmc[,keep_cols]
      }else{
        U = object$gibbs_mcmc
      }
      eta = etaCalc(fD_out$X, fD_out$Z, beta = fixef(object), U = U)
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
      # In future, create code for option else(!is.null(re.form))
  }
  return(pred)
}

#' @importFrom stats predict
#' @export
residuals.pglmmObj = function(object, type = "response"){
  # Add more type options?
  Y = object$y
  mu = fitted(object)
  res = Y - mu
  return(res)
}

################################################################### 
# Print functions:
# cat_f = function(...){
#   cat(..., fill = T)
# }

prt_family = function(object){
  f = object$family
  if(is.character(f)){
    fam = get(f, mode = "function")
    family = fam()
  }else if(is.function(f)){
    family = f()
  }
  cat(" Family:", family$family, paste(" (", family$link, ")"), fill = T)
}

prt_call <- function(object) {
  # Copied some code from lme4 source code
  call = object$call
  if (!is.null(cc <- call$formula))
    cat("Formula:", deparse(cc), fill = T)
  if (!is.null(cc <- call$data))
    cat("   Data:", deparse(cc), fill = T)
  if (!is.null(cc <- call$weights))
    cat("Weights:", deparse(cc), fill = T)
  if (!is.null(cc <- call$offset))
    cat(" Offset:", deparse(cc), fill = T)
  # if (!is.null(cc <- call$subset))
  #   cat(" Subset:", deparse(cc), fill = T)
}

prt_fixef = function(object, digits){
  if(length(fef <- fixef(object)) > 0) {
    cat("Fixed Effects:\n")
    print.default(format(fef, digits = digits),
                  print.gap = 2L, quote = FALSE)
  } else cat("No fixed effect coefficients\n")
}

prt_ranef = function(object, digits = 4){
  cat("Random Effects:\n")
  # Create character matrix
  cnms = c("Group","Name","Variance")
  # Assume only one group
  group = object$group[[1]] 
  group_name = names(object$group)
  ref = lapply(ranef(object), function(x) colnames(x))
  ref_names = ref[[1]]
  sigma = object$sigma
  output = matrix(0, nrow = length(ref_names), ncol = length(cnms))
  output[,2] = ref_names
  output[,1] = group_name
  output[,3] = round(diag(sigma), digits = digits)
  colnames(output) = cnms
  rownames(output) = rep("", nrow(output))
  
  print(output, quote = F)
}

prt_nobsgrps = function(object){
  cat(sprintf("Number Observations: %i,  groups: %s, %i", 
              nobs(object), names(object$group), nlevels(object$group[[1]])))
}

#' importFrom stats print
#' @export 
print.pglmmObj = function(object, digits = c(4,4)){
  # ToDo: Add in (best) lambda values, BIC, logLik
  
  # Title
  cat("Penalized generalized linear mixed model fit by Monte Carlo Expectation Conditional Minimization (MCECM)",  
      "  algorithm", " (", object$sampling, ") ", " ['", class(out_glmmPen), "'] ", fill = T, sep = "")
  # Family
  prt_family(object)
  # Call information: formula, data, weights, offset, (subset?)
  prt_call(object)
  # Fixed effects information
  prt_fixef(object, digits = digits[[1]])
  # Random effects information
  prt_ranef(object, digits = digits[[2]])
  # Number obs and groups
  prt_nobsgrps(object)
  
  invisible(object)
}

# # @exportMethod show
# setMethod("show",  "pglmmObj", function(object) print.pglmmObj(object))

# Functions to add:
## offset() - will not have $offset option
## weights() - will not have $weights option
### Ideas for offset() and weights() - extract from model.frame ?

