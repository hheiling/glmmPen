
# Edited versions of lFormula and glFormula functions from lme4 package
# lme4 code source: https://github.com/lme4/lme4/blob/master/R/modular.R
# Additional helper functions: https://github.com/lme4/lme4/blob/master/R/utilities.R 
# Edits:
#   Removed some checks on number of random effects
#   Removed predvars attributes from model frame


# check that the number of levels of the grouping factor(s) are at least 2 and < number of 
#   observations in the data. Edited version of checkNLevels from lme4
checkGrpLevels = function(flist, n){
  stopifnot(is.numeric(n))
  nlevelVec <- vapply(flist, function(x) nlevels(factor(x, exclude=NA)), 1)
  
  if(any(nlevelVec < 2)){
    stop("The grouping factor must have >= 2 levels")
  }
  
  # For now, allow n levels
  if(any(nlevelVec > n)){
    stop("number of levels of grouping factor must be <= number observations")
  }
}

# check if bad formula or data input
# Borrowed from lme4 'utilities.R' file https://github.com/lme4/lme4/blob/master/R/utilities.R
checkFormulaData <- function(formula, data, checkLHS=TRUE,
                             checkData=TRUE, debug=FALSE) {
  wd <- tryCatch(force(data), error = identity)
  if (bad.data <- inherits(wd,"error")) {
    bad.data.msg <- wd$message
  }
  
  ## data not found (this *should* only happen with garbage input,
  ## OR when strings used as formulae -> drop1/update/etc.)
  ##
  if (bad.data || debug) {
    varex <- function(v, env) exists(v, envir=env, inherits=FALSE)
    allvars <- all.vars(as.formula(formula))
    allvarex <- function(env, vvec=allvars) all(vapply(vvec, varex, NA, env))
  }
  if (bad.data) { ## Choose helpful error message:
    if (allvarex(environment(formula))) {
      stop("bad 'data', but variables found in environment of formula: ",
           "try specifying 'formula' as a formula rather ",
           "than a string in the original model",call.=FALSE)
    } else {
      stop("bad 'data': ", bad.data.msg, call. = FALSE)
    }
  } else {
    denv <- ## The data as environment
      if (is.null(data)) {
        if (!is.null(ee <- environment(formula))) {
          ee ## use environment of formula
        } else {
          ## e.g. no environment, e.g. because formula is a character vector
          ## parent.frame(2L) works because [g]lFormula (our calling environment)
          ## has been called within [g]lmer with env=parent.frame(1L)
          ## If you call checkFormulaData in some other bizarre way such that
          ## parent.frame(2L) is *not* OK, you deserve what you get
          ## calling checkFormulaData directly from the global
          ## environment should be OK, since trying to go up beyond the global
          ## environment keeps bringing you back to the global environment ...
          parent.frame(2L)
        }
      } else ## data specified
        list2env(data)
  }
  ##
  ## FIXME: set enclosing environment of denv to environment(formula), or parent.frame(2L) ?
  if (debug) {
    cat("Debugging parent frames in checkFormulaData:\n")
    ## find global environment -- could do this with sys.nframe() ?
    glEnv <- 1L
    while (!identical(parent.frame(glEnv),.GlobalEnv)) {
      glEnv <- glEnv+1L
    }
    ## where are vars?
    for (i in 1:glEnv) {
      OK <- allvarex(parent.frame(i))
      cat("vars exist in parent frame ", i)
      if (i == glEnv) cat(" (global)")
      cat(" ",OK, "\n")
    }
    cat("vars exist in env of formula ", allvarex(denv), "\n")
  } ## if (debug)
  
  stopifnot(!checkLHS || length(as.formula(formula,env=denv)) == 3)  ## check for two-sided formula
  return(denv)
}

# checks for model.matrix X
# Note: these checks of X will also provide a check for Z because Z is later restricted to 
# contain variables that are a subset of X (see fD_adj function in glmmPen.R R script)
checkXmatrix = function(X){
  # For now, do not allow input of character variables into X model matrix
  if (typeof(X)=="character") stop("input variables must be numeric", call.=FALSE)
  # Make sure individuals did not input an additional intercept or a variable with only one value
  # for all observations
  col_vars = apply(X[,-1,drop=F], 2, var, na.rm=TRUE)
  if(any(col_vars == 0)){
    # If true, more than one column with zero variance
    # If only one column, class(col_vars) == "numeric"
    stop("Variable(s) in formula has zero variance (constant column).
         Remove this variable from formula")
  }
}

# Edited version of lme4 package function glFormula
#' @importFrom lme4 factorize mkReTrms nobars subbars findbars
glFormula_edit <- function(formula, data=NULL, family,
                           subset, weights, na.action, offset,
                           contrasts = NULL, ...) {
  
  # Edits by Hillary:
  # Remove/change some checks
  ## Allow p > n for the random effects
  ## Don't check scale of X (will perform scaling of X later in glmmPen function)
  # Removed predvars.fixed and predvars.random attributes to model frame
  
  mf <- mc <- match.call()
  ## Note: input family is the family function passed from glmmPen
  if (family$family %in% c("quasibinomial", "quasipoisson", "quasi"))
    stop('"quasi" families cannot be used in glmmPen package')
  
  denv <- checkFormulaData(formula, data, checkLHS = T, checkData = T, debug = F)
  mc$formula <- formula <- as.formula(formula, env = denv)    ## substitute evaluated version
  
  
  m <- match(c("data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  fr.form <- subbars(formula) # substitute "|" by "+"
  environment(fr.form) <- environment(formula)
  ## model.frame.default looks for these objects in the environment
  ## of the *formula* (see 'extras', which is anything passed in '...'),
  ## so they have to be put there:
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x=.(i)))))
      assign(i, get(i, parent.frame()), environment(fr.form))
  }
  mf$formula <- fr.form
  fr <- eval(mf, parent.frame())
  ## convert character vectors to factor (defensive)
  fr <- factorize(fr.form, fr, char.only = TRUE)
  ## store full, original formula & offset
  attr(fr,"formula") <- formula
  attr(fr,"offset") <- mf$offset
  
  n <- nrow(fr)
  ## random effects and terms modules
  reTrms <- mkReTrms(findbars(formula), fr)
  wmsgNlev <- checkGrpLevels(reTrms$flist, n = n)
  
  ## fixed-effects model matrix X - remove random effect parts from formula:
  fixedform <- nobars(formula)
  X <- model.matrix(fixedform, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
  checkXmatrix(X)
  
  list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
       wmsgs = c(Nlev = wmsgNlev))
}

#####################################################################################################