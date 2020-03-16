
#' Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
#' Minimization (MCECM)
#' 
#' \code{glmmPen} is used to fit penalized generalized mixed models via Monte Carlo Expectation 
#' Conditional Minimization (MCECM)
#' 
#' @inheritParams formulaData
#' @inheritParams fit_dat
#' @param offset an optional vector of an \emph{a priori} known component to be included in the 
#' linear predictor during fitting. One or more \code{\link[stats]{offset}} terms can be included in 
#' formula instead or as well. If more than one is specified, their sum is use. 
#' See \code{\link[stats]{model.offset}}.
#' @param weights an optional vector of 'prior weights' to be used in the fitting process.
#' Should be \code{NULL} or a numeric vector.
#' @param control a list (of correct class, resulting from lambdaControl() or selectControl()) 
#' containing control parameters. If the user wants to run the algorithm using one specific set of
#' penalty parameters \code{lambda0} and \code{lambda1}, then use \code{lambdaControl()}. 
#' If the user wants to run the algorithm over multiple possible \code{lambda0} and \code{lambda1},
#' then use \code{lambdaControl{}}. See the \code{\link{lambdaControl}} and \code{\link{selectControl}}
#' documentation for details
#' 
## #' @inheritSection fit_dat
#' 
#' @return An reference class object of class \code{\link{pglmmObj}} for which many methods are 
#' available (e.g. \code{methods(class = "pglmmObj")})
#'  
#' @importFrom stringr str_to_lower
#' @export
glmmPen = function(formula, data = NULL, family = "binomial", na.action = na.omit,
                  offset = NULL, weights = NULL, penalty = "grMCP",
                  nMC = 100, nMC_max = 2000, returnMC = T,
                  maxitEM = 100, trace = 0, conv = 0.001, 
                  alpha = 1, gibbs = NULL, control = lambdaControl()
                  ){
  # Things to address / Questions to answer:
  ## Add option for different penalties
  ## Specify what fit_dat output will be
  ## Provide option for offset, weights
  ## gibbs T / F: internally specify T or F depending on X and Z dimensions
  
  # Input modification and restriction for family
  if(is.character(family)){
    # library(stringr)
    family = str_to_lower(family)
  }
  if(is.function(family)){
    family = family$family
  }
  if(!(family %in% c("binomial"))){
    print(family)
    stop("'family' not recognized")
  }

  # Acceptable input types and input restrictions - vectors, integers, positive numbers ...
  if(class(data) != "data.frame"){
    stop("data must be of class 'data.frame'")
  }
  if(sum(c(nMC, nMC_max, maxitEM) %% 1) > 0 | sum(c(nMC, nMC_max, maxitEM) <= 0) > 0){
    stop("nMC, nMC_max, and maxitEM must be positive integers")
  }
  if(nMC_max < nMC){
    warning("nMC_max should not be smaller than nMC \n", immediate. = T)
  }

  # Convert formula and data to useful forms to plug into fit_dat
  fD_out = formulaData(formula, data, na.action)
  
  ## Convert group to numeric factor - for fit_dat
  ## Even if already numeric, convert to levels 1,2,3,... (consecutive integers)
  group_num = as.factor(as.numeric(fD_out$group))
  # if(any(is.character(fD_out$group))){
  #   group = as.factor(as.numeric(fD_out$group))
  # }else{
  #   group = fD_out$group
  # }
  
  # Standardize X and Z
  std_out = XZ_std(fD_out, group_num)
  
  data_input = list(y = fD_out$y, X = std_out$X_std, Z = std_out$Z_std, group = group_num)
  
  coef_names = list(fixed = colnames(fD_out$X), random = fD_out$cnms, group = fD_out$group_name)
  
  if(is.null(gibbs)){
    if(length(fD_out$cnms) > 20){
      gibbs = T
    }else{
      gibbs = F
    }
  }else if(!is.logical(gibbs)){
    stop("gibbs is a logical variable; must be TRUE or FALSE")
  }
  
  # Things that should be included in call:
  ## formula, data, any other items included in glmmPen function call
  call = match.call(expand.dots = F)
  
  # rho = environment()
  if(!is.list(control) | !inherits(control, "pglmmControl")){
    stop("control parameter must be a list of type lambdaControl() or selectControl()")
  }
  
  if(inherits(control, "lambdaControl")){
    lambda0 = control$lambda0
    lambda1 = control$lambda1
    if(lambda0 < 0 | lambda1 < 0){
      stop("lambda0 and lambda1 cannot be negative")
    }
    # Call fit_dat function
    # fit_dat function found in "/R/fit_dat.R" file
    fit = fit_dat(dat = data_input, lambda0 = lambda0, lambda1 = lambda1, nMC = nMC,
                  family = family, trace = trace, penalty = penalty,
                  conv = conv, nMC_max = nMC_max, returnMC = returnMC, gibbs = gibbs,
                  maxitEM = maxitEM, alpha = alpha)
    
  }else if(inherits(control, "selectControl")){
    stop("selectControl option not yet available")
    
    if(is.null(control$lambda0_seq)){
      lambda0_range = LambdaRange(X = data_input$X, y = data_input$y, nlambda = control$nlambda)
    }else{
      lambda0_range = control$lambda0_seq
      if(!is.numeric(lambda0_range) | any(lambda0_range < 0)){
        stop("lambda0_seq must be a positive numeric sequence")
      }
    }
    if(is.null(control$lambda1_seq)){
      lambda1_range = lambda0_range
    }else{
      lambda1_range = control$lambda1_seq
      if(!is.numeric(lambda1_range) | any(lambda1_range < 0)){
        stop("lambda1_seq must be a positive numeric sequence")
      }
    }
    
    fit = select_tune(dat = data_input, lambda0_range = lambda0_range, lambda1 = lambda1_range,
                      penalty = penalty, returnMC = returnMC,
                      nMC = nMC, nMC_max = nMC_max, family = family, trace = trace,
                      ufull = NULL, coeffull = NULL, gibbs = gibbs,
                      maxitEM = maxitEM, alpha = alpha)
    
  }
  
  # Things that should be included in fit_dat:
  ## (fill in later)
  
  if(gibbs){
    sampling = "Gibbs Sampling"
  }else{
    if(fit$rej_to_gibbs < 3){
      sampling = "Rejection Sampling"
    }else{
      sampling = "Gibbs Sampling"
    }
  }
  
  # Format Output - create pglmmObj object
  output = c(fit, list(call = call, formula = formula, data = data, Y = fD_out$y, 
                       X = fD_out$X, Z = fD_out$Z, group = fD_out$flist, 
                       coef_names = coef_names, family = family,
                       offset = offset, weights = weights, frame = fD_out$frame,
                       sampling = sampling, std_out = std_out))

  out_object = pglmmObj$new(output)
  return(out_object)

}

#' @name lambdaControl
#' @aliases selectControl
#' 
#' @title Control of Penalized Generalized Linear Mixed Model Fitting
#' 
#' Constructs control structures for penalized mixed model fitting.
#' 
#' @inheritParams fit_dat
#' @param lambda0_seq a range of non-negative numeric penalty parameter for the fixed 
#' effects parameters. If \code{NULL}, then a range will be calculated as described by (...).
#' @param lambda1_seq a range of non-negative numeric penalty parameter for the grouped random 
#' effects covariance parameters. If \code{NULL}, then a range will be calculated as described 
#' by (...).
#' 
#' @return The *Control functions return a list (inheriting from class "\code{pglmmControl}") 
#' containing penalization parameter values, presented either as an individual set or as a range of
#' possible values.
#' 
#' @export
lambdaControl = function(lambda0 = 0, lambda1 = 0){
  structure(list(lambda0 = lambda0, lambda1 = lambda1), 
            class = c("lambdaControl","pglmmControl"))
}

#' @rdname lambdaControl
#' @export
selectControl = function(lambda0_seq = NULL, lambda1_seq = NULL, nlambda = 40){
  structure(list(lambda0_seq = lambda0_seq,
                 lambda1_seq = lambda1_seq,
                 nlambda = nlambda),
            class = c("selectControl", "pglmmControl"))
}

#' @importFrom ncvreg std
#' @export
XZ_std = function(fD_out, group_num){
  # Standardize X - ncvreg::std method
  X = fD_out$X
  X_noInt_std = std(X[,-1])
  X_std = cbind(1, X_noInt_std)
  X_center = attr(X_noInt_std, "center")
  X_scale = attr(X_noInt_std, "scale")
  # Note: X_noInt_std = (X[,-1] - X_center) / X_scale
  
  var_subset = (colnames(X) %in% fD_out$cnms)
  Z_center = X_center[var_subset]
  Z_scale = X_scale[var_subset]
  
  # Standardize Z using X_std output
  Z_sparse = fD_out$Z
  d = nlevels(group_num)
  num_vars = ncol(Z_sparse) / d
  
  Z_std = Matrix(data = 0, nrow = nrow(Z_sparse), ncol = ncol(Z_sparse), sparse = T)
  
  for(v in 1:num_vars){ 
    if("(Intercept)" %in% cnms) next # Don't need to scale intercept
    cols = seq(from = (v - 1)*d + 1, to = v*d, by = 1)
    for(k in 1:nlevels(group_num)){
      ids = which(group_num == k)
      Z_std[ids, cols[k]] = (Z_sparse[ids, cols[k]] - Z_center[v-1]) / Z_scale[v-1]
    }
  }
  
  return(list(X_std = X_std, Z_std = Z_std, X_center = X_center, X_scale = X_scale,
              Z_center = Z_center, Z_scale = Z_scale))
}


#' @importFrom ncvreg setupLambda
#' @export
LambdaRange = function(X, y, family, alpha = 1, lambda.min = NULL, nlambda = 40,
                       penalty.factor = NULL){
  # Borrowed elements from `ncvreg` function
  n = nrow(X)
  p = ncol(X)
  
  if(family == "gaussian"){
    yy = y = mean(y)
  }else{
    yy = y
  }
  
  if(is.null(lambda.min)){
    lambda.min = ifelse(n>p, 0.001, 0.05)
  }
  
  if(is.null(penalty.factor)){
    penalty.factor = rep(1, p)
  }
  
  # lambda = calcLambda(X, yy, family, alpha, lambda.min, nlambda, penalty.factor)
  lambda = setupLambda(X, yy, family, alpha, lambda.min, nlambda, penalty.factor)
  
  return(lambda)
  
}