##' @export
glmmPen = function(formula, data, family = "binomial", na.action = "na.omit",
                  offset = NULL, weights = NULL, # penalty,
                  lambda0 = 0, lambda1 = 0, nMC = 100, nMC_max = 2000, returnMC = T, gibbs = T,
                  maxitEM = 100, trace = 0, vartol = 0.00, conv = 0.001, pnonzerovar = 0, 
                  alpha = 1){
  # Things to address / Questions to answer:
  ## dat$pnonzero and pnonzerovar should equal ... ?
  ## Add option for different penalties
  ## Specify what fit_dat output will be
  ## Provide option for offset, weights
  
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
  if(lambda0 < 0 | lambda1 < 0){
    stop("lambda0 and lambda1 cannot be negative")
  }

  # Deal with NAs
  if(na.action == "na.omit"){ # Need to use character? Is na.action a function? ...
    na.omit(data)
  }else{
    warning("This function not equipted to deal with NA values. \n
            Please check that data does not contain NA values. \n", immediate. = T)
  }

  # Convert formula and data to useful forms to plug into fit_dat
  ## substitute | for +
  mod_frame_full = subbars(formula)
  ## Add offset = offset, weights = weights in model.frame function?
  frame_full = model.frame(mod_frame_full, data = data)

  ## Identify random effects
  ### If no | (no random effects listed) then stop(suggestion: use glmnet or ncvreg package instead)
  reExprs = findbars(formula)
  reGrpList = mkReTrms_glmmPen(reExprs, frame_full)
  
  Z = t(as.matrix(reGrpList$Zt))
  
  # Change group condition below?
  if(length(reGrpList$flist) > 1){
    stop("procedure can only handle one group")
  }else{
    group = reGrpList$flist[[1]]
    group_name = names(reGrpList$flist)
  }
  
  ## Get fixed effects X matrix
  formula_nobars = nobars(formula)
  X = model.matrix(formula_nobars, data)
  # Problem if variable specified in formula is the intercept / a constant column 
  # and -1 not included in formula
  constant_cols = X[,apply(X, 2, var, na.rm=TRUE) == 0]
  ## Change to matrix (X not data.frame)
  if(class(constant_cols) == "data.frame"){ # If true, more than one column with zero variance
    stop("Variable(s) in formula has zero variance (constant column).
         Either remove this variable from formula or specify -1 in formula")
  } # If only one column, class(constant_cols) == "numeric"

  ## Make sure colnames random effects subset of colnames X
  cnms = reGrpList$cnms[[1]]
  if(sum(!(cnms %in% colnames(X))) > 0){
    stop("random effects must be a subset of fixed effects")
  }

  ## Creating response variable Y
  Y = model.response(frame_full)
  ### Other options if including offset and/or weights?
  ### Use mkRespMod if offset / weights used
  if(is.null(offset) && is.null(weights)){
    Y = model.response(frame_full)
  }else{
    # This currently doesn't work; need to experiment with and fix
    # Y = lme4::mkRespMod(frame_full, family = family, y = model.response(frame_full),
    #               offset = offset, weights = weights)
    stop("glmmPen function currently not enabled to use offset and weight parameters")
  }

  data_input = list(y = Y, X = X, Z = Z, group = group, pnonzero = ncol(X))
  
  coef_names = list(fixed = colnames(X), random = cnms, group = group_name)
  
  # Things that should be included in call:
  ## formula, data, any other items included in glmmPen function call
  call = match.call(expand.dots = F)

  # Call fit_dat function - adjust to use match.call object?
  # fit_dat object found in "/R/fit_dat.R" file
  fit = fit_dat(dat = data_input, lambda0_scad = lambda0, lambda1_scad = lambda1, nMC = nMC,
                   family = family, group = group, trace = trace, vartol = vartol,
                   conv = conv, nMC_max = nMC_max, returnMC = returnMC, gibbs = gibbs,
                   pnonzerovar = ncol(Z), maxitEM = maxitEM, alpha = alpha)
  
  # Things that should be included in fit_dat:
  ## (fill in later)

  # Format Output - create pglmmObj object
  output = c(fit, list(call = call, formula = formula, data = data, Y = Y, X = X, Z = Z,
                       group = group, coef_names = coef_names, family = family,
                       offset = offset, weights = weights))

  
  out_object = pglmmObj$new(output)
  return(out_object)

  }