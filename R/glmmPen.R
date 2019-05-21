##' @export

glmmPen = function(formula, data, family = "binomial", na.action = "na.omit", 
                  offset = NULL, weights = NULL,
                  lambda0 = 0, lambda1 = 0, nMC = 100, nMC_max = 2000, returnMC = T, gibbs = T,
                  maxitEM = 100, trace = 0, vartol = 0.00, conv = 0.001, pnonzerovar = 0, alpha = 1){
  
  # Input modification and restriction for family
  if(is.character(family)){
    family = get(family, mode = "function", envir = parent.frame())
    # Include number in parent.frame (2)?
  }
  if(is.function(family)){
    family = family()
  }
  if(!(family %in% c(binomial()))){
    stop(sprintf("'%s' family is not available in pglmer function", family))
  }
  
  # Acceptable input types and input restrictions - vectors, integers, positive numbers ...
  if(class(data) != "data.frame"){
    stop("data must be of class 'data.frame'")
  }
  # if(!(group %in% colnames(data))){
  #   stop("'group' must be a column in data")
  # }
  # group = data$groupID
  # if(!is.vector(group) | !is.numeric(group) | !is.factor(group)){
  #   stop("'group' must be a numeric factor vector") 
  #   # What exactly does the group vector need to be like?
  # }
  if(!is.integer(c(nMC, nMC_max, maxitEM)) | !is.positive(c(nMC, nMC_max, maxitEM))){
    stop("nMC, nMC_max, and maxitEM must be positive integers")
  }
  if(nMC_max < nMC){
    warning("nMC_max should not be smaller than nMC \n", immediate. = T)
  }
  if(is.negative(c(lambda0, lambda1))){
    stop("lambda0 and lambda1 cannot be negative")
  }
  
  # Other processing - match.call
  # call = match.call(expand.dots = F)
  # Needs work? What does this do, and how can it be used properly?
  
  # Deal with NAs
  if(na.action = "na.omit"){ # Need to use character? Is na.action a function? ...
    na.omit(data)
  }else{
    warning("This function not equipted to deal with NA values. \n
            Please check that data does not contain NA values. \n", immediate. = T)
  }
  
  # Convert formula and data to useful forms to plug into fit_dat
  ## substitute | for + 
  mod_frame_full = subbars(formula)
  frame_full = model.frame(mod_frame_full, offset = offset, weights = weights)
  
  ## Identify random effects
  ### If no | (no random effects listed) then stop(suggestion: use glmnet or ncvreg package instead)
  reExprs = findbars(formula)
  reGrpList = mkReTrms_new(reExprs, frame_full)
  ### Make sure only one group is specified
  if(length(reGrpList) > 1){
    stop("only 1 grouping variable can be specified")
  }else{group = reGrpList[[1]]$ff}
  Z = reGrpList[[1]]$Z
  rownames(Z) = rownames(data)
  cnms = reGrpList[[1]]$cnms
  ### Check regarding including intercept?
  
  ## Get fixed effects X matrix
  formula_nobars = nobars(formula)
  X = model.matrix(formula_nobars, data)
  
  ## Make sure colnames Z subset of colnames X
  if(!(cnms %in% colnames(X))){
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
    stop("glmmPen functino currently not enabled to use offset and weight parameters")
  }
  
  
  data_input = list(Y, X, Z)
  
  # Call fit_dat function - adjust to use match.call object?
  # fit_dat object found in "/R/fit_dat.R" file
  output = fit_dat(dat = data_input, lambda0_scad = lambda0, lambda1_scad = lambda1, nMC = nMC, 
                   family = family, group = group, trace = trace, vartol = vartol, 
                   conv = conv, nMC_max = nMC_max, returnMC = returnMC, gibbs = gibbs, 
                   pnonzerovar = pnonzerovar, maxitEM = maxitEM, alpha = alpha)
  
  # Format Output - MerMod object ...
  
  return(output)
  
  }

```