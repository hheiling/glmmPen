#' @importFrom stringr str_to_lower
#' @export
glmmPen = function(formula, data = NULL, family = "binomial", na.action = na.omit,
                  offset = NULL, weights = NULL, penalty = "grMCP",
                  lambda0 = 0, lambda1 = 0, nMC = 100, nMC_max = 2000, returnMC = T,
                  maxitEM = 100, trace = 0, conv = 0.001, 
                  alpha = 1){
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
  if(lambda0 < 0 | lambda1 < 0){
    stop("lambda0 and lambda1 cannot be negative")
  }

  # Convert formula and data to useful forms to plug into fit_dat
  fD_out = formulaData(formula, data, na.action)
  
  ## Convert group to numeric factor - for fit_dat
  if(any(is.character(fD_out$group))){
    group = as.factor(as.numeric(fD_out$group))
  }else{
    group = fD_out$group
  }

  data_input = list(y = fD_out$y, X = fD_out$X, Z = fD_out$Z, group = group)
  
  coef_names = list(fixed = colnames(fD_out$X), random = fD_out$cnms, group = fD_out$group_name)
  
  # Things that should be included in call:
  ## formula, data, any other items included in glmmPen function call
  call = match.call(expand.dots = F)
  
  # Call fit_dat function - adjust to use match.call object?
  # fit_dat object found in "/R/fit_dat.R" file
  fit = fit_dat(dat = data_input, lambda0 = lambda0, lambda1 = lambda1, nMC = nMC,
                   family = family, trace = trace, penalty = penalty,
                   conv = conv, nMC_max = nMC_max, returnMC = returnMC, gibbs = T,
                   maxitEM = maxitEM, alpha = alpha)
  
  # Things that should be included in fit_dat:
  ## (fill in later)

  # Format Output - create pglmmObj object
  output = c(fit, list(call = call, formula = formula, data = data, Y = fD_out$y, 
                       X = fD_out$X, Z = fD_out$Z, group = fD_out$flist, 
                       coef_names = coef_names, family = family,
                       offset = offset, weights = weights, frame = fD_out$frame))

  out_object = pglmmObj$new(output)
  return(out_object)

}

# ## substitute | for +
# form_full = subbars(formula)
# 
# # Deal with NAs
# if(na.action == na.omit){ # Need to use character? Is na.action a function? ...
#   data = na.omit(data[,colnames(model.frame(form_full, data = data))])
# }else{
#   warning("This function not equipted to deal with NA values. \n
#           Please check that data does not contain NA values. \n", immediate. = T)
# }
# 
# ## Add offset = offset, weights = weights in model.frame function?
# frame_full = model.frame(form_full, data = data)
# 
# ## Identify random effects
# ### If no | (no random effects listed) then stop - mkBlist called by mkReTrms give this error
# reExprs = findbars(formula)
# reTrms = mkReTrms(reExprs, frame_full)
# # t(Zt) from mkReTrms: columns organized by group level, then vars within group level
# Zt = reTrms$Zt
# 
# # Change group condition below?
# if(length(reTrms$flist) > 1){
#   stop("procedure can only handle one group")
# }else{
#   group = reTrms$flist[[1]]
#   group_name = names(reTrms$flist)
# }
# 
# d = nlevels(group[[1]])
# Z = Matrix(0, nrow = ncol(Zt), ncol = nrow(Zt), sparse = T)
# # Want Z columns organized by vars, then levels of group within vars
# for(lev in 1:d){
#   Z[,(d*(lev-1)+1):(d*lev)] = Matrix::t(Zt[seq(lev, by = d, length.out = nrow(Zt)/d),])
# }
# # Z_dense = Matrix::as.matrix(Z) # Convert Z to dense matrix
# 
# ## Get fixed effects X matrix
# formula_nobars = nobars(formula)
# X = model.matrix(formula_nobars, data)
# 
# # Problem if variable specified in formula is the intercept / a constant column 
# # and -1 not included in formula
# constant_cols = X[,apply(X, 2, var, na.rm=TRUE) == 0]
# ## Change to matrix (X not data.frame)
# if(class(constant_cols) == "matrix"){ # If true, more than one column with zero variance
#   stop("Variable(s) in formula has zero variance (constant column).
#        Either remove this variable from formula or specify -1 in formula")
# } # If only one column, class(constant_cols) == "numeric"
# 
# ## Make sure colnames random effects subset of colnames X
# cnms = reTrms$cnms[[1]]
# if(sum(!(cnms %in% colnames(X))) > 0){
#   stop("random effects must be a subset of fixed effects")
# }
# 
# ## Creating response variable Y
# Y = model.response(frame_full)
# ### Other options if including offset and/or weights?
# ### Use mkRespMod if offset / weights used
# if(is.null(offset) && is.null(weights)){
#   Y = model.response(frame_full)
# }else{
#   # This currently doesn't work; need to experiment with and fix
#   # Y = lme4::mkRespMod(frame_full, family = family, y = model.response(frame_full),
#   #               offset = offset, weights = weights)
#   stop("glmmPen function currently not enabled to use offset and weight parameters")
# }
# 
