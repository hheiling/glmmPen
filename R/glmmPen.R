#' @importFrom stringr str_to_lower
#' @export
glmmPen = function(formula, data = NULL, family = "binomial", na.action = na.omit,
                  offset = NULL, weights = NULL, penalty = "grMCP",
                  lambda0 = 0, lambda1 = 0, nMC = 100, nMC_max = 2000, returnMC = T,
                  maxitEM = 100, trace = 0, conv = 0.001, 
                  alpha = 1, gibbs = NULL){
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
  ## Even if already numeric, convert to levels 1,2,3,... (consecutive integers)
  group_num = as.factor(as.numeric(fD_out$group))
  # if(any(is.character(fD_out$group))){
  #   group = as.factor(as.numeric(fD_out$group))
  # }else{
  #   group = fD_out$group
  # }

  data_input = list(y = fD_out$y, X = fD_out$X, Z = fD_out$Z, group = group_num)
  
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
  
  # Call fit_dat function - adjust to use match.call object?
  # fit_dat object found in "/R/fit_dat.R" file
  fit = fit_dat(dat = data_input, lambda0 = lambda0, lambda1 = lambda1, nMC = nMC,
                   family = family, trace = trace, penalty = penalty,
                   conv = conv, nMC_max = nMC_max, returnMC = returnMC, gibbs = gibbs,
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

