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
  
  if(gibbs){
    sampling = "Gibbs Sampling"
  }else{
    sampling = "Rejection Sampling"
  }
  
  
  # Things that should be included in call:
  ## formula, data, any other items included in glmmPen function call
  call = match.call(expand.dots = F)
  
  rho = environment()
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
    # fit_dat object found in "/R/fit_dat.R" file
    fit = fit_dat(dat = data_input, lambda0 = lambda0, lambda1 = lambda1, nMC = nMC,
                  family = family, trace = trace, penalty = penalty,
                  conv = conv, nMC_max = nMC_max, returnMC = returnMC, gibbs = gibbs,
                  maxitEM = maxitEM, alpha = alpha)
  }else if(inherits(control, "selectControl")){
    stop("selectControl option not yet available")
    
    if(is.null(control$lambda0_seq)){
      # Calculate lambda0_range
      const = 10^-3
      lam_max = 2
      lam_min = const * lam_max
      lambda0_range = seq(from = log(lam_max), to = log(lam_min), length.out = 40)
    }else{
      lambda0_range = control$lambda0_seq
      if(!is.numeric(lambda0_range) | any(lambda0_range < 0)){
        stop("lambda0_seq must be a positive numeric sequence")
      }
    }
    if(is.null(control$lambda1_seq)){
      # Calculate lambda0_range
      const = 10^-3
      lam_max = 2
      lam_min = const * lam_max
      lambda1_range = seq(from = log(lam_max), to = log(lam_min), length.out = 40)
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
                      maxitEM = maxitEM, alpha = alpha, environ = rho)
    
  }
  
  # Things that should be included in fit_dat:
  ## (fill in later)

  # Format Output - create pglmmObj object
  output = c(fit, list(call = call, formula = formula, data = data, Y = fD_out$y, 
                       X = fD_out$X, Z = fD_out$Z, group = fD_out$flist, 
                       coef_names = coef_names, family = family,
                       offset = offset, weights = weights, frame = fD_out$frame,
                       sampling = sampling))

  out_object = pglmmObj$new(output)
  return(out_object)

}

#' @export
lambdaControl = function(lambda0 = 0, lambda1 = 0){
  structure(list(lambda0 = lambda0, lambda1 = lambda1), 
            class = c("lambdaConrol","pglmmControl"))
}

#' @export
selectControl = function(lambda0_seq = NULL, lambda1_seq = NULL){
  structure(list(lambda0_seq = lambda0_seq,
                 lambda1_seq = lambda1_seq),
            class = c("selectControl", "pglmmControl"))
}

