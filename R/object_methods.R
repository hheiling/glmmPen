# Functions to define object's S3 methods

# Object name: pglmmObj
# Properties: 
## Reference Class object


#' @importFrom lme4 fixef
#' @export
fixef.pglmmObj = function(object){
  object$fixef
}

#' @importFrom lme4 ranef
#' @export
ranef.pglmmObj = function(object){
  object$ranef
}

#' @importFrom stats sigma
#' @export
sigma.pglmmObj = function(object){
  
  fam_fun = object$family
  family = fam_fun$family
  # Return covariance matrix of the random effects and any scale parameters
  group = object$data$group
  grp = names(group)
  scale = object$scale[[1]]
  out = list()
  out[[grp]] = object$sigma
  if(family == "gaussian"){
    out[["Residual StdDev"]] = sqrt(scale)
  }else if(family == "negbin"){
    out[["phi"]] = scale
  }
  
  return(out)
}

#' @importFrom stats coef
#' @export
coef.pglmmObj = function(object){
  # Find combined coefficients
  ## Borrowed code elements from coef.merMod in lme4
  
  fixef = data.frame(rbind(object$fixef), check.names = F)
  ranef = object$ranef
  group = object$data$group
  
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
    for(nm in nmsi) {
      output[[i]][[nm]] <- output[[i]][[nm]] + refi[,nm]
    }
  }
  
  return(output)
  
}

#' @importFrom stats family
#' @export
family.pglmmObj = function(object){
  object$family
}

#' @export
nobs.pglmmObj = function(object){
  nrow(object$data$X)
}

# @export
# ngrps.pglmmObj = function(object){
#   nlevels(object$data$group[[1]])
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
  frame = object$data$frame
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
  frame = object$data$frame
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
#   frame_rand = lapply(bars, function(j) object$data$frame[as.factor(vars.random[[j]])])
#   return(frame_rand)
# }

#' @importFrom stats model.matrix
#' @export
model.matrix.pglmmObj = function(object, type = c("fixed", "random")) {
  
  switch(type[1],
         "fixed" = object$data$X,
         "random" = Matrix(object$data$Z_std, sparse = T))
}

etaCalc = function(X, Z, beta, U){

  if(class(U) == "numeric"){ # U = vector instead of matrix; highly unlikely
    gamma = mean(U)
  }else{ # U is a matrix
    gamma = colMeans(U)
  }
  
  eta = X %*% beta + Z %*% gamma
  
  return(eta)
}

invLink = function(family, eta){
  fam = family$family
  # Inverse link functions (canonical links only)
  if(fam == "binomial"){
    mu = exp(eta) / (1+exp(eta))
  }else if(fam == "poisson"){
    mu = exp(eta)
  }else if(fam == "gaussian"){
    mu = eta
  }else{
    stop("fitted means for given family not available")
  }
  
  return(mu)
}

#' @importFrom stats fitted
#' @export
fitted.pglmmObj = function(object, fixed.only = T){
  
  offset = object$data$offset
  
  if(fixed.only){
    eta = object$data$X %*% matrix(object$fixef, ncol = 1) + offset
  }else{
    eta = etaCalc(X = object$data$X, Z = object$data$Z_std, beta = object$fixef, U = object$posterior_draws) + offset
  }
  
  mu = invLink(family = object$family, eta)
  
  return(as.numeric(mu))
  
}

#' @importFrom stats predict drop.terms reformulate
#' @importFrom lme4 nobars
#' @export
predict.pglmmObj = function(object, newdata = NULL, type = c("link","response"),
                            fixed.only = T, na.action = na.pass){
  ## Other arguments used by lme4: re.form, random.only = F, allow.new.levels = F, newparams = NULL
  
  if((!is.null(newdata)) & (class(newdata) != "data.frame")){
    stop("newdata must be a dataframe")
  }
  
  if(is.null(newdata)){
    # Calculate prediction result for original data
    if(!fixed.only){ # fixed.only = F
      pred = switch(type[1], # if unspecified, default = link output (linear predictor)
                    link = etaCalc(X = object$data$X, Z = object$data$Z_std, beta = object$fixef, 
                                   U = object$posterior_draws),
                    response = fitted(object, fixed.only = fixed_only))
    }else{ # fixed.only = T
      eta = object$data$X %*% object$fixef
      pred = switch(type[1],
                    link = eta,
                    response = invLink(family = object$family, eta = eta))
    }
    
    data = object$data$X
    
  }else{ # Calculate prediction for newdata
    if(!fixed.only){
      stop("prediction using random effects not appropriate for new data")
    }
    
    fixef = object$fixef
    
    # Make sure newdata has relevant variables needed for prediction
    # fixed_vars: names of the fixed effects used in analysis (from get_all_vars(fixed_formula[-2],data))
    fixed_vars = object$fixed_vars
    # If NA terms in var_names, this indicates that formula was specified using a matrix
    #   (e.g. formula = y ~ X + (X | group)) instead of specified using either the names of column
    #   names of a data frame or named vectors
    if(!any(is.na(fixed_vars))){
      if(!(fixed_vars %in% colnames(newdata))){
        cat("newdata must have the following variables: \n",
            fixed_vars, "\n")
        stop("newdata must have the variables listed in 'fixed_vars' from pglmmObj object")
      }
      
      formula_fixed = nobars(object$formula)
      formula_fixed = formula_fixed[-2] # remove response
      
      # Create model.matrix using newdata and fixed effects formula with non-zero fixed effects
      X = model.matrix(formula_fixed, data = newdata)
      
    }else{ # if any(is.na(fixed_vars))
      # Check that number of columns of newdata matches number of fixed effects
      cat("it is assumed that the columns of newdata match the order of the fixef coefficients")
      # Ignoring intercept (model always contains intercept), check that number of variables in
      # newdata match number of variables in the fixed effects coefficients
      if(ncol(newdata) != (length(fixef)-1)){ 
        stop("number of newdata columns ", ncol(newdata),
             " does not match the number of fixed effects coefficients in fixef (ignoring intercept) ",
             (length(fixef)-1))
      }
    } # End if-else any(is.na(fixed_vars))
    
    X = cbind(1, as.matrix(newdata))
    eta = X %*% fixef
    pred = switch(type[1], 
                  link = eta,
                  response = invLink(family = object$family, eta = eta))
    
    
    data = newdata
    
  } # End if-else is.null(newdata)
  
  if(class(pred) %in% c("Matrix","matrix")){
    if(is.null(rownames(data))){
      rownames(pred) = seq_len(nrow(data))
    }else{
      rownames(pred) = rownames(data)
    }
  }else if(class(pred) == "numeric"){
    if(is.null(rownames(data))){
      names(pred) = seq_len(nrow(data))
    }else{
      names(pred) = rownames(data)
    }
  }
  
  return(pred)
}


# Function for residuals.pglmmObj
var_hat = function(family, mu, sig2 = NULL, phi = NULL){
  if(family == "binomial"){
    v_hat = mu*(1-mu)
  }else if(family == "poisson"){
    v_hat = mu
  }else if(family == "gaussian"){
    v_hat = sig2
  }else if(family == "negbin"){
    v_hat = mu + phi * mu^2 
  }
  return(v_hat)
}

#' @export
residuals.pglmmObj = function(object, type = c("deviance","pearson","response","working")){
  # What is working response?
  y = object$data$y
  mu = Matrix::as.matrix(fitted(object))
  type = match.arg(type)
  fam = family(object)
  family = fam$family
  sig2 = object$scale$Gaus_sig2 # Residual error variance only used if family = gaussian
  phi = object$scale$phi # Only used if family = negbin
  
  if(type == "deviance"){ # reference: mcemGLM package
    if(family == "binomial"){
      res = ifelse(y, 1, -1) * sqrt(-2*(y * log(mu) + (1 - y) * log(1 - mu)))
    }else if(family == "poisson"){
      res = sign(y - mu) * sqrt(2 * ifelse(y > 0, y * log(y/mu), 0) - 2 * (y - mu))
    }else if(family == "gaussian"){
      res = (y - mu) / sqrt(sig2)
    }else if(family == "negbin"){
      # Note: a0 from mcemGLM = 1/phi
      stop("family not available")
      res = sign(y - mu) * sqrt(2 * (ifelse(y > 0, y * log(y/mu), 0)) - 2 * (y + 1/phi) * log((y + 1/phi)/(mu + 1/phi)))
    }else{
      stop("family not available")
    }
  }else if(type == "pearson"){
    res = (y - mu) / sqrt(var_hat(family, mu, sig2, phi))
  }else if(type == "response"){
    res = y - mu
  }else{ # working residuals
    res = (y - mu) / mu
  }
  
  attr(res, "residual type") = type
  return(as.numeric(res))
}

################################################################### 
# Print functions:

prt_family = function(object){
  f = object$family
  cat(" Family:", f$family, paste(" (", f$link, ")"), fill = T)
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
  if(length(fef <- object$fixef) > 0) {
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
  group = object$data$group[[1]] 
  group_name = names(object$data$group)
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
  cat(sprintf("Number Observations: %i,  groups: %s, %i \n", 
              nobs(object), names(object$data$group), nlevels(object$data$group[[1]])))
}

#' @export 
print.pglmmObj = function(object, digits = c(fef = 4, ref = 4)){
  # ToDo: Add in (best) lambda values, BIC, logLik
  
  # Title
  cat("Penalized generalized linear mixed model fit by Monte Carlo Expectation Conditional Minimization (MCECM)",  
      "  algorithm", " (", object$sampling, ") ", " ['", class(object), "'] ", fill = T, sep = "")
  # Family
  prt_family(object)
  # Call information: formula, data, weights, offset, (subset?)
  prt_call(object)
  # Fixed effects information
  prt_fixef(object, digits = digits[1])
  # Random effects information
  prt_ranef(object, digits = digits[2])
  # Number obs and groups
  prt_nobsgrps(object)
  
  invisible(object)
}

summary_ranef = function(object, digits = 4){
  cat("Random Effects:\n")
  # Create character matrix
  cnms = c("Group","Name","Variance","Std.Dev.")
  # Assume only one group
  group = object$data$group[[1]] 
  group_name = names(object$data$group)
  ref = lapply(ranef(object), function(x) colnames(x))
  ref_names = ref[[1]]
  sigma = object$sigma
  output = matrix(0, nrow = length(ref_names), ncol = length(cnms))
  output[,2] = ref_names
  output[,1] = group_name
  output[,3] = round(diag(sigma), digits = digits)
  output[,4] = round(sqrt(diag(sigma)), digits = digits)
  colnames(output) = cnms
  rownames(output) = rep("", nrow(output))
  
  print(output, quote = F)
}

#' @importFrom stringr str_to_title
prt_resids = function(resids, type = "Pearson", digits) {
  cat(str_to_title(type), "residuals: ", "\n", sep = " ")
  
  resid_quant <- setNames(zapsmall(quantile(resids, na.rm=TRUE), digits + 1L),
                          c("Min", "1Q", "Median", "3Q", "Max"))
  print(resid_quant, digits = digits)
  cat("\n")
}

#' @export
summary.pglmmObj = function(object, digits = c(fef = 4, ref = 4), 
                            resid_type = switch(object$family$family, gaussian = "pearson", "deviance")){
  # ToDo: Add in (best) lambda values, BIC, logLik
  
  # Title
  cat("Penalized generalized linear mixed model fit by Monte Carlo Expectation Conditional Minimization (MCECM)",
      "  algorithm", " (", object$sampling, ") ", " ['", class(object), "'] ", fill = T, sep = "")
  # Family
  prt_family(object)
  # Call information: formula, data, weights, offset, (subset?)
  prt_call(object)
  cat("\n")
  # Fixed effects information
  prt_fixef(object, digits = digits[1])
  cat("\n")
  # Random Effects information
  summary_ranef(object, digits = digits[2])
  # Number obs and groups
  prt_nobsgrps(object)
  cat("\n")
  # Residuals Summary
  prt_resids(residuals(object, type = resid_type), type = resid_type, digits = digits[1])
  cat("\n")
}

# @export
# posterior_draws = function(object, sampler = c("stan","random_walk","independence"), nMC = 10^4,
#                            nMC_burnin = 250){
#   
#   if(length(sampler) > 1){
#     sampler = sampler[1]
#   }
#   if(!(sampler %in% c("stan","random_walk","independence"))){
#     stop("sampler must be 'stan', 'random_walk', or 'independence'")
#   }
#   
#   # Converge group factor into a numeric factor
#   group = as.factor(as.numeric(object$data$group[[1]]))
#   d = nlevels(group)
#   
#   sigma = object$sigma
#   # Specify which random effect variables are non-zero 
#   ranef_idx = which(diag(sigma) > 0)
#   
#   # Adjust Z to incorporate the cholesky composition of the sigma matrix
#   gamma = t(chol(sigma))
#   Z_std = object$data$Z_std
#   Znew2 = Z_std
#   for(i in 1:d){
#     Znew2[group == i,seq(i, ncol(Z), by = d)] = Z[group == i, seq(i, ncol(Z), by = d)] %*% gamma
#   }
#   
#   family = object$family$family
#   # if(family == "gaussian"){
#   #   sig_g =
#   # }
#     
#   # Needed arguments for E_step function
#   # coef, ranef_idx, y, X, Znew2, group, nMC, nMC_burnin, family, link, phi, sig_g,
#   # sampler, d, uold, proposal_SD, batch, batch_length,
#   # offset_increment, trace, num_cores
#   draws = E_step(coef = object$fixef, ranef_idx = ranef_idx, y = object$data$y, X = object$data$X,
#                  Znew2 = Znew2, group = group, nMC = nMC, nMC_burnin = nMC_burnin,
#                  family = family, link = object$family$link,
#                  sampler = sampler, d = d)
# }

#' @export
logLik.pglmmObj = function(object){ 
  
  results_optim = object$results_optim
  ll_elem = which(colnames(results_optim) == "LogLik")
  
  ll = results_optim[ll_elem]
  names(ll) = "logLik"
  
  return(ll)

}


# @method extractBIC pglmmObj
#' @importFrom stringr str_detect
#' @export
extractBIC = function(object){ # extractBIC.pglmmObj
  
  results_optim = object$results_optim
  BIC_elem = which(str_detect(colnames(results_optim), "BIC"))
  BIC_names = colnames(results_optim[,BIC_elem, drop = F])
  BIC_out = results_optim[,BIC_elem]
  names(BIC_out) = BIC_names
  
  return(BIC_out)
  
}

#' @title Plot Diagnostics for MCMC Posterior Draws of the Random Effects
#' 
#' @description Provides graphical diagnostics of the random effect posterior draws from the (best) model.
#' Availabile diagnostics include the sample path, histograms, cummulative sums, and autocorrelation.
#' 
#' @param object an object of class \code{pglmmObj} output from either \code{\link{glmmPen}} 
#' or \code{\link{glmmPen_FineSearch}}.
#' @param plots a character string or a vector of character strings specifying which graphical
#' diagnostics to provide. 
#' @param grps a character string or a vector of character strings specifying which groups 
#' should have diagnostics provided. The names of the groups match the input group factor levels.
#' Default is set to 'all' for all groups.
#' @param vars a character string or a vector of character strings specifying which variables
#' should have diagnostics provided. Tip: can find the names of the random effect variables in
#' the output sigma matrix found in the \code{pglmmObj} object. Default is set to
#' 'all', which picks all variables with non-zero random effects. 
#' @param numeric.grps if TRUE, specifies that the groups factor should be converted to numeric 
#' values. This option could be used to ensure that the organization of the groups is in the 
#' proper numeric order.
#' @param bin_width optional binwidth argument for \code{geom_histogram} from the \code{ggplot2} 
#' package. Default set to \code{NULL}, which specifies the default \code{geom_histogram} binwidth.
#' 
#' @return a list of ggplot graphics, each faceted by group and random effect variable. 
#' Type of plots specified in the \code{plots} argument.
#' 
#' @importFrom reshape2 melt
#' @importFrom stringr str_c str_detect str_sub str_remove str_locate
#' @import ggplot2 
#' @export 
plot_mcmc = function(object, plots = c("all","sample.path","histogram","cumsum","autocorr"),
                     grps = "all", vars = "all", numeric.grps = F, bin_width = NULL){ #plot_mcmc.pglmmObj
  
  if(class(object) != "pglmmObj"){
    stop("'object' must be an object of class pglmmObj output from the glmmPen function")
  }
  
  if(!is.vector(grps) | !is.vector(vars)){
    stop("specified grps and vars must be vectors")
  }
  if("all" %in% plots){
    type = c("sample.path","histogram","cumsum","autocorr")
  }else{
    if(any(!(plots %in% c("sample.path","histogram","cumsum","autocorr")))){
      stop("Specified plot option(s) not recognized. ",
           "Plots allowed: sample.path, histogram, cumsum, and/or autocorr")
    }
    type = plots
  }
  
  U = object$posterior_draws
  U_cols = colnames(U)
  # colnames organization = var_name:grp_name
  
  if(any(vars == "all")){
    d = nlevels(object$data$group[[1]])
    non0 = (diag(object$sigma) != 0)
    non0cols = rep(non0, each = d)
    U_non0 = U[,non0cols]
    if(ncol(U_non0) != ncol(U)){
      vars_non0 = str_sub(colnames(U_non0), start = 1, 
                          end = (str_locate(colnames(U_non0), ":")[,1]-1))
      vars_long = str_remove(str_remove(vars_non0, "[(]"), "[)]")
      vars_string = vars_long[seq(from = 1, by = d, length.out = length(vars_long) / d)]
      # Create reduced U (U_red)
      for(v in vars_string){
        if(vars_string[1] == v){
          U_red = U[,str_detect(U_cols, str_c("[(]?",v,"[)]?",":"))]
        }else{
          U_red = cbind(U_red, U[,str_detect(U_cols, str_c("[(]?",v,"[)]?",":"))])
        }
      }
      U = U_red
      U_cols = colnames(U)
      var_names = vars_string
    }else{
      vars_all = str_sub(U_cols, start = 1, end = (str_locate(U_cols, ":")[,1]-1))
      vars_all2 = str_remove(str_remove(vars_all, "[(]"), "[)]")
      var_names = vars_all2[seq(from = 1, by = d, length.out = length(vars_all2) / d)]
    }
    
  } # End if(any(vars) == "all")
  
  
  
  # If include intercept, convert to recognized form "Intercept"
  if(any(vars != "all")){
    v_int = c("Intercept","intercept","Int","int")
    v_intTRUE = vars %in% v_int[-1]
    vars[v_intTRUE] = "Intercept"
    var_names = vars
  }
  
  if((any(grps == "all")) & (any(vars == "all"))){
    d = nlevels(object$data$group[[1]])
    var_num = ncol(U) / d
    U_keep = U
    grp_names = levels(object$data$group[[1]])
    # var_names = colnames(ranef(object)[[1]])
    
  }else if((any(grps == "all")) & (any(vars != "all"))){
    d = nlevels(object$data$group[[1]])
    var_num = length(vars)
    # Concatenate [)]?, var_name, [)]?, and :
    var_cat = str_c("[(]?",vars,"[)]?",":")
    for(v in var_cat){
      if(v == var_cat[1]){
        U_keep = U[,str_detect(U_cols, v)]
      }else{
        U_keep = cbind(U_keep, U[,str_detect(U_cols, v)])
      }
    }
    grp_names = levels(object$data$group[[1]])
    # var_names = vars
    
  }else if((any(vars == "all")) & (any(grps != "all"))){
    d = length(grps)
    var_num = ncol(U) / nlevels(object$data$group[[1]])
    # Concatenate :, grp_name, and $ to :grp_name$
    grp_cat = str_c(":",grps,"$")
    cols_U_keep0 = NULL
    for(g in grp_cat){
      if(g == grp_cat[1]){
        U_keep0 = U[,str_detect(U_cols, g)]
        cols_U_keep0 = U_cols[str_detect(U_cols, g)]
      }else{
        U_keep0 = cbind(U_keep0, U[,str_detect(U_cols, g)])
        cols_U_keep0 = c(cols_U_keep0, U_cols[str_detect(U_cols, g)])
      }
    }
    U_keep0 = as.matrix(U_keep0)
    colnames(U_keep0) = cols_U_keep0
    grp_names = grps
    # var_names = colnames(ranef(object)[[1]])
    U_keep = matrix(0, nrow = nrow(U_keep0), ncol = ncol(U_keep0))
    cols_U_keep = character(length = ncol(U_keep0))
    for(v in 1:length(var_names)){
      U_keep[,(1+(v-1)*d):(v*d)] = U_keep0[,seq(from = v, by = var_num, length.out = d)]
      cols_U_keep[(1+(v-1)*d):(v*d)] = cols_U_keep0[seq(from = v, by = var_num, length.out = d)]
    }
    colnames(U_keep) = cols_U_keep
    
  }else{ # both vars and grps != "all"
    d = length(grps)
    var_num = length(vars)
    grp_cat = str_c(":",grps,"$")
    for(g in grp_cat){
      if(grp_cat[1] == g){
        U_keep1 = U[,str_detect(U_cols, g)]
      }else{
        U_keep1 = cbind(U_keep1, U[,str_detect(U_cols, g)])
      }
    }
    U_cols_keep1 = colnames(U_keep1)
    var_cat = str_c("[(]?",vars,"[)]?",":")
    for(v in var_cat){
      if(var_cat[1] == v){
        U_keep2 = U_keep1[,str_detect(U_cols_keep1, v)]
      }else{
        U_keep2 = cbind(U_keep2, U_keep1[,str_detect(U_cols_keep1, v)])
      }
    }
    U_keep = U_keep2
    grp_names = grps
    # var_names = vars
  }
  
  if(numeric.grps){
    grp_names = as.numeric(grp_names)
  }
  
  num_plots = d*var_num
  if(num_plots > 100){
    warning("Number plots specified will be > 100. Consider limiting the groups or variables \n",
            immediate. = T)
  }
  
  U_t = data.frame(U_keep, t = 1:nrow(U_keep))
  colnames(U_t) = c(colnames(U_keep), "t")
  U_long = melt(U_t, id = "t")
  U_plot = data.frame(U_long, var_names = rep(var_names, each = d*nrow(U_keep)),
                      grp_names = rep(rep(grp_names, each = nrow(U_keep)), times = var_num))
  
  plots_return = list()
  
  if("sample.path" %in% type){
    plot_sp = ggplot(U_plot, mapping = aes(x = t, y = value)) + geom_path() +
      facet_grid(var_names ~ grp_names) + xlab("iteration t") + ylab("draws")
    
    plots_return$sample_path = plot_sp
  }
  if("histogram" %in% type){
    hist_U = ggplot(U_plot) + geom_histogram(mapping = aes(x = value), binwidth = bin_width) + 
      facet_grid(var_names ~ grp_names) + xlab("draws")
    
    plots_return$histogram = hist_U
  }
  if("cumsum" %in% type){
    U_means = colMeans(U_keep)
    U_means = data.frame(rbind(U_means))[rep.int(1L, nrow(U_keep)), , drop = FALSE]
    U_tmeans = apply(U_keep, 2, cumsum) / 1:nrow(U_keep)
    U_tdiff = U_tmeans - U_means
    U_cumsum = apply(U_tdiff, 2, cumsum)
    U_t = data.frame(U_cumsum, t = 1:nrow(U_cumsum))
    colnames(U_t) = c(colnames(U_keep), "t")
    U_long = melt(U_t, id = "t")
    U_plot = data.frame(U_long, var_names = rep(var_names, each = d*nrow(U_keep)),
                        grp_names = rep(rep(grp_names, each = nrow(U_keep)), times = var_num)) 
    plot_cumsum = ggplot(U_plot) + 
      geom_path(mapping = aes(x = t, y = value), color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_grid(var_names ~ grp_names) + xlab("iteration t") + ylab("Cumulative Sum")
    
    plots_return$cumsum = plot_cumsum
  }
  if("autocorr" %in% type){
    grp_index = rep(grp_names, times = var_num)
    var_index = rep(var_names, each = d)
    for(j in 1:ncol(U_keep)){
      ACF = acf(U_keep[,j], plot=F, lag.max = 40)
      ACF_df = with(ACF, data.frame(lag,acf))
      ACF_df$grp_names = grp_index[j]
      ACF_df$var_names = var_index[j]
      if(j == 1){
        ACF_all = ACF_df
      }else{
        ACF_all = rbind(ACF_all, ACF_df)
      }
    }
    
    plot_acf = ggplot(data = ACF_all, mapping = aes(x = lag, y = acf)) +
      geom_hline(mapping = aes(yintercept = 0)) + 
      geom_segment(mapping = aes(xend = lag, yend = 0)) +
      facet_grid(var_names ~ grp_names)
    
    plots_return$autocorr = plot_acf
  }
  
  return(plots_return)
  
}

#' @importFrom stringr str_c str_to_title
#' @import ggplot2 
#' @method plot pglmmObj
#' @export
plot.pglmmObj = function(object, x = fitted(object, fixed.only = F), 
                         y = switch(object$family$family,
                                    gaussian = residuals(object, type = "pearson"),
                                    residuals(object, type = "deviance")),
                         x_label = "Fitted Values", y_label = NULL){
  
  if(length(x) != length(y)){
    stop("x and y need to be of same length")
  }
  
  x_mat = Matrix::as.matrix(x)
  y_mat = Matrix::as.matrix(y)
  
  if(is.null(y_label)){
    if(!is.null(attr(y, "residual type"))){
      y_label = str_c(str_to_title(attr(y, "residual type")), "Residuals", sep = " ")
    }else{
      y_label = "specify y_label in args"
    }
  }
  data = data.frame(x = x_mat, y = y_mat)
  
  p = ggplot(data = data) + geom_point(mapping = aes(x = x, y = y), color = "blue") +
    ylab(y_label) + xlab(x_label) 
  if(!is.null(attr(y, "residual type"))){
    p = p + geom_hline(yintercept = 0, color = "black")
  }
  return(p)
}


# Functions to add:
## offset() - will not have $offset option
## weights() - will not have $weights option
### Ideas for offset() and weights() - extract from model.frame ?

