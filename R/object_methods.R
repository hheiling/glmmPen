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
                            fixed.only = F, na.action = na.pass){
  ## Other arguments used by lme4: re.form, random.only = F, allow.new.levels = F, newparams = NULL
  
  if(!is.null(newdata) && class(newdata) != "data.frame"){
    stop("newdata must be a dataframe")
  }
  
  if(is.null(newdata)){
    if(!fixed.only){
      pred = switch(type[1], # if unspecified, default = link output (linear predictor)
                    link = etaCalc(X = object$X, Z = object$Z, beta = fixef(object), 
                                   U = object$gibbs_mcmc),
                    response = fitted(object))
    }else{ # Fixed.only = T
      eta = object$X %*% fixef(object)
      pred = switch(type[1],
                    link = eta,
                    response = invLink(family = family(object), eta = eta))
    }
    
    data = object$X
    
  }else{
    
    fD_out = formulaData(formula = formula(object), data = newdata, na.action = na.action)
    if(!fixed.only){
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
      pred = switch(type[1], 
                    link = eta,
                    response = invLink(family = family(object), eta = eta))
      
    }else{ # fixed.only = T
      eta = fD_out$X %*% fixef(object)
      pred = switch(type[1], 
                    link = eta,
                    response = invLink(family = family(object), eta = eta))
    }
    
    data = newdata
    
    # In future, create code for option else(!is.null(re.form))
  }
  
  if(class(pred) %in% c("dgeMatrix","matrix")){
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

var_hat = function(family, mu, sig2 = NULL){
  if(family == "binomial"){
    v_hat = mu*(1-mu)
  }else if(family == "poisson"){
    v_hat = mu
  }else if(family == "gaussian"){
    v_hat = sig2
  }
  return(v_hat)
}

ll = function(family, mu, y, v = NULL){
  if(family == "binomial"){
    llik = dbinom(y, size = 1, prob = mu, log = T)
  }else if(family == "poisson"){
    llik = dpois(y, lambda = mu, log = T)
  }else if(family == "gaussian"){
    llik = dnorm(y, mean = mu, sd = sqrt(v), log = T)
  }
  return(llik)
}

#' @importFrom stats predict
#' @export
residuals.pglmmObj = function(object, type = c("deviance","pearson","response","working")){
  # What is working response?
  Y = object$y
  mu = Matrix::as.matrix(fitted(object))
  type = match.arg(type)
  
  if(type == "deviance"){
    res = sign(Y - mu)*sqrt(-2*(ll(family(object), mu, Y)))
  }else if(type == "pearson"){
    res = (Y - mu) / sqrt(var_hat(family(object), mu))
  }else if(type == "response"){
    res = Y - mu
  }else{
    res = (Y - mu) / mu
  }
  
  attr(res, "residual type") = type
  return(res)
}

################################################################### 
# Print functions:

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
  cat(sprintf("Number Observations: %i,  groups: %s, %i \n", 
              nobs(object), names(object$group), nlevels(object$group[[1]])))
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
  group = object$group[[1]] 
  group_name = names(object$group)
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
summary.pglmmObj = function(object, digits = c(fef = 4, ref = 4), resid_type = "deviance"){
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

#' @importFrom reshape2 melt
#' @importFrom stringr str_c str_detect str_sub str_remove str_locate
# @method plot_mcmc pglmmObj
#' @export 
plot_mcmc = function(object, plots = c("all","sample.path","histogram","cumsum","autocorr"),
                     grps = "all", vars = "all", numeric.grps = F){ #plot_mcmc.pglmmObj
  # ToDo: Remove cols from U associated with vars with zero variance?
  
  if(object$sampling != "Gibbs Sampling"){
    stop("The plots in plot_mcmc are only relevant when Gibbs sampling is used; \n
         Rejection sampling was used in this case")
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
  
  U = object$gibbs_mcmc
  U_cols = colnames(U)
  # colnames organization = var_name:grp_name
  
  if(vars == "all"){
    d = nlevels(object$group[[1]])
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
    }
  }
  
  # If include intercept, convert to recognized form "Intercept"
  if(vars != "all"){
    v_int = c("Intercept","intercept","Int","int")
    v_intTRUE = vars %in% v_int[-1]
    vars[v_intTRUE] = "Intercept"
  }
  
  if(grps == "all" && vars == "all"){
    d = nlevels(object$group[[1]])
    var_num = ncol(U) / d
    U_keep = U
    grp_names = levels(object$group[[1]])
    var_names = colnames(ranef(object)[[1]])
    
  }else if(grps == "all" && vars != "all"){
    d = nlevels(object$group[[1]])
    var_num = length(vars)
    # Concatenate [)]?, var_name, [)]?, and :
    var_cat = str_c("[(]?",vars,"[)]?",":")
    for(v in var_cat){
      if(var_cat[1] == v){
        U_keep = U[,str_detect(U_cols, v)]
      }else{
        U_keep = cbind(U_keep, U[,str_detect(U_cols, v)])
      }
    }
    grp_names = levels(object$group[[1]])
    var_names = vars
    
  }else if(vars == "all" && grps != "all"){
    d = length(grps)
    var_num = ncol(U) / nlevels(object$group[[1]])
    # Concatenate :, grp_name, and $ to :grp_name$
    grp_cat = str_c(":",grps,"$")
    for(g in grp_cat){
      if(grp_cat[1] == g){
        U_keep0 = U[,str_detect(U_cols, g)]
      }else{
        U_keep0 = cbind(U_keep0, U[,str_detect(U_cols, g)])
      }
    }
    grp_names = grps
    var_names = colnames(ranef(object)[[1]])
    U_keep = matrix(0, nrow = nrow(U_keep0), ncol = ncol(U_keep0))
    cols_U_keep = numeric(length = ncol(U_keep0))
    for(v in 1:length(var_names)){
      U_keep[,(1+(v-1)*d):(v*d)] = U_keep0[,seq(from = v, by = var_num, length.out = d)]
      cols_U_keep[(1+(v-1)*d):(v*d)] = colnames(U_keep0)[seq(from = v, by = var_num, length.out = d)]
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
    var_names = vars
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
    hist_U = ggplot(U_plot) + geom_histogram(mapping = aes(x = value)) + 
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
    plot_cumsum = ggplot(U_plot) + geom_smooth(mapping = aes(x = t, y = value), color = "black") +
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

#' @method plot pglmmObj
#' @export
plot.pglmmObj = function(object, x = fitted(object), y = residuals(object, type = "deviance"),
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

