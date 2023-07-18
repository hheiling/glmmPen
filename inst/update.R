
# Updating Namespace and package

library(devtools)
library(remotes)

# Download the latest version of the glmmPen package:

devtools::load_all(path = "C:/Users/hheiling/Documents/GitHub/glmmPen")

devtools::document("C:/Users/hheiling/Documents/GitHub/glmmPen")

devtools::install(pkg = "C:/Users/hheiling/Documents/GitHub/glmmPen")


devtools::check(pkg = "C:/Users/hheiling/Documents/GitHub/glmmPen")

# Build package
devtools::build(pkg = "C:/Users/hheiling/Documents/GitHub/glmmPen",
                path = "C:/Users/hheiling/Documents/Longleaf")
# devtools::build(pkg = "C:/Users/hheiling/Documents/GitHub/glmmPen",
#                 path = "C:/Users/hheiling/Documents/GitHub/CRAN")
# devtools::build(pkg = "C:/Users/hheiling/Documents/GitHub/glmmPen",
#                 path = "C:/Users/hheiling/Documents/LaTeX/glmmPen JSS paper")

# Build manual
devtools::build_manual(pkg = "C:/Users/hheiling/Documents/GitHub/glmmPen", 
                       path = "C:/Users/hheiling/Documents/GitHub/glmmPen/inst")

# Update stan model files for E-step sampling
library(rstantools)
# use_rstan(pkgdir = "C:/Users/hheiling/Documents/GitHub/glmmPen")
# Creates or update package-specific system files to compile .stan model files found in inst/stan.
rstantools::rstan_config(pkgdir = "C:/Users/hheiling/Documents/GitHub/glmmPen")

# Can add examples and specify that they do not run or are not checked, see https://stuff.mit.edu/afs/athena/software/r/current/RStudio/resources/roxygen_help.html
# 
# \donotrun{}, \donttest{}


# Stan: look-up available functions
# lookup("distribution")

# Versions:
## 1.5.2.6: Added thresholding checks for Poisson family variable selection
## 1.5.2.5: Removed Z_standardized option, added standardization = T/F option (X and Z
##    either both standardized or both non-standardized)
## 1.5.2.4: Added alternative convergence options
## 1.5.2.3: Implemented proximal line search algorithm

# other branch
# install_github("hheiling/glmmPen", ref = "alt_branch", force = TRUE)
# master branch
# install_github("hheiling/glmmPen", force = TRUE)
# M_Step branch
# install_github("hheiling/glmmPen", ref = "M_Step", force = TRUE)

# unlink("C:/Users/hheiling/Documents/R/win-library/3.6/00LOCK-glmmPen", recursive = TRUE)

# Update Namespace and man files:

# Package working directory: "C:/Users/hheiling/Documents/GitHub/glmmPen"
# Run document line when new @export functions added to package
# pkgbuild::compile_dll("C:/Users/hheiling/Documents/GitHub/glmmPen")
# devtools::document("C:/Users/hheiling/Documents/GitHub/glmmPen")
# roxygen2::roxygenize("C:/Users/hheiling/Documents/GitHub/glmmPen")




# Tips for building vignette:
# https://r-pkgs.org/vignettes.html
# usethis::use_vignette("glmmPen")
# https://kbroman.org/pkg_primer/pages/vignettes.html
devtools::build_vignettes()

# data("basal")
# dim(basal$X)
# names(basal)
# basal$group = factor(basal$group, 
#                      labels = c("UNC_PDAC","TCGA_PDAC","TCGA_Bladder","UNC_Breast"))
# save(basal, file = "data/basal.RData")

# Rcpp.package.skeleton(name = "C:/Users/hheiling/Documents/GitHub/Rcpp_pkg")
###########################################################################################

# Errors and fixes

# Error: (replace "stan_fit4lm_mod" with appropriate model)
## Error in Rcpp::loadModule(module = “stan_fit4lm_mod”, what = TRUE, env = ns, : *
## Unable to load module “stan_fit4lm_mod”: Failed to initialize module pointer: Error in FUN(X[[i]], …): no such symbol _rcpp_module_boot_stan_fit4lm_mod in package rstanlmdemo*
# Fix:
## Add back @useDynLib glmmPen to fit_dat_B documentation
## Add the following line manually to the NAMESPACE file: 
## useDynLib(glmmPen, .registration = TRUE)

# Error (R 4.1.1):
## Error in Rcpp::loadModule(module = "stan_fit4binomial_logit_model_mod",  : 
##    Unable to load module "stan_fit4binomial_logit_model_mod": cannot allocate vector of size 9840.8 Gb
# Fix:
## Add line to Documents.R/Makevars.win (hit "Enter" at end of line; without, didn't
##    register line end)
# Error (R 4.1.1):
## Errors of ".0 file format not recognized" ... 
# Fix:
## See "stan_commands.R" for Version 4.1.1
# Extras: 
# write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)
# writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
###########################################################################################


###########################################################################################
###########################################################################################
# @export
# extractBICq = function(object, BICq_posterior = NULL,
#                        tuning_options = NULL){
#   
#   ###################################################################################
#   # Extraction of data/arguments from pglmmObj
#   # Input argument checks
#   ###################################################################################
#   
#   
#   ###################################################################################
#   # Fit 'full' model, save posterior samples for BIC-ICQ calculation
#   ###################################################################################
#   
#   
#   ###################################################################################
#   # BIC-ICQ calculation
#   ###################################################################################
#   
# }



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

# @importFrom ncvreg std
# @export
# pglmmObj_mod = function(object, trace = 0){
#   
#   X = object$data$X
#   X_std = cbind(1, std(X[,-1, drop = F]))
#   
#   y = object$data$y
#   Z = Matrix::as.matrix(object$data$Z_std)
#   group = object$data$group
#   group_num = as.factor(as.numeric(group[[1]]))
#   
#   dat = list(y = y, X = X_std, Z = Z, group = group_num)
#   
#   fam_fun = object$family
#   
#   fit = c(object$Estep_material, list(J = object$J, sigma = object$sigma, 
#              sigma_gaus = ifelse(fam_fun == "gaussian", sqrt(object$scale$Gaus_sig2), 1.0),
#              object$Gibbs_info))
#   
#   Estep_out = E_step_final(dat = dat, fit = fit, optim_options = object$control_info$optim_options, 
#                            fam_fun = fam_fun, extra_calc = T, 
#                            adapt_RW_options = object$control_info$adapt_RW_options, trace = trace)
#   
#   
#   object$posterior_draws = Estep_out$post_out
#   
#   # Random effects coefficients
#   rand = Estep_out$post_modes
#   d = nlevels(group_num)
#   q = ncol(Z) /d
#   
#   ## Organization of rand: Var1 group levels 1, 2, ... Var2 group levels 1, 2, ...
#   ref = as.data.frame(matrix(rand, nrow = d, ncol = q, byrow = F) )
#   rownames(ref) = rownames(object$ranef[[1]])
#   colnames(ref) = colnames(object$ranef[[1]])
#   ranef = lapply(group, function(j) ref)
#   object$ranef = ranef
#   
#   optim_results = object$results_optim
#   idx = which(colnames(optim_results) %in% c("BICh","BIC","BICNgrp","LogLik"))
#   optim_results[,idx] = c(Estep_out$BICh, Estep_out$BIC, Estep_out$BICNgrp, Estep_out$ll)
#   object$results_optim = optim_results
#   
#   return(object)
#   
# }
