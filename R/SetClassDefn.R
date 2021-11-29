# Setting class definitions
# Example from merMod in lme4: https://github.com/lme4/lme4/blob/master/R/AllClass.R

#' @title Class \code{pglmmObj} of Fitted Penalized Generalized Mixed-Effects Models for 
#' package \code{glmmPen}
#' 
#' @description The functions \code{\link{glmm}}, \code{\link{glmmPen}}, and 
#' \code{\link{glmmPen_FineSearch}} from the package \code{glmmPen} 
#' output the reference class object of type \code{pglmmObj}.
#' 
#' @name pglmmObj-class 
#' @aliases pglmmObj-class pglmmObj
#' show, pglmmObj-method, coef.pglmmObj, fitted.pglmmObj, fixef.pglmmObj, formula.pglmmObj,
#' logLik.pglmmObj, model.frame.pglmmObj, model.matrix.pglmmObj, plot.pglmmObj, 
#' predict.pglmmObj, print.pglmmObj, ranef.pglmmObj, residuals.pglmmObj,
#' sigma.pglmmObj, summary.pglmmObj
#' 
#' @docType class
#' 
#' @keywords classes
#' 
#' @return The pglmmObj object returns the following items:
#' \item{fixef}{vector of fixed effects coefficients}
#' \item{ranef}{matrix of random effects coefficients for each 
#' explanatory variable for each level of the grouping factor}
#' \item{sigma}{random effects covariance matrix}
#' \item{scale}{if family is Gaussian, returns the residual error variance}
#' \item{posterior_samples}{Samples from the posterior distribution of the random effects,
#' taken at the end of the model fit (after convergence or after maximum iterations allowed).
#' Can be used for diagnositics purposes. Note: These posterior samples are from a single chain.}
#' \item{sampling}{character string for type of sampling used to calculate the posterior samples
#' in the E-step of the algorithm}
#' \item{results_all}{matrix of results from all model fits during variable selection (if selection
#' performed). Output for each model includes: penalty parameters for fixed (lambda0) and random 
#' (lambda1) effects, BIC-derived quantities and the log-likelihood 
#' (note: the arguments \code{BIC_option} and \code{logLik_calc} in \code{\link{selectControl}}
#' determine which of these quantities are calculated for each model), 
#' the number of non-zero fixed and random effects (includes intercept),
#' number of EM iterations used for model fit, whether or not the 
#' model converged (0 for no vs 1 for yes), and the fixed and random effects coefficients}
#' \item{results_optim}{results from the 'best' model fit; see results_all for details. 
#' BICh, BIC, BICNgrp, and LogLik computed for this best model if not previously calculated.}
#' \item{family}{Family}
#' \item{penalty_info}{list of penalty information}
#' \item{call}{arguments plugged into \code{glmm}, \code{glmmPen}, or \code{glmmPen_FineSearch}}
#' \item{formula}{formula}
#' \item{fixed_vars}{names of fixed effects variables}
#' \item{data}{list of data used in model fit, including the response y, the fixed effects 
#' covariates matrix X, the random effects model matrix Z (which is composed of values from the 
#' standardized fixed effects model matrix), 
#' the grouping factor, offset, model frame,
#' and standarization information used to standardize the fixed effects covariates}
#' \item{optinfo}{Information about the optimization of the 'best' model}
#' \item{control_info}{optimization parameters used for the model fit}
#' \item{Estep_init}{materials that can be used to initialize another E-step (for 
#' use in \code{glmmPen_FineSearch})}
#' \item{Gibbs_info}{list of materials to perform diagnositics on the Metropolis-within-Gibbs
#' sample chains, including the Gibbs acceptance rates (included for both the independence
#' and adaptive random walk samplers) and the final proposal standard deviations 
#' (included for the adaptive random walk sampler only))}
#' 
#' showClass("pglmmObj")
#' methods(class = "pglmmObj")
#' 
#' @importFrom stringr str_c
#' @export
pglmmObj = setRefClass("pglmmObj",
            fields = list(
              fixef = "numeric", # fixed effects coefficients
              ranef = "list", # random effects coefficients
              sigma = "matrix", 
              scale = "list",
              posterior_samples = "matrix",
              sampling = "character",
              results_all = "matrix",
              results_optim = "matrix",
              family = "family",
              penalty_info = "list",
              call = "call",
              formula = "formula",
              fixed_vars = "character",
              data = "list",
              optinfo = "list",
              control_info = "list",
              Estep_init = "list",
              Gibbs_info = "list"
            ),
            methods = list(
              initialize = function(x){ # x = input list object
                # ToDo: make sure items are formatted properly and have proper names associated with them
                ## Will do above when finalize fit_dat output and put finishing touches on input checks
                
                # Group
                group = x$group
                ## For now, assume only one group designation allowed
                # d = lapply(group, function (j) nlevels(j))
                # levs = lapply(group, function(j) levels(j))
                d = nlevels(group[[1]])
                levs = levels(group[[1]])
                
                # y, X, Z, frame
                
                y = x$y
                X = x$X
                Z_std = x$Z_std
                  rand_vars = rep(x$coef_names$random, each = d)
                  grp_levs = rep(levs, times = length(x$coef_names$random))
                  colnames(Z_std) = noquote(paste(rand_vars, grp_levs, sep = ":"))
                  rownames(Z_std) = rownames(X)
                frame = x$frame
                offset = x$offset
                std_info = list(X_center = x$std_out$X_center,
                                  X_scale = x$std_out$X_scale)
                data <<- list(y = y, X = X, Z_std = Z_std, group = group, 
                              offset = offset, frame = frame, std_info = std_info)
                
                # Fixed effects coefficients - need to unstandardize
                
                  p = ncol(X)
                  beta = x$coef[1:p]
                  center = std_info$X_center
                  scale_std = std_info$X_scale
                  beta[1] = beta[1] - sum(center*beta[-1]/scale_std)
                  beta[-1] = beta[-1] / scale_std
                names(beta) = x$coef_names$fixed
                fixef <<- beta
                
                # Covariance matrix of random effects
                sigma <<- x$sigma
                colnames(sigma) <<- x$coef_names$random
                rownames(sigma) <<- x$coef_names$random
                
                # coef <<- x$coef
                # names(coef) = c(x$coef_names$fixed, str_c("Gamma",0:length(coef[-c(1:p)])))
                
                # Return MCMC results - potentially for MCMC diagnostics if desired
                posterior_samples <<- x$Estep_out$post_out
                colnames(posterior_samples) <<- colnames(Z_std)
                
                # Random effects coefficients
                rand = x$Estep_out$post_modes
                q = ncol(Z_std) / d
                
                ## Organization of rand: Var1 group levels 1, 2, ... Var2 group levels 1, 2, ...
                ref = as.data.frame(matrix(rand, nrow = d, ncol = q, byrow = F) )
                rownames(ref) = levs
                colnames(ref) = x$coef_names$random
                ranef <<- lapply(group, function(j) ref)
                
                family <<- x$family
                # If Gaussian family, return gaussian residual error variance estimate
                if(family$family == "gaussian"){
                  scale <<- list(Gaus_sig2 = x$sigma_gaus^2)
                }else if(family$family == "negbin"){
                  scale <<- list(phi = x$phi)
                }else{
                  scale <<- list(scale = NULL)
                }
                
                formula <<- x$formula
                fixed_vars <<- x$fixed_vars
                call <<- x$call
                sampling <<- x$sampling
                results_all <<- x$selection_results
                results_optim <<- x$optim_results
                
                if(nrow(results_all) == 1){ # glmm, not glmmPen
                  prescreen_ranef = NULL
                }else{
                  prescreen_ranef = x$ranef_keep
                  names(prescreen_ranef) = x$coef_names$random
                }
                
                penalty_info <<- list(penalty = x$penalty, gamma_penalty = x$gamma_penalty, 
                                      alpha = x$alpha, fixef_noPen = x$fixef_noPen,
                                      prescreen_ranef = prescreen_ranef)
                optinfo <<- list(iter = x$EM_iter, conv = x$EM_conv, warnings = x$warnings,
                                 control_options = x$control_options$optim_options)
                control_info <<- x$control_options
                Estep_init <<- list(u_init = x$Estep_out$u_init, coef = x$coef,
                                        updated_batch = x$updated_batch)
                Gibbs_info <<- list(gibbs_accept_rate = x$gibbs_accept_rate, proposal_SD = x$proposal_SD)
                
              },
              show = function(){
                print(.self)
              }
            ))