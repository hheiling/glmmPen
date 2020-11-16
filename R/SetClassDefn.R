# Setting class definitions

#' @export
pglmmObj = setRefClass("pglmmObj",
            fields = list(
              fixef = "numeric", # fixed effects coefficients
              ranef = "list", # random effects coefficients
              sigma = "matrix", 
                # sigma: Turn in to "dgCMatrix"? - sparse for large dimensions matrix
                # alternatively, if large dims, assume diagonal, could report matrix with ncol = 1
              posterior_draws = "matrix",
              sampling = "character",
              results_all = "matrix",
              results_optim = "matrix",
              family = "family",
              J = "dgCMatrix", # sparse matrix
              penalty_info = "list",
              nonzeroFE = "numeric", # non-zero fixed effects coefficients (colnames)
              # nonzerovarRE = "character", # coefficients with non-zero variance (colnames)
              call = "call",
              formula = "formula",
              data = "list",
              optinfo = "list",
              control_info = "list",
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
                  center = x$std_out$X_center
                  scale = x$std_out$X_scale
                  beta[1] = beta[1] - sum(center*beta[-1]/scale)
                  beta[-1] = beta[-1] / scale
                names(beta) = x$coef_names$fixed
                fixef <<- beta
                nonzeroFE <<- fixef[which(fixef != 0)]
                
                # Covariance matrix of random effects
                sigma <<- x$sigma
                colnames(sigma) <<- x$coef_names$random
                rownames(sigma) <<- x$coef_names$random
                
                # Return MCMC results - potentially for MCMC diagnostics if desired
                  U = x$u
                  Gam = matrix(x$J%*%matrix(x$coef[-c(1:ncol(X))], ncol = 1), ncol = ncol(Z_std)/d)
                  q = ncol(Z_std) / d # Number random variables
                  post_U = matrix(0, nrow = nrow(U), ncol = ncol(U))
                  for(g in 1:d){
                    idx = seq(from = g, to = (d*q), by = d)
                    post_U[,idx] = U[,idx] %*% t(Gam)
                  }
                posterior_draws <<- post_U
                colnames(posterior_draws) <<- colnames(Z_std)
                
                # Random effects coefficients
                rand = colMeans(posterior_draws) 
                
                ## Organization of rand: Var1 group levels 1, 2, ... Var2 group levels 1, 2, ...
                ref = as.data.frame(matrix(rand, nrow = d, ncol = q, byrow = F) )
                rownames(ref) = levs
                colnames(ref) = x$coef_names$random
                ranef <<- lapply(group, function(j) ref)
                
                family <<- x$family
                J <<- x$J
                formula <<- x$formula
                call <<- x$call
                sampling <<- x$sampling
                results_all <<- x$selection_results
                results_optim <<- x$optim_results
                
                penalty_info <<- list(penalty = x$penalty, gamma_penalty = x$gamma_penalty, 
                                      alpha = x$alpha, fixef_noPen = x$fixef_noPen)
                optinfo <<- list(control_options = x$control_options$optim_options,
                                 warnings = x$warnings, iter = x$EM_iter, conv_val = x$EM_conv)
                control_info <<- x$control_options
                Gibbs_info <<- list(gibbs_accept_rate = x$gibbs_accept_rate, proposal_SD = x$proposal_SD)
                
              },
              show = function(){
                print(.self)
              }
            ))