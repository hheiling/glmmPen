# Setting class definitions

##' @export
pglmmObj = setRefClass("pglmmObj",
            fields = list(
              fixef = "numeric", # fixed effects coefficients
              ranef = "data.frame", # random effects coefficients
              # coefficients = "numeric", # sum of fixed and random effects coefficients
              group = "factor",
              covar = "matrix",
              gibbs_mcmc = "matrix",
              family = "character",
              J = "dgCMatrix", # sparse matrix
              # BIC_ICQ = "numeric",
              # penalty = "list", # type, chosen optimal values for fixed and random effects
              nonzeroFE = "numeric", # non-zero fixed effects coefficients (colnames)
              # nonzerovarRE = "character", # coefficients with non-zero variance (colnames)
              call = "call",
              formula = "formula",
              y = "numeric",
              X = "matrix",
              Z = "matrix"
            ),
            methods = list(
              initialize = function(x){ # x = input list object
                # ToDo: make sure items are formatted properly and have proper names associated with them
                ## Will do above when finalize fit_dat output and put finishing touches on input checks
                
                # Group
                group <<- x$group
                d = nlevels(group)
                
                # y, X, Z
                y <<- x$Y
                X <<- x$X
                Z <<- x$Z
                
                # Fixed effects coefficients
                p = ncol(X)
                fixef <<- x$coef[1:p]
                names(fixef) <<- x$coef_names$fixed
                nonzeroFE <<- fixef[which(fixef != 0)]
                
                q = length(x$coef_names$random) # Number random effects
                # Note: q > 0 always; error in glmmPen function if no random effect specified
                Gam_vec = x$coef[-c(1:p)]
                # if(Gamma_diag){ # If, in high dimensional case, Gamma assumed diagonal
                #   Gamma = diag(Gam_vec)
                # }else{ # In low dimensional case, have lower-triangular Gamma
                #   Gamma = matrix(0, nrow = q, ncol = q)
                #   Gamma[lower.tri(Gamma, diag=T)] = Gam_vec
                # }
                Gamma = matrix(0, nrow = q, ncol = q)
                Gamma[lower.tri(Gamma, diag=T)] = Gam_vec
                # Covariance matrix of random effects: from Rashid and Li paper, Gamma * t(Gamma)
                covar <<- Gamma %*% t(Gamma)
                # nonzervarRE <<- sum(rowSums(covar) != 0)
                # Add names of random effects to matrix (either rownames or colnames)
                
                # Return MCMC results - potentially for MCMC diagnostics if desired
                gibbs_mcmc <<- x$u 
                # colnames(gibbs_mcmc) = 
                
                # Random effects
                rand = colMeans(gibbs_mcmc)
                ## Organization of rand: Var1 group levels 1, 2, ... Var2 group levels 1, 2, ...
                ranef <<- as.data.frame(matrix(rand, nrow = d, ncol = q, byrow = F) )
                rownames(ranef) <<- levels(x$group)
                colnames(ranef) <<- x$coef_names$random
                
                family <<- x$family
                J <<- x$J
                formula <<- x$formula
                call <<- x$call
                # BIC_ICQ <<- x$BIC
                # penalty <<- x$penalty
                
              }
            ))