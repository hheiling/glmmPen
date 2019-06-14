# Setting class definitions

#' @export
pglmmObj = setRefClass("pglmmObj",
            fields = list(
              fixef = "numeric", # fixed effects coefficients
              ranef = "list", # random effects coefficients
              # coefficients = "numeric", # sum of fixed and random effects coefficients
              group = "list",
              sigma = "matrix", 
                # sigma: Turn in to "dgCMatrix"? - sparse for large dimensions matrix
                # alternatively, if large dims, assume diagonal, report matrix with ncol = 1
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
              Z = "dgCMatrix",
              frame = "data.frame",
              sampling = "character"
            ),
            methods = list(
              initialize = function(x){ # x = input list object
                # ToDo: make sure items are formatted properly and have proper names associated with them
                ## Will do above when finalize fit_dat output and put finishing touches on input checks
                
                # Group
                group <<- x$group
                ## For now, assume only one group designation allowed
                # d = lapply(group, function (j) nlevels(j))
                # levs = lapply(group, function(j) levels(j))
                d = nlevels(group[[1]])
                levs = levels(group[[1]])
                
                # y, X, Z, frame
                y <<- x$Y
                X <<- x$X
                Z <<- x$Z
                  rand_vars = rep(x$coef_names$random, each = d)
                  grp_levs = rep(levs, times = length(x$coef_names$random))
                  colnames(Z) <<- noquote(paste(rand_vars, grp_levs, sep = ":"))
                  rownames(Z) <<- rownames(X)
                frame <<- x$frame
                
                # Fixed effects coefficients
                p = ncol(X)
                fixef <<- x$coef[1:p]
                names(fixef) <<- x$coef_names$fixed
                nonzeroFE <<- fixef[which(fixef != 0)]
                
                # Covariance matrix of random effects
                sigma <<- x$sigma
                # nonzervarRE <<- sum(rowSums(covar) != 0)
                # Add names of random effects to matrix (either rownames or colnames)
                
                # Return MCMC results - potentially for MCMC diagnostics if desired
                gibbs_mcmc <<- x$u 
                colnames(gibbs_mcmc) <<- colnames(Z)
                
                # Random effects coefficients
                rand = colMeans(gibbs_mcmc) 
                q = ncol(Z) / d # Number random variables
                ## Organization of rand: Var1 group levels 1, 2, ... Var2 group levels 1, 2, ...
                ref = as.data.frame(matrix(rand, nrow = d, ncol = q, byrow = F) )
                rownames(ref) = levs
                colnames(ref) = x$coef_names$random
                ranef <<- lapply(group, function(j) ref)
                
                family <<- x$family
                J <<- x$J
                formula <<- x$formula
                call <<- x$call
                # BIC_ICQ <<- x$BIC
                # penalty <<- x$penalty
                sampling <<- x$sampling
              },
              show = function(){
                print(.self)
              }
            ))