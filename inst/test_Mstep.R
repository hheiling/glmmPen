# Test M-step IRLS and regular coordinate descent

library(ncvreg)
library(glmmPen)

########################################################################################################

# Create dataset to test IRLS (no random effects)
set.seed(2020)
N = 100
X = matrix(rnorm(5*N), nrow = N, ncol = 5)
X = cbind(1, X)
beta = rep(1, times = ncol(X))

eta = X %*% beta
y = rbinom(n = N, size = 1, prob = exp(eta)/(1+exp(eta)))

# Test IRLS

## The gold standard result 
fit1 = glm(formula = y ~ X[,-1], family = binomial())

int_only = glm(formula = y ~ 1, family = binomial())

## Matrix inversion IRLS (can only use for low-dimensional cases)
fit2 = M_step(y, X, family = binomial(),
              coef = c(int_only$coefficients, rep(0, times = ncol(X)-1)),
              maxit = 1000, conv = 0.0001, fit_type = 1, lambda = 0)

## Coordinate descent with lambda == 0 (coord descent version of IRLS)
fit3 = M_step(y, X, family = binomial(),
              coef = c(int_only$coefficients, rep(0, times = ncol(X)-1)),
              penalty = "lasso", maxit_CD = 1000, conv = 0.0001, fit_type = 2, lambda = 0)

## Grouped coordinate descent with all gropus comprised on one covariate, lambda == 0
fit4 = M_step(y, X, family = binomial(),
              coef = c(int_only$coefficients, rep(0, times = ncol(X)-1)),
              penalty = "lasso", maxit_CD = 1000, conv = 0.0001, fit_type = 3, lambda = 0)

cbind(fit1$coefficients, fit2$coef, fit3$coef, fit4$coef)

########################################################################################################

# Create dataset to test Coordinate Descent (no random effects)

set.seed(2020)
N = 100
X = matrix(rnorm(5*N), nrow = N, ncol = 10)
## Standardize the X matrix
X_std = std(X)
X_toUse = cbind(1, X_std)
beta = c(0, 2, 2, rep(0, times = 8))

eta = X_toUse %*% beta
y = rbinom(n = N, size = 1, prob = exp(eta)/(1+exp(eta)))

# Test coordinate descent

fitA = ncvreg(X_std, y, family = "binomial", penalty = "MCP", nlambda = 100,
              max.iter = 1000, alpha = 1.0)

## Fit intercept-only model
int_only = glm(formula = y ~ 1, family = binomial())
## Use same lambda values as ncvreg
lambda_vec = fitA$lambda

beta = matrix(NA, nrow = 11, ncol = length(lambda_vec)) # Same set-up as ncvreg$beta

for(i in 1:length(lambda_vec)){
  if(i == 1){
    beta_init = c(int_only$coefficients, rep(0, times = ncol(X_toUse)-1))
  }else{
    beta_init = beta[,i-1]
  }
  
  fitB = M_step(y, X_toUse, family = binomial(),
                coef = beta_init, penalty = "MCP", alpha = 1.0,
                maxit_CD = 1000, conv = 0.0001, fit_type = 2, lambda = lambda_vec[i])
  beta[,i] = fitB$coef
}

colnames(beta) = round(lambda_vec, 5)

# Side-by-side comparison of beta output
# for(i in 1:5){
#   print(cbind(beta[,i], fitA$beta[,i]))
# }


beta2 = matrix(NA, nrow = 11, ncol = length(lambda_vec)) # Same set-up as ncvreg$beta

for(i in 1:length(lambda_vec)){
  if(i == 1){
    beta_init = c(int_only$coefficients, rep(0, times = ncol(X_toUse)-1))
  }else{
    beta_init = beta2[,i-1]
  }
  
  fitC = M_step(y, X_toUse, family = binomial(),
                coef = beta_init, penalty = "MCP", alpha = 1.0,
                maxit_CD = 1000, conv = 0.0001, fit_type = 3, lambda = lambda_vec[i])
  beta2[,i] = fitC$coef
}

colnames(beta) = round(lambda_vec, 5)

# Side-by-side comparison of beta output
for(i in 1:5){
  print(cbind(beta[,i], beta2[,i], fitA$beta[,i]))
}


########################################################################################################

# Create dataset to test Grouped Coordinate Descent 
## (no random effects, linear predictor = X*Beta only)
library(mvtnorm)
library(ncvreg)

set.seed(2020)
N = 100
x1 = std(rmvnorm(N, sigma = matrix(c(1,0.7,0.7,1), nrow = 2)))
x1_SVD = svd(x1)
x2 = std(rmvnorm(N, sigma = matrix(c(1,-0.5,-0.5,1), nrow = 2, byrow = T)))
x2_SVD = svd(x2)
x3 = std(rmvnorm(N, sigma = matrix(c(1,0.5,0.1,0.5,1,0.5,0.1,0.5,1), nrow = 3, byrow = T)))
x3_SVD = svd(x3)

# Use orthonormal groups
X = cbind(x1%*%x1_SVD$v%*%diag(sqrt(N)/x1_SVD$d), 
          x2%*%x2_SVD$v%*%diag(sqrt(N)/x2_SVD$d), 
          x3%*%x3_SVD$v%*%diag(sqrt(N)/x3_SVD$d))

# Check: The following should all equal the Identity matrix
t(X[,1:2])%*%X[,1:2] / N
t(X[,3:4])%*%X[,3:4] / N
t(X[,5:7])%*%X[,5:7] / N

X_toUse = cbind(1, X)
group_X = c(0, rep(1,2), rep(2,2), rep(3,3))

beta = c(0, 2, 2, 1, 1, rep(0, times = 3))

eta = X_toUse %*% beta
y = rbinom(n = N, size = 1, prob = exp(eta)/(1+exp(eta)))

# Test coordinate descent

fit_grpA = grpreg(X, y, group = group_X[-1], family = "binomial", penalty = "grLasso", nlambda = 20,
                  max.iter = 1000, alpha = 1.0)

## Fit intercept-only model
int_only = glm(formula = y ~ 1, family = binomial())
## Use same lambda values as ncvreg
lambda_vec = fit_grpA$lambda

beta = matrix(NA, nrow = ncol(X_toUse), ncol = length(lambda_vec)) # Same set-up as ncvreg$beta

for(i in 1:length(lambda_vec)){
  if(i == 1){
    beta_init = c(int_only$coefficients, rep(0, times = ncol(X_toUse)-1))
  }else{
    beta_init = beta[,i-1]
  }
  
  fit_grpB = M_step(y, X_toUse, group_X = group_X, family = binomial(),
                    coef = beta_init, penalty = "lasso", alpha = 1.0,
                    maxit_CD = 1000, conv = 0.0001, fit_type = 3, lambda = lambda_vec[i])
  beta[,i] = fit_grpB$coef
}

colnames(beta) = round(lambda_vec, 5)

# Side-by-side comparison of beta output
for(i in c(1:5, length(lambda_vec))){
  print(cbind(beta[,i], fit_grpA$beta[,i], fit_grpA$beta[,i] - beta[,i]))
}

########################################################################################################

# # @importFrom ncvreg std
# # @export
# new_XZ = function(X, Z, var_subset = colnames(X), group_num, group_X){
#   # may eventually need to change var_subset = which(colnames(X) %in% fD_out$cnms) 
#   ## or just var_subset = fD_out$cnms ?
#   
#   # Assumptions:
#   ## X = model matrix with intercept present (Z also includs intercept)
#   ## group_X = vector composed of consecutive integers from 0 to J = total number groups,
#   ##    and first element corresponds to the intercept and equals 0
#   ## group_X = 0 means the covariate will not be penalized
#   
#   if(!is.factor(group_num)){
#     group_num <- factor(group_num)
#   }
#   
#   # Standardize X - ncvreg::std method
#   ## Assume input model matrix (with intercept)
#   X_std = std(X[,-1])
#   X_center = attr(X_std, "center")
#   X_scale = attr(X_std, "scale")
#   # Note: X_std = (X[,-1] - X_center) / X_scale
#   X_std = cbind(1, X_std)
#   
#   var_idx = which(colnames(X) %in% var_subset)
#   Z_center = X_center[var_idx]
#   Z_scale = X_scale[var_idx]
#   group_Z = group_X[var_idx]
#   
#   # Orthogonalize X
#   ortho_out = orthogonalize(X_std, group_X)
#   X_ortho = ortho_out$XX
#   SVD = ortho_out$SVD_lst
#   
#   # Standardize Z using X_std and X_orthog output
#   Z = Matrix::as.matrix(Z) ## Convert to non-sparse form if necessary
#   d = nlevels(group_num)
#   num_vars = ncol(Z) / d
#   
#   Z_std = matrix(0, nrow = nrow(Z), ncol = ncol(Z))
#   Z_ortho = matrix(0, nrow = nrow(Z), ncol = ncol(Z))
#   
#   # Standardize Z cols
#   for(v in 1:num_vars){ 
#     if("(Intercept)" %in% cnms) next # Don't need to scale/orthogonalize intercept
#     cols = seq(from = (v - 1)*d + 1, to = v*d, by = 1)
#     for(k in 1:nlevels(group_num)){
#       ids = which(group_num == k)
#       # Standardize Z cols
#       Z_std[ids, cols[k]] = (Z_sparse[ids, cols[k]] - Z_center[v-1]) / Z_scale[v-1]
#     }
#   }
#   
#   # Orthogonalize Z cols
#   
#   
#   return(list(X_std = X_std, Z_std = Z_std, X_center = X_center, X_scale = X_scale,
#               Z_center = Z_center, Z_scale = Z_scale))
# }

# orthogonalize = function(X, group_X){
#   
#   n = nrow(X)
#   J = max(group_X)
#   SVD_lst = list()
#   
#   XX = matrix(NA, nrow = nrow(X), ncol = ncol(X))
#   
#   # No adjustments for group_X == 0
#   XX[,which(group_X == 0)] == X[,which(group_X == 0)]
#   
#   for(j in 1:J){
#     
#     idx = which(group == j)
#     Xj = X[,idx]
#     SVD = svd(Xj, nu = 0)
#     SVD_lst[[j]][["v"]] = SVD$v
#     SVD_lst[[j]][["d"]] = SVD$d
#     
#     r = which(SVD$d > 1e-10)
#     
#     if(length(r) < length(idx)){
#       stop("Algorithm currently can't accomodate groups with rank < size of group. Error in group ", j, "\n")
#     }
#     
#     XX[,idx[r]] = Xj %*% SVD$v[,r] %*% diag(sqrt(n) / SVD$d[r])
#   }
#   
#   return(list(XX = XX, SVD_lst = SVD_lst))
#   
# }

library(Matrix)
d = 5
q = 5
Z = matrix(0, nrow = 1, ncol = d*q)
covar = "unstructured"

if(covar == "unstructured"){ # Originally: ncol(Z)/d <= 15 
  # create J, q2 x q*(q+1)/2
  J = Matrix(0, (ncol(Z)/d)^2, (ncol(Z)/d)*((ncol(Z)/d)+1)/2, sparse = T) #matrix(0, (ncol(Z)/d)^2, (ncol(Z)/d)*((ncol(Z)/d)+1)/2)
  index = 0
  indexc = 0
  sumy = 0
  sumx = 0
  zeros = 0
  covgroup = NULL
  for(i in 1:(ncol(Z)/d)){
    J[ sumx + zeros + 1:(ncol(Z)/d - (i-1)), sumy + 1:(ncol(Z)/d - (i-1))] = diag((ncol(Z)/d - (i-1)))
    sumy = sumy + (ncol(Z)/d - (i-1))
    sumx = sumy 
    zeros = zeros + i
    covgroup = rbind(covgroup, rep(i, (ncol(Z)/d)))
  }
  covgroup = covgroup[lower.tri(covgroup, diag = T)]
}else{ # covar == "independent". Originally: ncol(Z)/d > 15
  J = Matrix(0,(ncol(Z)/d)^2, (ncol(Z)/d), sparse = T) #matrix(0, (ncol(Z)/d)^2, (ncol(Z)/d))
  index = 0
  indexc = 0
  sumy = 0
  sumx = 0
  zeros = 0
  covgroup = NULL
  for(i in 1:(ncol(Z)/d)){
    J[ sumx + zeros + 1, i] = 1
    sumy = sumy + (ncol(Z)/d - (i-1))
    sumx = sumy 
    zeros = zeros + i
  }
  covgroup = rep(1:(ncol(Z)/d))
}

########################################################################################################