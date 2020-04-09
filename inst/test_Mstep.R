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

cbind(fit1$coefficients, fit2$coef, fit3$coef)

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

fitA = ncvreg(X_std, y, family = "binomial", penalty = "lasso", nlambda = 100,
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
                coef = beta_init, penalty = "lasso", alpha = 1.0,
                maxit_CD = 1000, conv = 0.0001, fit_type = 2, lambda = lambda_vec[i])
  beta[,i] = fitB$coef
}

colnames(beta) = round(lambda_vec, 5)

# beta[,1:10]
# fitA$beta[,1:10]
for(i in 1:5){
  print(cbind(beta[,i], fitA$beta[,i]))
}



########################################################################################################