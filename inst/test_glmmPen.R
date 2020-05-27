# Test glmmPen function using following settings:

# Setting 1
## fit_dat_B function
## binomial example
## Input formula: y ~ X + (X | group) where input just vectors and matrices
## Check both covar = "unstructured" and covar = "independent"

dat = sim.data2(n = 500, ptot = 2, pnonzero = 2, nstudies = 5, sd_raneff = 1.0, family = 'binomial',
                      slopes = T, seed = 2020, imbalance = 1, pnonzerovar = 0, beta = c(0, 2, 2))

y = dat$y
X = dat$X[,-1]
grp = dat$group
out1 = glmmPen(formula = y ~ X + (X | grp), family = "binomial", 
               optim_options = optimControl(nMC_start = 200, nMC_max = 1000, maxitEM = 10, fit_type=2), 
               tuning_options = lambdaControl())

out1B = glmmPen(formula = y ~ X + (X | grp), family = "binomial", 
               optim_options = optimControl(nMC_start = 200, nMC_max = 1000, maxitEM = 10, covar = "independent", fit_type=2), 
               tuning_options = lambdaControl())

head(predict(out1B, type = "response"))

# Setting 2
## fit_dat_B function
## gaussian example
## Input formula: y ~ X + (X | group) where input just vectors and matrices

# dat = sim.data2(n = 500, ptot = 2, pnonzero = 2, nstudies = 5, sd_raneff = 1.0, family = 'gaussian',
#                 slopes = T, seed = 2020, imbalance = 1, pnonzerovar = 0, beta = c(0, 2, 2))
# 
# y = dat$y
# X = dat$X[,-1]
# grp = dat$group
# out2A = glmmPen(formula = y ~ X + (X | grp), family = "gaussian", 
#                optim_options = optimControl(nMC_start = 200, nMC_max = 1000, maxitEM = 10, fit_type=2), 
#                tuning_options = lambdaControl())

# Setting 3
## select_tune function (nlambda = 3, 5 total covariates, covar = "unstructured")
## binomial example
## Input formula: y ~ X + (X | group) where input just vectors and matrices

dat = sim.data2(n = 500, ptot = 2, pnonzero = 2, nstudies = 5, sd_raneff = 1.0, family = 'binomial',
                slopes = T, seed = 2020, imbalance = 1, pnonzerovar = 0, beta = c(0, 2, 2))

y = dat$y
X = dat$X[,-1]
grp = dat$group
out3 = glmmPen(formula = y ~ X + (X | grp), family = "binomial", 
               optim_options = optimControl(nMC_start = 200, nMC_max = 1000, maxitEM = 10, fit_type=2), 
               tuning_options = selectControl(nlambda=2))

##################################################################################################