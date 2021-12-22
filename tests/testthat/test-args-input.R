
####################################################################################################
# Test input arguments
####################################################################################################


####################################################################################################
# Test glmmPen input arguments
####################################################################################################

test_that("test mis-specification of glmmPen arguments",{
  
  dat = sim.data(n = 500, ptot = 4, pnonzero = 4, nstudies = 5,
                 sd_raneff = 1.0, family = 'binomial',
                 seed = 1618, imbalance = 1, 
                 pnonzerovar = 0, beta = c(0,rep(1,4)))
  
  y = dat$y
  X = dat$X[,-1]
  group = dat$group
  
  expect_error(glmmPen(y ~ X + (X | group), family = "binomial",
                       penalty = "grMCP"),
               regexp = "not available")
  
  expect_error(glmmPen(y ~ X + (X | group), family = "binomial",
                       covar = "independence"),
               regexp = "covariance structure 'covar' must be")
  
  expect_error(glmmPen(y ~ X + (X | group), family = "binomial",
                       penalty = "MCP", gamma_penalty = 0.5),
               regexp = "gamma_penalty must be > 1 when using MCP penalty")
  
  expect_error(glmmPen(y ~ X + (X | group), family = "binomial",
                       penalty = "SCAD", gamma_penalty = 2),
               regexp = "gamma_penalty must be > 2 when using SCAD penalty")
  
  expect_error(glmmPen(y ~ X + (X | group), family = "binomial",
                       fixef_noPen = c(1,0)),
               regexp = "length of fixef_noPen")
  
})

####################################################################################################
# Test optimControl input arguments (sampling of arguments)
####################################################################################################

test_that("test mis-specification of optimControl arguments",{
  
  dat = sim.data(n = 500, ptot = 4, pnonzero = 4, nstudies = 5,
                 sd_raneff = 1.0, family = 'binomial',
                 seed = 1618, imbalance = 1, 
                 pnonzerovar = 0, beta = c(0,rep(1,4)))
  
  y = dat$y
  X = dat$X[,-1]
  group = dat$group
  
  expect_error(glmmPen(y ~ X + (X | group), family = "binomial",
                       optim_options = optimControl(sampler = "STAN")),
               regexp = "sampler must be specified")
  
  expect_error(glmmPen(y ~ X + (X | group), family = "binomial",
                       optim_options = optimControl(var_start = "rec")),
               regexp = "var_start must either be 'recommend'")
  
  expect_error(glmmPen(y ~ X + (X | group), family = "binomial",
                       optim_options = optimControl(var_start = -1)),
               regexp = "positive numeric value")
  
  expect_error(glmmPen(y ~ X + (X | group), family = "binomial",
                       optim_options = optimControl(nMC_start = 0)),
               regexp = "must be positive integers")
  
})

####################################################################################################
# Test selectControl input arguments
####################################################################################################

test_that("test mis-specification of selectControl arguments",{
  
  dat = sim.data(n = 500, ptot = 4, pnonzero = 4, nstudies = 5,
                 sd_raneff = 1.0, family = 'binomial',
                 slopes = T, seed = 1618, imbalance = 1, 
                 pnonzerovar = 0, beta = c(0,rep(1,4)))
  
  y = dat$y
  X = dat$X[,-1]
  group = dat$group
  
  expect_error(glmmPen(y ~ X + (X | group), family = "binomial",
                       tuning_options = selectControl(lambda0_seq = seq(from = 1, to = 0.1, by = -0.1))),
               regexp = "must be a sequence of ascending order")
  
  expect_error(glmmPen(y ~ X + (X | group), family = "binomial",
                       tuning_options = selectControl(lambda1_seq = seq(from = 1, to = 0.1, by = -0.1))),
               regexp = "must be a sequence of ascending order")
  
  expect_error(glmmPen(y ~ X + (X | group), family = "binomial",
                       tuning_options = selectControl(BIC_option = "BICH")),
               regexp = "BIC_option must be")
})

####################################################################################################
####################################################################################################
