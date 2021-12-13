
####################################################################################################
# Test data input
####################################################################################################

####################################################################################################
# Test response matches family
####################################################################################################

test_that("test response input appropriate for family", {
  
  set.seed(1618)
  n = 100
  p = 4
  X = matrix(rnorm(n*p), nrow = n, ncol = p)
  group = rep(1:5, times = 20)
  
  # Binomial family
  ## response is coded as categorical variable with too many response types (not binary)
  y_cat = sample(letters[1:3], size = 100, replace = T)
  ## response is coded as numeric values that are not 0 or 1
  y_num = sample(c(10,11), size = 100, replace = T)
  
  expect_error(glmmPen(y_cat ~ X + (X | group), family = "binomial"), 
               regexp = "response must be binary")
  
  expect_error(glmmPen(y_num ~ X + (X | group), family = "binomial"), 
               regexp = "binary variable with levels 0 vs 1")
  
  # Poisson family
  ## response with non-integers
  y_nonint = sample(seq(from = 0, to = 3, by = 0.5), size = 100, replace = T)
  ## response with negative values
  y_neg = sample(c(-1,0,1), size = 100, replace = T)
  
  expect_error(glmmPen(y_nonint ~ X + (X | group), family = "poisson"), 
               regexp = "must be integer")
  
  expect_error(glmmPen(y_neg ~ X + (X | group), family = "poisson"), 
               regexp = "must be non-negative integer")
  
})


####################################################################################################
# Test random effects
####################################################################################################

test_that("test random effects specifications meet model assumptions", {
  
  set.seed(1618)
  n = 100
  p = 4
  X = matrix(rnorm(n*p), nrow = n, ncol = p)
  colnames(X) = c("var1","var2","var3","var4")
  group = rep(1:5, times = 20)
  y = sample(c(0,1), size = 100, replace = T)

  data = data.frame(y = y, X, group = group)
  
  # Random effects must be subset of fixed effects
  expect_error(glmmPen(data = data, formula = y ~ var1 + var2 + (var3 + var4 | group), 
                       family = "binomial"),
               regexp = "random effects must be a subset of fixed effects")
  
  # Random effects must include random intercept
  expect_error(glmmPen(data = data, formula = y ~ var1 + var2 + (-1 + var1 | group), 
                       family = "binomial"),
               regexp = "Model requires a random intercept term")
  
  
})


####################################################################################################
# Test restriction of single group factor
####################################################################################################

test_that("test only single grouping factor allowed", {
  
  set.seed(1618)
  n = 100
  p = 4
  X = matrix(rnorm(n*p), nrow = n, ncol = p)
  colnames(X) = c("var1","var2","var3","var4")
  group = rep(1:5, times = 20)
  group2 = rep(1:2, times = 50)
  y = sample(c(0,1), size = 100, replace = T)
  
  data = data.frame(y = y, X, group = group, group2 = group2)
  
  expect_error(glmmPen(data = data, 
                       formula = y ~ var1 + var2 + var3 + var4 + (var1 | group) + (var4 | group2),
                       family = "binomial"),
               regexp = "can only handle one group")
  
})


