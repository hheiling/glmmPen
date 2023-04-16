dat = sim.data(n = 500, ptot = 5, pnonzero = 5, nstudies = 5, sd_raneff = 0.5,
               family = "binomial", seed = 1618, beta = c(0, rep(1,5)))
(lam_seq = LambdaSeq(dat$X[,-1], dat$y, family = "binomial", lambda.min = 0.05, nlambda = 10))

y = dat$y
X = dat$X[,-1]
group = dat$group

fit = glmm(y ~ X + (X | group), family = "binomial", covar = "independent",
           tuning_options = lambdaControl(lambda0 = max(lam_seq), lambda1 = max(lam_seq)))
