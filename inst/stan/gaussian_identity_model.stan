data {
  int<lower=0> N; // Number of observations in group k (N_k)
  int<lower=0> q; // Number of random effect variables
  vector[N] eta_fef; // fixed effects portion of linear predictor for individuals in group k
  real y[N]; // Outcome values in group k
  matrix[N,q] Z; // Portion of Z * Gamma matrix for group k
  real<lower=0> sigma; // standard deviation (sqrt of variance of residual error)
}

parameters {
  vector[q] alpha; // Random effects to sample
}

model {
  alpha ~ normal(0,1); // Prior of random effects (Gamma matrix factored out into Z matrix)
  // Distribution of y based on fixed and random effects
  y ~ normal(eta_fef + Z * alpha, sigma); 
}
