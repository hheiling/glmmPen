data {
  int<lower=0> N; // Number of observations in group k (N_k)
  int<lower=0> q; // Number of latent variables (glmmPen: number random effects; glmmPen_FA: number latent common factors)
  vector[N] eta_fef; // fixed effects portion of linear predictor for individuals in group k
  array[N] int<lower=0, upper=1> y; // Outcome values in group k
  matrix[N,q] Z; // Portion of Z * Gamma matrix for group k
}

parameters {
  vector[q] alpha; // Random effects to sample
}

model {
  alpha ~ normal(0,1); // Prior of random effects (Gamma matrix factored out into Z matrix)
  y ~ bernoulli_logit(eta_fef + Z * alpha); // Distribution of y based on fixed and random effects
}
