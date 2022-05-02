// Fit a penalized generalized linear mixed model (no random effects)

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utility_glm.h"
#include "utility_CD.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec pglm_fit(arma::vec y, arma::mat X, arma::vec dims,
                  arma::vec beta, arma::vec offset,
                  const char* family, int link,
                  const char* penalty, arma::vec params,
                  arma::vec penalty_factor, int trace) {
  
  int p = dims(0); // number fixed effects (including intercept)
  int N = dims(1); // total number observations (length(y))
  int i=0;
  
  // penalty parameters (MCP, SCAD, elastic net)
  double lambda = params(0);
  double gamma = params(1);
  double alpha = params(2);
  
  double epsilon = 1e-8; // Other threshold
  double wsum = 0.0; // Sum of weights
  
  arma::vec mu(N); // expected means 
  arma::vec eta(N); // linear predictors
  arma::vec deriv(N); // (d_eta / d_mu) = g'(mu) where g = link function
  arma::vec Vmu(N); // variance = b''(theta) where b(theta) from exponential family
  arma::vec weights(N);
  arma::vec resid(N); // IRLS residuals
  arma::vec constant(N);
  arma::vec mu_check(N);
  arma::vec beta_new(p); beta_new.zeros();
  
  // Recall link coding (see "fit_dat_B.R" for logic):
  // 10: logit, 11: probit, 12: cloglog
  // 20: log, 30: identity, 40: inverse
  
  // Initialize mu and eta
  mu = initial_mu(family, y, N);
  mu_check = muvalid(family, mu);
  eta = X * beta + offset; // Will input initial beta
  
  // Initialize residuals and weights
  deriv = dlink(link, mu);
  Vmu = varfun(family, mu, 1.0); // Generic place-holder for phi (negbin family not fit in initial coefficient estimate)
  resid = deriv % (y - mu);
  weights = constant.ones() / (deriv % deriv % Vmu);
  for(i=0; i<N; i++){
    if((weights(i) <= 1e-8) || (mu_check(i) == 0)){
      resid(i) = 0.0;
      weights(i) = 0.0;
    }
  }
  
  // Checks
  wsum = sum(weights);
  if(wsum < epsilon){
    stop("sum of weights too small");
  }
  
  //-------------------------------------------------------------------------------//
  // Coordinate Descent (ungrouped)
  //-------------------------------------------------------------------------------//
  
  beta_new = coord_desc(y, X, weights, resid, eta, dims, beta, penalty, lambda, gamma, alpha, family, link, penalty_factor, trace);
  
  return beta_new;
}

