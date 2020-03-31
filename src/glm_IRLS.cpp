// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utility_glm.h"
#include "fit_calls.h"

using namespace Rcpp;

// // declare internal functions
// arma::vec initial_mu(const char* family, arma::vec y, int N);
// arma::vec dlink(int link, arma::vec mu);
// arma::vec linkfun(int link, arma::vec mu);
// arma::vec invlink(int link, arma::vec eta);
// arma::vec varfun(const char* family, arma::vec mu);

//' @export
// [[Rcpp::export]]
arma::vec glm_fit(const arma::vec& y, const arma::mat& X, const arma::vec& dims,
                   arma::vec beta,
                   const char* family, int link, int fit_type,
                   const char* penalty, double lambda, arma::vec params) {
  
  int p = dims(0); // number covariates (ncol(X))
  int N = dims(1); // total number observations (length(y))
  double conv = dims(2); // Convergence threshold
  int maxit = dims(3);
  int maxit_CD = dims(4);
  
  int i=0;
  int iter=0;
  int converged = 0;
  
  // penalty parameters (MCP, SCAD, elastic net)
  double gamma = params(0);
  double alpha = params(1);
  
  double epsilon = 1e-8; // Other threshold
  double wsum = 0.0; // Sum of weights
  
  // arma::vec beta(p); // coefficients
  arma::vec mu(N); // expected means 
  arma::vec eta(N); // linear predictors
  arma::vec deriv(N); // (d_eta / d_mu) = g'(mu) where g = link function
  arma::vec Vmu(N); // variance = b''(theta) where b(theta) from exponential family
  arma::vec weights(N);
  arma::vec resid(N); // IRLS residuals
  arma::vec z(N); // 'working response'
  arma::vec constant(N);
  // arma::vec fitted(N); // Estimate of mu given initial beta input
  
  // Recall link coding (see "M_step.R" for logic):
  // 10: logit, 11: probit, 12: cloglog
  // 20: log, 30: identity, 40: inverse
  
  // Initialize mu and eta
  mu = initial_mu(family, y, N);
  eta = X*beta; // Will input initial beta. linkfun(link, fitted) instead?
  
  // Initialize residuals and weights
  deriv = dlink(link, mu);
  Vmu = varfun(family, mu);
  resid = deriv % (y - mu);
  weights = constant.ones() / (deriv % deriv % Vmu);
  
  // Checks
  wsum = sum(weights);
  if(wsum < epsilon){
    stop("sum of weights too small");
  }
  
  //-------------------------------------------------------------------------------//
  // Fit Algorithms
  //-------------------------------------------------------------------------------//
  
  if(fit_type == 1){ // Regular IRLS
    
    //-------------------------------------------------------------------------------//
    // IRLS Algorithm
    //-------------------------------------------------------------------------------//
    
    if(p > 10){
      stop("IRLS method not appropriate for given dimension of X");
    }
    
    // Initialize variables
    arma::mat wX(N,p); // Element-wise multiplication of each column of X by sqrt(weights)
    arma::vec wZ(N); // Element-wise multiplication of sqrt(weights) * z
    double euclid_dist = 1000.0; // Euclidean distance between past beta and updated beta
    
    // Caluclate initial 'working response'
    z = eta + resid;
    
    // Save initial beta
    arma::vec beta0 = beta;
    
    iter = 0;
    converged = 0;
    
    while(iter < maxit & converged == 0){
      
      // Update beta given latest weights and working responses
      wX = sqrt(weights) % X.each_col();
      wZ = sqrt(weights) % z;
      beta = (wX.t() * wX).i() * wX.t() * wZ;
      
      // Update mu and eta
      eta = X * beta;
      mu = invlink(link, eta);
      
      // Update residuals and weights
      deriv = dlink(link, mu);
      Vmu = varfun(family, mu);
      resid = deriv % (y - mu);
      weights = constant.ones() / (deriv % deriv % Vmu);
      
      // Update z (working response)
      z = eta + resid;
      
      // Compare old vs new beta (test convergence criteria)
      euclid_dist = as_scalar(sqrt((beta0-beta)%(beta0-beta)));
      if(euclid_dist < conv){
        converged = 1;
      }
      
      // Store new beta
      beta0 = beta;
    }
    
  }else if(fit_type == 2){ // Coordinate descent (regular)
    
    //-------------------------------------------------------------------------------//
    // Coordinate Descent (ungrouped)
    //-------------------------------------------------------------------------------//
    
    beta = coord_desc(y, X, weights, resid, dims, beta, penalty, lambda, gamma, family, link);
    
  }else if(fit_type == 3){ // Grouped coordinate descent
    
    //-------------------------------------------------------------------------------//
    // Grouped Coordinate Descent
    //-------------------------------------------------------------------------------//
    
    stop("grouped coordinate descent not yet available \n");
    
  }
  

  
  return beta;
}

