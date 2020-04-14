// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utility_glm.h"
#include "fit_calls.h"

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
arma::vec glm_fit(arma::vec y, arma::mat X, arma::vec dims,
                  arma::vec beta, arma::vec offset,
                  const char* family, int link, int fit_type,
                  arma::vec group_X, arma::vec K_X, // K_X = vector of size of groups in X
                  const char* penalty, double lambda, arma::vec params) {
  
  int p = dims(0); // number covariates (ncol(X))
  int N = dims(1); // total number observations (length(y))
  double conv = dims(2); // Convergence threshold
  int maxit = dims(3);
  
  int i=0;
  int iter=0;
  int converged = 0;
  int r = 0; // rank of X.t() * W * X
  
  // penalty parameters (MCP, SCAD, elastic net)
  double gamma = params(0);
  double alpha = params(1);
  
  double epsilon = 1e-8; // Other threshold
  double wsum = 0.0; // Sum of weights
  double v0 = 0.0; // sum(residuals^2) / N
  double v0_last = 0.0; // v0 from last iteration
  
  arma::vec mu(N); // expected means 
  arma::vec eta(N); // linear predictors
  arma::vec deriv(N); // (d_eta / d_mu) = g'(mu) where g = link function
  arma::vec Vmu(N); // variance = b''(theta) where b(theta) from exponential family
  arma::vec weights(N);
  arma::vec resid(N); // IRLS residuals
  arma::vec z(N); // 'working response'
  arma::vec constant(N);
  arma::vec mu_check(N);
  
  // Recall link coding (see "M_step.R" for logic):
  // 10: logit, 11: probit, 12: cloglog
  // 20: log, 30: identity, 40: inverse
  
  // Initialize mu and eta
  mu = initial_mu(family, y, N);
  mu_check = muvalid(family, mu);
  eta = X * beta + offset; // Will input initial beta. linkfun(link, fitted) instead?
  
  // Initialize residuals and weights
  deriv = dlink(link, mu);
  Vmu = varfun(family, mu);
  resid = deriv % (y - mu);
  weights = constant.ones() / (deriv % deriv % Vmu);
  for(i=0; i<N; i++){
    if(mu_check(i)==0){
      weights(i) = 0.0;
      resid(i) = 0.0;
    }
  }
  
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
    
    // Initialize variables
    arma::mat wX(N,p); // Element-wise multiplication of each column of X by sqrt(weights)
    arma::vec wZ(N); // Element-wise multiplication of sqrt(weights) * z
    // double euclid_dist = 1000.0; // Euclidean distance between past beta and updated beta
    
    // Check if X.t() * W * X is singular
    wX = sqrt(weights) % X.each_col();
    r = arma::rank(wX.t() * wX);
    // If rank < p, small penalty > 0?
    if((p > 10) | (r < p)){
      // If dimension too high or rank of X.t() * W * X < p, 
      // use lasso coordinate descent with lambda = 0.0 and alpha = 1.0
      const char* lasso = "lasso";
      beta = coord_desc(y, X, weights, resid, eta, offset, dims, beta, lasso, 0.0, gamma, 1.0, family, link);
      
    }else{
      
      // Caluclate initial 'working response'
      z = eta + resid;
      
      // Save initial beta
      arma::vec beta0 = beta;
      
      iter = 0;
      converged = 0;
      v0_last = sum(resid % resid) / N;
      
      while((iter < maxit) & (converged == 0)){
        
        // Update beta given latest weights and working responses
        wX = sqrt(weights) % X.each_col();
        wZ = sqrt(weights) % z;
        beta = (wX.t() * wX).i() * wX.t() * wZ;
        
        // Update mu and eta
        eta = X * beta + offset;
        mu = invlink(link, eta);
        mu_check = muvalid(family, mu);
        
        // Update residuals and weights
        deriv = dlink(link, mu);
        Vmu = varfun(family, mu);
        resid = deriv % (y - mu);
        weights = constant.ones() / (deriv % deriv % Vmu);
        for(i=0; i<N; i++){
          if(mu_check(i) == 0){
            weights(i) = 0.0;
            resid(i) = 0.0;
          }
        }
        
        // Update z (working response)
        z = eta + resid;
        
        // Check convergence criteria
        v0 = sum(resid % resid) / N;
        if(fabs(v0 - v0_last)/(v0_last+0.1) < conv*0.001){
          converged = 1;
        }
        
        // Alternative convergence criteria
        // Compare old vs new beta (test convergence criteria)
        // euclid_dist = sqrt(sum((beta0-beta)%(beta0-beta)));
        // if(euclid_dist < conv){
        //   converged = 1;
        // }
        
        // Store new beta
        beta0 = beta;
        // Record latest v0
        v0_last = v0;
        
        // Update counter
        iter = iter + 1;
        
        if((iter == maxit) & (converged == 0)){
          warning("iteration limit reached for IRLS w/o convergence \n");
        }
        
      }
      
    }
    
  }else if(fit_type == 2){ // Coordinate descent (regular)
    
    //-------------------------------------------------------------------------------//
    // Coordinate Descent (ungrouped)
    //-------------------------------------------------------------------------------//
    
    beta = coord_desc(y, X, weights, resid, eta, offset, dims, beta, penalty, lambda, gamma, alpha, family, link);
    
  }else if(fit_type == 3){ // Grouped coordinate descent
    
    //-------------------------------------------------------------------------------//
    // Grouped Coordinate Descent
    //-------------------------------------------------------------------------------//
    
    beta = grp_CD(y, X, weights, resid, eta, dims, beta, group_X, K_X, penalty, lambda, gamma, alpha, family, link);
    
  }
  
  return beta;
}

