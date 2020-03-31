// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utility_glm.h"

using namespace Rcpp;

// declare internal functions
arma::vec initial_mu(const char* family, arma::vec y, int N);
arma::vec dlink(int link, arma::vec mu);
arma::vec linkfun(int link, arma::vec mu);
arma::vec invlink(int link, arma::vec eta);
arma::vec varfun(const char* family, arma::vec mu);

// zeta, nu, gamma, and lambda defined in coord_desc function

// [[Rcpp::export]]
double soft_thresh( double zeta, double lambda ) {
  
  // # Soft-thresholding function - returns a scalar;
  // Input;
  // z: continuous value;
  // lambda: continuous value;
  
  arma::mat Z2(1,1);
  Z2 = zeta * zeta;
  double abs_z = as_scalar( sqrt(Z2) );
  
  double val;
  
  if( zeta > 0 && lambda < abs_z ){
    val = zeta - lambda;
  }
  if( zeta < 0 && lambda < abs_z ){
    val = zeta + lambda;
  }
  if( lambda >= abs_z ){
    val = 0;
  }
  
  return val;
  
}

// [[Rcpp::export]]
double MCP_soln(double zeta, double nu, double lambda, double gamma){ // gamma > 1/nu
  
  double val;
  arma::mat Z2(1,1);
  Z2 = zeta * zeta;
  double abs_z = as_scalar( sqrt(Z2) );
  
  if(abs_z <= nu*gamma*lambda){
    val = soft_thresh(zeta, lambda) / (nu - 1.0 / gamma);
  }else{
    val = zeta / nu;
  }
  
  return val;
  
}

// [[Rcpp::export]]
double SCAD_soln(double zeta, double nu, double lambda, double gamma){ // gamma > 1 + 1/nu
  
  double val;
  arma::mat Z2(1,1);
  Z2 = zeta * zeta;
  double abs_z = as_scalar( sqrt(Z2) );
  
  if(abs_z <= gamma*(nu+1)){ // |zeta| <= lambda*(nu+1)
    val = soft_thresh(zeta, lambda) / nu;
  }else if(abs_z <= nu*gamma*lambda){ // lambda*(nu+1) < |zeta| <= nu*gamma*lambda
    val = soft_thresh(zeta, gamma*lambda/(gamma-1.0)) / (nu-1.0/(gamma-1.0));
  }else{ // |zeta| > nu*lambda*gamma
    val = zeta / nu;
  }
  
  return val;
  
}


// [[Rcpp::export]]
arma::vec coord_desc(arma::vec y, arma::mat X, arma::vec weights,
                     arma::vec resid, arma::vec dims, arma::vec beta, 
                     const char* penalty, double lambda, double gamma, // penalty type and parameters
                     const char* family, int link){
  
  // Define available penalties
  const char* lasso = "lasso";
  const char* mcp = "MCP";
  const char* scad = "SCAD";
  
  int p = dims(0); // number covariates (ncol(X))
  int N = dims(1); // total number observations (length(y))
  int j=0;
  int iter=0;
  int converged=0; // 0 if not converged, 1 if converged
  
  double conv = dims(2);
  int maxit_CD = dims(4);
  
  
  // Notation: W = diag(weights), beta0 = beta from past iteration
  double nuj=0; // 1/N * t(X_j) %*% W %*% X_j
  double zetaj=0; // 1/N * t(X_j) %*% W %*% resid + nuj * beta0_j
  double euclid_dist=0;
  
  arma::vec beta0 = beta; // Input initial beta / beta from last round of updates
  
  arma::vec wXj(N); // sqrt(weights) * X_j (jth column of X)
  arma::vec wR(N); // sqrt(weights) * resid
  // arma::vec mu0 = mu; // Initial input mu / mu from last round of updates
  arma::vec eta(N); // Linear predictor (X %*% Beta)
  arma::vec mu(N); // expected mean for observations (g(eta) where g = link function)
  arma::vec deriv(N); // (d_eta / d_mu) = g'(mu) where g = link function
  arma::vec Vmu(N); // variance = b''(theta) where b(theta) from exponential family
  arma::vec constant(N);
  
  
  // Need to do something different for intercept?
  
  while(iter<maxit_CD & converged==0){
    
    for(j=0; j<p; j++){
      
      if(iter>=1 & beta(j) == 0){
        // If beta penalized to zero in past round, will stay zero in further rounds
        // Therefore, skip to next covariate
        continue; 
      }
      
      // Update zeta and nu given most recent weights and residuals
      wXj = sqrt(weights) % X.col(j);
      wR = sqrt(weights) % resid;
      nuj = as_scalar(wXj.t() * wXj) / (double)N;
      zetaj = as_scalar(wXj.t() * wR) / (double)N + nuj * beta(j);
      
      // Update beta
      if(j==0){
        // No penalization for the intercept
        beta(j) = zetaj / nuj;
      }else{
        if(std::strcmp(penalty, lasso) == 0){
          beta(j) = soft_thresh(zetaj, lambda);
        }else if(std::strcmp(penalty, mcp) == 0){
          beta(j) = MCP_soln(zetaj, nuj, lambda, gamma);
        }else if(std::strcmp(penalty, scad) == 0){
          beta(j) = SCAD_soln(zetaj, nuj, lambda, gamma);
        }
      }
      
      // Update eta and mu
      eta = X * beta;
      mu = invlink(link, eta);
      
      // Update residuals and weights
      deriv = dlink(link, mu);
      Vmu = varfun(family, mu);
      resid = deriv % (y - mu);
      weights = constant.ones() / (deriv % deriv % Vmu);
      
    } // End j for loop
    
    // Compare old vs new beta, check convergence criteria
    // Change convergence criteria?
    euclid_dist = as_scalar(sqrt((beta0 - beta) % (beta0 - beta)));
    if(euclid_dist < conv){ 
      converged = 1;
    }
    
    // Update counter
    iter += 1;
    
    if(iter == maxit_CD & converged == 0){
      warning("coord_desc algorithm did not converge \n");
    }
    
    // Update beta0
    beta0 = beta;
    
  } // End while loop
  
  return(beta);
  
}