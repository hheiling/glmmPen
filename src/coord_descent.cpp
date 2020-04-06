// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utility_glm.h"

using namespace Rcpp;

// // declare internal functions
// arma::vec initial_mu(const char* family, arma::vec y, int N);
// arma::vec dlink(int link, arma::vec mu);
// arma::vec linkfun(int link, arma::vec mu);
// arma::vec invlink(int link, arma::vec eta);
// arma::vec varfun(const char* family, arma::vec mu);

// zeta, nu, gamma, and lambda defined in coord_desc function

// [[Rcpp::export]]
double soft_thresh( double zeta, double lambda ) {
  
  // # Soft-thresholding function - returns a scalar;
  double abs_z = sqrt(zeta*zeta);
  double val = 0;
  
  if((zeta > 0) & (lambda < abs_z)){
    val = zeta - lambda;
  }else if((zeta < 0) & (lambda < abs_z)){
    val = zeta + lambda;
  } else if(lambda >= abs_z ){
    val = 0;
  }
  
  return val;
  
}

// [[Rcpp::export]]
double MCP_soln(double zeta, double nu, double lambda, double gamma){ // gamma > 1/nu
  
  double val = 0;
  double st = 0;
  double abs_z = sqrt(zeta*zeta);
  
  if(abs_z <= nu*gamma*lambda){
    st = soft_thresh(zeta, lambda);
    val = st / (nu - 1.0 / gamma);
  }else{
    val = zeta / nu;
  }
  
  return val;
  
}

// [[Rcpp::export]]
double SCAD_soln(double zeta, double nu, double lambda, double gamma){ // gamma > 1 + 1/nu
  
  double val = 0;
  double st = 0;
  double abs_z = abs(zeta);
  
  if(abs_z <= gamma*(nu+1)){ // |zeta| <= lambda*(nu+1)
    st = soft_thresh(zeta, lambda);
    val = st / nu;
  }else if(abs_z <= nu*gamma*lambda){ // lambda*(nu+1) < |zeta| <= nu*gamma*lambda
    st = soft_thresh(zeta, gamma*lambda/(gamma-1.0));
    val = st / (nu-1.0/(gamma-1.0));
  }else{ // |zeta| > nu*lambda*gamma
    val = zeta / nu;
  }
  
  return val;
  
}


// [[Rcpp::export]]
arma::vec coord_desc(arma::vec y, arma::mat X, arma::vec weights,
                     arma::vec resid, arma::vec eta, arma::vec dims, arma::vec beta, 
                     const char* penalty, double lambda, double gamma, // penalty type and parameters
                     const char* family, int link){
  
  // Define available penalties
  const char* lasso = "lasso";
  const char* mcp = "MCP";
  const char* scad = "SCAD";
  
  int p = dims(0); // number covariates (ncol(X))
  int N = dims(1); // total number observations (length(y))
  double conv = dims(2); // convergence criteria
  int maxit_CD = dims(4); // max number iterations for coordinate descent
  int i=0;
  int j=0;
  int iter=0;
  int converged=0; // 0 if not converged, 1 if converged
  
  // Notation: W = diag(weights), beta0 = beta from past iteration
  double nuj=0.0; // 1/N * t(X_j) %*% W %*% X_j
  double zetaj=0.0; // 1/N * t(X_j) %*% W %*% resid + nuj * beta0_j
  // double st=0.0; // soft_thresh result
  double euclid_dist=0.0;
  
  arma::vec beta0(p); // Input initial beta / beta from last round of updates
  
  arma::vec mu_pre(N); // Initial estimate of mu
  arma::vec mu(N); // expected mean for observations (g(eta) where g = link function)
  arma::vec deriv(N); // (d_eta / d_mu) = g'(mu) where g = link function
  arma::vec Vmu(N); // variance = b''(theta) where b(theta) from exponential family
  arma::vec constant(N);
  arma::vec mu_check(N); // 1 if mu a valid number for given family, 0 otherwise
  
  
  // Need to do something different for intercept?
  
  while((iter<maxit_CD) & (converged==0)){
    
    // Add to iter
    iter = iter + 1;
    // Save latest value of beta
    beta0 = beta;
    
    // Element-wise update of beta
    for(j=0; j<p; j++){
      
      if((iter>=5) & (beta(j) == 0)){
        // If beta penalized to zero in past round, will stay zero in further rounds
        // Therefore, skip to next covariate
        continue; 
      }
      
      // Update zeta and nu
      nuj = sum(X.col(j) % weights % X.col(j)) / N;
      zetaj = sum(X.col(j) % weights % resid) / N + nuj * beta(j);
      
      // Update beta
      if(j==0){
        // No penalization for the intercept
        beta(j) = zetaj / nuj;
      }else{
        if(std::strcmp(penalty, lasso) == 0){
          beta(j) = soft_thresh(zetaj, lambda) / nuj;
        }else if(std::strcmp(penalty, mcp) == 0){
          beta(j) = MCP_soln(zetaj, nuj, lambda, gamma);
        }else if(std::strcmp(penalty, scad) == 0){
          beta(j) = SCAD_soln(zetaj, nuj, lambda, gamma);
        }
      }
      
      // Update eta and mu
      eta = X * beta;
      mu = invlink(link,eta);
      // mu_check = muvalid(family, mu);
      // if(sum(mu_check) == 0.0){
      //   Rprintf("all mu invalid at iter %i \n", iter);
      // }
      
      // Update residuals and weights
      deriv = dlink(link, mu);
      Vmu = varfun(family, mu);
      resid = deriv % (y - mu);
      weights = constant.ones() / (deriv % deriv % Vmu);
      
      for(i=0; i<N; i++){
        if(weights(i) < 0.0001){
          weights(i) = 0.00001; // minimum weight allowed (same as ncvreg)
          resid(i) = (y(i) - mu(i)) / weights(i);
        }
      }
      // resid = (y - mu) / weights;
      
      // for(i=0; i<N; i++){
      //   if(mu_check(i) == 0){
      //     resid(i) = 0.0;
      //     weights(i) = 0.0;
      //   }
      // }
      
    } // End j for loop
      
    
    // Diagnose problems: print statements
    // Rcout << "beta" << std::endl << beta.t();
    // Rcout << "mu" << std::endl << mu.t();
    // Rcout << "weights" << std::endl << weights.t();
    // Rcout << "resids" << std::endl << resid.t();
    
    // Compare old vs new beta, check convergence criteria
    // Change convergence criteria?
    euclid_dist = sqrt(sum((beta0 - beta) % (beta0 - beta)));
    if(euclid_dist < conv){ 
      converged = 1;
    }
    
    if((iter == maxit_CD) & (converged == 0)){
      warning("coord_desc algorithm did not converge \n");
    }
    
    if(beta.has_nan()){
      Rprintf("beta has nan in iter %i \n", iter);
    }
    
  } // End while loop
  
  Rprintf("Number of iterations needed: %i \n", iter);
  
  return(beta);
  
}