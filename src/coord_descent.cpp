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

// Soft-thresholding function (used in Lasso, MCP, and SCAD penalties)
// [[Rcpp::export]]
double soft_thresh(double zeta, double lambda){
  
  // # Soft-thresholding function - returns a scalar;
  double abs_z = fabs(zeta);
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
double MCP_soln(double zeta, double nu, double lambda, double gamma, double alpha){ 
  // gamma > 1
  
  double val = 0;
  double abs_z = fabs(zeta);
  
  if(abs_z <= lambda*alpha){
    val = 0;
  }else if(abs_z <= gamma*(lambda*alpha)*(1.0+lambda*(1.0-alpha))){
    val = soft_thresh(zeta, lambda*alpha) / (nu*(1.0 - (1.0/gamma) + lambda*(1-alpha)));
  }else{
    val = zeta / (nu*(1.0+lambda*(1.0-alpha)));
  }
  
  return val;
  
}

// [[Rcpp::export]]
double SCAD_soln(double zeta, double nu, double lambda, double gamma, double alpha){ 
  // gamma > 2 
  
  double val = 0;
  double abs_z = fabs(zeta);
  
  if(abs_z <= lambda*alpha){
    val = 0.0;
  }else if(abs_z <= (lambda*alpha)*(2.0+lambda*(1.0-alpha))){
    val = soft_thresh(zeta, lambda*alpha) / (nu*(1.0 +lambda*(1.0-alpha)));
  }else if(abs_z <= gamma*(lambda*alpha)*(1.0+lambda*(1.0-alpha))){
    val = soft_thresh(zeta, gamma*lambda*alpha/(gamma-1.0)) / (nu * (1.0 - (1.0/(gamma-1.0)) + lambda*(1.0-alpha)));
  }else{
    val = zeta / (nu*(1.0+lambda*(1.0-alpha)));
  }
  
  return val;
  
}


// [[Rcpp::export]]
arma::vec coord_desc(arma::vec y, arma::mat X, arma::vec weights, arma::vec resid, 
                     arma::vec eta, arma::vec offset, arma::vec dims, arma::vec beta, 
                     const char* penalty, double lambda, double gamma, double alpha, // penalty type and parameters
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
  double euclid_dist=0.0;
  double v0=0.0; // sum(residuals^2) / N
  double v0_last=0.0; // value of v0 from past iteration
  
  arma::vec beta0(p); // Input initial beta / beta from last round of updates
  
  arma::vec mu(N); // expected mean for observations (g(eta) where g = link function)
  arma::vec deriv(N); // (d_eta / d_mu) = g'(mu) where g = link function
  arma::vec Vmu(N); // variance = b''(theta) where b(theta) from exponential family
  arma::vec constant(N);
  arma::vec mu_check(N); // 1 if mu a valid number for given family, 0 otherwise
  arma::vec change_metric(p);
  
  
  v0 = sum(resid % resid) / N;
  
  while((iter<maxit_CD) & (converged==0)){
    
    // Add to iter
    iter = iter + 1;
    // Save latest value of beta
    beta0 = beta;
    // Save latest v0
    v0_last = v0;
    
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
          // st = soft_thresh(zetaj, lambda);
          // beta(j) = st / nuj;
          beta(j) = soft_thresh(zetaj, lambda*alpha) / (nuj * (1.0 + lambda*(1.0 - alpha)));
        }else if(std::strcmp(penalty, mcp) == 0){
          beta(j) = MCP_soln(zetaj, nuj, lambda, gamma, alpha);
        }else if(std::strcmp(penalty, scad) == 0){
          beta(j) = SCAD_soln(zetaj, nuj, lambda, gamma, alpha);
        }
      }
      
      // Update eta and mu
      eta = X * beta + offset;
      mu = invlink(link,eta);
      mu_check = muvalid(family, mu);
      
      // Update residuals and weights
      deriv = dlink(link, mu);
      Vmu = varfun(family, mu);
      resid = deriv % (y - mu);
      weights = constant.ones() / (deriv % deriv % Vmu);
      for(i=0; i<N; i++){
        if(mu_check(i) == 0){
          resid(i) = 0.0;
          weights(i) = 0.0;
        }
      }
      
      // ncvreg check of binomial weights
      // for(i=0; i<N; i++){
      //   if(weights(i) < 0.0001){
      //     weights(i) = 0.00001; // minimum weight allowed (same as ncvreg)
      //     resid(i) = (y(i) - mu(i)) / weights(i);
      //   }
      // }
      
    } // End j for loop
    
    // Check convergence criteria
    v0 = sum(resid % resid) / N;
    if(fabs(v0 - v0_last)/(v0_last+0.1) < conv*0.001){
      converged = 1;
    }
    
    // Alternative convergence criteria
    // Compare old vs new beta, check convergence criteria
    // euclid_dist = sqrt(sum((beta0 - beta) % (beta0 - beta)));
    // if(euclid_dist < conv){ 
    //   converged = 1;
    // }
    
    if((iter == maxit_CD) & (converged == 0)){
      warning("coord_desc algorithm did not converge \n");
    }
    
  } // End while loop
  
  // Rprintf("Number of iterations needed: %i \n", iter);
  
  return(beta);
  
}