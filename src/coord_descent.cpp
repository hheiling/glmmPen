// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utility_glm.h"
#include "utility_CD.h"

using namespace Rcpp;
using namespace arma;

arma::vec coord_desc(arma::vec y, arma::mat X, arma::vec weights, arma::vec resid, 
                     arma::vec eta, arma::vec dims, arma::vec beta, 
                     const char* penalty, double lambda, double gamma, double alpha, // penalty type and parameters
                     const char* family, int link, arma::vec penalty_factor, int trace){
  
  // Define available penalties
  const char* lasso = "lasso";
  const char* mcp = "MCP";
  const char* scad = "SCAD";
  
  // const char* pois = "poisson";
  
  int p = dims(0); // number covariates (ncol(X))
  int N = dims(1); // total number observations (length(y))
  double conv = dims(2); // convergence criteria
  int maxit_CD = dims(3); // max number iterations for coordinate descent
  int i=0;
  int j=0;
  int iter=0;
  int converged=0; // 0 if not converged, 1 if converged
  
  // Notation: W = diag(weights), beta0 = beta from past iteration
  double nuj=0.0; // 1/N * t(X_j) %*% W %*% X_j
  double zetaj=0.0; // 1/N * t(X_j) %*% W %*% resid + nuj * beta0_j
  double v0=0.0; // sum(residuals^2) / N
  double v0_last=0.0; // value of v0 from past iteration
  
  arma::vec beta0(p); // Input initial beta / beta from last round of updates
  
  arma::vec mu(N); // expected mean for observations (g(eta) where g = link function)
  arma::vec deriv(N); // (d_eta / d_mu) = g'(mu) where g = link function
  arma::vec Vmu(N); // variance = b''(theta) where b(theta) from exponential family
  arma::vec constant(N); // Arbitrary vector used to create vector of ones
  arma::vec mu_check(N); // 1 if mu a valid number for given family, 0 otherwise
  
  v0 = sum(resid % resid) / N;
  
  while((iter<maxit_CD) && (converged==0)){
    
    // Add to iter
    iter = iter + 1;
    // Save latest value of beta
    beta0 = beta;
    // Save latest v0
    v0_last = v0;
    
    // Element-wise update of beta
    for(j=0; j<p; j++){
      
      if((iter>=5) && (beta(j) == 0)){
        // If beta penalized to zero in past round, will stay zero in further rounds
        // Therefore, skip to next covariate
        continue; 
      }
      
      // Update zeta and nu
      nuj = sum(X.col(j) % weights % X.col(j)) / N;
      zetaj = sum(X.col(j) % weights % resid) / N + nuj * beta(j);
      
      // Update beta
      if((j==0) || (penalty_factor(j)==0)){
        // No penalization for the intercept or variables specifically requested not to have penalization
        beta(j) = zetaj / nuj;
      }else{
        if(std::strcmp(penalty, lasso) == 0){
          beta(j) = soft_thresh(zetaj, lambda*alpha) / (nuj * (1.0 + lambda*(1.0 - alpha)));
        }else if(std::strcmp(penalty, mcp) == 0){
          beta(j) = MCP_soln(zetaj, nuj, lambda, gamma, alpha);
        }else if(std::strcmp(penalty, scad) == 0){
          beta(j) = SCAD_soln(zetaj, nuj, lambda, gamma, alpha);
        }
      }
      
      // Update eta and mu
      eta = eta + X.col(j) * (beta(j) - beta0(j));
      
      // Update residuals using fast update
      resid = resid - X.col(j) * (beta(j) - beta0(j));
      
    } // End j for loop
    
    // Update mu
    mu = invlink(link,eta);
    mu_check = muvalid(family, mu);
    
    // Update residuals and weights
    deriv = dlink(link, mu);
    Vmu = varfun(family, mu, 1.0); // Generic place-holder for phi (negbin family not fit in initial coefficient estimate)
    resid = deriv % (y - mu);
    weights = constant.ones() / (deriv % deriv % Vmu);
    
    // if((std::strcmp(family, pois) == 0) & (link == 20)){
    //   weights = Vmu;
    //   resid = (y - mu) / Vmu;
    // }else{
    //   
    //   
    // }
    
    for(i=0; i<N; i++){
      if((weights(i) <= 1e-8) || (mu_check(i) == 0)){
        resid(i) = 0.0;
        weights(i) = 0.0;
      }
    } 
    
    // Check convergence criteria
    v0 = sum(resid % resid) / N;
    if(fabs(v0 - v0_last)/(v0_last+0.1) < conv*0.001){
      converged = 1;
    }
    
    
    
  } // End while loop
  
  if((converged == 0) && (trace >= 1)){
    Rcpp::Rcout << "initial coordinate descent algorithm did not converge in " << maxit_CD << " iterations" << std::endl;
  }
  
  return beta;
  
}

