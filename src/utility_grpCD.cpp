// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utility_glm.h"


using namespace Rcpp;

// Calculates the M residuals for an individual given eta, updates nu and resid
// Canonical link option only
arma::vec resid_nu_i(double yi, arma::vec eta, const char* family, int link, double nu, double phi){
  
  int M = eta.n_elem;
  int m = 0;
  
  arma::vec mu(M);
  arma::vec mu_check(M);
  arma::vec resid(M);
  arma::vec weights(M);
  arma::vec const_ones(M); const_ones.ones();
  
  // const char* bin = "binomial";
  // const char* gaus = "gaussian";
  const char* pois = "poisson";
  
  double nu_tmp = 0.0;
  
  arma::vec output(M+1);
  
  // Update mu, resid, weights
  mu = invlink(link, eta);
  mu_check = muvalid(family, mu);
  
  // If not binomial or gaussian family with canonical link, calculate nu = max(weights)
  // non-canonical link options not yet available
  if((std::strcmp(family,pois) == 0) && (link == 20)){

    weights = varfun(family, mu, phi);

    for(m=0;m<M;m++){
      if(mu_check(m)==0){
        weights(m) = 0.0;
      }
    }
    nu_tmp = max(weights);
    if(nu_tmp > nu){
      nu = nu_tmp;
    }

    if(nu < 0.0001){
      nu = 0.0001;
    }
  } // End family if-else
  
  resid = ((yi*const_ones) - mu); // Will divide total zetaj by nu later
  for(m=0;m<M;m++){
    if(mu_check(m)==0){
      resid(m) = 0.0; // Ignore invalid mu
    }
  }
  
  
  output.subvec(0,M-1) = resid;
  output(M) = nu;
  
  return(output);
}

// Calculates the M residuals for an individual given eta
// Assumes update to nu occurs elsewhere
// Used in "grp_CD_XZ_step.cpp" code
arma::vec resid_i(double yi, arma::vec eta, const char* family, int link){
  
  int M = eta.n_elem;
  int m = 0;
  
  arma::vec mu(M);
  arma::vec mu_check(M);
  arma::vec resid(M);
  arma::vec const_ones(M); const_ones.ones();
  
  // Update mu and resid
  mu = invlink(link, eta);
  mu_check = muvalid(family, mu);
  
  resid = ((yi*const_ones) - mu); // Will divide by nu later
  for(m=0;m<M;m++){
    if(mu_check(m)==0){
      resid(m) = 0.0; // Ignore invalid mu
    }
  }
  
  return(resid);
}


// Calculates fixed effects zetaj using residuals matrix
arma::vec zeta_fixef_calc(arma::mat X, arma::mat resid, arma::uvec idxj){
  
  arma::vec zetaj(idxj.n_elem); zetaj.zeros();
  int N = resid.n_cols;
  
  int i = 0;
  
  arma::mat Xj = X.cols(idxj);
  arma::vec Xji(idxj.n_elem); Xji.zeros();
  
  for(i=0;i<N;i++){
    // Update zetaj sum
    Xji = trans(Xj.row(i));
    zetaj = zetaj + Xji * sum(resid.col(i));
    
  } // End i for loop
  
  return(zetaj);
}

// Quadratic approximation to Q function estimate (based on Taylor series expansion about 
//   the previous beta0 coefficient estimates)
// Note: in M-step, have already calculated Q evaluated at past coefficient value (Q0)
// Note: nu = 1 / step_size
// [[Rcpp::export]]
double Qfun_quad_beta(double Q0, double step_size, const arma::mat& diff0,
                      const arma::mat& eta, const arma::mat& eta0,
                      const arma::vec& beta, const arma::vec& beta0){
  
  int N = eta.n_cols;
  int M = eta.n_rows;
  int p_tot = beta.n_elem;
  
  int i = 0;
  
  double Q_quad = 0;
  double term1 = 0;
  double term2 = 0;
  
  arma::vec eta_diff(M);
  arma::vec beta_diff(p_tot);
  
  for(i=0;i<N;i++){
    eta_diff = eta.col(i) - eta0.col(i);
    term1 = term1 + sum(diff0.col(i) % eta_diff);
  }
  
  beta_diff = beta - beta0;
  term2 = sum(beta_diff % beta_diff);
  
  Q_quad = Q0 - term1 / M + (0.5 * N / step_size) * term2 ; 
  
  return(Q_quad);
  
}






/////////////////////////////////////////////////////////////////////////////////////////////////
