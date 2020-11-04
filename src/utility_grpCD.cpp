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
  // arma::vec deriv(M);
  // arma::vec Vmu(M);
  arma::vec resid(M);
  arma::vec weights(M);
  arma::vec const_ones(M); const_ones.ones();
  
  const char* bin = "binomial";
  const char* gaus = "gaussian";
  const char* pois = "poisson";
  
  double nu_tmp = 0.0;
  
  arma::vec output(M+1);
  
  // Update mu, resid, weights
  mu = invlink(link, eta);
  mu_check = muvalid(family, mu);
  
  // If not binomial or gaussian family with canonical link, calculate nu = max(weights)
  // non-canonical link options not yet available
  if((std::strcmp(family,pois) == 0) & (link == 20)){
    
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
    
    // if(nu < 0.0001){
    //   nu = 0.0001;
    // }
  } // End family if-else
  
  resid = ((yi*const_ones) - mu); // Will divide zetaj by nu later
  for(m=0;m<M;m++){
    if(mu_check(m)==0){
      resid(m) = 0.0; // Ignore invalid mu
    }
  }
  
  
  output.subvec(0,M-1) = resid;
  output(M) = nu;
  
  return(output);
}

// Non-canconical link options
// arma::vec resid_nu_i(double yi, arma::vec eta, const char* family, int link, double nu, double phi){
//   
//   int M = eta.n_elem;
//   int m = 0;
//   
//   arma::vec mu(M);
//   arma::vec mu_check(M);
//   arma::vec deriv(M);
//   arma::vec Vmu(M);
//   arma::vec resid(M);
//   arma::vec weights(M);
//   arma::vec const_ones(M); const_ones.ones();
//   
//   const char* bin = "binomial";
//   const char* gaus = "gaussian";
//   
//   double nu_tmp = 0.0;
//   
//   arma::vec output(M+1);
//   
//   // Update mu, resid, weights
//   mu = invlink(link, eta);
//   mu_check = muvalid(family, mu);
//   
//   // If not binomial or gaussian family with canonical link, calculate nu = max(weights)
//   if(!((std::strcmp(family, bin) == 0) & (link == 10)) & !((std::strcmp(family, gaus) == 0) & (link == 30))){
//     
//     deriv = dlink(link, mu);
//     Vmu = varfun(family, mu, phi);
//     
//     weights = const_ones / (deriv % deriv % Vmu);
//     for(m=0;m<M;m++){
//       if(mu_check(m)==0){
//         weights(m) = 0.0;
//       }
//     }
//     nu_tmp = max(weights);
//     if(nu_tmp > nu){
//       nu = nu_tmp;
//     }
//   } // End family if-else
//   
//   resid = ((yi*const_ones) - mu); // Will divide zetaj by nu later
//   for(m=0;m<M;m++){
//     if(mu_check(m)==0){
//       resid(m) = 0.0; // Ignore invalid mu
//     }
//   }
//   
//   
//   output.subvec(0,M-1) = resid;
//   output(M) = nu;
//   
//   return(output);
// }

// Calculates fixed effects zetaj using residuals matrix
arma::vec zeta_fixef_calc(arma::mat X, arma::mat resid, arma::uvec idxj){
  
  arma::vec zetaj(idxj.n_elem); zetaj.zeros();
  int N = resid.n_cols;
  
  int i = 0;
  
  arma::mat Xj = X.cols(idxj);
  
  for(i=0;i<N;i++){
    // Update zetaj sum
    arma::vec Xji = Xj.row(i);
    zetaj = zetaj + Xji * sum(resid.col(i));
    
  } // End i for loop
  
  return(zetaj);
}


/////////////////////////////////////////////////////////////////////////////////////////////////