// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utility_glm.h"
#include "penalties.h"

using namespace Rcpp;

// Calculates the M residuals for an individual given eta, updates nu and v0
arma::vec resid_nu_v0_i(double yi, arma::vec eta, const char* family, int link, double nu){
  
  int M = eta.n_elem;
  int m = 0;
  
  arma::vec mu(M);
  arma::vec mu_check(M);
  arma::vec deriv(M);
  arma::vec Vmu(M);
  arma::vec resid(M);
  arma::vec weights(M);
  arma::vec const_ones(M); const_ones.ones();
  
  const char* bin = "binomial";
  const char* gaus = "gaussian";
  
  double nu_tmp = 0.0;
  double v0 = 0.0;
  
  arma::vec output(M+2);
  
  // Update mu, resid, weights
  mu = invlink(link, eta);
  mu_check = muvalid(family, mu);
  
  deriv = dlink(link, mu);
  Vmu = varfun(family, mu);
  
  // If not binomial or gaussian family with canonical link, calculate nu = max(weights)
  if(!((std::strcmp(family, bin) == 0) & (link == 10)) & !((std::strcmp(family, gaus) == 0) & (link == 30))){
    
    weights = const_ones / (deriv % deriv % Vmu);
    for(m=0;m<M;m++){
      if(mu_check(m)==0){
        weights(m) = 0.0;
      }
    }
    nu_tmp = max(weights);
    if(nu_tmp > nu){
      nu = nu_tmp;
    }
  } // End family if-else
  
  resid = ((yi*const_ones) - mu); // Will divide zetaj by nu later
  for(m=0;m<M;m++){
    if(mu_check(m)==0){
      resid(m) = 0.0; // Ignore invalid mu
    }
  }
  
  v0 = sum(resid % resid); // Will divide v0 by nu^2 and M*N later
  
  output.subvec(0,M-1) = resid;
  output(M) = nu;
  output(M+1) = v0;
  
  return(output);
}

// Calculates zetaj for fixed effects starting with eta
arma::vec zeta_fixef(arma::vec y, arma::mat X, arma::mat eta, 
                     arma::uvec idxr, const char* family, int link, double nu){
  
  arma::vec zetar(idxr.n_elem); zetar.zeros();
  int M = eta.n_rows;
  int N = eta.n_cols;
  
  int i = 0;
  int m = 0;
  
  arma::vec mu(M);
  arma::vec mu_check(M);
  arma::vec deriv(M);
  arma::vec Vmu(M);
  arma::vec resid(M);
  arma::vec weights(M);
  arma::vec const_ones(M); const_ones.ones();
  arma::vec out(M+2);
  arma::vec output(idxr.n_elem + 2);
  
  const char* bin = "binomial";
  const char* gaus = "gaussian";
  
  double nu_tmp = 0.0;
  double v0 = 0.0;
  
  arma::mat Xr = X.cols(idxr);
  
  for(i=0;i<N;i++){
    
    out = resid_nu_v0_i(y(i), eta.col(i), family, link, nu);
    resid = out.subvec(0,M-1);
    nu = out(M);
    v0 = v0 + out(M+1);
    
    // Update zetar sum
    arma::vec Xri = Xr.row(i);
    zetar = zetar + Xri * sum(resid);
    
  } // End i for loop
  
  // List L = List::create(Named("zetar") = zetar, Named("nu") = nu, Named("v0") = v0);
  // 
  // return(L);
  
  output.subvec(0,idxr.n_elem-1) = zetar;
  output(idxr.n_elem) = nu;
  output(idxr.n_elem+1)= v0;
  
  return(output);
}

// Calculates residuals by group (used in calculation of zetaj for random effects covariate(s))
// Also updates nu and v0
arma::vec resid_nu_v0_k(arma::vec y, arma::vec eta, const char* family, int link, double nu){
  
  int Nk = eta.n_elem;
  int n = 0;
  
  arma::vec mu(Nk);
  arma::vec mu_check(Nk);
  arma::vec deriv(Nk);
  arma::vec Vmu(Nk);
  arma::vec resid(Nk);
  arma::vec weights(Nk);
  arma::vec const_ones(Nk); const_ones.ones();
  
  const char* bin = "binomial";
  const char* gaus = "gaussian";
  
  double nu_tmp = 0.0;
  double v0 = 0.0;
  
  arma::vec output(Nk+2);
  
  // Update mu, resid, weights
  mu = invlink(link, eta);
  mu_check = muvalid(family, mu);
  
  deriv = dlink(link, mu);
  Vmu = varfun(family, mu);
  // If not binomial or gaussian family with canonical link, calculate nu = max(weights)
  if(!((std::strcmp(family, bin) == 0) & (link == 10)) & !((std::strcmp(family, gaus) == 0) & (link == 30))){
    
    weights = const_ones / (deriv % deriv % Vmu);
    for(n=0;n<Nk;n++){
      if(mu_check(n)==0){
        weights(n) = 0.0;
      }
    }
    nu_tmp = max(weights);
    if(nu_tmp > nu){
      nu = nu_tmp;
    }
  }
  
  resid = (y - mu); // Will divide by nu later once at end of zetaj calculation (to make sure nu = max(all weights))
  for(n=0;n<Nk;n++){
    if(mu_check(n)==0){
      resid(n) = 0.0; // Ignore invalid mu
    }
  }
  
  v0 = sum(resid % resid); // Will divide by nu^2 and (N*M) later (to make sure nu = max(all weights))
  
  output.subvec(0,Nk-1) = resid;
  output(Nk) = nu;
  output(Nk+1) = v0;
  
  return(output);
}

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

