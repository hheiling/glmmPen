
#include <RcppArmadillo.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;

// [[Rcpp::export]]
double logLik_cpp(const arma::vec& U_means, const arma::mat& sigma, unsigned int M, 
                  const arma::vec& group, unsigned int n_levels, unsigned int df,
                  const arma::vec& y, const arma::vec& eta_fef, 
                  const arma::mat& Z, const arma::vec& coef_q,
                  const arma::vec& cols, arma::sp_mat& J, const char* family){
  
  unsigned int s = 0;
  unsigned int i = 0;
  unsigned int r = U_means.n_elem;
  unsigned int k = 0;
  unsigned int d = n_levels;
  unsigned int num_vars = r / d;
  unsigned int v = 0;
  unsigned int n = Z.n_rows;
  
  arma::mat samples(M, r);
  arma::uvec samp_index(num_vars);
  arma::vec mvt_means(num_vars);
  arma::mat samps(M,num_vars);
  arma::mat wts(M, d);
  arma::mat I = arma::eye(num_vars,num_vars);
  arma::vec z(num_vars); // Vector of zeros
  
  arma::mat Znew(n, J.n_cols);
  arma::vec eta(n);
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  
  double ll = 0;
  
  for(k = 0; k < d; k++){
    for(v = 0; v < num_vars; v++){
      mvt_means(v) = U_means(k+v*d);
      samp_index(v) = k*num_vars + v;
    }
    
    // Collect multivariate t samples
    samps = rmvt(M, mvt_means, sigma, df);
    samples.cols(samp_index) = samps;
    
    // Calculate unstandardized weights (w_star)
    wts.col(k) = dmvnorm(samps, z, I, false) / dmvt(samps, mvt_means, sigma, df, false);
    
    // Calculated standardized weights
    wts.col(k) = wts.col(k) / sum(wts.col(k));
  }
  
  Function Znew_mat("Znew_mat");
  
  for(s = 0; s < M; s++){
    Znew = as<arma::mat>(Znew_mat(samples.row(s), Z, group, cols, n, num_vars, d, J));
    eta = eta_fef + Znew * coef_q;
    if(std::strcmp(family, bin) == 0){
      for(i = 0; i < n; i++){
        ll = ll + R::dbinom(y(i), 1, (exp(eta(i)) / (1+exp(eta(i)))), true) * wts(s, group(i)-1);
        
      }
    }else if(std::strcmp(family, pois) == 0){
      for(i = 0; i < n; i++){
        ll = ll + R::dpois(y(i), exp(eta(i)), true) * wts(s, group(i)-1);
      }
    } 
    
  } 
  
  return(ll);
}



// [[Rcpp::export]]
arma::mat Znew_mat(const arma::vec& U, const arma::mat& Z, const arma::vec& g, const arma::vec& cols, 
                   unsigned int n, unsigned int q, unsigned int d, 
                   arma::sp_mat& J){ 
  
  unsigned int i = 0;
  // unsigned int j = 0;
  // unsigned int nMC = U.n_rows;
  unsigned int index = 0;
  int gr = 0;
  arma::mat Usub(1, q);
  arma::mat Zsub(1, q);
  arma::mat out(n, J.n_cols);
  // XPtr<BigMatrix> pMat(pBigMat); 
  // MatrixAccessor<double> out2(*pMat); 
  
  for(i = 0; i<n; i++){
    gr = g(i);
    for(index = 0;index<q; index++){
      Zsub(0,index) = Z(i,cols(index)-1 + gr - 1);
    }
    for(index = 0;index<q; index++){
      Usub(0,index) = U(cols(index)-1 + gr - 1);
    }
    out.row(i) = kron(Usub, Zsub) * J;
    
  }
  
  // for(i = 0; i<n*nMC; i++){
  //   for(j = 0; j<out.n_cols; j++){
  //     out2[j][i] = out(i, j);
  //   }
  // }
  
  return(out);
}
