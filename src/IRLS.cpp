// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec IRLS(arma::vec y, arma::mat X, arma::vec mu, arma::vec weights,
               arma::vec z, arma::vec dims){
  
  int p = dims(0); // number covariates (ncol(X))
  int N = dims(1); // total number observations (length(y))
  int j=0;
  
  arma::vec wXj(N); // sqrt(weights) * jth col of X
  
  arma::vec beta(p);
  
  for(j=0; j<p; j++){
    wXj = sqrt(weights) * X.col(j);
    // beta(j) = 
  }
  
  return(beta);
}