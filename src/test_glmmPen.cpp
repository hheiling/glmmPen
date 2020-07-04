// #define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

using namespace Rcpp ;

// [[Rcpp::export]]
void arma_test_value( arma::mat x){}

// [[Rcpp::export]]
void arma_test_ref( arma::mat& x){}

// [[Rcpp::export]]
void arma_test_const_ref( const arma::mat& x){}

//' @export
// [[Rcpp::export]]
void Znew_gen2( const arma::mat& U, const arma::mat& Z, const arma::vec& g, const arma::vec& cols, 
                unsigned int n, unsigned int q, unsigned int d, SEXP pBigMat, 
                arma::sp_mat& J){ 
  
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int nMC = U.n_rows;
  unsigned int index = 0;
  int gr = 0;
  arma::mat Usub(1, q);
  arma::mat Zsub(1, q);
  arma::mat out(n*nMC, J.n_cols);
  XPtr<BigMatrix> pMat(pBigMat); 
  MatrixAccessor<double> out2(*pMat); 
  
  for(i = 0; i<n; i++){
    gr = g(i);
    for(index = 0;index<q; index++){
      Zsub(0,index) = Z(i,cols(index)-1 + gr - 1);
    }
    for(j = 0; j < nMC; j++){
      for(index = 0;index<q; index++){
        Usub(0,index) = U(j,cols(index)-1 + gr - 1);
      }
      out.row((i*nMC + j)) = kron(Usub, Zsub) * J;
    }
  }
  
  for(i = 0; i<n*nMC; i++){
    for(j = 0; j<out.n_cols; j++){
      out2[j][i] = out(i, j);
    }
  }
  
}



// [[Rcpp::export]]
NumericMatrix orthog_inner( arma::mat& X, const arma::vec& g, int gmax, int gmin, int n){
  
  using namespace arma  ;
  // for svd
  arma::mat U;
  arma::mat V;
  arma::mat X2;
  arma::vec S;
  
  // other vars
  unsigned int j; //i removed
  
  // run svd
  svd(U, S, V, X);
  // transform X and save T
  for(j = 0; j < S.size(); j++){
    if(S(j) < pow(10, -10)){
      continue;
    }else{
      V.col(j) *= sqrt(n)/S(j); //fix this
    }
  }
  //  X2 = X.cols(gmin, gmax) * V;
  return wrap(V); 
}

/*
void standardize2(Rcpp::XPtr<BigMatrix> pBigMat, const arma::vec& c, const arma::vec& s, int n) {

// Declarations
MatrixAccessor<double> out2(*pBigMat);

int p = c.Size();

for (int j=0; j<p; j++) {
// Center
c[j] = 0;
for (int i=0; i<n; i++) {
c[j] += X[j][i];
}
c[j] = c[j] / n;
for (int i=0; i<n; i++) X[j][i] = X[j][i] - c[j];

// Scale
s[j] = 0;
for (int i=0; i<n; i++) {
s[j] += pow(X[j][i], 2);
}
s[j] = sqrt(s[j]/n);
for (int i=0; i<n; i++) X[j][i] = X[j][i]/s[j];
}
}
*/
