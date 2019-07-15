
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix sample_mc_inner(arma::mat f, // matrix
                              arma::mat z, // matrix
                              arma::vec y, // vector
                              arma::vec t, // vector
                              int NMC, //integer
                              int trace){ //integer
  arma::mat fitted=f;
  arma::mat Z=z;
  arma::vec Y=y;
  arma::vec tau =t;
  int nMC = NMC;
  
  
  int q = Z.n_cols;
  int n = Z.n_rows;
  
  int i = 0;
  int naccept = 0;
  int index = 0;
  double sum = 0;
  double stao = 0;
  double w = 0;
  double acc_rate = 0;
  arma::mat out(nMC, q);
  arma::mat error_out(nMC-3, q);
  arma::vec e(q);
  arma::vec etae(n);
  
  RNGScope scope;
  
  for(i = 0; i < n; i++) stao = stao + tau(i);
  
  while( (naccept < nMC) & (index < pow(10,8)*nMC)){
    
    if((naccept == 0) & (index > pow(10,8))){
      break;
    }
    
    if(index == round(0.25*pow(10,8))){
      acc_rate = naccept / index;
      if(acc_rate < pow(10,-3)){
        break;
      }
    }
    
    //populate e
    for(i = 0; i < q; i++) e(i) = R::rnorm(0.0, 1.0);
    
    //calculate etae
    etae = fitted + Z*e;
    
    //generate w
    w = log(R::runif(0.0,1.0));
    
    // calculate sum
    sum = 0;
    
    for(i = 0; i < n; i++){
      sum = sum + R::dbinom(Y(i), 1.0, exp(etae(i))/(1+exp(etae(i))), 1);
    }
    
    // check to save
    if(w < sum - stao){
      out.row(naccept) = e.t();
      naccept++;
    }
    index++;
    
  } 
  
  // Info for acceptance rate:
  if(trace > 1){
    acc_rate = naccept / index;
    Rprintf("index: %d, naccept: %d, accept. rate: %f  \n", index, naccept, acc_rate);
  }
  
  // acc_rate = naccept / index;
  Rprintf("index: %d, naccept: %d, accept. rate: %f  \n", index, naccept, naccept / index);
  
  if(naccept == nMC){
    return(wrap(out));
  }else{
    return(wrap(error_out));
  }
  
}
