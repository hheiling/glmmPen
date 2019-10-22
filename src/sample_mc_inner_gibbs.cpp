// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// Original: NumericMatrix sample_mc_inner_gibbs(...)

// [[Rcpp::export]]
List sample_mc_inner_gibbs(arma::mat f, // matrix
                          arma::mat z, // matrix
                          arma::vec y, // vector
                          arma::vec t, // vector
                          int NMC, // integer
                          arma::vec u0, //matrix
                          int trace){ // integer
  
  arma::mat fitted = f;
  arma::mat Z=z;
  arma::vec Y=y;
  arma::vec tau =t;
  int nMC = NMC; 
  arma::vec uold=u0;
  
  int q = Z.n_cols;
  int n = Z.n_rows;
  
  int i = 0;
  int j = 0; 
  // int k = 0; 
  int l = 0; 
  int naccept = 0;
  int index = 0;
  // int batch = 1;
  double sum = 0;
  double sumn = 0;
  // double stao = 0;
  double w = 0;
  double ep = 0; 
  double e0 = 0;
  arma::mat out(nMC, q);
  arma::vec e(q);
  arma::vec rate(q);
  arma::vec etae(n);
  arma::vec etaen(n);
  arma::vec index2(q);
  arma::vec var(q);
  arma::vec acc_rate(q);
  
  RNGScope scope;
  
  //initialize e 
  for(i = 0; i < q; i++){
    e(i) = uold(i);
    index2(i) = 0.0;
    var(i) = 1.5;
  }
  
  //iteratively update e
  while(naccept < nMC){
    
    for(j = 0;j < q; j++){
      
      // calculate etae
      etae = fitted + Z*e;
      
      // save current value of e(i)
      e0 = e(j);
      
      // generate w
      w = log(R::runif(0.0,1.0));
      
      // generate proposal e
      ep = R::rnorm(0.0, var(j));
      e(j) = ep;
      
      // calculate updated etae
      etaen = fitted + Z*e;
      
      // calculate ratio
      sum = sumn = 0;
      for(l = 0; l < n; l++){
        sum = sum + R::dbinom(Y(l), 1.0, exp(etae(l))/(1+exp(etae(l))), 1);
        sumn = sumn + R::dbinom(Y(l), 1.0, exp(etaen(l))/(1+exp(etaen(l))), 1);
      }
      
      // check for acceptance
      if(w < sumn - sum){
        // ep left in e(i)
        index2(j) = index2(j) + 1.0;
      }else{
        e(j) = e0;
      }
      //  }
    }
    
    if(index == 1+5*naccept){
      out.row(naccept) = e.t();
      naccept++;
    }
    index++;
    
  }
  
  // Info for acceptance rate:
  for(i = 0; i < q; i++){
    acc_rate(i) = index2(i) / index;
  }

  // if(trace == 2){
  //   Rprintf("index: %d, naccept: %d, accept. rate: %f  \n", index, naccept, acc_rate);
  // }
  
  return(List::create(Named("u") = wrap(out), Named("acc_rate") = acc_rate));
  // return(wrap(out));
  
}
