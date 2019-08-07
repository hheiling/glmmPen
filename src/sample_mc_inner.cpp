
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

//' Obtaining Rejection Sampling Monte Carlo Draws
//' 
//' Description
//' 
//' @param f a one-column matrix of the fixed effects component of the linear predictor for 
//' the individuals/observations in a specific group (i.e. group k)
//' @param z a sub-matrix of the full Z matrix (see output information of \code{\link{formulaData}} 
//' for details); the columns pertain to the columns of Z associated 
//' with a particular group, the rows peratin to the rows of Z that are associated with
//' the group of interest
//' @param y the sub-vector of the original response vector that is associated with the particular
//' group of interest
//' @param NMC positive integer, number of Monte Carlo draws to obtain
//' @param trace an integer specifying print output to include as function runs. If trace = 2, then 
//' will output acceptance rate information
//' 
//' @return A matrix of Rejection sampling Monte Carlo draws. nrows: \code{NMC}, columns: 
//' each column pertains to one of the q random effects variables that are associated with the 
//' group of interest (ncols = q). If the algorithm errored-out due to an unacceptably low 
//' acceptance rate, then a matrix of all zeros of dimension (NMC-3)xq is output.
//' 
//' @export
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
    
    if(index == pow(10,6)){
      acc_rate = ((double)naccept) / index;
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
  if(trace == 2){
    acc_rate = ((double)naccept) / index;
    Rprintf("index: %d, naccept: %d, accept. rate: %f  \n", index, naccept, acc_rate);
  }
  
  if(naccept == nMC){
    return(wrap(out));
  }else{
    return(wrap(error_out));
  }
  
}
