#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Helper functions for IRLS

// [[Rcpp::export]]
arma::vec initial_mu(const char* family, arma::vec y, int N) {
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  const char* gaus = "gaussian";
  const char* Gamma = "Gamma";
  
  int i;
  
  arma::vec mu(N);
  
  if(std::strcmp(family, bin) == 0){
    for (i=0; i<N; i++) {
      if (y[i] < 0) {
        stop("negative values not allowed for the Binomial family");
      }
      if (y[i] > 1.0) {
        stop("# of success is larger than 1");
      }
      mu[i] = (y[i] + 0.5)/(1.0+1.0);
    }
  }else if(std::strcmp(family, pois) == 0){
    for (i=0; i<N; i++) {
      if (y[i] < 0) {
        stop("negative values not allowed for the Poisson family");
      }      
      mu[i] = y[i] + 0.1;
    }
  }else if(std::strcmp(family, gaus) == 0){
    for (i=0; i<N; i++) {
      mu[i] = y[i];
    }
  }else if(std::strcmp(family, Gamma) == 0){
    for (i=0; i<N; i++) {
      if (y[i] <= 0) {
        stop("non-poistive values not allowed for the Gamma family");
      }      
      mu[i] = y[i] + 0.1;
    }
  }
  
  return mu;
}

// [[Rcpp::export]]
arma::vec dlink(int link, arma::vec mu){
  
  int N = mu.n_elem;
  arma::vec out(N);
  
  switch(link){
  case 10: return(-1.0/(mu % (1.0-mu))); // logit
  case 11: return(out.zeros()); // probit not yet available
  case 12: return(1.0/(log(1.0-mu)%(1.0-mu))); // cloglog
  case 20: return(1.0/mu); // log
  case 30: return(out.ones()); // identity
  case 40: return(-1.0/(mu%mu)); // inverse
  default: return(out.zeros()); 
  }
  
}

// [[Rcpp::export]]
arma::vec linkfun(int link, arma::vec mu){
  
  int N = mu.n_elem;
  arma::vec out(N);
  
  switch(link){
  case 10: return(log(mu/(1.0-mu))); // logit
  case 11: return(out.zeros()); // probit not yet available
  case 12: return(log(-1.0*log(1.0-mu))); // cloglog
  case 20: return(log(mu)); // log
  case 30: return(mu); // identity
  case 40: return(1.0/mu); // inverse
  default: return(out.zeros());
  }
  
}

// [[Rcpp::export]]
arma::vec invlink(int link, arma::vec eta){
  
  int N = eta.n_elem;
  arma::vec out(N);
  
  switch(link){
  case 10: return(exp(eta)/(1.0+exp(eta))); // logit
  case 11: return(out.zeros()); // probit not yet available
  case 12: return(1.0-exp(-1.0*exp(eta))); // cloglog
  case 20: return(exp(eta)); // log
  case 30: return(eta); // identity
  case 40: return(-1.0/eta); // inverse
  default: return(out.zeros());
  }
  
}

arma::vec varfun(const char* family, arma::vec mu){ // double phi
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  const char* gaus = "gaussian";
  const char* Gamma = "Gamma";
  // const char* NB = "Negative Binomial(4)"; // MASS negative.binomial() family function
  
  int N = mu.n_elem;
  arma::vec V(N);
  
  if(std::strcmp(family, bin) == 0){
    V = mu%(1.0-mu);
  }else if(std::strcmp(family, pois) == 0){
    V = mu;
  }else if(std::strcmp(family, gaus) == 0){
    V = V.ones();
  }else if(std::strcmp(family, Gamma) == 0){
    V = mu%mu;
  }else{
    stop("invalid family \n");
  }
  
  // else if(std::strcmp(family, NB) == 0){
  //   V = mu + mu*mu*phi;
  // }
  
  return(V);
  
}


////////////////////////////////////////////////////////////////////////////////////////////////
