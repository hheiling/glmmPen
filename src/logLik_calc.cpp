
#include <RcppArmadillo.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;

// [[Rcpp::export]]
double logLik_cpp(const arma::vec& U_means, const arma::mat& sigma, unsigned int M, 
                  const arma::vec& group, unsigned int n_levels, unsigned int df,
                  const arma::vec& y, const arma::vec& eta_fef, 
                  const arma::mat& Z, const arma::mat& Gamma,
                  const char* family){
  
  unsigned int s = 0;
  unsigned int i = 0;
  unsigned int r = U_means.n_elem;
  unsigned int k = 0;
  unsigned int d = n_levels;
  unsigned int num_vars = r / d;
  unsigned int v = 0;
  unsigned int n = Z.n_rows;
  unsigned int grp = 0;
  
  arma::mat samples(M, r);
  arma::uvec samp_index(num_vars);
  arma::uvec samp_index2(num_vars);
  arma::vec mvt_means(num_vars);
  arma::uvec index(num_vars);
  arma::mat samps(M,num_vars);
  arma::mat wts(M, d);
  arma::mat I = arma::eye(num_vars,num_vars);
  arma::vec z(num_vars); // Vector of zeros
  
  arma::vec Z_i(r);
  arma::vec eta(M);
  NumericVector prob(M);
  NumericVector prob_cap1(M);
  NumericVector prob_cap2(M);
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  
  double ll = 0;
  
  // Draw samples and calculate standardized weights
  
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
  
  // Calculated log-likelihood
  
  for(i = 0; i < n; i++){
    grp = group(i)-1;
    
    // Choose relevant columns of Z and samples
    for(v = 0; v < num_vars; v++){
      index(v) = grp + v*d;
      samp_index2(v) = grp*num_vars + v;
    }
    
    // Calculate Linear Predictor
    Z_i = Z.row(i).t();
    eta = eta_fef(i) + samples.cols(samp_index2) * Gamma * Z_i(index);
    
    // Calculate log-density, multiply by standardized weights, and sum
    if(std::strcmp(family, bin) == 0){
      // Calculate prob, and cap so that no element < 10^-16 or > 1 - 10^-16
      // Reasoning: avoids dbinom result equalling -Inf
      prob = exp(eta) / (1+exp(eta));
      prob_cap1 = ifelse((prob < pow(10,-16)), pow(10,-16), prob);
      prob_cap2 = ifelse((prob_cap1 > (1-pow(10,-16))), (1-pow(10,-16)), prob_cap1);
      for(s = 0; s < M; s++){ 
        ll = ll + R::dbinom(y(i), 1, prob_cap2(s), true) * wts(s, grp);
      }
    }else if(std::strcmp(family, pois) == 0){
      for(s = 0; s < M; s++){
        ll = ll + R::dpois(y(i), exp(eta(s)), true) * wts(s, grp);
      }
    }
    
  }
  
  return(ll);
}