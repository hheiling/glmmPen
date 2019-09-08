#include <RcppDist.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;

// Harmonic mean estimate function

// [[Rcpp::export]]
double logLik_HM(const arma::mat& samples, unsigned int M, unsigned int n_levels, 
                  const arma::vec& group, const arma::vec& y, const arma::vec& eta_fef,
                  const arma::mat& Z, const arma::mat& Gamma, const char* family){
  unsigned int s = 0;
  unsigned int i = 0;
  unsigned int r = Z.n_cols ;
  unsigned int k = 0;
  unsigned int d = n_levels;
  unsigned int num_vars = r / d;
  unsigned int v = 0;
  
  arma::uvec index(num_vars);
  arma::mat samps(M,num_vars);
  
  arma::vec prob(M);
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  
  double l_i = 0;
  double ll_i = 0;
  double ll = 0;
  
  // Calculate Harmonic Mean
  
  for(k = 0; k < d; k++){
    for(v = 0; v < num_vars; v++){
      index(v) = k+v*d;
    }
    
    // Specify cols of MCMC samples associated with Group k
    samps = samples.cols(index);
    
    // Identify the n_k individuals corresponding to group k
    arma::uvec ids = find((group-1) == k);
    arma::mat Z_k = Z.rows(ids); 
    arma::rowvec y_k = y(ids).t();
    
    // Calculate Linear Predictor
    arma::rowvec eta_fef_k = eta_fef(ids).t(); // 1 x n_k, fixed effects
    arma::mat eta_ref = samps * Gamma * Z_k.cols(index).t() ; // M x n_k, random effects
    arma::mat eta_mat = eta_fef_k + eta_ref.each_row(); // M x n_k, combined
      // Result: each column = full eta values for one individual 
    
    // Calculate HM estimate of log-likelihood 
    
    if(std::strcmp(family, bin) == 0){
      arma::mat prob_mat = exp(eta_mat) / (1+exp(eta_mat)); // M x n_k
      for(i = 0; i < eta_mat.n_cols; i++){
        // Calculate individual contribution to likelihood
        l_i = 0;
        prob = prob_mat.col(i);
        for(s = 0; s < M; s++){
          l_i = l_i + (1.00 / (R::dbinom(y_k(i), 1, prob(s), false))) ;
        }
        l_i = l_i / ((double)M);
        ll_i = log(1.00 / l_i);
        // Add individual contribution to overall log lik
        ll = ll_i;
      }
      
    }else if(std::strcmp(family, pois) == 0){
      for(i = 0; i < eta_mat.n_cols; i++){
        l_i = 0;
        for(s = 0; s < M; s++){
          l_i = l_i + (1.00 / (R::dpois(y_k(i), exp(eta_mat(s,i)), false))) ;
        }
        l_i = l_i / ((double)M);
        ll_i = log(1.00 / l_i);
        // Add individual contribution to overall log lik
        ll = ll + ll_i;
      }
      
    }
    
  }
  
  return(ll);
  
}

// Modified Harmonic Mean Estimate Function

// [[Rcpp::export]]
double logLik_mHM(const arma::mat& samples, const arma::vec& modes, const arma::mat& sigma,
                  unsigned int M, unsigned int n_levels, double chisq_cutoff,
                  const arma::vec& group, const arma::vec& y, const arma::vec& eta_fef,
                  const arma::mat& Z, const arma::mat& Gamma, const char* family){
  unsigned int s = 0;
  unsigned int i = 0;
  unsigned int r = modes.n_elem;
  unsigned int k = 0;
  unsigned int d = n_levels;
  unsigned int num_vars = r / d;
  unsigned int v = 0;
  
  arma::vec mvnorm_means(num_vars);
  arma::uvec index(num_vars);
  arma::mat samps(M,num_vars);
  arma::vec trunc_norm(M);
  arma::vec trunc_region(M);
  // arma::vec trunc_decision(M);
  arma::vec prior(M);
  arma::mat I = arma::eye(num_vars,num_vars);
  arma::vec z(num_vars); // Vector of zeros
  
  arma::vec prob(M);
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  
  double l_i = 0;
  double ll_i = 0;
  double ll = 0;
  
  // Calculate Modified Harmonic Mean (Gelfand and Dey, 1994)
  
  for(k = 0; k < d; k++){
    for(v = 0; v < num_vars; v++){
      index(v) = k+v*d;
      mvnorm_means(v) = modes(index(v));
    }
    
    // Specify cols of MCMC samples associated with Group k
    samps = samples.cols(index);
    
    // Calculate Truncated Normal Density
    // Based on recommendation from Chan and Grant (2015)
    trunc_norm = dmvnorm(samps, mvnorm_means, sigma, false);
    for(s = 0; s < M; s++){
      trunc_region(s) = as_scalar((samps.row(s) - mvnorm_means.t()) * sigma.i() * trans(samps.row(s) - mvnorm_means.t()));
    } 
    arma::uvec trunc_cutoff = find(trunc_region > chisq_cutoff);
    trunc_norm(trunc_cutoff).zeros();
    
    Rprintf("Number above cutoff: %i \n", trunc_cutoff.n_elem);
    
    // Calculate Prior Density for MCMC Samples
    prior = dmvnorm(samps, z, I, false); 
    
    // Identify the n_k individuals corresponding to group k
    arma::uvec ids = find((group-1) == k);
    arma::mat Z_k = Z.rows(ids); 
    arma::rowvec y_k = y(ids).t();
    
    // Calculate Linear Predictor
    arma::rowvec eta_fef_k = eta_fef(ids).t(); // 1 x n_k, fixed effects
    arma::mat eta_ref = samps * Gamma * Z_k.cols(index).t() ; // M x n_k, random effects
    arma::mat eta_mat = eta_fef_k + eta_ref.each_row(); // M x n_k, combined
    
    // Calculate modified HM estimate of log-likelihood 
    
    if(std::strcmp(family, bin) == 0){
      arma::mat prob_mat = exp(eta_mat) / (1+exp(eta_mat)); // M x n_k
      for(i = 0; i < eta_mat.n_cols; i++){
        // Calculate individual contribution to likelihood
        l_i = 0;
        prob = prob_mat.col(i);
        for(s = 0; s < M; s++){
          l_i = l_i + trunc_norm(s) / (R::dbinom(y_k(i), 1, prob(s), false) * prior(s)) ;
        }
        l_i = l_i / ((double)M);
        ll_i = log(1.00 / l_i);
        // Add individual contribution to overall log lik
        ll = ll + ll_i;
      }
      
    }else if(std::strcmp(family, pois) == 0){
      for(i = 0; i < eta_mat.n_cols; i++){
        l_i = 0;
        for(s = 0; s < M; s++){
          l_i = l_i + trunc_norm(s) / (R::dpois(y_k(i), exp(eta_mat(s,i)), false) * prior(s)) ;
        }
        l_i = l_i / ((double)M);
        ll_i = log(1.00 / l_i);
        // Add individual contribution to overall log lik
        ll = ll + ll_i;
      }
      
    }
    
  }
  
  return(ll);
  
}