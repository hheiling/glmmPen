
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
  
  double ll_i = 0;
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
    
    ll_i = 0;
    
    // Calculate log-density, multiply by standardized weights, and sum
    if(std::strcmp(family, bin) == 0){
      
      prob = exp(eta) / (1+exp(eta));
      for(s = 0; s < M; s++){ 
        ll_i = ll_i + R::dbinom(y(i), 1, prob(s), false) * wts(s, grp);
      }
      
    }else if(std::strcmp(family, pois) == 0){
      for(s = 0; s < M; s++){
        ll_i = ll_i + R::dpois(y(i), exp(eta(s)), false) * wts(s, grp);
      }
      
    }
    
    ll = ll + log(ll_i / M);
    
  }
  
  return(ll);
}

// [[Rcpp::export]]
double logLik_modif(const arma::vec& U_means, const arma::mat& sigma, unsigned int M, 
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
  
  arma::vec mvt_means(num_vars);
  arma::uvec index(num_vars);
  arma::mat samps(M,num_vars);
  arma::vec wts(M);
  arma::mat I = arma::eye(num_vars,num_vars);
  arma::vec z(num_vars); // Vector of zeros
  
  arma::vec prob(M);
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  
  double l_i = 0;
  double ll_k = 0;
  double ll = 0;
  
  // Draw samples and calculate standardized weights
  
  for(k = 0; k < d; k++){
    for(v = 0; v < num_vars; v++){
      index(v) = k+v*d;
      mvt_means(v) = U_means(index(v));
    }
    
    // Collect multivariate t samples
    samps = rmvt(M, mvt_means, sigma, df);
    
    // Calculate unstandardized weights (w_star)
    wts = dmvnorm(samps, z, I, false) / dmvt(samps, mvt_means, sigma, df, false);
    
    // Identify the n_k individuals corresponding to group k
    arma::uvec ids = find((group-1) == k);
    arma::mat Z_k = Z.rows(ids);
    arma::rowvec y_k = y(ids).t();
    
    // Calculate Linear Predictor
    arma::rowvec eta_fef_k = eta_fef(ids).t(); // 1 x n_k
    arma::mat eta_ref =  samps * Gamma * Z_k.cols(index).t() ; // M x n_k
    arma::mat eta_mat = eta_fef_k + eta_ref.each_row(); // M x n_k
    
    // Calculate log-likelihood - group contribution
    // Function dbinom_mat("dbinom_mat");
    ll_k = 0;
    
    if(std::strcmp(family, bin) == 0){
      // Calculate density - vectorized version:
      // arma::mat prob_mat = exp(eta_mat) / (1+exp(eta_mat)); // M x n_k
      // arma::mat prob_log = log(prob_mat);
      // arma::mat logdens = y_k % prob_log.each_row() + (1-y_k) % prob_log.each_row();
      // arma::mat dens = exp(logdens);
      // // // arma::mat dens = as<arma::mat>(dbinom_mat(y(ids),1,prob_mat));
      // arma::mat dens_wt = wts % dens.each_col();
      // arma::rowvec l_i = sum(dens_wt,0) / M;
      // ll_k = sum(log(l_i));
      // Unvectorized version: 
      arma::mat prob_mat = exp(eta_mat) / (1+exp(eta_mat)); // M x n_k
      for(i = 0; i < eta_mat.n_cols; i++){
        // Calculate individual contribution to likelihood
        l_i = 0;
        prob = prob_mat.col(i);
        for(s = 0; s < M; s++){
          l_i = l_i + R::dbinom(y_k(i), 1, prob(s), false) * wts(s);
        }
        // Add individual lik contribution to group log-lik contribution
        ll_k = ll_k + log(l_i / M);
      }
      
    }else if(std::strcmp(family, pois) == 0){
      for(i = 0; i < eta_mat.n_cols; i++){
        l_i = 0;
        for(s = 0; s < M; s++){
          l_i = l_i + R::dpois(y_k(i), exp(eta_mat(s,i)), false) * wts(s);
        }
        ll_k = ll_k + log(l_i / M);
      }
      
    }
    
    // Add up group contribution
    ll = ll + ll_k;
  }
  
  return(ll);
}

// [[Rcpp::export]]
double logLik_MCI(unsigned int M, const arma::vec& group, unsigned int n_levels,
                  const arma::vec& y, const arma::vec& eta_fef, 
                  const arma::mat& Z, const arma::mat& Gamma,
                  const char* family){
  
  unsigned int s = 0;
  unsigned int i = 0;
  unsigned int r = Z.n_cols;
  unsigned int k = 0;
  unsigned int d = n_levels;
  unsigned int num_vars = r / d;
  unsigned int v = 0;
  
  arma::uvec index(num_vars);
  arma::mat samps(M,num_vars);
  arma::mat I = arma::eye(num_vars,num_vars);
  arma::vec z(num_vars); // Vector of zeros
  
  arma::vec prob(M);
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  
  double l_i = 0;
  double ll_k = 0;
  double ll = 0;
  
  // Draw samples and calculate log-lik
  
  for(k = 0; k < d; k++){
    for(v = 0; v < num_vars; v++){
      index(v) = k+v*d;
    }
    
    // Collect stnadard MVN samples
    samps = rmvnorm(M, z, I);
    
    // Identify individuals corresponding to group k
    arma::uvec ids = find((group-1) == k);
    arma::mat Z_k = Z.rows(ids);
    arma::rowvec y_k = y(ids).t();
    
    // Calculate Linear Predictor
    arma::rowvec eta_fef_k = eta_fef(ids).t(); // 1 x n_k
    arma::mat eta_ref =  samps * Gamma * Z_k.cols(index).t() ; // M x n_k
    arma::mat eta_mat = eta_fef_k + eta_ref.each_row(); // M x n_k
    
    // Calculate log-likelihood - group contribution
    
    ll_k = 0;
    
    if(std::strcmp(family, bin) == 0){
      arma::mat prob_mat = exp(eta_mat) / (1+exp(eta_mat)); // M x n_k
      for(i = 0; i < eta_mat.n_cols; i++){
        // Calculate individual contribution to log-likelihood
        l_i = 0;
        prob = prob_mat.col(i);
        for(s = 0; s < M; s++){ 
          l_i = l_i + R::dbinom(y_k(i), 1, prob(s), false) ;
        }
        ll_k = ll_k + log(l_i / M);
      }
      
    }else if(std::strcmp(family, pois) == 0){
      for(i = 0; i < eta_mat.n_cols; i++){
        l_i = 0;
        for(s = 0; s < M; s++){
          l_i = l_i + R::dpois(y_k(i), exp(eta_mat(s,i)), false);
        }
        ll_k = ll_k + log(l_i / M);
      }
      
    }
    
    // Add up group contributions
    ll = ll + ll_k;
    
  }
  
  return(ll);
}