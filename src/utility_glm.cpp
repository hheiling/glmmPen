#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>

using namespace Rcpp;

// Helper functions for IRLS, Q-function calcuation (glmm() and glmmPen() implementation),
// and Gaussian sigma estimate (for E-step)


arma::vec initial_mu(const char* family, arma::vec y, int N) {
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  const char* negbin = "negbin";
  const char* gaus = "gaussian";
  const char* Gamma = "Gamma";
  
  int i;
  
  arma::vec mu(N);
  
  if(std::strcmp(family, bin) == 0){
    for (i=0; i<N; i++) {
      if (y(i) < 0) {
        stop("negative values not allowed for the Binomial family");
      }
      if (y(i) > 1.0) {
        stop("# of success is larger than 1");
      }
      mu(i) = (y(i) + 0.5)/(1.0+1.0);
    }
  }else if((std::strcmp(family, pois) == 0) || (std::strcmp(family, negbin) == 0)){
    for (i=0; i<N; i++) {
      if (y(i) < 0) {
        stop("negative values not allowed for the Poisson family");
      }      
      mu(i) = y(i) + 0.1;
    }
  }else if(std::strcmp(family, gaus) == 0){
    for (i=0; i<N; i++) {
      mu(i) = y(i);
    }
  }else if(std::strcmp(family, Gamma) == 0){
    for (i=0; i<N; i++) {
      if (y(i) <= 0) {
        stop("non-poistive values not allowed for the Gamma family");
      }      
      mu(i) = y(i) + 0.1;
    }
  }
  
  return mu;
}


arma::vec muvalid(const char* family, arma::vec mu) {
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  const char* negbin = "negbin";
  const char* gaus = "gaussian";
  const char* Gamma = "Gamma";
  
  double minb = 0.0001; // minimum allowed binomial mu value
  double maxb = 0.9999; // maximum allowed binomial mu value
  double minp = 0.0001; // minimum allowed poisson and negbin mu value
  double gammaMin = 0.001; // miminum allowed gamma mu value
  
  int i = 0;
  int N = mu.n_elem;
  arma::vec valid(N);
  
  if(std::strcmp(family, bin) == 0){
    for(i=0; i<N; i++){
      valid(i) = (mu(i) > minb && mu(i) < maxb);
    }
  }else if((std::strcmp(family, pois) == 0) || (std::strcmp(family, negbin) == 0)){
    for(i=0; i<N; i++){
      valid(i) = (mu(i) > minp);
    }
  }else if(std::strcmp(family, gaus) == 0){
    for(i=0; i<N; i++){
      valid(i) = 1;
    }
  }else if(std::strcmp(family, Gamma) == 0){
    for(i=0; i<N; i++){
      valid(i) = (mu(i) > gammaMin);
    }
  }else{
    stop("invalid family \n");
  }
  
  return(valid);
  
}


arma::vec mu_adjust(const char* family, arma::vec mu) {
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  const char* negbin = "negbin";
  const char* gaus = "gaussian";
  const char* Gamma = "Gamma";
  
  double minb = 0.001; // minimum allowed binomial mu value
  double maxb = 0.999; // maximum allowed binomial mu value
  double minp = 0.001; // minimum allowed poisson and negbin mu value
  double gammaMin = 0.001; // miminum allowed gamma mu value
  
  int i = 0;
  int N = mu.n_elem;
  arma::vec mu_new = mu;
  
  if(std::strcmp(family, bin) == 0){
    for(i=0; i<N; i++){
      if(mu(i) < minb){
        mu_new(i) = minb;
      }else if(mu(i) > maxb){
        mu_new(i) = maxb;
      }
    }
  }else if((std::strcmp(family, pois) == 0) || (std::strcmp(family, negbin) == 0)){
    for(i=0; i<N; i++){
      if(mu(i) < minp){
        mu_new(i) = minp;
      }
    }
  }else if(std::strcmp(family, Gamma) == 0){
    for(i=0; i<N; i++){
      if(mu(i) < gammaMin){
        mu_new(i) = gammaMin;
      }
    }
  }else if(std::strcmp(family, gaus) == 0){ 
    // No invalid mu
    mu_new = mu;
  }else{
    stop("invalid family \n");
  }
  
  // Note: gaussian family does not have any invalid mu
  
  return(mu_new);
  
}

// general equation:
// g(mu) = eta <---> mu = g^(-1)(eta)

// dlink: take the derivative of the link function with respect to mu 
// (d eta / d mu) = d (g(mu)) / d mu
arma::vec dlink(int link, arma::vec mu){
  
  int N = mu.n_elem;
  int i = 0;
  arma::vec out(N);
  arma::vec ones_vec = out.ones();
  arma::vec zeros_vec = out.zeros();
  double tmp;
  
  // If probit link
  if(link == 11){
    out.zeros();
    for(i=0;i<N;i++){
      tmp = R::qnorm5(mu(i), 0.0, 1.0, 1, 0); // lower.tail = T, log = F
      out(i) = 1.0 / R::dnorm4(tmp, 0.0, 1.0, 0); // log = F
    }
  }
  
  switch(link){
  case 10: return(ones_vec / (mu % (ones_vec - mu))); // logit
  case 11: return(out); // probit
  case 12: return(ones_vec/(log(ones_vec-mu)%(ones_vec-mu))); // cloglog
  case 20: return(ones_vec/mu); // log
  case 21: return((0.5 * ones_vec)/sqrt(mu)); // sqrt
  case 30: return(ones_vec); // identity
  case 31: return(zeros_vec - ones_vec/(mu%mu)); // inverse
  default: return(zeros_vec); 
  }
  
}

// linkfun: return g(mu) to get eta estimate
arma::vec linkfun(int link, arma::vec mu){
  
  int N = mu.n_elem;
  int i=0;
  arma::vec out(N);
  arma::vec ones_vec = out.ones();
  // double tmp;
  
  if(link == 11){ // probit
    out.zeros();
    for(i=0;i<N;i++){
      out(i) = R::qnorm5(mu(i), 0.0, 1.0, 1, 0); // lower.tail = T, log = F
    }
  }
  
  switch(link){
  case 10: return(log(mu/(ones_vec-mu))); // logit
  case 11: return(out); // probit
  case 12: return(log(-ones_vec % log(ones_vec-mu))); // cloglog
  case 20: return(log(mu)); // log
  case 21: return(sqrt(mu)); // sqrt
  case 30: return(mu); // identity
  case 31: return(ones_vec / mu); // inverse
  default: return(out.zeros());
  }
  
}

// invlink: return g^(-1)(eta) to get mu estimate
// [[Rcpp::export]]
arma::vec invlink(int link, arma::vec eta){
  
  int N = eta.n_elem;
  int i=0;
  arma::vec out(N);
  arma::vec ones_vec = out.ones();
  arma::vec zeros_vec = out.zeros();
  
  if(link == 11){ // probit
    out.zeros();
    for(i=0;i<N;i++){
      out(i) = R::pnorm5(eta(i), 0.0, 1.0, 1, 0); // lower.tail = T, log = F
    }
  }
  
  switch(link){
  case 10: return(exp(eta) / (ones_vec + exp(eta))); // logit
  case 11: return(out); // probit 
  case 12: return(ones_vec - exp(zeros_vec - ones_vec % exp(eta))); // cloglog
  case 20: return(exp(eta)); // log
  case 21: return(eta%eta); // sqrt
  case 30: return(eta); // identity
  case 31: return(zeros_vec - ones_vec / eta); // inverse
  default: return(zeros_vec);
  }
  
}

// varfun: return the variance function of the family
arma::vec varfun(const char* family, arma::vec mu, double phi){
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  const char* negbin = "negbin"; 
  const char* gaus = "gaussian";
  const char* Gamma = "Gamma";
  
  int N = mu.n_elem;
  arma::vec V(N);
  arma::vec empty(N);
  arma::vec ones_vec = empty.ones();
  
  if(std::strcmp(family, bin) == 0){
    V = mu%(ones_vec-mu);
  }else if(std::strcmp(family, pois) == 0){
    V = mu;
  }else if(std::strcmp(family, negbin) == 0){
    V = mu + phi * (mu % mu);
  }else if(std::strcmp(family, gaus) == 0){
    V = ones_vec;
  }else if(std::strcmp(family, Gamma) == 0){
    V = mu%mu;
  }else{
    stop("invalid family \n");
  }
  
  return(V);
  
}

// Q function estimate (positive version of estimate)
// [[Rcpp::export]]
double Qfun(const arma::vec& y, const arma::mat& X, const arma::mat& Z, SEXP pBigMat, 
            const arma::vec& group, const arma::sp_mat& J_q,
            const arma::vec& beta, const arma::vec offset, arma::vec dims,
            const char* family, int link, double sig_g, double phi){
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  const char* negbin = "negbin";
  const char* gaus = "gaussian";
  const char* Gamma = "Gamma";
  
  // Provide access to the big.matrix of posterior draws
  XPtr<BigMatrix> pMat(pBigMat);
  arma::Mat<double> post((double*) pMat->matrix(), pMat->nrow(), pMat->ncol(),false);
  
  int M = pMat->nrow();
  int N = y.n_elem;
  int p = X.n_cols; // Number fixed effects covariates
  int p2 = beta.n_elem; // Total number of coefficients 
  
  int m=0;
  arma::uvec m0(1);
  int i=0;
  int k=0;
  int f=0;
  int n_k=0;
  
  int d = dims(0); // Number groups
  int q = dims(1); // Number random effect variables
  
  arma::uvec col_idx(q);
  arma::vec const_ones(N); const_ones.ones();
  
  arma::mat eta(M,N);
  arma::vec mu(M);
  
  arma::mat u(M,q);
  arma::rowvec Zki(q);
  arma::mat A(M,J_q.n_cols);
  
  // parameters needed for gaussian distribution
  // double s2 = sig_g * sig_g; // estimate of sigma^2 
  
  double llQ = 0.0;
  
  // Calculate MxN eta matrix
  for(k=0;k<d;k++){
    
    // Rows of X and Z to use
    arma::uvec ids = find(group == (k+1));
    arma::mat Xk = X.rows(ids);
    // Index of appropriate columns of Z and u
    for(f=0;f<q;f++){
      col_idx(f) = k + f*d;
    }
    
    arma::mat Zk = Z.submat(ids, col_idx);
    n_k = ids.n_elem;
    
    u = post.cols(col_idx);
    
    for(i=0;i<n_k;i++){
      Zki = Zk.row(i);
      A = kron(u, Zki) * J_q;
      eta.col(ids(i)) = as_scalar(Xk.row(i) * beta.subvec(0,p-1)) + A * beta.subvec(p,p2-1) + offset(ids(i));
    }
    
  } // End k for loop
  
  
  // Calculate llQ (The Q1 portion of the Q function specified in Section 4 of Rashid et al paper)
  
  llQ = 0.0;
  
  for(i=0;i<N;i++){
    
    mu = invlink(link, eta.col(i));
    mu = mu_adjust(family, mu);
    
    if(std::strcmp(family, bin) == 0){
      for(m=0;m<M;m++){
        llQ = llQ + R::dbinom(y(i), 1.0, mu(m), 1); // log value
      }
    }else if(std::strcmp(family, pois) == 0){
      for(m=0;m<M;m++){
        llQ = llQ + R::dpois(y(i), mu(m), 1); // log value
      }
    }else if(std::strcmp(family, negbin) == 0){
      for(m=0;m<M;m++){
        llQ = llQ + R::dnbinom(y(i), (1.0/phi), mu(m), 1); // log value
      }
    }else if(std::strcmp(family, gaus) == 0){
      for(m=0;m<M;m++){
        llQ = llQ + R::dnorm(y(i), mu(m), sig_g, 1); // log value
      }
    }else if(std::strcmp(family, Gamma) == 0){
      stop("Gamma not yet available");
    }else{
      stop("invalid family \n");
    }
    
  }
  
  llQ = llQ / M;
  
  return(llQ);

}

// estimate of standard deviation for linear regression error term
// [[Rcpp::export]]
double sig_gaus(const arma::vec& y, const arma::mat& X, const arma::mat& Z, SEXP pBigMat, 
                const arma::vec& group, const arma::sp_mat& J_q,
                const arma::vec& beta, const arma::vec offset, arma::vec dims, int link){
  
  // Provide access to the big.matrix of posterior draws
  XPtr<BigMatrix> pMat(pBigMat);
  arma::Mat<double> post((double*) pMat->matrix(), pMat->nrow(), pMat->ncol(),false);
  
  int M = pMat->nrow();
  int N = y.n_elem;
  int p = X.n_cols; // Number fixed effects covariates
  int p2 = beta.n_elem; // Total number of coefficients (p+H)
  int H = J_q.n_cols; // From paper, J_q
  
  int i=0;
  // int m=0;
  arma::uvec m0(1);
  int k=0;
  int f=0;
  
  int d = dims(0); // Number groups
  int q = dims(1); // Number random effect variables
  int n_k=0; // Number of observations in group k
  
  arma::uvec col_idx(q);
  arma::vec const_ones(M); const_ones.ones();
  
  arma::mat eta(M,N);
  arma::vec mu(M);
  arma::mat u(M,q);
  arma::rowvec Zki(q);
  arma::mat A(M,H);
  
  // parameters needed for gaussian distribution
  double s2 = 0.0; // estimate of sigma^2 
  double sig = 0.0; // estimate of standard deviation
  
  for(k=0;k<d;k++){
    
    // Rows of X and Z to use
    arma::uvec ids = find(group == (k+1));
    arma::mat Xk = X.rows(ids);
    // Index of appropriate columns of Z and u
    for(f=0;f<q;f++){
      col_idx(f) = k + f*d;
    }
    
    arma::mat Zk = Z.submat(ids, col_idx);
    n_k = ids.n_elem;
    
    u = post.cols(col_idx);
    
    for(i=0;i<n_k;i++){
      Zki = Zk.row(i);
      A = kron(u, Zki) * J_q;
      eta.col(ids(i)) = as_scalar(Xk.row(i) * beta.subvec(0,p-1)) + A * beta.subvec(p,p2-1) + offset(ids(i));
    }
    
  } // End k for loop
  
  // If gaussian family, calculate sigma^2 
  s2 = 0.0;
  for(i=0;i<N;i++){
    mu = invlink(link, eta.col(i));
    s2 += sum((y(i)*const_ones - mu) % (y(i)*const_ones - mu));
  }
  
  sig = sqrt(s2 / (M*N));
  
  return(sig);
  
}

/**********************************************************************
 *
 * Generalized from MASS/negbin.R 
 * (generalized from all N independent observations to N*M total observations)
 * 
 
 Package: MASS
Priority: recommended
Version: 7.3-5
Date: 2010-01-03
Depends: R (>= 2.10.1), grDevices, graphics, stats, utils
Suggests: lattice, nlme, survival
Author: S original by Venables & Ripley. R port by Brian Ripley
<ripley@stats.ox.ac.uk>, following earlier work by Kurt Hornik
and Albrecht Gebhardt.
Maintainer: Brian Ripley <ripley@stats.ox.ac.uk>
Description: Functions and datasets to support Venables and Ripley,
'Modern Applied Statistics with S' (4th edition).
Title: Main Package of Venables and Ripley's MASS
License: GPL-2 | GPL-3
URL: http://www.stats.ox.ac.uk/pub/MASS4/
LazyLoad: yes
LazyData: yes
Packaged: 2010-01-03 10:50:27 UTC; ripley
Repository: CRAN
Date/Publication: 2010-01-03 14:05:40

**********************************************************************/

/**********************************************************************
*
* score_info
*
* score and Fisher information, i.e., the first and negative second 
* derivative of likelihood of phi 
*
**********************************************************************/

void score_info(double theta, arma::mat eta, arma::vec y, int link,
                double* score, double* info){
  
  int N = eta.n_cols;
  int M = eta.n_rows;
  int i=0;
  int m=0;
  arma::vec mu(M);
  double score1=0.0, info1=0.0;
  double mui, yi, scorei, infoi, thMui;
  
  for (i=0; i<N; i++) {
    
    mu = invlink(link, eta.col(i));
    if(mu.has_inf()){
      Rprintf("mu has at least one inf value in individual %i \n", i);
    }
    yi = y(i);
    
    for(m=0; m<M; m++){
      
      mui = mu(m);
      
      thMui   = theta + mui;
      scorei  = R::digamma(yi + theta) - R::digamma(theta) - (theta + yi)/thMui;
      score1 += (scorei - log(thMui) + 1.0 + log(theta));
      
      infoi   = R::trigamma(theta) - R::trigamma(yi + theta) + (mui - yi)/(thMui*thMui);
      info1  += (infoi + 1/thMui - 1/theta);
    }
    
  }
  
  *score = score1;
  *info  = info1;
  
} // End score_info

void score_info_init(double theta, arma::vec mu, arma::vec y, int link,
                     double* score, double* info){
  
  int N = y.n_elem;
  int i=0;
  double score1=0.0, info1=0.0;
  double mui, yi, scorei, infoi, thMui;
  
  for (i=0; i<N; i++) {
    
    yi = y(i);
    mui = mu(i);
    
    thMui   = theta + mui;
    scorei  = R::digamma(yi + theta) - R::digamma(theta) - (theta + yi)/thMui;
    score1 += (scorei - log(thMui) + 1 + log(theta));
    
    infoi   = R::trigamma(theta) - R::trigamma(yi + theta) + (mui - yi)/(thMui*thMui);
    info1  += (infoi + 1/thMui - 1/theta);
    
  }
  
  *score = score1;
  *info  = info1;
  
} // End score_info_init

// score_info function from hmmcov package
// void score_info(int N, double theta, double* mu, double* y, 
//                 double* score, double* info, double *prior)
// {
//   int i;
//   double score1=0.0, info1=0.0;
//   double mui, yi, scorei, infoi, thMui;
//   
//   for (i=0; i<N; i++) {
//     yi  = y[i];
//     mui = mu[i];
//     
//     thMui   = theta + mui;
//     scorei  = digamma(yi + theta) - digamma(theta) - (theta + yi)/thMui;
//     score1 += prior[i]*(scorei - log(thMui) + 1 + log(theta));
//     
//     infoi   = trigamma(theta) - trigamma(yi + theta) + (mui - yi)/(thMui*thMui);
//     info1  += prior[i]*(infoi + 1/thMui - 1/theta);
//   }
//   
//   *score = score1;
//   *info  = info1;
// }


/**********************************************************************
*
* phi_ml
*
* MLE of phi (over-dispersion parameter), given mu 
*
* Actually we find MLE of 1/phi here and then take inverse
*
**********************************************************************/

// phi_ml used within M step after phi initialized at start of algorithm
// [[Rcpp::export]]
double phi_ml(arma::vec y, arma::mat eta, int link, int limit, double eps, double phi){
  
  // int N = y.n_elem;
  
  double theta0, del;
  double score=0.0;
  double info=0.0;
  // double n=0;
  int it=0;
  double minTheta = 1e-5;
  double maxTheta = 1.0/minTheta;
  // int fail = 0;
  
  double phi_new;
  
  if((phi > minTheta) && (phi < maxTheta)){
    theta0 = 1.0/(phi);
  }else if(phi < minTheta){
    theta0 = maxTheta;
  }else{
    theta0 = minTheta;
  }
  
  it  = 0;
  del = 1.0;
  
  while((it < limit) && (fabs(del) > eps)) {
    score_info(theta0, eta, y, link, &score, &info);
    Rcout << "score: " << score << std::endl;
    Rcout << "info: " << info << std::endl;
    del     = score/info;
    theta0 += del;
    it     += 1;
    
    if (theta0 > maxTheta) {
      // phi truncated at 1/maxTheta, no overdispersion
      theta0 = maxTheta;
    }
    
    if(theta0 < minTheta) {
      // phi truncated at 1/minTheta, too much overdispersion
      theta0 = minTheta;
    }
    
  }
  
  Rcout << "theta0: " << theta0 << std::endl;
  
  if(it == limit) {
    Rprintf("  phi.ml: iteration limit reached in phi_ml\n");
  }
  
  phi_new = 1.0/theta0;
  
  return(phi_new);
  
} // End phi_ml function

// phi_ml_init used at start of algorithm to produce initial phi estimate using fixed effects only
// [[Rcpp::export]]
double phi_ml_init(arma::vec y, arma::vec eta, int link, int limit, double eps){
  
  int N = y.n_elem;
  
  double theta0, del, tmp;
  double score=0.0;
  double info=0.0;
  // double n=0;
  int i, it=0;
  double minTheta = 1e-5;
  double maxTheta = 1.0/minTheta;
  // int fail = 0;
  
  arma::vec mu(N);
  
  double phi=0.0;
  
  mu = invlink(link, eta);
  
  theta0 = 0.0;
  
  for (i=0; i<N; i++) {
    tmp = y(i)/mu(i) - 1.0;
    theta0 += tmp*tmp;
  }
  
  // From theta.ml function of MASS package, n = number of observations 
  theta0 = N/theta0;
  
  it  = 0;
  del = 1.0;
  
  while(it < limit && fabs(del) > eps) {
    score_info_init(theta0, mu, y, link, &score, &info);
    del     = score/info;
    theta0 += del;
    it     += 1;
    
    if (theta0 > maxTheta) {
      // phi truncated at 1/maxTheta, no overdispersion
      theta0 = maxTheta;
    }
    
    if(theta0 < minTheta) {
      // phi truncated at 1/minTheta, too much overdispersion
      theta0 = minTheta;
    }
    
  }
  
  if(it == limit) {
    Rprintf("  phi.ml: iteration limit reached in phi_ml \n");
  }
  
  phi = 1/theta0;
  
  return(phi);
  
} // End phi_ml_init function

// phi_ml function for N independent observations
// double phi_ml(arma::vec y, arma::vec mu, int N, int limit, double eps, 
//               double phi, int initPhi, int trace){
//   
//   double theta0, del, tmp;
//   double score=0.0;
//   double info=0.0;
//   double n=0;
//   int i, it=0;
//   double minTheta = 1e-5;
//   double maxTheta = 1.0/minTheta;
//   int tryPoisson = 0;
//   int tryZINB = 0;
//   int fail = 0;
//   double phi_new;
//   
//   //To try later:  If theta previous hit boundary, then set init=0
//   //This is to help when theta is too large on initial penalizations
//   //if(initPhi & *phi != minTheta & *phi != maxTheta){
//   if(initPhi){
//     theta0 = 1.0/(phi);
//   }else{
//     theta0 = 0.0;
//     for (i=0; i<N; i++) {
//       tmp = y[i]/mu[i] - 1.0;
//       theta0 += tmp*tmp;
//     }
//     
//     // From theta.ml function of MASS package, n = number of observations or the sum of the weights
//     // (not coordinate descent weights)
//     theta0 = N/theta0; 
//   }
//   
//   it  = 0;
//   del = 1.0;
//   
//   if(trace > 5) Rprintf("  phi.ml: initial phi = %.2e\n", 1/theta0);
//   
//   while(it < limit && fabs(del) > eps) {
//     score_info(N, theta0, mu, y, &score, &info);
//     del     = score/info;
//     theta0 += del;
//     it     += 1;
//     
//     if (theta0 > maxTheta) {
//       // phi truncated at 1/maxTheta, no overdispersion
//       theta0 = maxTheta;
//       
//     }
//     
//     if(theta0 < minTheta) {
//       // phi truncated at 1/minTheta, too much overdispersion
//       theta0 = minTheta;
//       
//     }
//     
//   }
//   
//   if(it == limit) {
//     if(trace > 3)
//       Rprintf("  phi.ml: iteration limit reached in phi_ml\n");
//   }
//   
//   phi_new = 1/theta0;
//   
//   return(phi_new);
//   
// } // End phi_ml function

// original phi_ml from hmmcov package
// int phi_ml(double* y, double* mu, int N, int limit, double eps, 
//            double* phi, int initPhi, int trace, double *prior)
// {
//   double theta0, del, tmp;
//   double score=0.0;
//   double info=0.0;
//   double n=0;
//   int i, it=0;
//   double minTheta = 1e-5;
//   double maxTheta = 1.0/minTheta;
//   int tryPoisson = 0;
//   int tryZINB = 0;
//   int fail = 0;
//   
//   //To try later:  If theta previous hit boundary, then set init=0
//   //This is to help when theta is too large on initial penalizations
//   //if(initPhi & *phi != minTheta & *phi != maxTheta){
//   if(initPhi){
//     theta0 = 1.0/(*phi);
//   }else{
//     theta0 = 0.0;
//     for (i=0; i<N; i++) {
//       tmp = y[i]/mu[i] - 1.0;
//       theta0 += tmp*tmp;
//     }
//     
//     for(i=0;i<N;i++) n=n+prior[i];
//     theta0 = n/theta0;
//   }
//   
//   it  = 0;
//   del = 1.0;
//   
//   if(trace > 5) Rprintf("  phi.ml: initial phi = %.2e\n", 1/theta0);
//   
//   while(it < limit && fabs(del) > eps) {
//     score_info(N, theta0, mu, y, &score, &info, prior);
//     del     = score/info;
//     theta0 += del;
//     it     += 1;
//     
//     if(trace > 5) Rprintf("  phi.ml: iter %d, phi=%.2e, score=%.2e, info=%.2e\n", 
//        it,  1/theta0, score, info);
//     
//     if (theta0 > maxTheta) {
//       theta0 = maxTheta;
//       if(trace > 3)
//         Rprintf("    phi is truncated at %.2e, no overDispersion?\n", 1/maxTheta);
//       
//       //tryPoisson = 1;
//       //break;
//     }
//     
//     if(theta0 < minTheta) {
//       theta0 = minTheta;
//       if(trace > 3)
//         Rprintf("    phi is truncated at %.2e, too much overDispersion?\n", 1/minTheta);
//       
//       //tryZINB = 1;
//       //break;
//     }
//     
//   }
//   
//   if(it == limit) {
//     //fail = 1;
//     if(trace > 3)
//       Rprintf("  phi.ml: iteration limit reached in phi_ml\n");
//   }
//   
//   *phi = 1/theta0;
//   
//   return(tryPoisson + 2*tryZINB + 4*fail);
// }


////////////////////////////////////////////////////////////////////////////////////////////////
