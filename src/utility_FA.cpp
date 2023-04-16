#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>
#include "utility_glm.h"

using namespace Rcpp;
using namespace arma;

// Helper functions for FA (factor analysis) implementation used in glmm_FA and glmmPen_FA



// Calculate Q function estimate (assuming random effects follow a factor analysis method)
// Q function used in BIC-ICQ calculation, and possibly M-step (Poisson, eventually non-canonical links)
// For families that are NOT Cox Proportional Hazards
// [[Rcpp::export]]
double Qfun_FA(const arma::vec& y, const arma::mat& X, const arma::mat& Z, SEXP pBigMat, 
                const arma::vec& group, const arma::sp_mat& J_f,
                const arma::vec& beta, const arma::vec offset, arma::vec dims,
                const char* family, int link, double sig_g, double phi){
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  // const char* negbin = "negbin";
  const char* gaus = "gaussian";
  // const char* Gamma = "Gamma";
  
  // Provide access to the big.matrix of posterior draws
  XPtr<BigMatrix> pMat(pBigMat);
  arma::Mat<double> post((double*) pMat->matrix(), pMat->nrow(), pMat->ncol(),false);
  
  int M = pMat->nrow();
  int N = y.n_elem;
  int p = X.n_cols; // Number fixed effects covariates
  
  int m=0;
  int i=0;
  int k=0;
  int f=0;
  int n_k;
  
  int d = dims(0); // Number groups
  int q = dims(1); // Number random effect variables
  int r = dims(2); // Number common effects
  
  arma::uvec col_idxz(q);
  arma::uvec col_idxu(r);
  
  arma::mat eta(M,N);
  arma::vec mu(M);
  
  arma::uvec fef_cols = as<arma::uvec>(wrap(seq(0,p-1)));
  arma::uvec ref_cols = as<arma::uvec>(wrap(seq(p, p + q*r - 1)));
  
  arma::mat u(M,r);
  arma::rowvec Zki(q);
  arma::mat A(M,q*r);
  
  double llQ = 0.0;
  
  // Calculate MxN eta matrix
  for(k=0;k<d;k++){
    
    // Rows of X and Z to use
    arma::uvec ids = find(group == (k+1));
    arma::mat Xk = X.rows(ids);
    // Index of appropriate columns of Z and u
    for(f=0;f<q;f++){
      col_idxz(f) = k + f*d;
    }
    for(f=0;f<r;f++){
      col_idxu(f) = k + f*d;
    }
    
    arma::mat Zk = Z.submat(ids, col_idxz);
    n_k = ids.n_elem;
    
    u = post.cols(col_idxu);
    
    for(i=0;i<n_k;i++){
      Zki = Zk.row(i);
      A = kron(u, Zki) * J_f;
      eta.col(ids(i)) = as_scalar(Xk.row(i) * beta.elem(fef_cols)) + A * beta.elem(ref_cols) + offset(ids(i));
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
    }else if(std::strcmp(family, gaus) == 0){
      for(m=0;m<M;m++){
        llQ = llQ + R::dnorm(y(i), mu(m), sig_g, 1); // log value
      }
    }else{
      stop("invalid family \n");
    }
    
    // To add later when package can handle other families
    // else if(std::strcmp(family, negbin) == 0){
    //   for(m=0;m<M;m++){
    //     llQ = llQ + R::dnbinom(y(i), (1.0/phi), mu(m), 1); // log value
    //   }
    // }else if(std::strcmp(family, Gamma) == 0){
    //   stop("Gamma not yet available");
    // }else{
    //   stop("invalid family \n");
    // }
    
  }
  
  llQ = llQ / M;
  
  return(llQ);
  
}

// estimate of standard deviation for linear regression error term
// [[Rcpp::export]]
double sig_gaus_FA(const arma::vec& y, const arma::mat& X, const arma::mat& Z, SEXP pBigMat,
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
  int r = dims(2); // Number of latent common factors
  int n_k=0; // Number of observations in group k

  arma::uvec col_idxz(q);
  arma::uvec col_idxu(r);
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
      col_idxz(f) = k + f*d;
    }
    for(f=0;f<r;f++){
      col_idxu(f) = k + f*d;
    }
    
    arma::mat Zk = Z.submat(ids, col_idxz);
    n_k = ids.n_elem;
    
    u = post.cols(col_idxu);

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


//////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////

// Calculate initial eta, an M x N matrix of linear predictor values
//    N = number total subjects, M = number of posterior draws taken during the previous E-step
// arma::mat eta_zeta_calc_FA(arma::mat& X, arma::mat& Z, arma::vec& group, arma::vec& offset,
//                           arma::vec& beta, arma::vec& beta0, arma::mat& eta0,
//                           const arma::sp_mat& J, SEXP pBigMat,
//                           int N, int M, int d, int q, int r, int p,
//                           int init){
//   
//   // init: 1 if initializing eta at start of M-step algorithm,
//   //    0 if updating random effects portion of eta (linear predictor) throughout M-step
//   // beta: if init = 1, initialized beta values; if init = 0, most recent beta values in algorithm
//   // beta0: if init = 0, beta values before most recent beta update
//   
//   XPtr<BigMatrix> pMat(pBigMat);
//   arma::Mat<double> post((double*) pMat->matrix(), pMat->nrow(), pMat->ncol(),false);
//   
//   int k=0; // Index for number of groups
//   int i=0; // Index for number of subjects
//   int f=0; // Index for columns of Z matrix or matrix of posterior draws
//   
//   int n_k=0; // Number of observations in group k
//   
//   // Indices for coefficients: fixed effects coefficients, random effects covariance matrix coefficients
//   arma::uvec fef_cols = as<arma::uvec>(wrap(seq(0,p-1)));
//   arma::uvec ref_cols = as<arma::uvec>(wrap(seq(p, p + q*r - 1)));
//   
//   arma::vec col_idxz(q); // Column indices for Z matrix
//   arma::vec col_idxu(r); // Column indices for posterior draws matrix u
//   
//   int H = q*r;
//   arma::mat A(M,H);
//   arma::mat u(M,q);
//   arma::rowvec Zki(q);
//   
//   arma::mat eta(M,N);
//   
//   for(k=0;k<d;k++){
//     
//     // Rows of X and Z to use
//     arma::uvec ids = find(group == (k+1));
//     arma::mat Xk = X.rows(ids);
//     // Index of appropriate columns of Z and u
//     for(f=0;f<q;f++){
//       col_idxz(f) = k + f*d;
//     }
//     for(f=0;f<r;f++){
//       col_idxu(f) = k + f*d;
//     }
//     
//     arma::mat Zk = Z.submat(ids, col_idxz);
//     n_k = ids.n_elem;
//     
//     u = post.cols(col_idxu);
//     
//     for(i=0;i<n_k;i++){
//       Zki = Zk.row(i);
//       A = kron(u, Zki) * J;
//       if(init == 1){ // If initializing eta matrix at start of M-step
//         eta.col(ids(i)) = as_scalar(Xk.row(i) * beta.elem(fef_cols)) + A * beta.elem(ref_cols) + offset(ids(i));
//       }else if(init == 0){ // If updating random effects portion of eta matrix throughout M-step
//         // Update random effects component of eta for each individual
//         eta.col(ids(i)) = eta0.col(ids(i)) + A * (beta.elem(ref_cols) - beta0.elem(ref_cols));
//       }
//      
//     }
//     
//   } // End k for loop
//   
//   return(eta);
// }
