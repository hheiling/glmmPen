#ifndef _COORDES_H_
#define _COORDES_H_

#include <Rcpp.h>
#include <RcppArmadillo.h>

// Coordinate Descent Penalties: see "coord_descent.cpp" file for function details
/* penalties: lasso, MCP, SCAD */
/* penalty parameters: lambda, gamma, alpha */
double soft_thresh(double zeta, double lambda);

double MCP_soln(double zeta, double nu, double lambda, double gamma, double alpha);

double SCAD_soln(double zeta, double nu, double lambda, double gamma, double alpha);

// Other helper functions for grp_CD_XZ function
// See "utility_grpCD.cpp" file for function details

arma::vec zeta_fixef(arma::vec y, arma::mat X, arma::mat eta, 
                     arma::uvec idxr, const char* family, int link, double nu);

arma::vec resid_nu_v0_k(arma::vec y, arma::vec eta, const char* family, int link, double nu);

// Call to Metropolis-within-Gibbs Adaptive Random Walk 
// Creates draws from posterior, see"sample_mc_random_walk.cpp" file for function details

// arma::mat posterior(arma::vec yk, arma::mat Xk, arma::mat Zk, arma::vec beta_fef,
//                     int q, int M, arma::mat sigma, arma::vec uold, 
//                     arma::rowvec proposal_SD_k, int batch, int batch_length,
//                     int offset, int burnin_batchnum, int trace);

// Rcpp::NumericMatrix sample_mc_gibbs_adapt_rw(arma::mat f, arma::mat z, arma::vec y, 
//                                              arma::vec t, int NMC, arma::vec u0, 
//                                              arma::rowvec proposal_SD, double batch,
//                                              double batch_length, double offset,
//                                              double burnin_batchnum, int trace); 

#endif
