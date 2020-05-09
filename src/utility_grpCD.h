#ifndef _COORDES_H_
#define _COORDES_H_

#include <Rcpp.h>
#include <RcppArmadillo.h>

//----------------------------------------------------------------------------------------//
// Coordinate Descent Penalties: see "coord_descent.cpp" file for function details
//----------------------------------------------------------------------------------------//

/* penalties: lasso, MCP, SCAD */
/* penalty parameters: lambda, gamma, alpha */
double soft_thresh(double zeta, double lambda);

double MCP_soln(double zeta, double nu, double lambda, double gamma, double alpha);

double SCAD_soln(double zeta, double nu, double lambda, double gamma, double alpha);

//----------------------------------------------------------------------------------------//
// Other helper functions for grp_CD_XZ function
// See "utility_grpCD.cpp" file for function details
//----------------------------------------------------------------------------------------//

// Calculates zetaj for fixed effects starting with eta
arma::vec zeta_fixef(arma::vec y, arma::mat X, arma::mat eta, 
                     arma::uvec idxr, const char* family, int link, double nu);

// Calculating residuals by group (used for zetaj calc for random effects covariates)
arma::vec resid_nu_v0_k(arma::vec y, arma::vec eta, const char* family, int link, double nu);

// Calculating residuals by individual using eta as input
arma::vec resid_nu_v0_i(double yi, arma::vec eta, const char* family, int link, double nu);

// Calculates fixed effects zetaj using residuals matrix
arma::vec zeta_fixef_calc(arma::mat X, arma::mat resid, arma::uvec idxj);

#endif
