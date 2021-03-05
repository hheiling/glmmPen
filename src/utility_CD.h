#ifndef _COORDES_H_
#define _COORDES_H_

#include <Rcpp.h>
#include <RcppArmadillo.h>

//----------------------------------------------------------------------------------------//
// Utility functions for regular and grouped Coordinate Descent
//----------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------//
// Coordinate Descent Penalties: see "penalties.cpp" file for function details
//----------------------------------------------------------------------------------------//

/* penalties: lasso, MCP, SCAD */
/* penalty parameters: lambda, gamma, alpha */
/* See penalties.cpp for details  */
double soft_thresh(double zeta, double lambda);

double MCP_soln(double zeta, double nu, double lambda, double gamma, double alpha);

double SCAD_soln(double zeta, double nu, double lambda, double gamma, double alpha);

//----------------------------------------------------------------------------------------//
// Coordinate Descent (ungrouped case): see "coord_descent.cpp" file for function details
//----------------------------------------------------------------------------------------//

/* glm coordinate descent: lasso, MCP, SCAD */
/* penalty parameters: lambda, gamma, alpha */
arma::vec coord_desc(arma::vec y, arma::mat X, arma::vec weights, arma::vec resid, 
                     arma::vec eta, arma::vec dims, arma::vec beta, 
                     const char* penalty, double lambda, double gamma, double alpha, 
                     const char* family, int link, arma::vec penalty_factor, int trace);

//----------------------------------------------------------------------------------------//
// Helper functions for grp_CD_XZ (grouped coordinate descent)function
// See "utility_grpCD.cpp" file for function details
//----------------------------------------------------------------------------------------//

// Calculating residuals by individual using eta as input
arma::vec resid_nu_i(double yi, arma::vec eta, const char* family, int link, double nu, double phi);

// Calculates fixed effects zetaj using residuals matrix
arma::vec zeta_fixef_calc(arma::mat X, arma::mat resid, arma::uvec idxj);

#endif
