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
// Helper functions for grp_CD_XZ (grouped coordinate descent) function
// See "utility_grpCD.cpp" file for function details
//----------------------------------------------------------------------------------------//

// Calculating residuals by individual using eta as input - grp_CD_XZ_fast function
arma::vec resid_nu_i(double yi, arma::vec eta, const char* family, int link, double nu, double phi);

// Calculating residuals by individual using eta as input - grp_CD_XZ_step function
arma::vec resid_i(double yi, arma::vec eta, const char* family, int link);

// Calculates fixed effects zetaj using residuals matrix
arma::vec zeta_fixef_calc(arma::mat X, arma::mat resid, arma::uvec idxj);

//----------------------------------------------------------------------------------------//
// Q-function calculations: needed for step-size adjustment calculations
// See "utility_glm.cpp" file for function details
//----------------------------------------------------------------------------------------//
// Calculates Q-function (positive version of "Q1" expression in Rashid et al. (2020))
// for original glmm() and glmmPen() implementation.
// See function in "utility_glm.cpp"
double Qfun(const arma::vec& y, const arma::mat& X, const arma::mat& Z, SEXP pBigMat, 
            const arma::vec& group, const arma::sp_mat& J_q,
            const arma::vec& beta, const arma::vec offset, arma::vec dims,
            const char* family, int link, double sig_g, double phi);

// Calculates Q-function (positive version of "Q1" expression in Rashid et al. (2020))
// for glmm_FA() and glmmPen_FA() implementation
// See function in "utility_FA.cpp"
double Qfun_FA(const arma::vec& y, const arma::mat& X, const arma::mat& Z, SEXP pBigMat, 
            const arma::vec& group, const arma::sp_mat& J_q,
            const arma::vec& beta, const arma::vec offset, arma::vec dims,
            const char* family, int link, double sig_g, double phi);


// Calculates quadratic approximation to the Q-function for step-size calculation purposes
// About beta
double Qfun_quad_beta(double Q0, double step_size, const arma::mat& diff0,
                      const arma::mat& eta, const arma::mat& eta0,
                      const arma::vec& beta, const arma::vec& beta0);

#endif
