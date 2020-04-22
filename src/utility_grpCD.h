#ifndef _COORDES_H_
#define _COORDES_H_

#include <RcppArmadillo.h>

/* penalties: lasso, MCP, SCAD */
/* penalty parameters: lambda, gamma, alpha */
double soft_thresh(double zeta, double lambda);

double MCP_soln(double zeta, double nu, double lambda, double gamma, double alpha);

double SCAD_soln(double zeta, double nu, double lambda, double gamma, double alpha);

// Other helper functions for grp_CD_XZ function

arma::vec zeta_fixef(arma::vec y, arma::mat X, arma::mat eta, 
                     arma::uvec idxr, const char* family, int link, double nu);

arma::vec resid_nu_v0(double yi, arma::vec eta, const char* family, int link, double nu);

#endif
