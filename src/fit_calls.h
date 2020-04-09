#ifndef _FIT_H_
#define _FIT_H_

#include <RcppArmadillo.h>

/* glm coordinate descent: lasso, MCP, SCAD */
/* penalty parameters: lambda, gamma, alpha */
arma::vec coord_desc(arma::vec y, arma::mat X, arma::vec weights, arma::vec resid, 
                     arma::vec eta, arma::vec offset, arma::vec dims, arma::vec beta, 
                     const char* penalty, double lambda, double gamma, double alpha, // penalty type and parameters
                     const char* family, int link);

#endif
