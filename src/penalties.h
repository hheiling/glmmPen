#ifndef _COORDES_H_
#define _COORDES_H_

#include <RcppArmadillo.h>

/* penalties: lasso, MCP, SCAD */
/* penalty parameters: lambda, gamma, alpha */
double soft_thresh(double zeta, double lambda);

double MCP_soln(double zeta, double nu, double lambda, double gamma, double alpha);

double SCAD_soln(double zeta, double nu, double lambda, double gamma, double alpha);

#endif
