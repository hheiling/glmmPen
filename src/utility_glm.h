#ifndef _GLM_H_
#define _GLM_H_

#include <RcppArmadillo.h>

//----------------------------------------------------------------------------------------//
// Utility functions for IRLS (used in grouped and regular Coordinate Descent) and
// functions related to glm fit (see "utility_glm.cpp" for further details)
//----------------------------------------------------------------------------------------//

/* Link */

// #define LOGIT     10
// #define PROBIT    11
// #define CLOGLOG   12
// #define LOG       20
// #define IDENTITY  30
// #define INVERSE   40

/* GLM definition functions */

arma::vec initial_mu(const char* family, arma::vec y, int N);
arma::vec muvalid(const char* family, arma::vec mu);
arma::vec mu_adjust(const char* family, arma::vec mu);
arma::vec dlink(int link, arma::vec mu);
arma::vec linkfun(int link, arma::vec mu);
arma::vec invlink(int link, arma::vec eta);
arma::vec varfun(const char* family, arma::vec mu, double phi);

/* Negative Binomial helper functions */

/* score and info for solving MLE of phi (called within phi_ml) */

void score_info(double theta, arma::mat eta, arma::vec y, int link,
                double* score, double* info);
 
/* MLE of phi */

double phi_ml(arma::vec y, arma::mat eta, int link, int limit, double eps, 
              double phi);



#endif
