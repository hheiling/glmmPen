#ifndef _GLM_H_
#define _GLM_H_

#include <RcppArmadillo.h>

// #define BINOMIAL  1
// #define POISSON   2
// #define GAUSSIAN  3
// #define GAMMA     4
// #define NB        5
// #define ZINB      6

/* Link */

#define LOGIT     10
#define PROBIT    11
#define CLOGLOG   12
#define LOG       20
#define IDENTITY  30
#define INVERSE   40

/* GLM definition functions */

arma::vec initial_mu(const char* family, arma::vec y, int N);
arma::vec muvalid(const char* family, arma::vec mu);
arma::vec mu_adjust(const char* family, arma::vec mu);
arma::vec dlink(int link, arma::vec mu);
arma::vec linkfun(int link, arma::vec mu);
arma::vec invlink(int link, arma::vec eta);
arma::vec varfun(const char* family, arma::vec mu);

// Q-function calculation
double Qfun(const char* family, double yi, arma::vec mu);


/* Fit a base model */

// int glmFit(int* familyR, int* linkR, int* dims, int* nIter,
//            double *y, double *offset, double *z,
//            double *X, double *nTotal_binom, double *convR, 
//            int *rank, double *Xb, double *fitted, double *resid, 
//            double *weights, double *phi, int* trace, 
//            double *scale, int *df_resid, double* beta);

// /*  log likelihood of Poisson */
// double loglik_Poisson(int N, double* mu, double* y);
// 
// /* log likelihood of negative binomial */
// 
// double loglik_NB(int N, double phi, double* mu, double* y,double *prior);
// 
// /* Score test for additional terms */
// 
// void glm_score_test(int* dims, double *Z, double *resid, 
//                     double *weights, double *Xb, double* scaleR,
//                     double *chi2, int *df);
// 
// /* score and infor for solving MLE of phi */
// 
// void score_info(int N, double theta, double* mu, double *y, 
//                 double* score, double* info, double *prior);
// 
// /* MLE of phi */
// 
// int phi_ml(double* y, double* mu, int N,  
//            int limit, double eps, double* phi, int initPhi, int trace, double *prior);
// 
// /* glmNB */
// 
// int glmNB(int *dims, int *nIter, double *y, double *z, 
//           int *linkR, double *offset, double *X, double *convR, 
//           int *rank, double *Xb, double *fitted, double *resid, 
//           double *weights, double *phi, double *scale, 
//           int *df_resid, int* family, double *twologlik, 
//           double *scoreTestP, int *trace, double *beta);
// 
// int glmNBlog(int *dimsNew, int *nIter, double *pY, double *z, 
//              int *linkR, double *offset, double *pX, double *conv, 
//              double *convGLM, int *rank, double *Xb, double *fitted, 
//              double *resid, double *weights, double *phi, double *scale, 
//              int *df_resid, int* family, double *twoLL_trec, 
//              double *scoreTestP, int *trace, double *beta,
//              double *fitted2, double *offsetN);

/* glmIAL*/


#endif
