
#include <Rcpp.h>
#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

// crossprod and maxprod code copied and adapted from ncvreg code from
// src/ncvreg_init.c and src/maxprod.c, respectively

double crossprod(double *X, double *y, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i];
  return(val);
}

// [[Rcpp::export]]
double maxprod(SEXP X_, SEXP y_, SEXP v_, SEXP m_,
               int n, int p) {
  
  // Declarations
  // SEXP z = PROTECT(z = allocVector(REALSXP, 1));
  // REAL(z)[0] = 0;
  double z=0;
  double zz;
  double *X = REAL(X_);
  double *y = REAL(y_);
  double *m = REAL(m_);
  int *v = INTEGER(v_);
  
  for (int j=0; j<p; j++) {
    zz = crossprod(X, y, n, v[j]-1) / m[v[j]-1];
    if (fabs(zz) > z) z = fabs(zz);
    // if (fabs(zz) > REAL(z)[0]) REAL(z)[0] = fabs(zz);
  } 
  
  // Return list
  return(z);
}
