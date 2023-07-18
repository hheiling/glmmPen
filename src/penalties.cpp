// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// zeta, nu, gamma, and lambda defined in coord_desc and grp_CD functions
// contained in files "coord_descent.cpp" and "grp_CD_XZ_fast.cpp"

// Soft-thresholding function (used in Lasso, MCP, and SCAD penalties)
// [[Rcpp::export]]
double soft_thresh(double zeta, double lambda){
  
  // # Soft-thresholding function - returns a scalar;
  double abs_z = fabs(zeta);
  double val = 0;
  
  if((zeta > 0) && (lambda < abs_z)){
    val = zeta - lambda;
  }else if((zeta < 0) && (lambda < abs_z)){
    val = zeta + lambda;
  } else if(lambda >= abs_z ){
    val = 0;
  }
  
  return val;
  
}

// [[Rcpp::export]]
double MCP_soln(double zeta, double nu, double lambda, double gamma, double alpha){ 
  // gamma > 1
  
  double val = 0;
  double abs_z = fabs(zeta);
  
  double lam1 = lambda*alpha;
  double lam2 = lambda*(1.0-alpha);
  
  if(abs_z <= lam1){
    val = 0;
  }else if(abs_z <= gamma*lam1*(1.0+lam2)){
    val = soft_thresh(zeta, lam1) / (nu*(1.0 - (1.0/gamma) + lam2));
  }else{
    val = zeta / (nu*(1.0+lam2));
  }
  
  return val;
  
}

// [[Rcpp::export]]
double SCAD_soln(double zeta, double nu, double lambda, double gamma, double alpha){ 
  // gamma > 2 
  
  double val = 0;
  double abs_z = fabs(zeta);
  
  double lam1 = lambda*alpha;
  double lam2 = lambda*(1.0-alpha);
  
  if(abs_z <= lam1){
    val = 0.0;
  }else if(abs_z <= lam1*(2.0+lam2)){
    val = soft_thresh(zeta, lam1) / (nu*(1.0 +lam2));
  }else if(abs_z <= gamma*lam1*(1.0+lam2)){
    val = soft_thresh(zeta, gamma*lam1/(gamma-1.0)) / (nu * (1.0 - (1.0/(gamma-1.0)) + lam2));
  }else{
    val = zeta / (nu*(1.0+lam2));
  }
  
  return val;
  
}

////////////////////////////////////////////////////////////////////////////////////////////////
