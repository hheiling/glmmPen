#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>
#include "utility_glm.h"
#include "utility_CD.h"

using namespace Rcpp;
using namespace arma;

// IRLS coordinate descent update for the factor analysis assumption model, where
//    random effects gamma (of dim q) can be expressed as B * alpha (B matrix qxr, 
//    alpha of dim r where r = number of latent common factors)
// Update the residuals as specified in the grpreg and ncvreg paper (to speed things up)
// Assume no orthogonalization/standardization needed for random effects

// Grouped coordinate descent with both fixed (X) and random (Z) effects

// Majorization-minimization approximation for the second derivative of the loss function
//    (a.k.a. max weight from IRLS weights), represented by "nu" in code
// Binomial with canonical (logit) link: nu = 0.25 
// Gaussian family with canonical (identity) link: nu = 1.0
// All other cases: nu = 1 / step_size, where step_size is determined by a line search
// algorithm. Reference: Book "Proximal Algorithms" by Parihk and Boyd (2013),
// see Section 4.2

// [[Rcpp::export]]
arma::vec grp_CD_XZ_FA_step(const arma::vec& y, const arma::mat& X, const arma::mat& Z,
                    const arma::vec& group, 
                    SEXP pBigMat, const arma::sp_mat& J_f, arma::vec dims,
                    arma::vec beta, const arma::vec& offset,
                    double step_size, double sig_g,
                    const char* family, int link, int init, double phi,
                    const arma::uvec& X_group, arma::uvec K, // covariate group index and size of covariate groups
                    const char* penalty, arma::vec params, int trace) {
  
  
  // y = response vector
  // X = fixed effecs covariate matrix
  // Z = random effects matrix
  // group = index of observation groups
  // pBigMat = big.matrix of MCMC draws (cols organized same as Z)
  // J = sparse matrix such that vec(Gamma) = J * gamma (converts lower triangular matrix to a vector)
  // dims = collection of useful values
  // beta = coefficient vector (for both fixed and random effects)
  // init = 1 if initializing beta from a previous model (using previous lambda in sequence)
  //    and = 0 if initializing beta from previous fit with same lambda
  // X_group = grouping of fixed effect covariates
  // K = size of X_groups
  // penalty = type of penalty (lasso, MCP, SCAD)
  // params = additional penalty parameters (gamma and alpha)
  
  // Provide access to the big.matrix of posterior draws
  XPtr<BigMatrix> pMat(pBigMat);
  arma::Mat<double> post((double*) pMat->matrix(), pMat->nrow(), pMat->ncol(),false);
  
  unsigned int p = dims(0); // number fixed effect covariates (ncol(X))
  int N = dims(1); // total number observations (length(y))
  int d = dims(2); // number of groups within observations
  int q = dims(3); // number of random effect covariates
  int M = dims(4); // number of MCMC draws (nrow(u))
  int J_X = dims(5); // number of covariate groups in fixed effects
  double conv = dims(6); // Convergence threshold
  int maxit = dims(7); // maximum number of iterations
  
  int Kj = 0; // Size of covariate group j
  
  int r = post.n_cols / d; // number common factors
  
  arma::uvec n(1); // other observation counter
  int i=0;
  int j=0; // covariate group counter
  int f=0;
  int k=0; // index of observation group (from 0 to d-1)
  int iter=0; // while loop iteration
  int converged = 0; // 1 if algorithm converged, 0 otherwise
  int n_k=0; // Number of individuals in group k
  
  int H = q*r;
  
  // Define available penalties
  const char* lasso = "lasso";
  const char* mcp = "MCP";
  const char* scad = "SCAD";
  
  // Families with specific majorization-minimization nu values (theoretically determined)
  // definition of nu provided below
  const char* bin = "binomial";
  const char* gaus = "gaussian";
  const char* negbin = "negbin";
  
  // penalty parameters (for lasso, MCP, SCAD, elastic net)
  double lambda0 = params(0);
  double lambda1 = params(1);
  double gamma = params(2);
  double alpha = params(3);
  
  // Values used to update step size in Poisson case, or generic cases outside of 
  // Binomial or Gaussian families with canonical links
  double Q0=0; // Q-function approximation estimate (evaluated at previous beta0), gives "Q1" function listed in Rashid et al. (2020) paper
  double Q=0; // Q-function approximation estimate (evaluated at current beta), gives "Q1" function listed in Rashid et al. (2020) paper
  double Q_quad=0; // Quadratic approximation to the Q-function based on Taylor series expansion
  arma::mat eta0(M,N); // linear predictors from previous iterations
  arma::mat diff0(M,N); // (y - mu) residuals from previous iterations
  int check_conv = 0;
  int update = 0;
  
  arma::mat eta(M,N); // linear predictors
  arma::mat resid(M,N); // residuals
  arma::vec tmp_out(M+1);
  
  double nu=0.0; // max second deriv of loss function (max possible weight)
  double zetaj_L2=0.0; // L2 norm of zetaj vector
  
  arma::vec beta0 = beta; // Input initial beta / beta from last round of updates
  double bj=0.0; // place-holder for value of element of beta
  arma::vec vec_ones(1); vec_ones.ones();
  arma::uvec fef_cols = as<arma::uvec>(wrap(seq(0,p-1)));
  arma::uvec ref_cols = as<arma::uvec>(wrap(seq(p, p + q*r - 1)));
  
  arma::uvec col_idxz(q); // columns of Z to choose at a time
  arma::uvec col_idxu(r); // columns of u to choose at a time
  arma::uvec m0(1);
  // int H = J_f.n_cols; // From paper, J_f
  arma::mat u(M,q);
  arma::rowvec Zki(q);
  
  arma::mat A(M,H);
  // arma::mat A_c(M,H);
  arma::rowvec shift_fef(N);
  arma::vec shift_ref(M);
  
  arma::vec zeta_ref(H);
  arma::uvec idx_ref(H); // index values for all random effects
  for(f=0;f<H;f++){
    idx_ref(f) = p+f;
  }
  
  arma::uvec idx_bt(r); // index for coefficients associated with t^th row of B matrix
  
  // arma::uvec covgroup = XZ_group.subvec(p,p+H-1) - XZ_group(p); // Group index for random effects
  // Make sure K is a vector only corresponding to fixed effects
  arma::vec K_vec = arma::conv_to<arma::vec>::from(K); // Size of groups
  arma::vec lam0(J_X); // lambda to use for fixed effects; lambda_j = lambda * sqrt(Kj) where Kj = size of jth group
  arma::vec lam1(q); // lambda to use for random effects; lambda_j = lambda * sqrt(r) where r = length of row of B matrix (number common factors)
  lam0 = lambda0 * sqrt(K_vec);
  for(j=0; j<q; j++){
    lam1(j) = lambda1 * sqrt(r); // group size will always be r (row of B matrix)
  }
  
  // Recall link coding:
  // 10: logit, 11: probit, 12: cloglog
  // 20: log, 30: identity, 40: inverse
  
  //-------------------------------------------------------------------------------//
  // Initialize nu = max possible weight
  //-------------------------------------------------------------------------------//
  
  if((std::strcmp(family, bin) == 0) && (link == 10)){
    // nu always 0.25 for binomial family with logit link
    nu = 0.25;
  }else if((std::strcmp(family, gaus) == 0) && (link == 30)){
    // nu always 1 for gaussian family with identity link
    nu = 1.0;
  }else{
    // Will search for appropriate step size. Value nu = 1 / step_size,
    // initialize with input step_size value
    nu = 1 / step_size;
  }
  
  //-------------------------------------------------------------------------------//
  // Initialize eta and residual matrices
  //-------------------------------------------------------------------------------//
  
  // Initialize eta
  for(k=0;k<d;k++){
    
    // Rows of X and Z to use
    arma::uvec ids = find(group == (k+1));
    arma::mat Xk = X.rows(ids);
    // Index of appropriate columns of Z and u
    for(f=0;f<q;f++){
      col_idxz(f) = k + f*d;
    }
    for(f=0;f<r;f++){
      col_idxu(f) = k + f*d;
    }
    
    arma::mat Zk = Z.submat(ids, col_idxz);
    n_k = ids.n_elem;
    
    u = post.cols(col_idxu);
    
    for(i=0;i<n_k;i++){
      Zki = Zk.row(i);
      A = kron(u, Zki) * J_f;
      eta.col(ids(i)) = as_scalar(Xk.row(i) * beta.elem(fef_cols)) + A * beta.elem(ref_cols) + offset(ids(i));
    }
    
  } // End k for loop
  
  // Calculate initial residuals, determine initial nu
  for(i=0;i<N;i++){
    resid.col(i) = resid_i(y(i), eta.col(i), family, link);
  }
  
  // Save (y - mu) residuals
  diff0 = resid;
  
  // Note: resid in the above calculation is (y - mu), but really want (y-mu)/nu
  // Therefore, divide entire residual matrix by nu
  resid = resid / nu;
  
  
  // Calculation of Q-function using initialized coefficients
  // See Qfun_FA() in "utility_FA.cpp"
  arma::vec dims_Qcalc = arma::vec(3);
  dims_Qcalc(0) = d;
  dims_Qcalc(1) = q;
  dims_Qcalc(2) = r;
  Q0 = (-1.0) * Qfun_FA(y, X, Z, pBigMat, group, J_f, beta, offset, dims_Qcalc, family, link, sig_g, phi);
  // Save most recent eta (this terms will be updated throughout algorithm)
  eta0 = eta;
  
  //-------------------------------------------------------------------------------//
  // Grouped Coordinate Descent
  //-------------------------------------------------------------------------------//
  
  while((iter<maxit) & (converged==0)){
    
    // Add to iter
    iter = iter + 1;
    // Save last beta
    beta0 = beta;
    
    // ----------------------------------------------------------------------------------//
    // Element-wise update of fixed effects betaj first
    // ----------------------------------------------------------------------------------//
    
    for(j=0; j<J_X; j++){
      
      // Identify covariates belonging to group j
      arma::uvec idxj = find(X_group == j);
      Kj = idxj.n_elem;
      
      if((init == 1) && (iter>=3) && (sum(beta.elem(idxj)) == 0.0)){
        // If initializing beta using beta calculated from previous lambda in sequence of lambdas,
        // update all beta a few times regardless of whether betaj == 0.
        // If beta penalized to zero in past round for same lambda, will stay zero in further rounds
        // Therefore, skip to next covariate grouping
        continue;
      }else if((init == 0) && (sum(beta.elem(idxj)) == 0.0)){
        // If beta penalized to zero in past round of EM iteration for same lambda, 
        // will stay zero in further rounds of EM algorithm
        // Therefore, skip to next covariate grouping
        continue;
      }
      
      //----------------------------------------------------------------------------------//
      // Calculate zetaj for current group j given most recent eta update
      //----------------------------------------------------------------------------------//
      
      // Calculate zetaj vector (initialize as zero, then sum components of interest)
      arma::vec zetaj(Kj); zetaj.zeros();
      
      // Calculate zetaj for fixed effects covariates in group j
      // zeta_fixef_calc() function in 'utility_grpCD.cpp'
      zetaj = zeta_fixef_calc(X, resid, idxj);
      
      // Finish calc of zetaj
      zetaj = zetaj / (N*M) + beta.elem(idxj);
      
      // L2 norm of zetaj
      if(idxj.n_elem == 1){
        zetaj_L2 = as_scalar(zetaj);
      }else{
        zetaj_L2 = sqrt(sum(zetaj % zetaj));
      }
      
      // Update betaj
      if(j==0){
        // j = 0: fixed effects intercept, other fixed effects covariates given XZ_group values of 0
        // No penalization for above scenarios
        beta.elem(idxj) = zetaj;
      }else{
        if(std::strcmp(penalty, lasso) == 0){
          bj = soft_thresh(zetaj_L2*nu, lam0(j)*alpha) / (nu * (1.0 + lam0(j)*(1.0 - alpha)));
        }else if(std::strcmp(penalty, mcp) == 0){
          bj = MCP_soln(zetaj_L2*nu, nu, lam0(j), gamma, alpha);
        }else if(std::strcmp(penalty, scad) == 0){
          bj = SCAD_soln(zetaj_L2*nu, nu, lam0(j), gamma, alpha);
        }
        
        if(idxj.n_elem == 1){
          beta.elem(idxj) = bj * vec_ones;
        }else{
          beta.elem(idxj) = bj * (zetaj / zetaj_L2);
        }
        
      } // End if-else update to beta for group j
      
      
      // Update residuals and eta matrices
      arma::mat Xj = X.cols(idxj); 
      shift_fef = trans(Xj * (beta.elem(idxj) - beta0.elem(idxj)));
      resid = resid.each_row() - shift_fef;
      eta = eta.each_row() + shift_fef;
      
    } // End j for loop for fixed effects
    
    // ----------------------------------------------------------------------------------//
    // Calculate zetaj for random effect coefficients
    // Will not penalize coefficients corresponding to first row of B matrix 
    //    (always want to have a non-zero variance for the random intercept)
    // ----------------------------------------------------------------------------------//
    
    // Calculate zeta_ref vector (initialize as zero, then sum components of interest)
    // This will be used to identify zetaj for all random effects
    zeta_ref.zeros();
    
    for(k=0;k<d;k++){
      // Rows of Z to use
      arma::uvec ids = find(group == (k+1));
      // Index of appropriate columns for Z and u
      for(f=0;f<q;f++){
        col_idxz(f) = k + f*d;
      }
      for(f=0;f<r;f++){
        col_idxu(f) = k + f*d;
      }
      
      arma::mat Zk = Z.submat(ids, col_idxz);
      n_k = ids.n_elem;
      arma::mat Xj(n_k, Kj); Xj.zeros();
      
      u = post.cols(col_idxu);
      
      for(i=0;i<n_k;i++){
        Zki = Zk.row(i);
        A = kron(u, Zki) * J_f;
        zeta_ref = zeta_ref + A.t() * resid.col(ids(i));
      } // End i for loop
      
    } // End k for loop
    
    // Finish calc of zetaj
    zeta_ref = zeta_ref / (N*M) + beta.elem(idx_ref);
    
    // ----------------------------------------------------------------------------------//
    // Element-wise update of random effects effects betaj
    // First: update random intercept (no penalization on random intercept)
    // Then update all other random effects.
    // ----------------------------------------------------------------------------------//
    
    // Identify covariates belonging to random intercept 
    for(f=0;f<r;f++){
      idx_bt(f) = f;
    }
    arma::vec zetaj = zeta_ref(idx_bt);
    
    // Update beta (no penalization)
    beta.elem(idx_bt+p) = zetaj;
    
    // Identify index and zetaj for each of the remaining random effects
    // Update betaj
    
    for(j=1; j<q; j++){ // Start at j=1 because we already worked on j=0 (first row of B matrix)
      
      // Identify covariates belonging to group j
      // arma::uvec idxj = find(XZ_group == j);
      for(f=0;f<r;f++){
        idx_bt(f) = f + r*j;
      }
      
      if((sum(beta.elem(idx_bt+p)) == 0.0)){
        // Within EM algorithm for single penalty parameter combination: 
        // If coefficients penalized to zero in past round of EM iteration for same lambda
        // or set to zero in a past iteration within a single M-step,
        // will stay zero in further rounds of EM algorithm and all further M-step iterations.
        // Therefore, skip to next covariate grouping
        
        // Start of EM algorithm:
        // The initialized coefficients for each row of the B matrix will be set to 0 if 
        // there is reason to ignore these random effects (see "select_tune_FA.R", 
        // could be due to pre-screening step or past model fit ...) and set to non-zero
        // values if there is reason to consider these random effects in the model.
        
        // Consequently, if row of B is set to 0 in any case, keep this row as 0 in all
        // further iterations of the M-step
        continue;
      }
      
      // Identify correct elements of zeta_ref for zetaj
      arma::vec zetaj = zeta_ref.elem(idx_bt);
      
      zetaj_L2 = sqrt(sum(zetaj % zetaj));
      
      // Update betaj
      if(std::strcmp(penalty, lasso) == 0){
        bj = soft_thresh(zetaj_L2*nu, lam1(j)*alpha) / (nu * (1.0 + lam1(j)*(1.0 - alpha)));
      }else if(std::strcmp(penalty, mcp) == 0){
        bj = MCP_soln(zetaj_L2*nu, nu, lam1(j), gamma, alpha);
      }else if(std::strcmp(penalty, scad) == 0){
        bj = SCAD_soln(zetaj_L2*nu, nu, lam1(j), gamma, alpha);
      }
      
      beta.elem(idx_bt+p) = bj * (zetaj / zetaj_L2);
      
    } // End j for loop (end update to random effects coefficients)
    
    
    // ----------------------------------------------------------------------------------//
    // Update eta with last betaj random effect update 
    // Recalculate residuals using updated eta
    // Recalculate nu (if necessary)
    // ----------------------------------------------------------------------------------//
    
    // Update eta
    for(k=0;k<d;k++){
      // Rows of Z to use
      arma::uvec ids = find(group == (k+1));
      // Index of appropriate columns for Z and u
      for(f=0;f<q;f++){
        col_idxz(f) = k + f*d;
      }
      for(f=0;f<r;f++){
        col_idxu(f) = k + f*d;
      }
      
      arma::mat Zk = Z.submat(ids, col_idxz);
      n_k = ids.n_elem;
      
      u = post.cols(col_idxu);
      
      for(i=0;i<n_k;i++){
        Zki = Zk.row(i);
        A = kron(u, Zki) * J_f;
        // Update random effects component of eta for each individual
        eta.col(ids(i)) += A * (beta.elem(idx_ref) - beta0.elem(idx_ref));
      } // End i for loop
      
    } // End k for loop
    
    // Re-set nu if necessary
    // If exponential family NOT Binomial or Gaussian with canonical link,
    // update step size (aka update nu) and determine if appropriate to 
    // update parameters and check convergence.
    if((std::strcmp(family, bin) == 0) && (link == 10)){
      // nu always 0.25 for binomial family with logit link
      nu = 0.25;
      check_conv = 1;
      update = 1;
    }else if((std::strcmp(family, gaus) == 0) && (link == 30)){
      // nu always 1 for gaussian family with identity link
      nu = 1.0;
      check_conv = 1;
      update = 1;
    }else{
      // Will search for appropriate step size, then calculate nu = 1 / step_size
      // 1. (a) Calculate Q-function using updated coefficients (see Qfun_FA in "utility_FA.cpp")
      Q = (-1.0) * Qfun_FA(y, X, Z, pBigMat, group, J_f, beta, offset, dims_Qcalc, family, link, sig_g, phi);
      // 1. (b) Calculate quadratic approximation to the Q-function based on Taylor series
      // expansion around the beta0 coefficients, see Qfun_quad_beta in "utility_glm.cpp"
      Q_quad = Qfun_quad_beta(Q0, step_size, diff0, eta, eta0, beta, beta0);
      // 2. If Q-function with new coefficients greater than quadratic approximation, reduce step size
      if(Q > Q_quad){ 
        // Reduce step size
        step_size = step_size * 0.95;
        // Do no update coefficients, start loop over with original results
        update = 0;
        // Do not evaluate convergence
        check_conv = 0;
      }else{
        update = 1;
        check_conv = 1;
      }
      // Update nu = 1 / step_size
      nu = 1.0 / (step_size);
    }
    
    if(update){
      // Re-calculate residuals from updated eta
      for(i=0;i<N;i++){
        resid.col(i) = resid_i(y(i), eta.col(i), family, link);
      }
      // Save the (y-mu) residuals for this set of updated coefficients
      diff0 = resid;
      // Note: resid in the above calculation is (y - mu), but really want (y-mu)/nu
      // Therefore, divide entire residual matrix by nu
      resid = resid / nu;
      // Save Q0 and eta0 evaluated at most recent beta update
      Q0 = Q;
      eta0 = eta;
    }else{
      // Discard newly updated beta coefficients and corresponding eta and residuals, start again
      beta = beta0;
      resid = diff0 / nu; 
      eta = eta0;
    }
    
    // ----------------------------------------------------------------------------------//
    // Check Convergence Criteria
    // ----------------------------------------------------------------------------------//
    
    // Convergence criteria after every full beta update
    if(check_conv){
      if(arma::max(abs(beta0 - beta)) < conv){
        converged = 1;
      }
    }
    
  } // End while loop
  
  if(trace >= 1){
    Rprintf("Number iterations in grouped coord descent: %i \n", iter);
    Rcout << "final nu: " << nu << std::endl;
  }
  
  // If negative binomial family, update phi (dispersion) estimate
  // negbin variance: mu + phi * mu^2
  // double phi_ml(arma::vec y, arma::mat eta, int link, int limit, double eps, double phi);
  if(std::strcmp(family, negbin) == 0){
    phi = phi_ml(y, eta, link, 150, conv, phi);
    Rcout << "phi: " << phi << std::endl;
  }
  
  if((converged == 0) && (trace >= 1)){
    Rcout << "grouped coordinate descent algorithm did not converge" << std::endl;
  }
  
  // Return output vector: beta, step_size, and phi 
  arma::vec output(p+H+2);
  output.subvec(0,p+H-1) = beta;
  output(p+H) = step_size;
  output(p+H+1) = phi;
  
  return output;
  
}

