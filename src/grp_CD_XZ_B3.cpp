// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utility_glm.h"
#include "utility_grpCD.h"

using namespace Rcpp;
using namespace arma;

// Change from original grp_CD_XZ_B: 
// No centering or scaling of random effects covariates
// No orthonormalization
// Treat random effect covariates as if no adjustments are needed.

// Grouped coordinate descent with both fixed (X) and random (Z) effects
// Comparied to grp_CD_XZ_A: store eta in MxN matrix

// Note: implement a check that if U matrix cols = 0 for a variable, skip that ranef coef update (?)

// [[Rcpp::export]]
arma::vec grp_CD_XZ_B3(const arma::vec& y, const arma::mat& X, const arma::mat& Z,
                      const arma::vec& group,
                      const arma::mat& u, const arma::sp_mat& J_q, arma::vec dims,
                      arma::vec beta, const arma::vec& offset,
                      const char* family, int link, int init,
                      const arma::uvec& XZ_group, arma::uvec K, // covariate group index and size of covariate groups
                      const char* penalty, arma::vec params) {

  // y = response vector
  // X = fixed effecs covariate matrix
  // Z = random effects matrix
  // group = index of observation groups
  // u = matrix of MCMC draws (cols organized same as Z)
  // J = sparse matrix such that vec(Gamma) = J * gamma (converts lower triangular matrix to a vector)
  // dims = collection of useful values
  // beta = coefficient vector (for both fixed and random effects)
  // init = 1 if initializing beta from a previous model (using previous lambda in sequence)
  //    and = 0 if initializing beta from previous fit with same lambda
  // XZ_group = grouping of fixed and random effect covariates
  // K = size of XZ_groups
  // penalty = type of penalty (lasso, MCP, SCAD)
  // params = additional penalty parameters (gamma and alpha)

  unsigned int p = dims(0); // number fixed effect covariates (ncol(X))
  int N = dims(1); // total number observations (length(y))
  int d = dims(2); // number of groups within observations
  int q = dims(3); // number of random effect covariates
  int M = dims(4); // number of MCMC draws (nrow(u))
  int J_XZ = dims(5); // number of covariate groups (in fixed and random effects)
  double conv = dims(6); // Convergence threshold
  int maxit = dims(7); // maximum number of iterations
  int J_X = XZ_group(p); // Covariate group corresponding to random intercept

  int Kj = 0; // Size of covariate group j
  
  arma::uvec n(1); // other observation counter
  int j=0; // covariate group counter
  int m=0; // MCMC draw counter
  int f=0;
  int r=0;
  int k=0; // index of observation group (from 0 to d-1)
  int iter=0; // while loop iteration
  int converged = 0; // 1 if algorithm converged, 0 otherwise

  // Define available penalties
  const char* lasso = "lasso";
  const char* mcp = "MCP";
  const char* scad = "SCAD";

  // Families with specific majorization-minimization nu values (theoretically determined)
  // definition of nu provided below
  const char* bin = "binomial";
  const char* gaus = "gaussian";

  // penalty parameters (for lasso, MCP, SCAD, elastic net)
  double lambda0 = params(0);
  double lambda1 = params(1);
  double gamma = params(2);
  double alpha = params(3);

  double v0 = 0.0; // sum(residuals^2) / (N*M)
  double v0_last = 0.0; // v0 from last iteration

  arma::mat eta(M,N); // linear predictors
  arma::vec tmp_out(M+2);
  
  double nu=0.0; // max second deriv of loss function (max possible weight)
  double zetaj_L2=0.0; // L2 norm of zetaj vector

  arma::vec beta0 = beta; // Input initial beta / beta from last round of updates
  double bj=0.0; // place-holder for value of element of beta
  arma::vec vec_ones(1); vec_ones.ones();
  arma::uvec fef_cols = as<arma::uvec>(wrap(seq(0,p-1)));
  arma::uvec ref_cols = as<arma::uvec>(wrap(seq(p, p + J_q.n_cols - 1)));

  arma::uvec col_idx(q); // columns of Z and u to choose at a time
  arma::uvec m0(1);
  int H = J_q.n_cols; // From paper, J_q
  
  arma::uvec covgroup = XZ_group.subvec(p,p+H-1) - XZ_group(p); // Group index for random effects
  arma::vec K_vec = arma::conv_to<arma::vec>::from(K); // Size of groups
  arma::vec lam(J_XZ); // lambda to use; lambda_j = lambda * sqrt(Kj) where Kj = size of jth group
  lam.subvec(0,J_X-1) = lambda0 * sqrt(K_vec.subvec(0,J_X-1));
  lam.subvec(J_X,J_XZ-1) = lambda1 * sqrt(K_vec.subvec(J_X,J_XZ-1));


  // Recall link coding (see "M_step.R" for logic):
  // 10: logit, 11: probit, 12: cloglog
  // 20: log, 30: identity, 40: inverse
  
  //-------------------------------------------------------------------------------//
  // Initialize nu = max possible weight
  //-------------------------------------------------------------------------------//

  if((std::strcmp(family, bin) == 0) & (link == 10)){
    // nu always 0.25 for binomial family with logit link
    nu = 0.25;
  }else if((std::strcmp(family, gaus) == 0) & (link == 30)){
    // nu always 1 for gaussian family with identity link
    nu = 1.0;
  }else{
    // Will use max(weights) for other family and link combinations
    // Initialize to 0
    nu = 0.0;
  }

  //-------------------------------------------------------------------------------//
  // Initialize eta matrix
  //-------------------------------------------------------------------------------//

  for(k=0;k<d;k++){

    // Rows of X and Z to use
    arma::uvec ids = find(group == (k+1));
    arma::mat Xk = X.rows(ids);
    // Index of appropriate columns of Z and u
    for(f=0;f<q;f++){
      col_idx(f) = k + f*d;
    }

    arma::mat Zk = Z.submat(ids, col_idx);
    
    for(m=0;m<M;m++){
      m0 = m;
      arma::mat A = kron(u.submat(m0,col_idx), Zk) * J_q;
      
      eta.submat(m0,ids) = trans(Xk * beta.elem(fef_cols) + A * beta.elem(ref_cols) + offset(ids));
      
    }
  }
  
  v0 = R_PosInf;

  //-------------------------------------------------------------------------------//
  // Grouped Coordinate Descent
  //-------------------------------------------------------------------------------//
  
  while((iter<maxit) & (converged==0)){

    // Add to iter
    iter = iter + 1;
    // Save last beta
    beta0 = beta;
    // Save latest v0
    v0_last = v0;

    // ----------------------------------------------------------------------------------//
    // Element-wise update of fixed effects betaj first
    // ----------------------------------------------------------------------------------//

    for(j=0; j<J_X; j++){

      // Identify covariates belonging to group j
      arma::uvec idxj = find(XZ_group == j);
      Kj = idxj.n_elem;

      if((init == 1) & (iter>=3) & (sum(beta.elem(idxj)) == 0.0)){
        // If initializing beta using beta calculated from previous lambda in sequence of lambdas,
        // update all beta a few times regardless of whether betaj == 0.
        // If beta penalized to zero in past round for same lambda, will stay zero in further rounds
        // Therefore, skip to next covariate grouping
        continue;
      }else if((init == 0) & (sum(beta.elem(idxj)) == 0.0)){
        // If beta penalized to zero in past round for same lambda, will stay zero in further rounds
        // Therefore, skip to next covariate grouping
        continue;
      }

      //----------------------------------------------------------------------------------//
      // Calculate zetaj for current group j given most recent eta update
      //----------------------------------------------------------------------------------//

      // Calculate zetaj vector (initialize as zero, then sum components of interest)
      arma::vec zetaj(Kj); zetaj.zeros();

      // Initialize mean(resid^2) as 0, then sum components of interest
      v0 = 0.0;

      // Calculate zetaj for fixed effects covariates in group j
        // zeta_fixef() function in 'utility_grpCD.cpp'
      arma::vec out = zeta_fixef(y, X, eta, idxj, family, link, nu);
      zetaj = out.subvec(0,Kj-1);
      nu = out(Kj);
      v0 = out(Kj+1);

      // resid used in above = (y-mu), need resid = (y-mu)/nu
      // Therefore, incorporate (1/nu) into zetaj and (1/nu^2) into v0 = mean(resid^2)
      // Finish calc of zetaj
      zetaj = zetaj / (N*M*nu) + beta.elem(idxj);
      v0 = v0 / (nu*nu*N*M);

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
          bj = soft_thresh(zetaj_L2*nu, lam(j)*alpha) / (nu * (1.0 + lam(j)*(1.0 - alpha)));
        }else if(std::strcmp(penalty, mcp) == 0){
          bj = MCP_soln(zetaj_L2*nu, nu, lam(j), gamma, alpha);
        }else if(std::strcmp(penalty, scad) == 0){
          bj = SCAD_soln(zetaj_L2*nu, nu, lam(j), gamma, alpha);
        }

        if(idxj.n_elem == 1){
          beta.elem(idxj) = bj * vec_ones;
        }else{
          beta.elem(idxj) = bj * (zetaj / zetaj_L2);
        }

      } // End if-else update to beta for group j


      // Update fixed effects component of eta
      arma::mat Xj = X.cols(idxj);
      arma::rowvec shift = trans(Xj * (beta.elem(idxj) - beta0.elem(idxj)));
      eta = eta.each_row() + shift;

      // Reset nu to 0 if not binomial or gaussian with canonical link
      if(!((std::strcmp(family, bin) == 0) & (link == 10)) & !((std::strcmp(family, gaus) == 0) & (link == 30))){
        nu = 0.0;
      }


    } // End j for loop for fixed effects
    
    // ----------------------------------------------------------------------------------//
    // Calculate zetaj for random intercept
    // Will not penalize random intercept
    // ----------------------------------------------------------------------------------//

    // Identify covariates belonging to random intercept (group j = J_X)
    arma::uvec idxj = find(XZ_group == J_X);
    Kj = idxj.n_elem;

    // Calculate zetaj vector (initialize as zero, then sum components of interest)
    arma::vec zetaj(Kj); zetaj.zeros();

    // Initialize mean(resid^2) as 0, then sum components of interest
    v0 = 0.0;
    
    for(k=0;k<d;k++){
      // Rows of Z to use
      arma::uvec ids = find(group == (k+1));
      // Index of appropriate columns for Z and u
      for(f=0;f<q;f++){
        col_idx(f) = k + f*d;
      }

      arma::mat Zk = Z.submat(ids, col_idx);
      arma::mat Xj(ids.n_elem, Kj); Xj.zeros();
      
      for(m=0;m<M;m++){
        m0 = m;
        arma::mat A = kron(u.submat(m0,col_idx),Zk) * J_q;
        Xj = A.cols(idxj - p);
        arma::vec tmp_out = resid_nu_v0_k(y(ids), trans(eta.submat(m0,ids)), family, link, nu);
        arma::vec resid = tmp_out.subvec(0,ids.n_elem-1);
        nu = tmp_out(ids.n_elem);
        v0 += tmp_out(ids.n_elem+1);
        zetaj = zetaj + Xj.t() * resid;

      } // End m for loop

    } // End k for loop
    
    // resid in above = (y-mu), need resid = (y-mu)/nu
    // Therefore, incorporate (1/nu) into zetaj and (1/nu^2) into v0 = mean(resid^2)
    // Finish calc of zetaj
    zetaj = zetaj / (N*M*nu) + beta.elem(idxj);
    v0 = v0 / (nu*nu*N*M);

    // L2 norm of zetaj
    if(idxj.n_elem == 1){
      zetaj_L2 = as_scalar(zetaj);
    }else{
      zetaj_L2 = sqrt(sum(zetaj % zetaj));
    }

    // Update beta (no penalization)
    beta.elem(idxj) = zetaj;

    // Reset nu to 0 if not binomial or gaussian with canonical link
    if(!((std::strcmp(family, bin) == 0) & (link == 10)) & !((std::strcmp(family, gaus) == 0) & (link == 30))){
      nu = 0.0;
    }

    r = J_X;
    
    // ----------------------------------------------------------------------------------//
    // Element-wise update of random effects effects betaj
    // ----------------------------------------------------------------------------------//

    for(j=(J_X+1); j<J_XZ; j++){

      // Identify covariates belonging to group j
      arma::uvec idxj = find(XZ_group == j);
      Kj = idxj.n_elem;

      if((init == 1) & (iter>=3) & (sum(beta.elem(idxj)) == 0.0)){
        // If initializing beta using beta calculated from previous lambda in sequence of lambdas,
        // update all beta a few times regardless of whether betaj == 0.
        // If beta penalized to zero in past round for same lambda, will stay zero in further rounds
        // Therefore, skip to next covariate grouping
        continue;
      }else if((init == 0) & (sum(beta.elem(idxj)) == 0.0)){
        // If beta penalized to zero in past round for same lambda, will stay zero in further rounds
        // Therefore, skip to next covariate grouping
        continue;
      }

      //----------------------------------------------------------------------------------//
      // Update eta from last updated beta (group r)
      // Calculate zetaj for current group j
      //----------------------------------------------------------------------------------//

      // Find index of covariates corresponding to previous j in sequence (r)
      arma::uvec idxr = find(XZ_group == r);
      
      // Calculate zetaj vector (initialize as zero, then sum components of interest)
      arma::vec zetaj(Kj); zetaj.zeros();

      // Initialize mean(resid^2) as 0, then sum components of interest
      v0 = 0.0;

      for(k=0;k<d;k++){
        // Rows of Z to use
        arma::uvec ids = find(group == (k+1));
        // Index of appropriate columns for Z and u
        for(f=0;f<q;f++){
          col_idx(f) = k + f*d;
        }

        arma::mat Zk = Z.submat(ids, col_idx);
        arma::mat Xj(ids.n_elem, Kj); Xj.zeros();

          if(H == q){ // sigma specified with independence structure (diagonal)

            for(m=0;m<M;m++){
              m0 = m;
              arma::mat A = kron(u.submat(m0, col_idx),Zk) * J_q;
              // Update random effects component of eta for each individual
              eta.submat(m0,ids) += trans(A.cols(idxr - p) * (beta.elem(idxr) - beta0.elem(idxr)));
              // Calculate zetaj
              Xj = A.cols(idxj - p);
              arma::vec tmp_out = resid_nu_v0_k(y(ids), trans(eta.submat(m0,ids)), family, link, nu);
              arma::vec resid = tmp_out.subvec(0,ids.n_elem-1);
              nu = tmp_out(ids.n_elem);
              v0 += tmp_out(ids.n_elem+1);
              zetaj = zetaj + Xj.t() * resid;
              
            }

          }else{ // H = q(q+1)/2, unstructured sigma
            
            for(m=0;m<M;m++){
              m0 = m;
              arma::mat A = kron(u.submat(m0, col_idx), Zk) * J_q;
              
              // Update random effects component of eta for each individual
              eta.submat(m0,ids) += trans(A.cols(idxr - p) * (beta.elem(idxr) - beta0.elem(idxr)));
              // Calculate zetaj
              Xj = A.cols(idxj - p);
              arma::vec tmp_out = resid_nu_v0_k(y(ids), trans(eta.submat(m0,ids)), family, link, nu);
              arma::vec resid = tmp_out.subvec(0,ids.n_elem-1);
              nu = tmp_out(ids.n_elem);
              v0 += tmp_out(ids.n_elem+1);
              zetaj = zetaj + Xj.t() * resid;
              
            } // End m for loop

          } // End if-else H == q

      } // End k for loop

      // resid in above = (y-mu), need resid = (y-mu)/nu
      // Therefore, incorporate (1/nu) into zetaj and (1/nu^2) into v0 = mean(resid^2)
      // Finish calc of zetaj
      zetaj = zetaj / (N*M*nu) + beta.elem(idxj);
      v0 = v0 / (nu*nu*N*M);

      // L2 norm of zetaj
      if(idxj.n_elem == 1){
        zetaj_L2 = as_scalar(zetaj);
      }else{
        zetaj_L2 = sqrt(sum(zetaj % zetaj));
      }

      // Update betaj
      if(std::strcmp(penalty, lasso) == 0){
        bj = soft_thresh(zetaj_L2*nu, lam(j)*alpha) / (nu * (1.0 + lam(j)*(1.0 - alpha)));
      }else if(std::strcmp(penalty, mcp) == 0){
        bj = MCP_soln(zetaj_L2*nu, nu, lam(j), gamma, alpha);
      }else if(std::strcmp(penalty, scad) == 0){
        bj = SCAD_soln(zetaj_L2*nu, nu, lam(j), gamma, alpha);
      }

      if(idxj.n_elem == 1){
        beta.elem(idxj) = bj * vec_ones;
      }else{
        beta.elem(idxj) = bj * (zetaj / zetaj_L2);
      }

      // Reset nu to 0 if not binomial or gaussian with canonical link
      if(!((std::strcmp(family, bin) == 0) & (link == 10)) & !((std::strcmp(family, gaus) == 0) & (link == 30))){
        nu = 0.0;
      }

      // Save r as the most recent updated j
      r = j;

    } // End j for loop

    // ----------------------------------------------------------------------------------//
    // Update eta with last betaj update (r = last j updated)
    // ----------------------------------------------------------------------------------//

    // Find index of covariates corresponding to previous j in sequence (r)
    arma::uvec idxr = find(XZ_group == r);
    
    // Initialize mean(resid^2) as 0, then sum components of interest
    v0 = 0.0;

    for(k=0;k<d;k++){
      // Rows of Z to use
      arma::uvec ids = find(group == (k+1));
      // Index of appropriate columns for Z and u
      for(f=0;f<q;f++){
        col_idx(f) = k + f*d;
      }

      arma::mat Zk = Z.submat(ids, col_idx);

      if(H == q){ // sigma specified with independence structure (diagonal)

        for(m=0;m<M;m++){
          m0 = m;
          arma::mat A = kron(u.submat(m0, col_idx),Zk) * J_q;
          // Update random effects component of eta for each individual
          eta.submat(m0,ids) += trans(A.cols(idxr - p) * (beta.elem(idxr) - beta0.elem(idxr)));
          // Update nu and v0
          arma::vec tmp_out = resid_nu_v0_k(y(ids), trans(eta.submat(m0,ids)), family, link, nu);
          nu = tmp_out(ids.n_elem);
          v0 += tmp_out(ids.n_elem+1);

        }

      }else{ // H = q(q+1)/2, unstructured sigma

        for(m=0;m<M;m++){
          m0 = m;
          arma::mat A = kron(u.submat(m0, col_idx), Zk) * J_q;
          
          // Update random effects component of eta for each individual
          eta.submat(m0,ids) += trans(A.cols(idxr - p) * (beta.elem(idxr) - beta0.elem(idxr)));
          // Update nu and v0
          arma::vec tmp_out = resid_nu_v0_k(y(ids), trans(eta.submat(m0,ids)), family, link, nu);
          nu = tmp_out(ids.n_elem);
          v0 += tmp_out(ids.n_elem+1);

        } // End m for loop

      } // End if-else H == q

    } // End k for loop

    // resid used above = (y-mu), true resid = (y-mu)/nu
    // In previous calculations, v0 = sum(resid^2), want = mean(resid^2)
    // should divide by (nu^2 * M * N)
    v0 = v0 / (nu*nu*M*N);

    // ----------------------------------------------------------------------------------//
    // Check Convergence Criteria
    // ----------------------------------------------------------------------------------//
   
    // Convergence criteria after every full beta update
    // if(arma::max(abs(beta0 - beta))*sqrt(nu) < conv){
    //   converged = 1;
    // }
    if(fabs(v0 - v0_last)/(v0_last+0.1) < conv*0.001){
      converged = 1;
    }

  } // End while loop
  
  Rprintf("Number iterations in grouped coord descent: %i \n", iter);

  if(converged == 0){
    // warning("grouped coordinate descent algorithm did not converge \n");
    Rcout << "grouped coordinate descent algorithm did not converge" << std::endl;
  }

  return beta;
}
