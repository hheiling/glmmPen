// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utility_glm.h"
#include "utility_grpCD.h"

using namespace Rcpp;

// Grouped coordinate descent with both fixed (X) and random (Z) effects
// Comparied to grp_CD_XZ_A: store eta in MxN matrix

// Note: implement a check that if U matrix cols = 0 for a variable, skip that ranef coef update (?)

//' @export
// [[Rcpp::export]]
arma::vec grp_CD_XZ_B(const arma::vec& y, const arma::mat& X, const arma::mat& Z, 
                      const arma::vec& group,
                      const arma::mat& u, const arma::sp_mat& J, const arma::vec& dims,
                      arma::vec beta, const arma::vec& offset,
                      const char* family, int link, int init,
                      const arma::uvec& XZ_group, arma::uvec K, // covariate group index and size of covariate groups
                      const char* penalty, double lambda0, double lambda1, arma::vec params) {
  
  // y = response vector
  // X = fixed effects covariate matrix
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
  
  int p = dims(0); // number fixed effect covariates (ncol(X))
  int N = dims(1); // total number observations (length(y))
  int d = dims(2); // number of groups within observations
  int q = dims(3); // number of random effect covariates
  int M = dims(4); // number of MCMC draws (nrow(u))
  int J_XZ = dims(5); // number of covariate groups (in fixed and random effects)
  double conv = dims(6); // Convergence threshold
  int maxit = dims(7); // maximum number of iterations
  
  int Kj = 0; // Size of covariate group j
  int Kr = 0; // Size of covariate group r
  
  int i=0; // observation counter
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
  
  // penalty parameters (MCP, SCAD, elastic net)
  double gamma = params(0);
  double alpha = params(1);
  
  double v0 = 0.0; // sum(residuals^2) / N
  double v0_last = 0.0; // v0 from last iteration
  
  arma::mat eta(M,N); // linear predictors
  arma::vec tmp_out(M+2);
  
  double nu=0.0; // max second deriv of loss function (max possible weight)
  double zj=0.0; // 1/N * t(X_j) %*% resid + beta0_j if can be made into scalar
  double zetaj_L2=0.0; // L2 norm of zetaj vector
  
  arma::vec beta0 = beta; // Input initial beta / beta from last round of updates
  double bj=0.0; // place-holder for value of element of beta
  arma::vec vec_ones(1); vec_ones.ones();
  arma::uvec fef_cols = seq(0,p-1);
  arma::uvec ref_cols = seq(p, p + J.n_cols - 1);
  
  int grp = 0; // observation group
  arma::uvec col_idx(q); // columns of Z and u to choose at a time
  arma::uvec m0(1);
  int H = J.n_cols; // From paper, J_q
  
  arma::uvec K_sub = K.subvec(p,p+H-1) - K(p); // Size of random effects groups
  // lambda_j = lambda * sqrt(Kj) where Kj = size of jth group
  arma::vec K_vec = arma::conv_to<arma::vec>::from(K); // Size of groups
  arma::vec lam(p+H); // lambda to use
  lam.subvec(0,p-1) = lambda0 * sqrt(K_vec.subvec(0,p-1));
  lam.subvec(p,p+H-1) = lambda1 * sqrt(K_vec.subvec(p,p+H-1));
  
  
  // Recall link coding (see "M_step.R" for logic):
  // 10: logit, 11: probit, 12: cloglog
  // 20: log, 30: identity, 40: inverse
  
  //-------------------------------------------------------------------------------//
  // Standardization and Orthogonalization
  //-------------------------------------------------------------------------------//
  
  // Standardize A matrix (from paper, matrix of (alpha_k kronecker Zki)T J_q)
  // Find eigenvectors and eigenvalues of t(A_j) * (A_j) for jth group
  // Note: Standardization sufficient if specified sigma matrix as diagonal
  
  arma::rowvec center(H); center.zeros();
  arma::rowvec scale(H); scale.zeros();
  
  for(k=0;k<d;k++){ // For each of d groups of individuals
    
    // Rows of Z to use
    arma::uvec ids = find(group == k);
    
    for(f=0;f<q;f++){
      col_idx(f) = k + f*d;
    }
    arma::mat Zk = Z.submat(ids, col_idx);
    
    for(m=0;m<M;m++){
      m0 = m;
      arma::mat A = kron(u.submat(m0,col_idx),Zk) * J; // n_k by H
      center = center + sum(A,0); // Add up column sums
      scale = scale + sum(A%A,0);
    }
    
  }
  
  // If ncols(J) > q then orthogonalize. Not necessary otherwise.
    // If sigma matrix diagonal, then coefficients relating to random effects not grouped
    // and orthogonalization not necessary. If sigma matrix diagonal, ncol(J) == q.
    // Otherwise, ncol(J) == q(q+1)/2
    
  List ortho_factor;
  
  if(H > q){
    
    // Initialize list for t(A_j) * A_j matrices
    List L;
    
    for(f=0;f<q;f++){
      arma::mat xTx(f+1,f+1); 
      if(f == 0){
        xTx.ones(); // Standardization sufficient for group of size 1, no additional orthogonalization necessary
      }else{
        xTx.zeros();
      }
      L[f] = xTx;
    }
    
    for(k=0;k<d;k++){
      // Rows of X and Z to use
      arma::uvec ids = find(group == k);
      
      // Index of columns of Z and u to use
      for(f=0;f<q;f++){
        col_idx(f) = k + f*d;
      }
      
      arma::mat Zk = Z.submat(ids,col_idx);
      
      for(m=0;m<M;m++){
        arma::mat A = kron(u.row(m), Zk) * J;
        arma::mat A_std = (A.each_col() - center) / scale;
        
        // Calculate t(A_j) * A_j for groups
        for(f=1;f<q;f++){
          arma::mat xTx = L[f];
          arma::uvec Kidx = find(K_sub == f);
          xTx = xTx + trans(A_std.cols(Kidx)) * A_std.cols(Kidx);
          L[f] = xTx;
        }
      }
    }
    
    // Divide t(Xj) * Xj by (N*M)
    for(f=1;f<q;f++){
      arma::mat xTx = L[f];
      xTx = xTx / (N*M);
      L[f] = xTx;
    }
    
    // Find eigenvectors and eigenvalues of the xTx matrices, create matrix used to orthonormalize
    // Note: Incorporate checks to make sure eigen values (evec) not equal to 0
    
    for(f=0;f<q;f++){
      if(f == 0){
        ortho_factor[f] = 1.0;
      }else{
        arma::mat xTx = L[f];
        arma::mat emat;
        arma::vec evec;
        arma::eig_sym(evec, emat, xTx);
        ortho_factor[f] = emat.each_col() % (1/sqrt(evec));
      }
    }
    
  }
  
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
    arma::uvec ids = find(group == k);
    arma::mat Xk = X.rows(ids);
    // Index of appropriate columns of Z and u
    for(f=0;f<q;f++){
      col_idx = k + f*d;
    }
    
    arma::mat Zk = Z.submat(ids, col_idx);
    
    for(m=0;m<M;m++){
      m0 = m;
      arma::mat A = kron(u.row(m), Zk) * J;
      arma::mat A_std = (A.each_col() - center) / scale;
      
      if(H == q){ // sigma specified with independence structure (diagonal)
        // Caculate eta for each individual using initial beta, weight by 1/M
        eta.submat(m0,ids) = trans(Xk * beta.elem(fef_cols) + A_std * beta.elem(ref_cols) + offset(ids)) / M;
      }else{ // H = q(q+1)/2, unstructured sigma
        // Orthonormalize A_std matrix
        arma::mat A_ortho(ids.n_elem,H);
        for(f=0;f<q;f++){
          arma::uvec Aidx = find(K_sub == f);
          if(f==0){
            A_ortho.cols(Aidx) = A_std.cols(Aidx);
          }else{
            arma::mat ortho_mat = ortho_factor[f];
            A_ortho.cols(Aidx) = A_std.cols(Aidx) * ortho_mat;
          }
        }
        eta.submat(m0,ids) = trans(Xk * beta.elem(fef_cols) + A_ortho * beta.elem(ref_cols) + offset(ids)) / M;
      }
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
    // Element-wise update of singular or grouped betaj
    // ----------------------------------------------------------------------------------//
    
    for(j=0; j<J_XZ; j++){
      
      // Identify covariates belonging to group j
      arma::uvec idxj = find(XZ_group == j);
      Kj = idxj.n_elem;
      
      if((init == 1) & (iter>=4) & (sum(beta.elem(idxj)) == 0.0)){
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
      Kr = idxr.n_elem;
      
      // Calculate zetaj vector (initialize as zero, then sum components of interest)
      arma::vec zetaj(Kj); zetaj.zeros(); 
      
      // Initialize mean(resid^2) as 0, then sum components of interest
      v0 = 0.0;
      
      if(max(idxr) < p){ // Last updated beta component was composed of fixed effects beta
        // Just need to update fixed effects portion of eta
        
        arma::mat Xr = X.cols(idxr);
        
        // Update fixed effects component of eta
        eta = eta.each_row() + trans(Xr * (beta.elem(idxr) - beta0.elem(idxr))) / M;
        
        if(max(idxj) < p){ // Calculate zetaj for fixed effect covariate(s)
          
          // Create vector including fixef zetaj and updated nu and v0
          arma::vec out = zeta_fixef(y, X, eta, idxj, family, link, nu); 
          zetaj = out.subvec(0,Kj-1);
          nu = out(Kj);
          v0 = out(Kj+1);
          
        }else{ // Calculate zetaj for random effect
          
          // If previously updated beta is part of fixed effects, next random effect
          // to update is random intercept.
          // Random intercept does not require orthogonalization
          
          for(k=0;k<d;k++){
            // Rows of Z to use
            arma::uvec ids = find(group == k);
            // Index of appropriate columns for Z and u
            for(f=0;f<q;f++){
              col_idx(f) = k + f*d;
            }
            
            arma::mat Zk = Z.submat(ids, col_idx);
            arma::mat Xj(ids.n_elem, Kj); Xj.zeros();
            
            for(m=0;m<M;m++){
              m0 = m;
              arma::mat A = kron(u.submat(m0,col_idx),Zk) * J;
              arma::mat A_std = (A.each_col() - center) / scale;
              Xj = A_std.cols(idxj - p);
              arma::vec tmp_out = resid_nu_v0_k(y, eta.submat(m0,ids), family, link, nu);
              arma::vec resid = tmp_out.subvec(0,ids.n_elem-1);
              nu = tmp_out(ids.n_elem);
              v0 += tmp_out(ids.n_elem+1);
              zetaj = zetaj + Xj.t() * resid;
              
            } // End m for loop
            
          } // End k for loop
          
        } // End if-else for fixed vs random effects zetaj
        
      }else{ // Need to update random effects portion of eta
        
        // Improvements in efficiency?
        
        // Note: if max(idxj) >= p, then calculate zetaj for a random effects group j
        //    else calculate zetaj for a fixed effects group j
        
        for(k=0;k<d;k++){
          // Rows of Z to use
          arma::uvec ids = find(group == k);
          // Index of appropriate columns for Z and u
          for(f=0;f<q;f++){
            col_idx(f) = k + f*d;
          }
          
          arma::mat Zk = Z.submat(ids, col_idx);
          arma::mat Xj(ids.n_elem, Kj); Xj.zeros();
          
            if(H == q){ // sigma specified with independence structure (diagonal)

              for(m=0;m<M;m++){
                m0 = m;
                arma::mat A = kron(u.submat(m0, col_idx),Zk) * J;
                arma::mat A_std = (A.each_col() - center) / scale;
                // Update random effects component of eta for each individual, weight by 1/M
                eta.submat(m0,ids) += A_std.cols(idxr - p) * (beta.elem(idxr) - beta0.elem(idxr)) / M;
                if(max(idxj) >= p){
                  Xj = A_std.cols(idxj - p);
                  arma::vec tmp_out = resid_nu_v0_k(y, eta.submat(m0,ids), family, link, nu);
                  arma::vec resid = tmp_out.subvec(0,ids.n_elem-1);
                  nu = tmp_out(ids.n_elem);
                  v0 += tmp_out(ids.n_elem+1);
                  zetaj = zetaj + Xj.t() * resid;
                }
              }

            }else{ // H = q(q+1)/2, unstructured sigma

              for(m=0;m<M;m++){
                m0 = m;
                arma::mat A = kron(u.submat(m0, col_idx), Zk) * J;
                arma::mat A_std = (A.each_col() - center) / scale;
                // Orthonormalize A_std matrix
                arma::mat A_ortho(ids.n_elem, H);
                for(f=0;f<q;f++){
                  arma::uvec Aidx = find(K_sub == f);
                  if(f==0){
                    A_ortho.cols(Aidx) = A_std.cols(Aidx);
                  }else{
                    arma::mat ortho_mat = ortho_factor[f];
                    A_ortho.cols(Aidx) = A_std.cols(Aidx) * ortho_mat;
                  }
                }
                // Update random effects component of eta for each individual, weight by 1/M
                eta.submat(m0,col_idx) = A_ortho.elem(idxr - p) * (beta.elem(idxr) - beta0.elem(idxr)) / M;
                if(max(idxj) >= p){
                  Xj = A_ortho.cols(idxj - p);
                  arma::vec tmp_out = resid_nu_v0_k(y, eta.submat(m0,ids), family, link, nu);
                  arma::vec resid = tmp_out.subvec(0,ids.n_elem-1);
                  nu = tmp_out(ids.n_elem);
                  v0 += tmp_out(ids.n_elem+1);
                  zetaj = zetaj + Xj.t() * resid;
                }
              } // End m for loop

            } // End if-else H == q
            
            if(max(idxj) < p){
              // Create list of fixef zetar and updated nu and v0
              arma::vec out = zeta_fixef(y(ids), X.rows(ids), eta.cols(ids), idxj, family, link, nu); 
              zetaj = zetaj + out.subvec(0,idxj.n_elem-1);
              nu = out(idxj.n_elem);
              v0 = out(idxj.n_elem+1);
            }
          
        } // End k for loop
        
        
      } // End if-else for fixed or random effects update to eta
      
      // resid in above = (y-mu), need resid = (y-mu)/nu 
      // Therefore, incorporate (1/nu) into zetaj and (1/nu^2) into v0 = mean(resid^2)
      // Finish calc of zetaj
      zetaj = zetaj / (N*M*nu) + beta.elem(idxj);
      v0 = v0 / (nu*nu);
      
      // L2 norm of zetaj
      if(idxj.n_elem == 1){
        zetaj_L2 = as_scalar(zetaj);
      }else{
        zetaj_L2 = sqrt(sum(zetaj % zetaj));
      }
      
      // Update betaj
      if((j==0) | (idxj(1) == p)){
        // j = 0: fixed effects intercept, other fixed effects covariates given XZ_group values of 0
        // idxj(1) == p: random effects intercept
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
      
      // Reset nu to 0 if not binomial or gaussian with canonical link
      if(!((std::strcmp(family, bin) == 0) & (link == 10)) & !((std::strcmp(family, gaus) == 0) & (link == 30))){
        nu = 0.0;
      }
      
      // Save r as the most recent updated j
      r = j;
      
    } // End j for loop
    
    // Convergence criteria after every full beta update
    if(fabs(v0 - v0_last)/(v0_last+0.1) < conv*0.001){
      converged = 1;
    }
    
  } // End while loop
  
  if(converged == 0){
    warning("grouped coordinate descent algorithm did not converge \n");
  }
  
  // Need to change scale of random effects componenet of beta back to original scale
  if(H == q){ // sigma matrix diagonal (independent)
    // Only need to undo standardization (centering and scaling)
    arma::vec beta_ranef = beta.subvec(p,p+H-1);
    beta(0) = beta(0) - sum((beta_ranef - center) / scale);
    beta.subvec(p,p+H-1) = beta_ranef / scale;
  }else{ // sigma matrix unstructured
    // Need to (a) unorthogonalize and (b) undo standardization (centering and scaling)
    arma::vec beta_ranef = beta.subvec(p,p+H-1);
    for(f=1;f<q;f++){
      arma::mat ortho_mat = ortho_factor[f];
      arma::uvec beta_idx = find(K_sub == f);
      beta_ranef.elem(beta_idx) = ortho_mat * beta_ranef.elem(beta_idx);
    }
    beta(0) = beta(0) - sum((beta_ranef - center) / scale);
    beta.subvec(p,p+H-1) = beta_ranef / scale;
  }
  
  return beta;
}

