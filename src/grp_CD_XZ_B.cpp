// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utility_glm.h"
#include "penalties.h"

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
                      const char* penalty, double lambda, arma::vec params) {
  
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
  
  int i=0; // observation counter
  int j=0; // covariate group counter
  int m=0; // MCMC draw counter
  int f=0; 
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
  
  arma::vec mu(M); // expected means 
  arma::mat eta(M,N); // linear predictors
  arma::vec deriv(M); // (d_eta / d_mu) = g'(mu) where g = link function
  arma::vec Vmu(M); // variance = b''(theta) where b(theta) from exponential family
  arma::vec weights(M);
  arma::vec resid(M); // IRLS residuals
  arma::vec const_ones(M); const_ones.ones(); // arbitrary vector of ones
  arma::vec mu_check(M);
  
  double nu=0.0; // max second deriv of loss function (max possible weight)
  double nu_tmp=0.0;
  double zj=0.0; // 1/N * t(X_j) %*% resid + beta0_j if can be made into scalar
  double zetaj_L2=0.0; // L2 norm of zetaj vector
  
  arma::vec beta0 = beta; // Input initial beta / beta from last round of updates
  double bj=0.0; // place-holder for value of element of beta
  arma::uvec fef_cols = seq(0,p-1);
  arma::uvec ref_cols = seq(p, p + J.n_cols - 1);
  
  int grp = 0; // observation group
  arma::uvec col_idx(q); // columns of Z and u to choose at a time
  arma::rowvec Zki(d);
  arma::rowvec Xki(p);
  arma::uvec id(1);
  int H = J.n_cols; // From paper, J_q
  arma::rowvec A(H); // From paper, (alpha_k kronecker Zki)^T * J_q
  arma::rowvec A_std(H); // Standardized B (center of 0, sum of square column values = 1)
  arma::rowvec A_ortho(H); // Orthogonalized (by group) version of A_std
  
  arma::vec K_vec = arma::conv_to<arma::vec>::from(K);
  arma::vec lam = lambda * sqrt(K_vec); // lambda_j = lambda * sqrt(Kj) where Kj = size of jth group
  
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
  
  for(i=0;i<N;i++){
    
    // Rows of X and Z to use
    id = i;
    Xki = X.row(i);
    // Identify observation group, find appropriate columns of Z and u
    grp = group(i);
    for(f=0;f<q;f++){
      col_idx(f) = grp + f*d;
    }
    
    Zki = Z.submat(id,col_idx);
    
    for(m=0;m<M;m++){
      A = kron(u.row(m),Zki) * J;
      center = center + A;
      scale = scale + A % A;
    }
    
  }
  
  // If ncols(J) > q then orthogonalize. Not necessary otherwise.
    // If sigma matrix diagonal, then coefficients relating to random effects not grouped
    // and orthogonalization not necessary. If sigma matrix diagonal, ncol(J) == q.
    // Otherwise, ncol(J) == q(q+1)/2
    
  List ortho_factor;
  
  if(H > q){
    
    arma::uvec K_sub = K.subvec(p,p+H-1) - K(p); // Includes integers from 0 to (q-1)
    
    // Initialize list for t(A_j) * A_j matrices
    List L;
    
    for(j=0;j<q;j++){
      arma::mat xTx(j+1,j+1); 
      if(j == 0){
        xTx.ones(); // Standardization sufficient for group of size 1, no additional orthogonalization necessary
      }else{
        xTx.zeros();
      }
      L[j] = xTx;
    }
    
    for(i=0;i<N;i++){
      
      // Rows of X and Z to use
      id = i;
      
      // Identify observation group, find appropriate columns of Z and u
      grp = group(i);
      for(f=0;f<q;f++){
        col_idx(f) = grp + f*d;
      }
      
      Zki = Z.submat(id,col_idx);
      
      for(m=0;m<M;m++){
        A = kron(u.row(m),Zki) * J;
        A_std = (A - center) / scale;
        
        // Calculate t(A_j) * A_j for groups
        for(j=1;j<q;j++){ // Skip j = 0, no need to orthogonalize group of size 1
          arma::mat xTx = L[j];
          arma::uvec Kidx = find(K_sub == j);
          xTx = xTx + trans(A_std(Kidx)) * A_std(Kidx);
          L[j] = xTx;
        }
      }
      
    }
    
    // Divide t(Xj) * Xj by (N*M)
    for(j=1;j<q;j++){
      arma::mat xTx = L[j];
      xTx = xTx / (N*M);
      L[j] = xTx;
    }
    
    // Find eigenvectors and eigenvalues of the xTx matrices, create matrix used to orthonormalize
    // Note: Incorporate checks to make sure eigen values (evec) not equal to 0
    
    for(j=0;j<q;j++){
      if(j == 0){
        ortho_factor[j] = 1.0;
      }else{
        arma::mat xTx = L[j];
        arma::mat emat;
        arma::vec evec;
        arma::eig_sym(evec, emat, xTx);
        ortho_factor[j] = emat.each_col() % (1/sqrt(evec));
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
    nu = 0.0;
  }
  
  //-------------------------------------------------------------------------------//
  // Initialize eta matrix
  //-------------------------------------------------------------------------------//
  
  for(i=0;i<N;i++){
    
    // Rows of X and Z to use
    id = i;
    Xki = X.row(i);
    // Identify observation group, find appropriate columns of Z and u
    grp = group(i);
    for(f=0;f<q;f++){
      col_idx(f) = grp + f*d;
    }
    
    Zki = Z.submat(id,col_idx);
    
    for(m=0;m<M;m++){
      A = kron(u.row(m),Zki) * J;
      A_std = (A - center) / scale;
      
      if(H == q){ // sigma specified with independence structure (diagonal)
        // Caculate eta for each individual using initial beta, weight by 1/M
        eta(m,i) = as_scalar((Xki * beta.elem(fef_cols) + A_std * beta.elem(ref_cols) + offset(i))) / M;
      }else{ // H = q(q+1)/2, unstructured sigma
        // Orthonormalize A_std matrix
        arma::uvec K_sub = K.subvec(p,p+H-1) - K(p); // Includes integers from 0 to (q-1)
        for(f=0;f<q;f++){
          arma::uvec Aidx = find(K_sub == f);
          if(f==0){
            A_ortho.cols(Aidx) = A_std.cols(Aidx);
          }else{
            arma::mat ortho_mat = ortho_factor[f];
            A_ortho.cols(Aidx) = A_std.cols(Aidx) * ortho_mat;
          }
        }
        eta(m,i) = as_scalar((Xki * beta.elem(fef_cols) + A_ortho * beta.elem(ref_cols) + offset(i))) / M;
      }
    }
    
    // Calculate zetaj for j = 0
    
    // Find covariates in group j = 0
    arma::uvec idxj = find(XZ_group == 0);
    
    // Calculate zetaj vector (initialize as zero, then sum components of interest)
    arma::vec zetaj(idxj.n_elem); zetaj.zeros(); 
    
    v0 = 0.0;
    
    arma::mat Xj = X.cols(idxj);
    
    for(i=0;i<N;i++){
      // Update mu, resid, weights
      mu = invlink(link, eta.col(i));
      mu_check = muvalid(family, mu);
      
      deriv = dlink(link, mu);
      Vmu = varfun(family, mu);
      resid = deriv % ((y(i)*const_ones) - mu);
      for(m=0;m<M;m++){
        if(mu_check(m)==0){
          resid(m) = 0.0; // Ignore invalid mu
        }
      }
      
      v0 = v0 + mean(resid % resid);
      
      // If not binomial or gaussian family with canonical link, calculate nu = max(weights)
      if(!((std::strcmp(family, bin) == 0) & (link == 10)) & !((std::strcmp(family, gaus) == 0) & (link == 30))){
        
        weights = const_ones / (deriv % deriv % Vmu);
        for(m=0;m<M;m++){
          if(mu_check(m)==0){
            weights(m) = 0.0;
          }
        }
        nu_tmp = max(weights);
        if(nu_tmp > nu){
          nu = nu_tmp;
        }
      }
      
      // Update zetaj sum
      arma::vec Xji = Xj.row(i);
      zetaj = zetaj + Xji * sum(resid);
      
    }
    
    zetaj = zetaj / (N*M) + beta.elem(idxj);
    
  }
  
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
    // Set zetaj = zeta(j+1) from last
    
    // Element-wise update of beta or grouped betaj
    for(j=0; j<J_XZ; j++){
      
      // Identify covariates belonging to group j
      arma::uvec idxj = find(XZ_group == j);
      
      if((init == 1) & (iter>=4) & (sum(beta.elem(idxj)) == 0)){
        // If beta penalized to zero in past round, will stay zero in further rounds
        // Therefore, skip to next covariate grouping
        continue;
      }else if((init == 0) & (sum(beta.elem(idxj)) == 0)){
        continue;
      }
      
      // Calculate zetaj vector (initialize as zero, then sum components of interest)
      arma::vec zetaj(idxj.n_elem); zetaj.zeros(); 
      
      v0 = 0.0;
      
      if(max(idxj) > p){ // Just need to update fixed effects portion of eta
        
        arma::mat Xj = X.cols(idxj);
        // Update fixed effects component of eta
        eta = eta.each_row() + trans(Xj * (beta.elem(idxj) - beta0.elem(idxj))) / M;
        
        for(i=0;i<N;i++){
          // Update mu, resid, weights
          mu = invlink(link, eta.col(i));
          mu_check = muvalid(family, mu);
          
          deriv = dlink(link, mu);
          Vmu = varfun(family, mu);
          resid = deriv % ((y(i)*const_ones) - mu);
          for(m=0;m<M;m++){
            if(mu_check(m)==0){
              resid(m) = 0.0; // Ignore invalid mu
            }
          }
          
          v0 = v0 + mean(resid % resid);
          
          // If not binomial or gaussian family with canonical link, calculate nu = max(weights)
          if(!((std::strcmp(family, bin) == 0) & (link == 10)) & !((std::strcmp(family, gaus) == 0) & (link == 30))){
            
            weights = const_ones / (deriv % deriv % Vmu);
            for(m=0;m<M;m++){
              if(mu_check(m)==0){
                weights(m) = 0.0;
              }
            }
            nu_tmp = max(weights);
            if(nu_tmp > nu){
              nu = nu_tmp;
            }
          }
          
          // Update zetaj sum
          arma::vec Xji = Xj.row(i);
          zetaj = zetaj + Xji * sum(resid);
          
        }
        
      }else{ // Need to update random effects portion of eta
        
        // Improvements in efficiency?
        
        arma::mat Xj(M,idxj.n_elem); Xj.zeros(); // Max dimension M x q
        
        for(i=0;i<N;i++){
          
          // Rows of Z to use
          id = i;
          // Identify observation group, find appropriate columns of Z and u
          grp = group(i);
          for(f=0;f<q;f++){
            col_idx(f) = grp + f*d;
          }
          
          Zki = Z.submat(id,col_idx);
          
          if(H == q){ // sigma specified with independence structure (diagonal)
            
            for(m=0;m<M;m++){
              A = kron(u.row(m),Zki) * J;
              A_std = (A - center) / scale;
              // Update random effects component of eta for each individual, weight by 1/M
              eta(m,i) += as_scalar(A_std.elem(idxj - p) * (beta.elem(idxj) - beta0.elem(idxj))) / M;
              Xj.row(m) = A_std.elem(idxj - p);
            }
            
          }else{ // H = q(q+1)/2, unstructured sigma
            
            for(m=0;m<M;m++){
              A = kron(u.row(m),Zki) * J;
              A_std = (A - center) / scale;
              // Orthonormalize A_std matrix
              arma::uvec K_sub = K.subvec(p,p+H-1) - K(p); // Includes integers from 0 to (q-1)
              for(f=0;f<q;f++){
                arma::uvec Aidx = find(K_sub == f);
                if(f==0){
                  A_ortho.cols(Aidx) = A_std.cols(Aidx);
                }else{
                  arma::mat ortho_mat = ortho_factor[f];
                  A_ortho.cols(Aidx) = A_std.cols(Aidx) * ortho_mat;
                }
              }
              // Update random effects componenet of eta for each individual, weight by 1/M
              eta(m,i) = as_scalar(A_ortho.elem(idxj - p) * (beta.elem(idxj) - beta0.elem(idxj))) / M;
              Xj.row(m) = A_ortho.elem(idxj - p);
            }
            
          }
          
          // Update mu, resid, weights
          mu = invlink(link, eta.col(i));
          mu_check = muvalid(family, mu);
          
          deriv = dlink(link, mu);
          Vmu = varfun(family, mu);
          resid = deriv % ((y(i)*const_ones) - mu);
          for(m=0; m>M; m++){
            if(mu_check(m)==0){
              resid(m) = 0.0; // Ignore invalid mu
            }
          }
          
          v0 = v0 + mean(resid % resid);
          
          // If not binomial or gaussian family with canonical link, calculate nu = max(weights)
          if(!((std::strcmp(family, bin) == 0) & (link == 10)) & !((std::strcmp(family, gaus) == 0) & (link == 30))){
            
            weights = const_ones / (deriv % deriv % Vmu);
            for(m=0; m>M; m++){
              if(mu_check(m)==0){
                weights(m) = 0.0;
              }
            }
            nu_tmp = max(weights);
            if(nu_tmp > nu){
              nu = nu_tmp;
            }
          }
          
          // Calculate zetaj
          zetaj = zetaj + Xj.t() * resid;
          
        } // End i for loop
        
      } // End if-else for fixed or random effects
      
      zetaj = zetaj / (N*M) + beta.elem(idxj);
      
      // L2 norm of zetaj
      if(idxj.n_elem == 1){
        zetaj_L2 = as_scalar(zetaj);
      }else{
        zetaj_L2 = sqrt(sum(zetaj % zetaj));
      }
      
      // Save latest value of beta
      beta0 = beta;
      
      // Update beta
      if(j==0){
        // No penalization for the intercept and other covariates given XZ_group values of 0
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
          arma::vec vec(1); vec.ones();
          beta.elem(idxj) = bj * vec;
        }else{
          beta.elem(idxj) = bj * (zetaj / zetaj_L2);
        }
        
      } // End if-else update to beta for group j
      
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
    arma::uvec K_sub = K.subvec(p,p+H-1) - K(p);
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

