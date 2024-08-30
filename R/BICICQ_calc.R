# Calculate BIC-ICQ
## NOTE: for internal use within package only
## Used to calculate BIC-ICQ within glmm, glmm_FA, phmm, or phmm_FA (single model fit)

#' @importFrom bigmemory attach.big.matrix describe
#' @importFrom mvtnorm dmvnorm
#' @importFrom Matrix Matrix
BICICQ_calc = function(ufull_describe, fit_type, 
                       family, link_int, sig_g, phi,
                       y, X, Z_Matrix, group, d,
                       coef, J, offset_fit,
                       r){
  
  # Extract big matrix of relevant posterior draws for BIC-ICQ calculation
  ufull = attach.big.matrix(ufull_describe)
  
  # Convert Z from class Matrix to matrix
  Z = Matrix::as.matrix(Z_Matrix)
  
  if(fit_type == "glmmPen"){
    q1 = Qfun(y, X, Z, ufull@address, group, J, coef, offset_fit, c(d,ncol(Z)/d), family, link_int, sig_g, phi)
    q2 = 0
    for(k in 1:d){
      cols_idx = seq(from = k, to = ncol(Z), by = d)
      post_k = ufull[,cols_idx]
      q2 = q2 + sum(dmvnorm(post_k, log=TRUE)) / nrow(ufull)
    }
  }else if(fit_type == "glmmPen_FA"){
    # const arma::vec& y, const arma::mat& X, const arma::mat& Z, SEXP pBigMat, 
    # const arma::vec& group, const arma::sp_mat& J_f,
    # const arma::vec& beta, const arma::vec offset, arma::vec dims,
    # const char* family, int link, double sig_g, double phi
    q1 = Qfun_FA(y, X, Z, ufull@address, group, J, coef, offset_fit, c(d, ncol(Z)/d, r), family, link_int, sig_g, phi)
    q2 = 0
    for(k in 1:d){
      cols_idx = seq(from = k, to = ncol(ufull), by = d)
      post_k = ufull[,cols_idx]
      q2 = q2 + sum(dmvnorm(post_k, log=TRUE)) / nrow(ufull)
    }
  }
  
  llq = q1 + q2
  BICq = -2*llq + sum(coef != 0)*log(nrow(X))
  
  if(is.na(BICq)){
    warning("BICq value calculated as NA due to divergent coefficient values. 
              Consider checking correlations in the covariates and/or adjusting
              model fit parameters.", immediate. = TRUE)
    cat("Fixed effects (scaled X): \n")
    cat(round(coef[1:ncol(X)],3), "\n")
    
    if(fit_type == "glmmPen"){
      gamma = matrix(J%*%matrix(coef[-c(1:ncol(X))], ncol = 1), ncol = ncol(Z)/d)
      cov = gamma %*% t(gamma)
    }else if(fit_type == "glmmPen_FA"){
      B = t(matrix(coef[-c(1:ncol(X))], nrow = r, ncol = ncol(Z)/d))
      cov = B %*% t(B)
    }
    
    if(nrow(cov) <= 5){
      cat("random effect covariance matrix: \n")
      print(round(cov,3))
    }else{
      cat("random effect covariance matrix diagonal: \n")
      cat(round(diag(cov),3), "\n")
    }
  }
  
  return(BICq)
  
}


############################################################################################################