 



############################################################################
# Updated for adaptive random walk metropolis w/in gibbs with random scan
# Adapt proposal SD only during burn-in period
# Output a big.matrix for the posterior draws 
############################################################################

#' @name sample_mc_adapt_BigMat
#' @aliases sample_mc2_BigMat
#' 
#' @title Calculate Monte Carlo draws using Metropolis-within-Gibbs
#' 
#' @description Samples draws from the posterior distribution of the random effects using 
#' Metropolis-within-Gibbs. \code{sample_mc_adapt_BigMat} samples using adaptive 
#' random walk and \code{sample_mc2_BigMat} samples using independence sampling. See
#' \code{\link{adaptControl}} for possible adaptive random walk adjustments. Draws are automatically
#' thinned to only record every fifth sample.
#' 
#' @inheritParams fit_dat_B
#' @param coef a numeric vector of the fixed effects coefficients (from the latest EM iteration or
#' the final model of interest)
#' @param ranef_idx a vector of integers describing which random effects are non-zero (i.e. which
#' diagonal elements of the sigma matrix are non-zero)
#' @param y a numeric vector of the response variable
#' @param X a model matrix of the fixed effects covariates
#' @param Z a sparse model matrix of the random effects multiplied by the lower triangular cholesky
#' decomposition of the sigma matrix (from the latest EM iteration or the final model of interest)
#' @param group a factor vector of the grouping variable, converted to a factor of 
#' consecutive numeric integers
#' @param d integer, the number of groups present (number factors of the group vector)
#' @param okindex ?
#' @param uold a matrix with a single row comprised of the last Monte Carlo draw from either the 
#' most recent E step of the EM algorithm or the final model
#' @param proposal_SD a matrix of dimension (number of groups)x(number of random effects) that
#' gives the proposal standard deviation for the adaptive random walk. The \code{glmmPen} MCEM 
#' algorithm initializes the proposal standard deviation as 1.0 for all variables and groups
#' and is updated at the start of every E step when the Metropolis-within-Gibbs adaptive random walk
#' algorithm is used. 
#' @param batch integer specifying how many batches of posterior draws has already been sampled in
#' the MCEM algorithm. As the batch number increases, the proposal standard deviation for the adaptive 
#' random walk algorithm is adjusted by a smaller amount.
#' @param batch_length an integer specifying the number of posterior draws (ignoring thinning) used
#' to evaluate the acceptance rate of the random walk algorithm. This acceptance rate is then used
#' to adjust the proposal standard deviation if necessary.
#' @param offset
#' @param burnin_batchnum an integer specifying the number of batches of (unthinned) posterior draws
#' to allow adjustments to the proposal standard deviation. For the MCEM algorithm, it is recommended
#' to allow for the adjustment of the proposal standard deviation at the beginning of each E step.
#' 
#' @return a list made of the following components:
#' \item{u0}{list of the information needed to attached to a big.matrix object. Use 
#' \code{bigmemory::attach.big.matrix(u0)} to extract a big.matrix of the Monte Carlo draws of
#' the random effect posterior distribution. Number rows = nMC, number columns = 
#' (number random effects)*(number groups) = ncol(Z). Organization of columns: first by random effect 
#' variable, then by group within variable (i.e. Var1:Grp1 Var1:Grp2 ... Var1:GrpK Var2:Grp1 ... Varq:GrpK)}
#' \item{gibbs_accept_rate}{matrix of dimension (number groups)x(number random effects). Each entry
#' gives the proportion of accepted draws for the particular random effect variable in group k 
#' either in the entire \code{nMC} draws (for \code{sample_mc2_BigMat}) or for the last 
#' \code{batch_length} draws.}
#' \item{proposal_SD}{a matrix of dimension (number of groups)x(number of random effects) that is
#' only output from \code{sample_mc_adapt_BigMat}. It records the updated proposal standard deviation
#' for the particular random effect variable in group k.}
#' \item{updated_batch}{integer specifying the updated number of batches of posterior draws that 
#' have been sampled.}
#' 
#' @importFrom bigmemory big.matrix describe
#' @export
sample_mc_adapt_BigMat = function(coef, ranef_idx, y, X, Z, nMC, trace = 0, family, group, d, okindex,
                                  uold, proposal_SD, batch, batch_length = 100, 
                                  offset = 0, burnin_batchnum = 500){
  
  f = get(family, mode = "function", envir = parent.frame())
  
  uhat  = rep(0, ncol(Z))
  eta = X %*% matrix(coef, ncol=1) 
  
  # matrix to hold accepted samples
  u0 = big.matrix(nrow = nMC, ncol = ncol(Z), init = 0)
  # u0 = matrix(rnorm(nMC*ncol(Z)) , nMC, ncol(Z))
  
  # fitted
  fitted_mat = as.matrix(X %*% matrix(coef, ncol=1))
  #generate samples for each i
  
  q = ncol(Z) / d
  gibbs_accept_rate = matrix(NA, nrow = d, ncol = q)
  gibbs_output = NULL

  for(i in 1:d){
    select = group == i
    index = seq(i, ncol(Z), by = d)
    ## new code to limit to non-zero z, skipping elements of q where diag(sigma) are 0
    # index = index[which(diag(cov) != 0)]
    index = index[ranef_idx]
    if(length(index) == 0) next
    var_index = ranef_idx
    
    gibbs_output = sample_mc_gibbs_adapt_rw(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                                            matrix(Z[select,index],ncol = length(index), nrow = sum(select)),  
                                            y[select], nMC, uold[index], 
                                            matrix(proposal_SD[i,var_index], nrow = 1), batch, 
                                            batch_length, offset, burnin_batchnum, trace)
    
    u0[,index] = gibbs_output[1:nMC,]
    gibbs_accept_rate[i,var_index] = matrix(gibbs_output[(nMC+1),], nrow = 1)
    proposal_SD[i,var_index] = matrix(gibbs_output[(nMC+2),], nrow = 1)
    
  }
  
  
  return(list(u0 = describe(u0), gibbs_accept_rate = gibbs_accept_rate, 
              proposal_SD = proposal_SD, updated_batch = batch))
  
} # End sample_mc_adapt_BigMat()

#' @rdname sample_mc_adapt_BigMat
#'
#' @importFrom bigmemory big.matrix describe
#' @export
sample_mc2_BigMat = function(coef, ranef_idx, y, X, Z, nMC, trace = 0, family, group, d, okindex,
                             uold){
  
  uhat  = rep(0, ncol(Z))
  eta = X %*% matrix(coef, ncol=1) 
  
  # matrix to hold accepted samples
  u0 = big.matrix(nrow = nMC, ncol = ncol(Z), init = 0)
  # u0 = matrix(rnorm(nMC*ncol(Z)) , nMC, ncol(Z))
  
  # fitted
  fitted_mat = as.matrix(X %*% matrix(coef, ncol=1))
  #generate samples for each i
  
  error_out = F
  q = ncol(Z) / d
  gibbs_accept_rate = matrix(NA, nrow = d, ncol = q)
  
  
  for(i in 1:d){
    select = group == i
    index = seq(i, ncol(Z), by = d)
    ## new code to limit to non-zero z, skipping elements of q where diag(sigma) are 0
    # index = index[which(diag(cov) != 0)]
    index = index[ranef_idx]
    if(length(index) == 0) next
    
    gibbs_list = sample_mc_inner_gibbs(matrix(fitted_mat[select], ncol = 1, nrow = sum(select)), 
                                       matrix(Z[select,index],ncol = length(index), nrow = sum(select)),  
                                       y[select], uhat[index], nMC, uold[index], 
                                       trace)
    u0[,index] = gibbs_list$u
    gibbs_accept_rate[i,] = matrix(gibbs_list$acc_rate, nrow = 1)
  }
  
  return(list(u0 = describe(u0), gibbs_accept_rate = gibbs_accept_rate))
  
  
} # End sample_mc2_BigMat

# tau calculate for rejection sampling - discontinued
# if(family == "binomial" & gibbs == F){
#   tau = dbinom(y, size = 1, prob = exp(eta)/(1+exp(eta)), log = T)
# }else if(family == "poisson" & gibbs == F){
#   tau = dpois(y, lambda = exp(eta), log = T)
# }else if(family == "gaussian" & gibbs == F){
#   s2 = sum((y-eta)*(y-eta)) / (length(y) - length(coef))
#   tau = dnorm(y, mean = eta, sd = sqrt(s2), log = T)
# }


###################################################################################################