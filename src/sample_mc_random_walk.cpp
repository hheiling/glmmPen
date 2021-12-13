// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "utility_glm.h"

using namespace Rcpp;

// Adaptive (Burn-in Only) Random Walk Metropolis within Gibbs with Random Scan

// [[Rcpp::export]]
arma::mat sample_mc_gibbs_adapt_rw(arma::mat f, // matrix
                                   arma::mat z, // matrix
                                   arma::vec y, // vector
                                   int NMC, // integer
                                   arma::vec u0, //matrix
                                   arma::rowvec proposal_SD, // row vector
                                   int batch, int batch_length, int offset,
                                   int nMC_burnin,
                                   const char* family, int link, double phi, double sig_g){ 
  
  arma::mat fitted = f;
  arma::mat Z=z;
  arma::vec Y=y;
  int nMC = NMC; 
  arma::vec uold=u0;
  
  int q = Z.n_cols;
  int n = Z.n_rows;
  
  int i = 0;
  int j = 0;
  int l = 0; 
  int k = 0;
  int naccept = 0;
  int index = 0;
  int index2 = 0;
  double sum = 0;
  double sumn = 0;
  double w = 0;
  double ep = 0; 
  double e0 = 0;
  double delta = 0;
  double increment = 0;
  arma::mat out(nMC+3, q); // Last three lines = acceptance rates, updated proposal_SD, and updated batch number, respectively
  arma::vec e(q);
  arma::vec rate(q);
  arma::vec etae(n);
  arma::vec etaen(n);
  arma::vec mu_denom(n);
  arma::vec mu_num(n);
  arma::vec prob_nb(n);
  arma::vec accept_index(q);
  arma::rowvec SD = proposal_SD; // At beginning of EM algorithm, proposal_SD = 1.0 for each variable
  arma::rowvec acc_rate(q);
  arma::uvec var_index(q);
  arma::uvec samp(1);
  arma::rowvec final_acc_rate(q);
  arma::rowvec batch_vec(q); batch_vec.zeros();
  
  const char* bin = "binomial";
  const char* pois = "poisson";
  const char* negbin = "Negative Binomial(1)";
  const char* gaus = "gaussian";
  
  RNGScope scope;
  
  //initialize e 
  for(i = 0; i < q; i++){
    e(i) = uold(i);
    accept_index(i) = 0.0;
    acc_rate(i) = 0.0;
    var_index(i) = i;
  }
  
  //iteratively update e 
  while(naccept < nMC){
    
    // Random scan: randomize order of updated random effect
    samp = RcppArmadillo::sample(var_index, q, 0);
    
    for(k = 0; k < q; k++){
      
      j = samp(k);
      // scan_index(j) = scan_index(j) + 1.0;
      
      // calculate etae (using current value of e)
      etae = fitted + Z*e;
      
      // save current value of e(i)
      e0 = e(j);
      
      // generate w
      w = log(R::runif(0.0,1.0));
      
      // generate proposal e
      ep = e0 + R::rnorm(0.0, SD(j));
      e(j) = ep;
      
      // calculate updated etae (using new proposed value of e)
      etaen = fitted + Z*e;
      
      // calculate Metropolis-Hastings (MH) ratio
      sum = sumn = 0;
      // sum: log of denominator of MH ratio
      // sumn: log of numerator of MH ratio
      mu_denom = invlink(link, etae);
      mu_num = invlink(link, etaen);
      // Note: last argument 1 indicates request for log value of loglik estimate
      if(std::strcmp(family, bin) == 0){
        for(l = 0; l < n; l++){
          sum = sum + R::dbinom(Y(l), 1.0, mu_denom(l), 1);
          sumn = sumn + R::dbinom(Y(l), 1.0, mu_num(l), 1);
        }
      }else if(std::strcmp(family, pois) == 0){
        for(l = 0; l < n; l++){
          sum = sum + R::dpois(Y(l), mu_denom(l), 1);
          sumn = sumn + R::dpois(Y(l), mu_num(l), 1);
        }
      }else if(std::strcmp(family, negbin) == 0){
        for(l = 0; l < n; l++){
          sum = sum + R::dnbinom(Y(l), (1.0/phi), (1.0 / (1.0 + mu_denom(l) * phi)), 1);
          sumn = sumn + R::dnbinom(Y(l), (1.0/phi), (1.0 / (1.0 + mu_num(l) * phi)), 1);
        }
      }else if(std::strcmp(family, gaus) == 0){
        for(l = 0; l < n; l++){
          // Need to replace 1.0 in below with appropriate sigma
          // sum = sum + R::dnorm4(Y(l), mu_denom(l), 1.0, 1);
          // sumn = sumn + R::dnorm4(Y(l), mu_num(l), 1.0, 1);
          sum = sum + R::dnorm(Y(l), mu_denom(l), sig_g, 1);
          sumn = sumn + R::dnorm(Y(l), mu_num(l), sig_g, 1);
        }
      }
      
      // for(l = 0; l < n; l++){
      //   sum = sum + R::dbinom(Y(l), 1.0, exp(etae(l))/(1+exp(etae(l))), 1);
      //   sumn = sumn + R::dbinom(Y(l), 1.0, exp(etaen(l))/(1+exp(etaen(l))), 1);
      // }
      
      sum = sum + R::dnorm(e0, 0.0, 1.0, 1) ;
      sumn = sumn + R::dnorm(ep, 0.0, 1.0, 1) ;
      // Assume same as R version of dnorm: use standard deviation
      
      // check for acceptance
      if(w < sumn - sum){
        // ep left in e(i); accept proposal
        accept_index(j) = accept_index(j) + 1.0;
      }else{
        // reject proposal, keep most recent accepted value
        e(j) = e0;
      }
      
    } // End k loop
    
    if(index == 1+5*naccept){
      out.row(naccept) = e.t();
      naccept++;
    }
    
    index++;
    index2++;
    
    
    // Update proposal variance
    // if((index % batch_length == 0) && (batch < burnin_batchnum)){ // if index = multiple of batch_length
    if((index % batch_length == 0) && (naccept < nMC_burnin)){ // if index = multiple of batch_length
      
      // Update batch information 
      batch = batch + 1;
      
      for(j = 0; j < q; j++){
        
        // Determine acceptance rate for latest batch
        acc_rate(j) = accept_index(j) / (double)batch_length;
        
        // Update proposal SD (separate for each variable)
        // delta = min(0.01, (T_b)^(-1/2))
        increment = 1 / sqrt(batch*batch_length + offset);
        if(increment < 0.01){
          delta = increment;
        }else{
          delta = 0.01;
        }
        
        // Update proposal standard deviation based on acceptance rate for last batch
        // log(std_dev) +/- delta --> std_dev * exp(+/- delta)
        // Based on papers (Roberts and Rosehthal 2009) and (Chimisov, Latuszynski, and Roberts 2018),
        // increase SD if acc_rate > 0.44, decrease SD if acc_rate < 0.44, which directly contradicts
        // the Givens and Hoeting textbook 
        if(acc_rate(j) > 0.5){
          SD(j) = SD(j) * exp(delta); 
        }else if(acc_rate(j) < 0.4){
          SD(j) = SD(j) * exp(-delta); 
        }
        
        // Set min and max cap of log(standard deviation)
        // No longer need if only adapting during burn-in?
        if(log(SD(j)) > log(5)){
          SD(j) = 5.0;
        }else if(log(SD(j)) < -log(5)){
          SD(j) = 0.20;
        }
        
        // Re-set accept_index - new acceptance rate for next batch
        accept_index(j) = 0.0;
        
        // Re-set index2 (only count index2 after burnin_batchnum limit reached)
        index2 = 0.0;
        
      } // End j for loop (w/in if statement)
      
    } // End if(index % (int)batch_lenth)
    
  } // End while loop
  
  for(j = 0; j < q; j++){
    final_acc_rate(j) = accept_index(j) / index2;
  }
  batch_vec(0) = (double)batch;
  
  // Rcout << "Last Batch Acceptance Rate" << std::endl << acc_rate;
  // Rcout << "Final Updated Proposal SD" << std::endl << SD;
  
  out.row(nMC) = final_acc_rate; // Third-to-last row
  out.row(nMC+1) = SD; // Second-to-last row
  out.row(nMC+2) = batch_vec; // Last row
  
  return(out); 
  
}
