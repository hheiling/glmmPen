// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// Original: NumericMatrix sample_mc_inner_gibbs(...)

// No adaptive change in proposal variation
// [[Rcpp::export]]
List sample_mc_inner_gibbs(arma::mat f, // matrix
                          arma::mat z, // matrix
                          arma::vec y, // vector
                          arma::vec t, // vector
                          int NMC, // integer
                          arma::vec u0, //matrix
                          int trace){ // integer
  
  arma::mat fitted = f;
  arma::mat Z=z;
  arma::vec Y=y;
  arma::vec tau =t;
  int nMC = NMC; 
  arma::vec uold=u0;
  
  int q = Z.n_cols;
  int n = Z.n_rows;
  
  int i = 0;
  int j = 0; 
  // int k = 0; 
  int l = 0; 
  int naccept = 0;
  int index = 0;
  // int batch = 1;
  double sum = 0;
  double sumn = 0;
  // double stao = 0;
  double w = 0;
  double ep = 0; 
  double e0 = 0;
  arma::mat out(nMC, q);
  arma::vec e(q);
  arma::vec rate(q);
  arma::vec etae(n);
  arma::vec etaen(n);
  arma::vec index2(q);
  arma::vec var(q);
  arma::vec acc_rate(q);
  
  RNGScope scope;
  
  //initialize e 
  for(i = 0; i < q; i++){
    e(i) = uold(i);
    index2(i) = 0.0;
    var(i) = 1.0;
  }
  
  //iteratively update e
  while(naccept < nMC){
    
    for(j = 0;j < q; j++){
      
      // calculate etae
      etae = fitted + Z*e;
      
      // save current value of e(i)
      e0 = e(j);
      
      // generate w
      w = log(R::runif(0.0,1.0));
      
      // generate proposal e
      ep = R::rnorm(0.0, var(j));
      e(j) = ep;
      
      // calculate updated etae
      etaen = fitted + Z*e;
      
      // calculate ratio
      sum = sumn = 0;
      for(l = 0; l < n; l++){
        sum = sum + R::dbinom(Y(l), 1.0, exp(etae(l))/(1+exp(etae(l))), 1);
        sumn = sumn + R::dbinom(Y(l), 1.0, exp(etaen(l))/(1+exp(etaen(l))), 1);
      }
      
      // check for acceptance
      if(w < sumn - sum){
        // ep left in e(i)
        index2(j) = index2(j) + 1.0;
      }else{
        e(j) = e0;
      }
      //  }
    }
    
    if(index == 1+5*naccept){
      out.row(naccept) = e.t();
      naccept++;
    }
    index++;
    
    
    
  }
  
  // Info for acceptance rate:
  for(i = 0; i < q; i++){
    acc_rate(i) = index2(i) / index;
  }

  // if(trace == 2){
  //   Rprintf("index: %d, naccept: %d, accept. rate: %f  \n", index, naccept, acc_rate);
  // }
  
  return(List::create(Named("u") = wrap(out), Named("acc_rate") = acc_rate));
  // return(wrap(out));
  
}

// Adaptive change in proposal variation
// [[Rcpp::export]]
NumericMatrix sample_mc_inner_gibbs2(arma::mat f, // matrix
                           arma::mat z, // matrix
                           arma::vec y, // vector
                           arma::vec t, // vector
                           int NMC, // integer
                           arma::vec u0, //matrix
                           arma::rowvec proposal_var, // row vector
                           double batch,
                           int trace){ // integer
  
  arma::mat fitted = f;
  arma::mat Z=z;
  arma::vec Y=y;
  arma::vec tau =t;
  int nMC = NMC; 
  arma::vec uold=u0;
  
  int q = Z.n_cols;
  int n = Z.n_rows;
  
  int i = 0;
  int j = 0; 
  // int k = 0; 
  int l = 0; 
  int naccept = 0;
  int index = 0;
  // int batch = 1;
  double sum = 0;
  double sumn = 0;
  // double stao = 0;
  double w = 0;
  double ep = 0; 
  double e0 = 0;
  double delta = 0;
  double increment = 0;
  double batch_length = 100.0;
  arma::mat out(nMC+2, q); // Last two lines = acceptance rates and updated proposal_var, respectively
  arma::vec e(q);
  arma::vec rate(q);
  arma::vec etae(n);
  arma::vec etaen(n);
  arma::vec index2(q);
  arma::rowvec var = proposal_var; // At beginning of EM algorithm, proposal_var = 1.0 for each variable
  arma::rowvec acc_rate(q);
  
  RNGScope scope;
  
  //initialize e 
  for(i = 0; i < q; i++){
    e(i) = uold(i);
    index2(i) = 0.0;
    acc_rate(i) = 0.0;
  }
  
  //iteratively update e
  while(naccept < nMC){
    
    for(j = 0; j < q; j++){
      
      // calculate etae
      etae = fitted + Z*e;
      
      // save current value of e(i)
      e0 = e(j);
      
      // generate w
      w = log(R::runif(0.0,1.0));
      
      // generate proposal e
      ep = R::rnorm(0.0, var(j));
      e(j) = ep;
      
      // calculate updated etae
      etaen = fitted + Z*e;
      
      // calculate ratio
      sum = sumn = 0;
      for(l = 0; l < n; l++){
        sum = sum + R::dbinom(Y(l), 1.0, exp(etae(l))/(1+exp(etae(l))), 1);
        sumn = sumn + R::dbinom(Y(l), 1.0, exp(etaen(l))/(1+exp(etaen(l))), 1);
      }
      
      sum = sum + R::dnorm4(e0, 0.0, 1.0, 1) + R::dnorm4(ep, 0.0, var(j), 1) ;
      sumn = sumn + R::dnorm4(ep, 0.0, 1.0, 1) + R::dnorm4(e0, 0.0, var(j), 1)  ;
        // Question: does dnorm4 use variance or standard deviation?
        // Assume same as R version of dnorm: use standard deviation
        // Should change "var" name to reflect this when this is confirmed
      
      // check for acceptance
      if(w < sumn - sum){
        // ep left in e(i)
        index2(j) = index2(j) + 1.0;
      }else{
        e(j) = e0;
      }
      
    } // End j for loop
    
    if(index == 1+5*naccept){
      out.row(naccept) = e.t();
      naccept++;
    }
    index++;
    
    // Update proposal variance
    if(index % (int)batch_length == 0){
      
      for(j = 0; j < q; j++){
        // Update batch information 
        batch = batch + batch_length;
        // Determine acceptance rate for latest batch
        acc_rate(j) = index2(j) / batch_length;
        
        // Update proposal variance (separate for each variable)
        // delta = min(0.01, (T_b)^(-1/2))
        increment = sqrt(1 / batch);
        if(increment < 0.01){
          delta = increment;
        }else{
          delta = 0.01;
        }
        
        // Update proposal standard deviation (variance?) based on acceptance rate for last batch
        // log(std_dev) +/- delta --> std_dev * exp(+/- delta)
        if(acc_rate(j) > 0.5){
          var(j) = var(j) * exp(-delta); // Or 2 delta - does dnorm4 use SD or var?
        }else if(acc_rate(j) < 0.4){
          var(j) = var(j) * exp(delta); // Or 2 delta - does dnorm4 use SD or var?
        }
        
        // Set min and max cap of log(standard deviation)
        if(log(var(j)) > 1.0){
          var(j) = exp(1.0);
        }else if(log(var(j)) < -1.0){
          var(j) = exp(-1.0);
        }
        
        // Re-set index2 - new acceptance rate for next batch
        index2(j) = 0.0;
        
      } // End j for loop (w/in if statement)
        
    } // End if(index % (int)batch_lenth)
      
  } // End while loop
  
  Rcout << "Final Acceptance Rate" << std::endl << acc_rate;
  Rcout << "Final Updated Proposal Var" << std::endl << var;
  
  out.row(nMC) = acc_rate; // Second-to-last row
  out.row(nMC+1) = var; // Last row
    
  return(wrap(out)); 
  
}
