// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// Adaptive Random Walk Metropolis within Gibbs 

//' @export
// [[Rcpp::export]]
NumericMatrix sample_mc_gibbs_rw(arma::mat f, // matrix
                                 arma::mat z, // matrix
                                 arma::vec y, // vector
                                 arma::vec t, // vector
                                 int NMC, // integer
                                 arma::vec u0, //matrix
                                 arma::rowvec proposal_SD, // row vector
                                 double batch,
                                 double batch_length,
                                 double offset,
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
  int l = 0; 
  int naccept = 0;
  int index = 0;
  double sum = 0;
  double sumn = 0;
  double w = 0;
  double ep = 0; 
  double e0 = 0;
  double delta = 0;
  double increment = 0;
  // double batch_length = 500.0; // Note: Total "retained" draws = total draws / 5
  arma::mat out(nMC+2, q); // Last two lines = acceptance rates and updated proposal_SD, respectively
  arma::vec e(q);
  arma::vec rate(q);
  arma::vec etae(n);
  arma::vec etaen(n);
  arma::vec accept_index(q);
  arma::rowvec SD = proposal_SD; // At beginning of EM algorithm, proposal_SD = 1.0 for each variable
  arma::rowvec acc_rate(q);
  
  RNGScope scope;
  
  //initialize e 
  for(i = 0; i < q; i++){
    e(i) = uold(i);
    accept_index(i) = 0.0;
    acc_rate(i) = 0.0;
  }
  
  //iteratively update e 
  while(naccept < nMC){
    
    for(j = 0; j < q; j++){
      
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
      for(l = 0; l < n; l++){
        sum = sum + R::dbinom(Y(l), 1.0, exp(etae(l))/(1+exp(etae(l))), 1);
        sumn = sumn + R::dbinom(Y(l), 1.0, exp(etaen(l))/(1+exp(etaen(l))), 1);
      }
      
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
      
    } // End j for loop
    
    if(index == 1+5*naccept){
      out.row(naccept) = e.t();
      naccept++;
    }
    index++;
    
    // Update proposal variance
    if(index % (int)batch_length == 0){ // if index = multiple of batch_length
      
      for(j = 0; j < q; j++){
        // Update batch information 
        batch = batch + batch_length;
        // Determine acceptance rate for latest batch
        acc_rate(j) = accept_index(j) / batch_length;
        
        // Update proposal SD (separate for each variable)
        // delta = min(0.01, (T_b)^(-1/2))
        increment = 1 / sqrt(batch + offset);
        if(increment < 0.01){
          delta = increment;
        }else{
          delta = 0.01;
        }
        
        // Update proposal standard deviationbased on acceptance rate for last batch
        // log(std_dev) +/- delta --> std_dev * exp(+/- delta)
        if(acc_rate(j) > 0.5){
          SD(j) = SD(j) * exp(-delta); 
        }else if(acc_rate(j) < 0.4){
          SD(j) = SD(j) * exp(delta); 
        }
        
        // if(batch < 2000.0){
        //   Rprintf("Intermediate Proposal SD(%u): %f \n", j, SD(j));
        // }
        
        // Set min and max cap of log(standard deviation)
        if(log(SD(j)) > 1.0){
          SD(j) = exp(1.0);
        }else if(log(SD(j)) < -1.0){
          SD(j) = exp(-1.0);
        }
        
        // Re-set accept_index - new acceptance rate for next batch
        accept_index(j) = 0.0;
        
      } // End j for loop (w/in if statement)
      
    } // End if(index % (int)batch_lenth)
    
  } // End while loop
  
  Rcout << "Final Acceptance Rate" << std::endl << acc_rate;
  Rcout << "Final Updated Proposal SD" << std::endl << SD;
  
  out.row(nMC) = acc_rate; // Second-to-last row
  out.row(nMC+1) = SD; // Last row
  
  return(wrap(out)); 
  
}

// Adaptive Random Walk Metropolis within Gibbs with Random Scan

//' @export
// [[Rcpp::export]]
NumericMatrix sample_mc_gibbs_rw_rs(arma::mat f, // matrix
                                 arma::mat z, // matrix
                                 arma::vec y, // vector
                                 arma::vec t, // vector
                                 int NMC, // integer
                                 arma::vec u0, //matrix
                                 arma::rowvec proposal_SD, // row vector
                                 double batch,
                                 double batch_length,
                                 double offset,
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
  int l = 0; 
  int k = 0;
  int naccept = 0;
  int index = 0;
  double sum = 0;
  double sumn = 0;
  double w = 0;
  double ep = 0; 
  double e0 = 0;
  double delta = 0;
  double increment = 0;
  // double batch_length = 500.0; // Note: Total "retained" draws = total draws / 5
  arma::mat out(nMC+2, q); // Last two lines = acceptance rates and updated proposal_SD, respectively
  arma::vec e(q);
  arma::vec rate(q);
  arma::vec etae(n);
  arma::vec etaen(n);
  arma::vec accept_index(q);
  arma::rowvec SD = proposal_SD; // At beginning of EM algorithm, proposal_SD = 1.0 for each variable
  arma::rowvec acc_rate(q);
  arma::uvec var_index(q);
  arma::uvec samp(1);
  // arma::vec scan_index(q); // Keep track of how many times each variable is randomly selected
  
  RNGScope scope;
  
  //initialize e 
  for(i = 0; i < q; i++){
    e(i) = uold(i);
    accept_index(i) = 0.0;
    acc_rate(i) = 0.0;
    var_index(i) = i;
    // scan_index(i) = 0.0;
  }
  
  //iteratively update e 
  while(naccept < nMC){
    
    // var_index = randperm(q); 
    // j = var_index(1);
    samp = RcppArmadillo::sample(var_index, q, 0);
    
    // if(index < 10.0){
    //   Rcout << "Random Scan Order" << std::endl << samp;
    // }
    
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
      for(l = 0; l < n; l++){
        sum = sum + R::dbinom(Y(l), 1.0, exp(etae(l))/(1+exp(etae(l))), 1);
        sumn = sumn + R::dbinom(Y(l), 1.0, exp(etaen(l))/(1+exp(etaen(l))), 1);
      }
      
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
    
    
    // Update proposal variance
    if(index % (int)batch_length == 0){ // if index = multiple of batch_length
      
      for(j = 0; j < q; j++){
        // Update batch information 
        batch = batch + batch_length;
        // Determine acceptance rate for latest batch
        acc_rate(j) = accept_index(j) / batch_length;
        
        // Update proposal SD (separate for each variable)
        // delta = min(0.01, (T_b)^(-1/2))
        increment = 1 / sqrt(batch + offset);
        if(increment < 0.01){
          delta = increment;
        }else{
          delta = 0.01;
        }
        
        // Update proposal standard deviation based on acceptance rate for last batch
        // log(std_dev) +/- delta --> std_dev * exp(+/- delta)
        if(acc_rate(j) > 0.5){
          SD(j) = SD(j) * exp(-delta); 
        }else if(acc_rate(j) < 0.4){
          SD(j) = SD(j) * exp(delta); 
        }
        
        // if(batch < 2000.0){
        //   Rprintf("Intermediate Proposal SD(%u): %f \n", j, SD(j));
        // }
        
        // Set min and max cap of log(standard deviation)
        if(log(SD(j)) > 1.0){
          SD(j) = exp(1.0);
        }else if(log(SD(j)) < -1.0){
          SD(j) = exp(-1.0);
        }
        
        // Re-set accept_index and scan_index - new acceptance rate for next batch
        accept_index(j) = 0.0;
        // scan_index(j) = 0.0;
        
      } // End j for loop (w/in if statement)
      
    } // End if(index % (int)batch_lenth)
    
  } // End while loop
  
  Rcout << "Final Acceptance Rate" << std::endl << acc_rate;
  Rcout << "Final Updated Proposal SD" << std::endl << SD;
  
  out.row(nMC) = acc_rate; // Second-to-last row
  out.row(nMC+1) = SD; // Last row
  
  return(wrap(out)); 
  
}

// Adaptive (Burn-in Only) Random Walk Metropolis within Gibbs with Random Scan

//' @export
// [[Rcpp::export]]
NumericMatrix sample_mc_gibbs_adapt_rw(arma::mat f, // matrix
                                    arma::mat z, // matrix
                                    arma::vec y, // vector
                                    arma::vec t, // vector
                                    int NMC, // integer
                                    arma::vec u0, //matrix
                                    arma::rowvec proposal_SD, // row vector
                                    double batch,
                                    double batch_length,
                                    double offset,
                                    double burnin,
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
  int l = 0; 
  int k = 0;
  int naccept = 0;
  int index = 0;
  double sum = 0;
  double sumn = 0;
  double w = 0;
  double ep = 0; 
  double e0 = 0;
  double delta = 0;
  double increment = 0;
  // double batch_length = 500.0; // Note: Total "retained" draws = total draws / 5
  arma::mat out(nMC+2, q); // Last two lines = acceptance rates and updated proposal_SD, respectively
  arma::vec e(q);
  arma::vec rate(q);
  arma::vec etae(n);
  arma::vec etaen(n);
  arma::vec accept_index(q);
  arma::rowvec SD = proposal_SD; // At beginning of EM algorithm, proposal_SD = 1.0 for each variable
  arma::rowvec acc_rate(q);
  arma::uvec var_index(q);
  arma::uvec samp(1);
  arma::rowvec final_acc_rate(q);
  
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
    
    // var_index = randperm(q); 
    // j = var_index(1);
    
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
      for(l = 0; l < n; l++){
        sum = sum + R::dbinom(Y(l), 1.0, exp(etae(l))/(1+exp(etae(l))), 1);
        sumn = sumn + R::dbinom(Y(l), 1.0, exp(etaen(l))/(1+exp(etaen(l))), 1);
      }
      
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
    
    
    // Update proposal variance
    if((index % (int)batch_length == 0) && (batch < burnin)){ // if index = multiple of batch_length
      
      for(j = 0; j < q; j++){
        // Update batch information 
        batch = batch + 1.0;
        // Determine acceptance rate for latest batch
        acc_rate(j) = accept_index(j) / batch_length;
        
        // Update proposal SD (separate for each variable)
        // delta = min(0.01, (T_b)^(-1/2))
        increment = 1 / sqrt(batch*batch_length + offset);
        if(increment < 0.005){
          delta = increment;
        }else{
          delta = 0.005;
        }
        
        // Update proposal standard deviation based on acceptance rate for last batch
        // log(std_dev) +/- delta --> std_dev * exp(+/- delta)
        if(acc_rate(j) > 0.5){
          SD(j) = SD(j) * exp(-delta); 
        }else if(acc_rate(j) < 0.4){
          SD(j) = SD(j) * exp(delta); 
        }
        
        // Set min and max cap of log(standard deviation)
          // No longer need if only adapting during burn-in
        // if(log(SD(j)) > log(1.5)){
        //   SD(j) = 1.5;
        // }else if(log(SD(j)) < -log(1.5)){
        //   SD(j) = 1/1.5;
        // }
        
        // Re-set accept_index - new acceptance rate for next batch
        accept_index(j) = 0.0;
        
      } // End j for loop (w/in if statement)
      
    } // End if(index % (int)batch_lenth)
    
  } // End while loop
  
  for(j = 0; j < q; j++){
    final_acc_rate(j) = accept_index(j) / (index - burnin*batch_length);
  }
  
  Rcout << "Final Acceptance Rate" << std::endl << acc_rate;
  Rcout << "Final Updated Proposal SD" << std::endl << SD;
  
  out.row(nMC) = final_acc_rate; // Second-to-last row
  out.row(nMC+1) = SD; // Last row
  
  return(wrap(out)); 
  
}


