
#' @export
select_tune = function(dat, nMC, lambda0_range,lambda1_range, family, 
                       penalty, returnMC = T,
                       conv, nMC_max, trace = 0, ufull = NULL, coeffull = NULL, 
                       gibbs = T, maxitEM = 100, alpha = 1, t = 10,
                       M = 10^4, MwG_sampler = c("random_walk","independence"), 
                       adapt_RW_options = adaptControl()){
  
  n1 = length(lambda0_range)
  n2 = length(lambda1_range)
  BICold = BIC = BIC2old2 = BIC2= Inf
  BIChold = BICh = Inf
  res = matrix(0, n1*n2, 9)
  # res2 = matrix(0, n2, 4) 
  coef = NULL
  outl = list()
  for(i in 1:n1){
    for(j in 1:n2){
      if(i == 1 & j == 1){
        # start from the full model values
        coeffull0 = coeffull
        ufullinit = ufull
      }else if(i != 1 & j == 1){
        # take the previous i value for j = 1
        coeffull0 = coeffulli0
        ufullinit = ui0
        rm(out)
      }else{
        # take the previous j values
        coeffull0 = out$coef
        ufullinit = out$u
        rm(out)
      }
      
      gc()
      print(sprintf("lambda0 %f lambda1 %f", lambda0_range[i], lambda1_range[j]))
      out = try(fit_dat(dat, lambda0 = lambda0_range[i], lambda1 = lambda1_range[j], 
                        nMC = nMC, family = family, penalty = penalty,
                        trace = trace, conv = conv, nMC_max = nMC_max, 
                        ufull = ufull, coeffull = coeffull0, 
                        gibbs = gibbs, maxitEM = maxitEM + maxitEM*(j==1), 
                        returnMC = returnMC, ufullinit = ufullinit, alpha = alpha, t = t,
                        M = M, MwG_sampler = MwG_sampler, adapt_RW_options = adapt_RW_options))
      
      
      if(is.character(out)) next
      
      # for restarting the next i for j = 1
      if(j == 1){
        coeffulli0 = out$coef
        ui0 = out$u 
      }
      
      fout = list()
      
      BICh = out$BICh
      print(BICh)
      if(!is.numeric(BICh)) BICh = Inf
      if(length(BICh) != 1) BICh = Inf
      
      if(BICh < BIChold){
        fout[["BICh"]] = out
        BIChold = BICh
      }
      
      BIC = out$BIC
      print(BIC)
      if(!is.numeric(BIC)) BIC = Inf
      if(length(BIC) != 1) BIC = Inf
      
      if(BIC < BICold){
        fout[["BIC"]] = out
        BICold = BIC
      }
      
      res[(i-1)*(n2)+j,1] = lambda0_range[i]
      res[(i-1)*(n2)+j,2] = lambda1_range[j]
      res[(i-1)*(n2)+j,3] = BICh
      # res[(i-1)*(n2)+j,4] = out$llb
      res[(i-1)*(n2)+j,4] = out$ll
      res[(i-1)*(n2)+j,5] = sum(out$coef[2:ncol(dat$X)] !=0)
      res[(i-1)*(n2)+j,6] = sum(diag(out$sigma) !=0)
      res[(i-1)*(n2)+j,7] = sum(out$coef !=0)
      # if(maxitEM > 1) res[(i-1)*(n2)+j,8] = loglik(dat = dat, coef = out$coef, u0 = out$u, nMC = 100000, J = out$J)
      if(maxitEM > 1) res[(i-1)*(n2)+j,8] = out$ll
      # res[(i-1)*(n2)+j,9] = -2*res[(i-1)*(n2)+j,8] + log(length(dat$y))*sum(out$coef!=0)
      res[(i-1)*(n2)+j,9] = BIC
      outl[[(i-1)*(n2)+j]] = 1
      coef = rbind(coef, out$coef)
      print(res[(i-1)*(n2)+j,])
      colnames(res) = c("lambda0","lambda1","BICh","LogLik","Non_0_fef","Non_0_ref", "Non_0_coef",
                        "ll_Pajor","BIC")
      
    }
  }
  
  # print("Start of ncvreg code for internal comparison")
  # 
  # library(ncvreg)
  # if(n2 > 1){ 
  #   outcv = try(cv.ncvreg(y = dat$y, X = dat$X[,-1], family = family, nlambda = n2, alpha = alpha))
  #   if(is.character(outcv)){
  #     out2 = NULL
  #   }else{
  #     bicglm = 2*outcv$fit$loss + 
  #       log(outcv$fit$n)*as.numeric(apply(outcv$fit$beta, 2, FUN = function(x) sum(x[-1]!=0)))
  #     lambda.min = outcv$fit$lambda[which.min(bicglm)]
  #     out2 = ncvreg(y = dat$y, X = dat$X[,-1], family = family, lambda = lambda.min, alpha = alpha)
  #   }
  #   outind = list()
  #   for(jj in 1:max(as.numeric(dat$group))){
  #     outcv = try(cv.ncvreg(y = dat$y[dat$group == jj], X = dat$X[dat$group == jj,-1], family = family, nlambda = n2, alpha = alpha))
  #     
  #     if(is.character(outcv)){
  #       outind[[jj]] = NULL
  #     }else{	
  #       bicglm = 2*outcv$fit$loss + log(outcv$fit$n)*as.numeric(apply(outcv$fit$beta, 2, FUN = function(x) sum(x[-1]!=0)))
  #       lambda.min = outcv$fit$lambda[which.min(bicglm)]
  #       outind[[jj]] = try(ncvreg(y = dat$y[dat$group == jj], X = dat$X[dat$group == jj,-1], family = family, lambda = lambda.min, alpha = alpha))
  #       if(is.character(outind[[jj]])){
  #         outind[[jj]] = list(beta = outcv$fit$beta[,which.min(bicglm)])
  #       }
  #     }
  #   }
  #   
  # }else{
  #   lambda.min = 0.000001
  #   outcv = NULL
  #   outind = NULL
  #   out2 = ncvreg(y = dat$y, X = dat$X[,-1], family = family, lambda = lambda.min, alpha = alpha)
  # }
  
  # return(list(result = res, out = fout, res2 = outcv, out2 = out2, coef = coef, outind = outind))
  
  return(list(result = res, out = fout, coef = coef))
}

