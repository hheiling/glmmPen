
## Update glmmPen package - chose to update none of the packages listed
library(devtools)
library(remotes)
install_github("hheiling/glmmPen", force = TRUE)

## Needed libraries
library(glmmPen)
library(stringr)
library(ggplot2)
library(reshape2)

## Load simulation results
load("glmer_Pajor_diagnostics.RData")
output_1 = output[c(101:200)] # sd_ranef = 1.0 cases

K = 5
q = 3

## Code to develop example data - no need to run, have .RData
# fit_dat_out = function(dat){
#   
#   set.seed(dat$seed)
#   # Set gibbs = T
#   fit_glmmPen = fit_dat(dat, lambda0 = 0, lambda1 = 0, conv = 0.001, nMC = 100, 
#                         family = "binomial", trace = 0, penalty = "grMCP",
#                         alpha = 1, nMC_max = 4000, 
#                         returnMC = T, ufull = NULL, coeffull = NULL, gibbs = T, maxitEM = 150, 
#                         ufullinit = NULL) 
#   
#   return(fit_glmmPen)
# }

# test_fit = fit_dat_out(output_1[[1]]$dat)
# save(test_fit, file = "test_fit.RData")

## Testing Adaptive Metropolis-within-Gibbs functions
AMCMC_test = function(fit_glmmPen, dat){
  
  # Set variables
  K = 5
  q = 3
  
  # Test acceptance rates for new Adaptive Metropolis-within-gibbs
  post_list = sample.mc3(fit = fit_glmmPen$fit, cov = fit_glmmPen$sigma, y = dat$y, X = dat$X,
                         Z = fit_glmmPen$extra$Znew2, nMC = 10^4, family = "binomial", 
                         group = dat$group, d = nlevels(dat$group), okindex = fit_glmmPen$extra$ok_index,
                         nZ = fit_glmmPen$extra$Znew2, gibbs = T, uold = fit_glmmPen$u, trace = 2,
                         proposal_var = matrix(1.0, nrow = K, ncol = q), batch = 0.0)
  
  post_U = post_list$u0
  
  gibbs_accept_rate = post_list$gibbs_accept_rate

  proposal_var = post_list$proposal_var
  
  return(list(dat = dat, post_U = post_U, gibbs_accept_rate = gibbs_accept_rate, proposal_var = proposal_var))
}

## Examine performance of Adaptive Metropolis-within-Gibbs
load("test_fit.RData") # gives test_fit object
test_obj1 = AMCMC_test(fit_glmmPen = test_fit, dat = output_1[[1]]$dat)

## Evaluate the performance of the chain

mcmc_diagnostics = function(post_U){ #plot_mcmc.pglmmObj
  
  d = 5
  var_num = 3
  
  type = c("sample.path","histogram","cumsum","autocorr")
  
  U_keep = post_U
  var_names = c("Intercept","X1","X2")
  grp_names = str_c("grp", 1:5)
  var_str = rep(var_names, each = d)
  grp_str = rep(grp_names, times = q)
  U_cols = str_c(var_str, ":", grp_str)
  
  vars = "all"
  grps = "all"
  
  U_t = data.frame(U_keep, t = 1:nrow(U_keep))
  colnames(U_t) = c(U_cols, "t")
  U_long = melt(U_t, id = "t")
  U_plot = data.frame(U_long, var_names = rep(var_names, each = d*nrow(U_keep)),
                      grp_names = rep(rep(grp_names, each = nrow(U_keep)), times = var_num))
  
  plots_return = list()
  
  if("sample.path" %in% type){
    plot_sp = ggplot(U_plot, mapping = aes(x = t, y = value)) + geom_path() +
      facet_grid(var_names ~ grp_names) + xlab("iteration t") + ylab("draws")
    
    plots_return$sample_path = plot_sp
  }
  if("histogram" %in% type){
    hist_U = ggplot(U_plot) + geom_histogram(mapping = aes(x = value)) + 
      facet_grid(var_names ~ grp_names) + xlab("draws")
    
    plots_return$histogram = hist_U
  }
  if("cumsum" %in% type){
    U_means = colMeans(U_keep)
    U_means = data.frame(rbind(U_means))[rep.int(1L, nrow(U_keep)), , drop = FALSE]
    U_tmeans = apply(U_keep, 2, cumsum) / 1:nrow(U_keep)
    U_tdiff = U_tmeans - U_means
    U_cumsum = apply(U_tdiff, 2, cumsum)
    U_t = data.frame(U_cumsum, t = 1:nrow(U_cumsum))
    colnames(U_t) = c(colnames(U_keep), "t")
    U_long = melt(U_t, id = "t")
    U_plot = data.frame(U_long, var_names = rep(var_names, each = d*nrow(U_keep)),
                        grp_names = rep(rep(grp_names, each = nrow(U_keep)), times = var_num)) 
    plot_cumsum = ggplot(U_plot) + geom_smooth(mapping = aes(x = t, y = value), color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_grid(var_names ~ grp_names) + xlab("iteration t") + ylab("Cumulative Sum")
    
    plots_return$cumsum = plot_cumsum
  }
  if("autocorr" %in% type){
    grp_index = rep(grp_names, times = var_num)
    var_index = rep(var_names, each = d)
    for(j in 1:ncol(U_keep)){
      ACF = acf(U_keep[,j], plot=F, lag.max = 40)
      ACF_df = with(ACF, data.frame(lag,acf))
      ACF_df$grp_names = grp_index[j]
      ACF_df$var_names = var_index[j]
      if(j == 1){
        ACF_all = ACF_df
      }else{
        ACF_all = rbind(ACF_all, ACF_df)
      }
    }
    
    plot_acf = ggplot(data = ACF_all, mapping = aes(x = lag, y = acf)) +
      geom_hline(mapping = aes(yintercept = 0)) + 
      geom_segment(mapping = aes(xend = lag, yend = 0)) +
      facet_grid(var_names ~ grp_names)
    
    plots_return$autocorr = plot_acf
  }
  
  return(plots_return)
  
}

## Plots: sample_path, histogram, cumsum, autocorr
test_plots = mcmc_diagnostics(post_U = test_obj1$post_U)

