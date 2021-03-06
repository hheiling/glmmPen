% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control_options.R
\name{optimControl}
\alias{optimControl}
\title{Control of Penalized Generalized Linear Mixed Model Fitting

Constructs the control structure for the optimization of the penalized mixed model fit algorithm.}
\usage{
optimControl(
  conv_EM = 0.0015,
  conv_CD = 5e-04,
  nMC_burnin = NULL,
  nMC_start = NULL,
  nMC_max = NULL,
  nMC_report = 5000,
  maxitEM = NULL,
  maxit_CD = 50,
  M = 10000,
  t = 2,
  mcc = 2,
  sampler = c("stan", "random_walk", "independence"),
  var_start = "recommend",
  max_cores = 1
)
}
\arguments{
\item{conv_EM}{a non-negative numeric convergence criteria for the convergence of the 
EM algorithm. Default is 0.001. EM algorithm is considered to have converge if the average Euclidean 
distance between the current coefficient estimates and the coefficient estimates from 
\code{t} EM iterations back is less than \code{conv_EM} \code{mcc} times in a row.
See \code{t} and \code{mcc} for more details.}

\item{conv_CD}{a non-negative numeric convergence criteria for the convergence of the 
grouped coordinate descent loop within the M step of the EM algorithm. Default 0.0001.}

\item{nMC_burnin}{positive integer specifying the number of posterior draws to use as
burnin for each E step in the EM algorithm. If set to \code{NULL}, the algorithm inputs
the following defaults: Default 250 when the number of random effects 
predictors is less than or equal to 10; default 100 otherwise. Function will not allow \code{nMC_burnin}
to be less than 100.}

\item{nMC_start}{a positive integer for the initial number of Monte Carlo draws. If set to
\code{NULL}, the algorithm inputs the following defaults: Default 250 when 
the number of random effects predictors is less than or equal to 10; default 100 otherwise.}

\item{nMC_max}{a positive integer for the maximum number of allowed Monte Carlo draws used
in each step of the EM algorithm. If set to \code{NULL}, the algorithm inputs the following 
defaults: When the number of random effect predictors is 10 or less, 
Default is set to 5000 when no selection is performed and 2500 when selection is performed.
Default is set to 1000 when the number of random effect predictors is greater than 10.}

\item{nMC_report}{a positive integer for the number of posterior draws to save from the final
model. These posterior draws can be used for diagnostic purposes, see \code{\link{plot_mcmc}}}

\item{maxitEM}{a positive integer for the maximum number of allowed EM iterations. 
If set to \code{NULL}, then the algorithm inputs the following defaults:
Default equals 50 for the Binomial and Poisson families, 100 for the Gaussian family.}

\item{maxit_CD}{a positive integer for the maximum number of allowed interations for the
coordinate descent algorithms used within the M-step of each EM iteration. Default equals 50.}

\item{M}{positive integer specifying the number of posterior draws to use within the 
Pajor log-likelihood calculation. Default is 10^4; minimum allowed value is 5000.}

\item{t}{the convergence criteria is based on the average Euclidean distance between 
the most recent coefficient estimates and the coefficient estimates from t EM iterations back.
Positive integer, default equals 2.}

\item{mcc}{the number of times the convergence critera must be met before the algorithm is
seen as having converged (mcc for 'meet condition counter'). Default set to 2. Value retricted 
to be no less than 2.}

\item{sampler}{character string specifying whether the posterior draws of the random effects
should be drawn using Stan (default, from package rstan) or the Metropolis-within-Gibbs procedure 
incorporating an adaptive random walk sampler ("random_walk") or an
independence sampler ("independence"). If using the random walk sampler, see \code{\link{adaptControl}}
for some additional control structure parameters.}

\item{var_start}{either the character string "recommend" or a positive number specifying the 
starting values to initialize the variance of the covariance matrix. Default "recommend" first
fits a simple model with a fixed and random intercept only using a Laplace estimate. The 
random intercept variance estimate from this model is then multiplied by 2 and used as the 
starting variance.}

\item{max_cores}{integer describing the number of cores available for computation. If 
\code{max_cores} is specified to be greater than 1 and the sampler is specified as "stan", then 
parallel computation using multiple cores is used to calculate the Stan MCMC samples within
each E step. The package authors do not advise using parallelization in the E step
unless \code{nMC_max} is set to a large number, such as 10^4 or more.}
}
\value{
Function returns a list (inheriting from class "\code{optimControl}") 
containing fit and optimization criteria values used in optimization routine.
}
\description{
Control of Penalized Generalized Linear Mixed Model Fitting

Constructs the control structure for the optimization of the penalized mixed model fit algorithm.
}
\details{
When the \code{optim_options} arugment in \code{\link{glmm}} and \code{\link{glmmPen}}
is set to "recommend", the default settings discussed in the given \code{optimControl} arguments are
used. These default settings depend on both the family of the data structure and the number 
of random effects predictors specified for use. If the user specifies 
\code{optim_options = optimControl()} with any argument specifications, no additional
adjustments will be performed on the arguments based on family or random effect predictors.
}
