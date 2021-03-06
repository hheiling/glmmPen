% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glmmPen_FineSearch.R
\name{glmmPen_FineSearch}
\alias{glmmPen_FineSearch}
\title{Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
Minimization (MCECM) using a finer penalty grid search

\code{glmmPen_FineSearch} finds the best model from the selection results of a \code{pglmmObj} object 
created by \code{glmmPen}, identifies a more targeted grid search around the optimum lambda penalty
values, and performs model selection on this finer grid search.}
\usage{
glmmPen_FineSearch(
  object,
  tuning_options = selectControl(),
  idx_range = 2,
  optim_options = NULL,
  adapt_RW_options = NULL,
  trace = 0,
  BICq_posterior = NULL
)
}
\arguments{
\item{object}{an object of class \code{pglmmObj} created by \code{glmmPen}. This object must 
contain model selection results.}

\item{tuning_options}{a list of class selectControl resulting from \code{\link{selectControl}} 
containing model selection control parameters. See the \code{\link{selectControl}}
documentation for details. The user can specify their own fine grid search, or if the 
lambda0_seq and lambda1_seq arguments are left as \code{NULL}, the algorithm will automatically
select a fine grid search based on the best model from the previous selection. See Details for 
more information. Default value set to 1.}

\item{idx_range}{a positive integer that determines what positions within the sequence of the 
fixed and random effect lambda penalty parameters used in the previous coarse grid search
will be used as the new fixed and random effect lambda penalty parameter ranges. See Details 
for more information.}

\item{optim_options}{an optional list of class "optimControl" created from function \code{\link{optimControl}}
that specifies optimization parameters. If set to the default \code{NULL}, will use the 
optimization parameters used for the previous round of selection stored within the 
\code{pglmmObj} object.}

\item{adapt_RW_options}{an optional list of class "adaptControl" from function \code{\link{adaptControl}} 
containing the control parameters for the adaptive random walk Metropolis-within-Gibbs procedure. 
Ignored if \code{\link{optimControl}} parameter \code{sampler} is set to "stan" or "independence".
If set to the default \code{NULL}, will use the adaptive random walk paraters used for the 
previous round of selection stored within the \code{pglmmObj} object.}

\item{trace}{an integer specifying print output to include as function runs. Default value is 0. 
See Details for more information about output provided when trace = 0, 1, or 2.}

\item{BICq_posterior}{an optional character string specifying the file containing the posterior
draws used to calculate the BIC-ICQ selection criterion if such a file was created in the
previous round of selection.}
}
\value{
A reference class object of class \code{\link{pglmmObj}} for which many methods are 
available (e.g. \code{methods(class = "pglmmObj")})
}
\description{
Fit a Penalized Generalized Mixed Model via Monte Carlo Expectation Conditional 
Minimization (MCECM) using a finer penalty grid search

\code{glmmPen_FineSearch} finds the best model from the selection results of a \code{pglmmObj} object 
created by \code{glmmPen}, identifies a more targeted grid search around the optimum lambda penalty
values, and performs model selection on this finer grid search.
}
\details{
The \code{glmmPen_FineSearch} function extracts the data, the penalty information (penalty type,
gamma_penalty, and alpha), and some other argument specifications from the \code{pglmmObj} object
created during a previous round of model selection. In this finer grid search, the user has
the ability to make the following adjustments: the user can change the BIC option used for selection,
any optimization control parameters, or any adaptive random walk parameters (if the sampler
specified in the optimization parameters is "random_walk"). The user could manually specify the 
lambda penalty grid to search over within the \code{\link{selectControl}} control parameters,
or the user could let the \code{glmmPen_FineSearch} algorithm calculate a finer grid search 
automatically (see next paragraph for details).

If the sequences of lambda penalty values are left unspecified in the \code{\link{selectControl}} 
tuning options, the \code{glmmPen_FineSearch} algorithm performs the following steps to find
the finer lambda grid search: (i) The lambda combination from the best model is identified 
from the earlier selection results saved in the \code{pglmmObj} object. (ii) For the fixed and
random effects separately, the new max and min lambda values are the lambda values \code{idx_range} 
positions away from the best lambda in the original lambda sequences for the fixed and random
effects. 

Trace details: The value of 0 outputs some general updates for each EM iteration (iteration number EM_iter,
number of MCMC draws nMC, average Euclidean distance between current coefficients and coefficients
from t iterations back EM_diff, and number of non-zero coefficients Non0 Coef). The value of 1
additionally outputs the updated coefficients, updated covariance matrix values, and the
number of coordinate descent iterations used for the M step for each
EM iteration. The value of 2 outputs all of the above plus gibbs acceptance rate information
for the adaptive random walk and independence samplers and the updated proposal standard deviation
for the adaptive random walk.
}
