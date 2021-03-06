% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/object_methods.R
\name{plot_mcmc}
\alias{plot_mcmc}
\title{Plot Diagnostics for MCMC Posterior Draws of the Random Effects}
\usage{
plot_mcmc(
  object,
  plots = c("all", "sample.path", "histogram", "cumsum", "autocorr"),
  grps = "all",
  vars = "all",
  numeric.grps = F,
  bin_width = NULL
)
}
\arguments{
\item{object}{an object of class \code{pglmmObj} output from either \code{\link{glmmPen}} 
or \code{\link{glmmPen_FineSearch}}.}

\item{plots}{a character string or a vector of character strings specifying which graphical
diagnostics to provide.}

\item{grps}{a character string or a vector of character strings specifying which groups 
should have diagnostics provided. The names of the groups match the input group factor levels.
Default is set to 'all' for all groups.}

\item{vars}{a character string or a vector of character strings specifying which variables
should have diagnostics provided. Tip: can find the names of the random effect variables in
the output sigma matrix found in the \code{pglmmObj} object. Default is set to
'all', which picks all variables with non-zero random effects.}

\item{numeric.grps}{if TRUE, specifies that the groups factor should be converted to numeric 
values. This option could be used to ensure that the organization of the groups is in the 
proper numeric order.}

\item{bin_width}{optional binwidth argument for \code{geom_histogram} from the \code{ggplot2} 
package. Default set to \code{NULL}, which specifies the default \code{geom_histogram} binwidth.}
}
\value{
a list of ggplot graphics, each faceted by group and random effect variable. 
Type of plots specified in the \code{plots} argument.
}
\description{
Provides graphical diagnostics of the random effect posterior draws from the (best) model.
Availabile diagnostics include the sample path, histograms, cummulative sums, and autocorrelation.
}
