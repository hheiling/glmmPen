# glmmPen
glmmPen package

This package fits a penalized generalized mixed model via Monte Carlo Expectation Conditional 
Minimization (MCECM). It can be downloaded using the following lines of code:

library(devtools)
library(remotes)
install_github("hheiling/glmmPen")

Family options: Currently, only the 'binomial' family option is available, and currently only the logit link for this family is fully functional.

The function 'glmmPen' is the user-friendly function in the package. While the documentation is not fully complete in the package, the available documentation should be sufficient for most user's purposes for the time being. 

The formula should be specified like glmer specifies their formulas. For instance, suppose we have a vector of response variables y, a matrix of covariates (not including an intercept) X, and a vector grp that specifies which group an observation belongs to (of same length as y). The formula would be specified as y ~ X + (X | grp). If the user has ideas about restricting the number of random effects, the X within the parentheses can be restricted.

If there are any covariates that the user thinks should always remain in the model (i.e. should not be penalized out of the fixed effects part of the model), then the argument fixef_noPen can be specified in the function (vector of same length as the number of columns of X, 0 if the variable should never be penalized, 1 otherwise).

The output object, of class pglmmObj, will have the following S3 methods of interest:
fixef.pglmmObj, ranef.glmmObj, fitted, predict, print, and summary, which act similar to the same S3 methods for the glmer output object. (There are a few others, but these are likely the most relevant for the time being).

There is also a plot_mcmc function that can perform some diagnostic plots on the posterior draw outputs if that is of interest. The documentation is not yet finished for this function, but hopefully the argument names make enough sensefor it to be usable

As this package is still in development, there are likely still some bugs to work out. 

Contact information: email hheiling@live.unc.edu