# glmmPen
glmmPen package

Generalized linear mixed models (GLMMs) are popular for their flexibility and their ability to estimate population-level effects while accounting for between-group heterogeneity. While GLMMs are very
versatile, the specification of fixed and random effects is a critical part of the modeling process. The package glmmPen simultaneously selects fixed and random effects from high dimensional penalized generalized linear mixed models (pGLMMs) using the funcion `glmmPen`. This function `glmmPen` fits a sequence of pGLMMs and chooses the best model using one of several Bayesian Information Criterion (BIC)-derived selection criteria. The package can also fit single GLMM models (with or without penalization) using the function `glmm`. Model parameters are estimated using a Monte Carlo Expectation Conditional Maximization (MCECM) algorithm, which leverages Stan and RcppArmadillo
to increase computational efficiency.

Our package supports the penalty functions MCP, SCAD, and LASSO, and the distributional families Binomial, Gaussian, and Poisson (with canonical links). The available BIC-derived selection criteria include the BIC-ICQ, the regular BIC, and the hybrid BICh (see documentation for further details). The user interface of the package was designed to be similar to the popular lme4 R package, including the specification of the model equation with fixed and random effects.Tools available in the package include
automated tuning parameter selection and automated initialization of the random effect variance. 

Windows users must first run the following lines of code in order to properly install the package:

```
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars.win")
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS=-O3 -march=corei7 -mtune=corei7",
    "CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y",
    "CXX11FLAGS=-O3 -march=corei7 -mtune=corei7",
    file = M, sep = "\n", append = TRUE)
```

The package can then be installed using the following lines of code:

```
#install.packages("devtools")
devtools::install_github("hheiling/glmmPen")
```

The manual is available in the 'inst/' folder of the package and gives more specifics about the required and optional arguments for the main functions `glmmPen` and `glmm`.

The output object, of class pglmmObj, includes the following S3 methods of interest:
fixef, ranef, sigma, fitted, predict, print, and summary, which act similar to the same S3 methods for the output object of lme4 functions. 

There is also a plot_mcmc function that can perform some diagnostic plots on the posterior draw outputs.

Contact information: email hheiling@live.unc.edu