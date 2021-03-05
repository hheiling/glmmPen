library(rstan)
library(rstantools)
rstan_create_package(path = "C:/Users/hheiling/Documents/GitHub/stanpkg")

rstan_create_package(path = "C:/Users/hheiling/Documents/GitHub/stanpkg2",
                     stan_files = "C:/Users/hheiling/Documents/GitHub/glmmPen/inst/stan/binomial_logit_model.stan")


# rstantools::rstan_package_skeleton(path = "C:/Users/hheiling/Documents/GitHub/stanpkg_test")

library(devtools)
pkgbuild::compile_dll("C:/Users/hheiling/Documents/GitHub/stanpkg")
devtools::document("C:/Users/hheiling/Documents/GitHub/stanpkg")

# Necessary step to make sure Makevars in packages depending on rstan work
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars.win")
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS=-O3 -march=corei7 -mtune=corei7",
    "CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y",
    "CXX11FLAGS=-O3 -march=corei7 -mtune=corei7",
    file = M, sep = "\n", append = TRUE)
