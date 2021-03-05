
# Updating Namespace and package

library(devtools)
library(remotes)

# Download the latest version of the glmmPen package:

devtools::load_all(path = "C:/Users/hheiling/Documents/GitHub/glmmPen")

devtools::install(pkg = "C:/Users/hheiling/Documents/GitHub/glmmPen")

# other branch
# install_github("hheiling/glmmPen", ref = "alt_branch", force = TRUE)
# master branch
# install_github("hheiling/glmmPen", force = TRUE)
# M_Step branch
# install_github("hheiling/glmmPen", ref = "M_Step", force = TRUE)

# unlink("C:/Users/hheiling/Documents/R/win-library/3.6/00LOCK-glmmPen", recursive = TRUE)

# Update Namespace and man files:

# Package working directory: "C:/Users/hheiling/Documents/GitHub/glmmPen"
# Run document line when new @export functions added to package
# pkgbuild::compile_dll("C:/Users/hheiling/Documents/GitHub/glmmPen")
devtools::document("C:/Users/hheiling/Documents/GitHub/glmmPen")
# roxygen2::roxygenize("C:/Users/hheiling/Documents/GitHub/glmmPen")

devtools::build(pkg = "C:/Users/hheiling/Documents/GitHub/glmmPen",
                path = "C:/Users/hheiling/Documents/Longleaf")

library(rstantools)
# use_rstan(pkgdir = "C:/Users/hheiling/Documents/GitHub/glmmPen")
rstan_config(pkgdir = "C:/Users/hheiling/Documents/GitHub/glmmPen")

# Build manual
devtools::build_manual()

# Tips for building vignette:
# https://r-pkgs.org/vignettes.html
# usethis::use_vignette("glmmPen")
# https://kbroman.org/pkg_primer/pages/vignettes.html
devtools::build_vignettes()

###########################################################################################

# Errors and fixes

# Error: (replace "stan_fit4lm_mod" with appropriate model)
## Error in Rcpp::loadModule(module = “stan_fit4lm_mod”, what = TRUE, env = ns, : *
## Unable to load module “stan_fit4lm_mod”: Failed to initialize module pointer: Error in FUN(X[[i]], …): no such symbol _rcpp_module_boot_stan_fit4lm_mod in package rstanlmdemo*
# Fix:
## Add back @useDynLib glmmPen to fit_dat_B documentation
## Add the following line manually to the NAMESPACE file: 
## useDynLib(glmmPen, .registration = TRUE)
###########################################################################################
