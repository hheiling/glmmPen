
# Updating Namespace and package

library(devtools)
library(remotes)

# Download the latest version of the glmmPen package:

# other branch
# install_github("hheiling/glmmPen", ref = "alt_branch", force = TRUE)
# master branch
# install_github("hheiling/glmmPen", force = TRUE)
# adaptive branch
install_github("hheiling/glmmPen", ref = "adaptive", force = TRUE)

# Update Namespace and man files:

# Package working directory: "C:/Users/hheiling/Documents/GitHub/glmmPen"
# Run document line when new @export functions added to package
devtools::document("C:/Users/hheiling/Documents/GitHub/glmmPen")

# devtools::build(pkg = "C:/Users/hheiling/Documents/GitHub/glmmPen",
#                 path = "C:/Users/hheiling/Documents/GitHub/glmmPen_extras")