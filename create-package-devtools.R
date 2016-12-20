
library(devtools)

# Package name
setwd("~/GitHub/")
pkgName <- "rp.flm.test"

# Install from GitHub - auth token should work for anyone
# install_github(paste("egarpor/", pkgName, sep = ""), username = NULL, ref = "master", subdir = NULL, auth_token = "5d85d98bdd4ad74f3897192760c4bb20e311173e", host = "api.github.com")

# Load package
load_all(pkg = pkgName, export_all = TRUE)

# Documentation
document(pkg = pkgName)#, roclets = c("rd", "collate"))

# Run examples
# run_examples(pkgName)

# Check
check(pkg = pkgName, build_args = "--resave-data=best", document = FALSE, manual = TRUE)

# Build
build(pkgName, manual = TRUE)

# Install locally
install(pkgName, args = c("--no-multiarch", "--no-test-load"))

# Load package
do.call(library, list(pkgName))

