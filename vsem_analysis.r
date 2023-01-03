#
# vsem_analysis.r
# Testing parameter calibration methods for the "Very Simple Ecosystem Model" (VSEM)
# implemented in the BayesianTools R package. 
#
# Andrew Roberts
#

# Installing BayesianTools:
#    - Older versions of the package were on CRAN, but BayesianTools has been removed from CRAN.
#    - Below code installs development version. 
library(devtools)
devtools::install_github(repo = "florianhartig/BayesianTools", subdir = "BayesianTools", dependencies = T, build_vignettes = T)

install.packages("BayesianTools")
