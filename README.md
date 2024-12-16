# Package Dependencies
Almost all of the packages required by this are available on CRAN, and thus 
can easily be installed with a call to `install.packages()`. However, there 
are a couple exceptions. 

## support
The R package `support` provides methods for efficiently subsampling a 
probability distribution (e.g., thinning MCMC output). While older versions 
of the package are available on CRAN, the newer versions no longer satisfy 
CRAN's requirements. The version of `support` I am using was emailed to me 
directly by the package author, Simon Mak. To install it, unzip the 
*tar.g.z* file to create the `support` directory, and then copy this directory 
to your R package directory. Then, from within the command line within the 
package directory, run the commands `R CMD build support` and then 
`R CMD INSTALL packagename`. Pay attention to the printed output to see where 
the package is installed, and move the package to the correct R package 
directory if needed. The package should now be available via 
`library(support)`, so long as the package directory is included in 
`.libPaths()`.