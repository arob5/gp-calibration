# parameter_calibration_1d_example.R
# Run MCMC for parameter calibration on toy example with simple scalar function.
# Note that the Stan function f() defined in parameter_calibration_1d_example.stan
# should be equivalent to the function f() defined in this file. This file 
# simulates a dataset and then performs parameter calibration using the 
# sample likelihood and priors used to simulate the dataset. It therefore 
# assumes a perfect model. 
#
# Andrew Roberts
# Working Directory: /projectnb2/dietzelab/arober/test_code

library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#
# Settings
#

# The function that is acting as the computer model. Should be equivalent to 
# the function of the same name in the Stan file. 
f <- function(u) {
  return(sin(u))
}

# Hyperparameters (shape and rate ) for Gamma hyperprior on tau (precision parameter). Using values 
# taken from an actual PEcAn run. 
a <- 7.5
b <- 1.0









