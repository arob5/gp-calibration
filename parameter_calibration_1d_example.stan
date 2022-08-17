/**
* parameter_calibration_1d_example.stan
* Run MCMC for parameter calibration on toy example with simple scalar function. 
* "Brute force" MCMC (no Gaussian Process approximation).
*
* Andrew Roberts
*
*/

functions {
#include gaussian_process_functions.stan

  // This defines the stand-in for the process-based model/computer simulation 
  // that is used in this test. 
  real f(real u) {
    return(u); 
  }
}

data {
  int<lower=1> n; // Number of observations
  vector[n] y;    // Observed data
  
  real<lower=0> a; // Shape parameter for Gamma hyperprior on tau
  real<lower=0> b; // Rate parameter for Gamma hyperprior on tau
  
  real u_mean;            // Mean for Gaussian prior on calibration parameter u
  real<lower=0> u_sigma; // Standard deviation for Gaussian prior on calibration parameter u
}

parameters {
  real u;        // Calibration parameters
  real<lower=0> tau; // Precision parameter for Gaussian likelihood
}

transformed parameters {
  real f_u; 
  real<lower=0> y_sigma; 
  
  f_u = f(u); 
  y_sigma = 1/sqrt(tau); 
}

model {
  u ~ normal(u_mean, u_sigma); 
  tau ~ gamma(a, b);
  // target += 0.5 * n * log(tau) - 0.5*tau*squared_distance_weighted(y', rep_row_vector(f_u, n), rep_vector(1.0, n)); 
  y ~ normal(f_u, y_sigma); 
}















