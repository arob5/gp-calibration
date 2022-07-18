/**
* parameter_calibration_1d_example.stan
* Run MCMC for parameter calibration on toy example with simple scalar function. 
*
* Andrew Roberts
*
*/

functions {
  // This defines the stand-in for the process-based model/computer simulation 
  // that is used in this test. 
  real f(real u) {
    return( sin(u) ); 
  }
}

data {
  int<lower=1> n; // Number of observations
  int<lower=1> k; // Number of calibration parameters
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
  y ~ normal(f_u, y_sigma); 
}





