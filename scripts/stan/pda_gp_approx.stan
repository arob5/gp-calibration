functions {
#include gaussian_process_functions.stan  
}

data {
  int<lower=1> n; // Number of observations for output variable
  int<lower=1> k; // Number of calibration parameters
  int<lower=1> N; // Number of design points (i.e. knots)
  matrix[N, k] X; // Design matrix
  vector[N] y;    // Process model evaluations at design points
  
  real<lower=0> a; // Shape parameter for Gamma hyperprior on tau
  real<lower=0> b; // Rate parameter for Gamma hyperprior on tau
  
  vector<lower=0>[k] gp_rho; // Vector of lengthscale parameters
  real<lower=0> gp_alpha; // Marginal standard deviation
  real<lower=0> gp_sigma; // Nugget (standard deviation)
  vector[k] gp_beta; // Vector of regression coefficients for linear mean function
}

transformed data {
  matrix[N, N] L; // Cholesky factor of covariance matrix
  vector[N] l;    // Pre-computed portion of predictive mean equation
  
  L = cholesky_decompose(cov_exp_quad_same(X, gp_rho, gp_alpha, gp_sigma)); 
  l = L' \ (L \ y); 
}

parameters {
  real<lower=0> tau; // Precision parameter in Gaussian model 
  vector[k] u; // Calibration parameters
}

model {
  tau ~ gamma(a, b); 
  // Priors on u are read in from R/PEcAn
  // TODO: custom likelihood function 
}


