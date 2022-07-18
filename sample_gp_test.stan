/**
* sample_gp_test.stan
* Testing user-defined functions in 'gaussian_process_functions.stan'. This script
* simply samples from a GP. 
*
* Andrew Roberts
*
*/

functions {
#include gaussian_process_functions.stan
}

data {
  int<lower=1> N; // Number of locations at which to sample
  int<lower=1> k; // Dimension of underlying space
  matrix[N, k] X; // Design matrix
  
  vector<lower=0>[k] gp_rho; // Vector of lengthscale parameters
  real<lower=0> gp_alpha; // Marginal standard deviation
  real<lower=0> gp_sigma; // Nugget (standard deviation)
  real gp_mean; // Vector of regression coefficients for linear mean function
  
  real<lower=0> gp_sigma_test; // Passed to cov_exp_quad() to test against user-defined functions.
}

transformed data {
  vector[N] mu; 
  matrix[N, N] K; 
  matrix[N, N] L; 
  matrix[N, N] K_baseline; // For testing 
  matrix[N, N] L_baseline; 
  row_vector[k] X_rowvec[N]; // cov_exp_quad() does not accept matrix as argument
  for(i in 1:N) {
    X_rowvec[i] = X[i]; 
  }
  
  mu = rep_vector(gp_mean, N); 
  K = cov_exp_quad_same(X, gp_rho, gp_alpha, gp_sigma); 
  // Can uncomment this and comment above line in order to test cov_exp_quad_cross()
  // K = cov_exp_quad_cross(X, X, gp_rho, gp_alpha) + diag_matrix(rep_vector(1e-10, N)); 
  K_baseline = cov_exp_quad(X_rowvec, gp_alpha, sqrt(gp_rho[1])); 
  
  L = cholesky_decompose(K); 
  L_baseline = cholesky_decompose(K_baseline + diag_matrix(rep_vector(1e-10 + gp_sigma_test^2, N)));
}

parameters {}

model {}

generated quantities {
  vector[N] y = multi_normal_cholesky_rng(mu, L);
  vector[N] y_baseline = multi_normal_cholesky_rng(mu, L_baseline);
}




