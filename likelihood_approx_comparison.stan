/**
* likelihood_approx_comparison.stan
* Utilized in the R script of the same name to investigate the effect of the 
* GP approximation on the likelihood function in a simple 1D toy example.
*
* Andrew Roberts
*
*/

functions {
#include gaussian_process_functions.stan
}

data {
  int<lower=1> N; // Number of locations at which to sample
  int<lower=1> n; // The number of observed data points (i.e. the dimension of the MVN likelihood)
  int<lower=1> m; // The number of prediction points (i.e. points at which to evaluate gp_approx())
  
  real<lower=0> tau; // Precision parameter to use for Gaussian likelihood
  matrix[N, 1] X; // Design matrix
  vector[N] y;    // Vector of sufficient statistics evaluated at design points
  vector[1] u_vals[m]; // Vector of 1D prediction points
  
  vector<lower=0>[1] gp_rho; // Vector of lengthscale parameters
  real<lower=0> gp_alpha; // Marginal standard deviation
  real<lower=0> gp_sigma; // Nugget (standard deviation)
  real gp_mean; // Vector of regression coefficients for linear mean function
}

transformed data {
  matrix[N, N] K; 
  matrix[N, N] L; 
  vector[N] K_inv_y; 
  
  K = cov_exp_quad_same(X, gp_rho, gp_alpha, gp_sigma);
  L = cholesky_decompose(K); 
  K_inv_y = L' \ (L \ (y - rep_vector(gp_mean, N)));
}

generated quantities {
  vector[m] u_vals_llik; 
  
  for(i in 1:m) {
    u_vals_llik[i] = gp_approx(L, K_inv_y, X, N, n, 1, gp_rho, gp_alpha, gp_sigma, gp_mean, tau, u_vals[i]); 
  }
  
}


