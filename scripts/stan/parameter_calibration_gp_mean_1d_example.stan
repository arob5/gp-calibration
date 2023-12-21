/**
* parameter_calibration_gp_mean_1d_example.stan
* Run MCMC for parameter calibration on toy example with simple scalar function. 
* Uses GP mean approximation of sufficient statistic in Gaussian likelihood. 
*
* Andrew Roberts
*
*/

functions {
#include gaussian_process_functions.stan
}

data {
  int<lower=1> N; // Number of design points (i.e. knots)
  int<lower=1> n; // The number of observed data points (i.e. the dimension of the MVN likelihood)
  
  vector[N] y;    // Vector of sufficient statistics evaluated at design points
  matrix[N, 1] X; // Design matrix

  real<lower=0> a; // Shape parameter for Gamma hyperprior on tau
  real<lower=0> b; // Rate parameter for Gamma hyperprior on tau
  
  real u_mean;            // Mean for Gaussian prior on calibration parameter u
  real<lower=0> u_sigma; // Standard deviation for Gaussian prior on calibration parameter u
  
  real u_lower; // Lower bound on u
  real u_upper; // Upper bound on u
  
  vector<lower=0>[1] gp_rho; // Vector of lengthscale parameters
  real<lower=0> gp_alpha;    // Marginal standard deviation
  real<lower=0> gp_sigma;    // Nugget (standard deviation)
  real gp_mean;              // Constant GP mean
}

transformed data {
  matrix[N, N] K = cov_exp_quad_same(X, gp_rho, gp_alpha, gp_sigma); 
  matrix[N, N] L = cholesky_decompose(K);
  vector[N] K_inv_y = L' \ (L \ (y - rep_vector(gp_mean, N))); 
}

parameters {
  vector<lower=u_lower, upper=u_upper>[1] u; // Calibration parameters
  real<lower=0> tau; // Precision parameter for Gaussian likelihood
}

model {
  u ~ normal(u_mean, u_sigma); 
  tau ~ gamma(a, b);
  
  // Increment log target density
  // target += 0.5 * n * log(tau) - 0.5 * tau * calc_gp_predictive_mean(to_vector(cov_exp_quad_cross(to_matrix(u, 1, 1, 1), X, gp_rho, gp_alpha)), K_inv_y, gp_mean);  
                                                                     
  
  target += gp_mean_gaussian_llik(u, X, K_inv_y, N, n, 1, gp_rho, gp_alpha, gp_mean, tau); 
}








