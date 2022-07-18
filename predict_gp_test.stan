/**
* predict_gp_test.stan
* Testing user-defined functions in 'gaussian_process_functions.stan'. This script
* tests the functions 'calc_gp_predictive_mean()' and 'calc_gp_predictive_var()'
* by computing pointwise predictions and GP variances at a set of test points.
* This Stan program does not have parameters or model blocks and doesn't actually
* do any sampling or optimization; it's simply a test to retrieve the values
* from the user-defined functions of interest. 
*
* Andrew Roberts
*
*/

functions{
#include gaussian_process_functions.stan
}

data {
  int<lower=1> k; // Number of calibration parameters
  int<lower=1> N; // Number of design points (i.e. knots)
  int<lower=1> m; // Number of test points
  matrix[N, k] X; // Design matrix
  vector[N] y;    // Process model evaluations at design points
  matrix[m, k] X_pred; // Test points
  
  vector<lower=0>[k] gp_rho; // Vector of lengthscale parameters
  real<lower=0> gp_alpha; // Marginal standard deviation
  real<lower=0> gp_sigma; // Nugget (standard deviation)
}

transformed data {
  matrix[N, N] K = cov_exp_quad_same(X, gp_rho, gp_alpha, gp_sigma); 
  matrix[N, N] L = cholesky_decompose(K);
  vector[N] K_inv_y = L' \ (L \ y); 
}

generated quantities {
  vector[m] mean_pred; 
  vector[m] var_pred; 
  matrix[N, N] K_out;
  matrix[N, m] K_cross_out; 
  
  for(i in 1:m) {
    matrix[1, k] x_mat = to_matrix(X_pred[i], 1, k, 1); 
    vector[N] k_xX = to_vector(cov_exp_quad_cross(x_mat, X, gp_rho, gp_alpha));
    mean_pred[i] = calc_gp_predictive_mean(k_xX, K_inv_y); 
    var_pred[i] = calc_gp_predictive_var(L, k_xX, x_mat, N, gp_rho, gp_alpha, gp_sigma); 
  }
  
  // Output covariances to test the kernel functions.
  K_out = K; 
  K_cross_out = cov_exp_quad_cross(X, X_pred, gp_rho, gp_alpha);
  
}





