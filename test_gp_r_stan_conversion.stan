/**
* test_gp_r_stan_conversion.stan
* Utilized in the R script 'test_gp_r_stan_conversion.R' to test the user-defined
* Stan kernel functions. 
* 
* Andrew Roberts
*
*/

functions {
#include gaussian_process_functions.stan
}

data {
  int<lower=1> N; // Number of design points
  int<lower=1> N_pred; // The number of prediction points
  int<lower=1> k;      // Dimension of input space
  
  matrix[N, k] X; // Design matrix
  vector[N] y;    // Vector of model evaluations at the design points
  matrix[N_pred, k] X_pred; // Matrix of prediction points
  
  vector<lower=0>[k] gp_rho; // Vector of lengthscale parameters
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
  matrix[N, N] K_out;
  matrix[N_pred, N_pred] K_pred; 
  matrix[N, N_pred] K_cross; 
  vector[N_pred] mean_pred; 
  vector[N_pred] var_pred; 
  vector[N] mean_design; 
  vector[N] var_design; 
  
  K_out = K; 
  K_pred = cov_exp_quad_same(X_pred, gp_rho, gp_alpha, gp_sigma);
  K_cross = cov_exp_quad_cross(X, X_pred, gp_rho, gp_alpha); 
  
  for(i in 1:N_pred) {
    matrix[1, k] x_mat = to_matrix(X_pred[i], 1, k, 1); 
    vector[N] k_xX = to_vector(cov_exp_quad_cross(x_mat, X, gp_rho, gp_alpha));
    mean_pred[i] = calc_gp_predictive_mean(k_xX, K_inv_y, gp_mean); 
    var_pred[i] = calc_gp_predictive_var(L, k_xX, x_mat, N, gp_rho, gp_alpha, gp_sigma);
  }
  
  for(i in 1:N) {
    matrix[1, k] x_mat = to_matrix(X[i], 1, k, 1); 
    vector[N] k_xX = to_vector(cov_exp_quad_cross(x_mat, X, gp_rho, gp_alpha));
    mean_design[i] = calc_gp_predictive_mean(k_xX, K_inv_y, gp_mean); 
    var_design[i] = calc_gp_predictive_var(L, k_xX, x_mat, N, gp_rho, gp_alpha, gp_sigma);
  }
  
}
