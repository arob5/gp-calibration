/**
 * A generalization of the built-in Stan function 'squared_distance()'.
 * Generalizes this function by allowing for scaling each dimension
 * by provided weights. This is only currently defined for row_vectors
 * given its use in 'cov_exp_quad_gen()'. 
 *
 * @param x First row_vector. 
 * @param y Second row_vector. Must be of same length as x. 
 * @param w Vector of weights. Must be of same length as x. 
 * 
 * @return The squared, weighted Euclidean distance between x and y. 
 *         The ith dimension is weighted as w[i] * (x[i] - y[i])^2
 */
real squared_distance_weighted(row_vector x, row_vector y, vector w) {
  real dist = 0.0; 
  int N = num_elements(x);
  
  for(i in 1:N) {
    dist += w[i] * (x[i] - y[i])^2; 
  }
  
  return(dist); 
}




/**
 * A generalization of the built-in Stan function 'cov_exp_quad()'.
 * Generalizes this function by allowing for separate lengthscale
 * parameters for each dimension and a (constant) nugget term. 
 * Defines the squared exponential kernel/covariance  function
 * and returns a square matrix formed by applying the covariance  
 * function to pairs of input vectors. To apply the covariance 
 * function to two separate matrices (e.g. train and test data)
 * see the function cov_exp_quad_cross(). 
 *
 * @param X Matrix of input vectors. Each input is a row of the matrix.
 * @param rho Vector of length scale parameters. Must have length equal to 
 *            number of columns of X. Dimension k is scaled by 1/rho[k]^2. 
 * @param alpha Marginal standard deviation of the GP. 
 * @param sigma Nugget parameter. sigma^2 is added to each element on the 
 *              diagonal of the resulting covariance matrix. 
 * 
 * @return The matrix formed by applying to covariance function 
 *         to each pair of inputs. Has dimension (rows(X), rows(X))
 *         and the (i, j) element of the returned matrix is the 
 *         covariance between X[i] and X[j]. 
 */
matrix cov_exp_quad_same(matrix X, vector rho, real alpha, real sigma) {
  int N = rows(X); 
  int k = num_elements(rho); 
  vector[k] weights = 1.0 ./ square(rho);
  
  real eps;
  matrix[N, N] K;
    
  if(sigma == 0.0) {
    eps = sqrt(machine_precision()); 
  } else{
    eps = sigma^2; 
  }
    
  for (i in 1:(N - 1)) {
    K[i, i] = alpha^2 + eps;
    for (j in (i + 1):N) {
      K[i, j] = alpha^2 * exp(-0.5 * squared_distance_weighted(X[i], X[j], weights));
      K[j, i] = K[i, j];
    }
  }
  K[N, N] = alpha^2 + eps;
  
  return(K); 
}



/**
 * A generalization of the built-in Stan function 'cov_exp_quad()'.
 * Generalizes this function by allowing for separate lengthscale
 * parameters for each dimension. 
 * Defines the squared exponential kernel/covariance  function
 * and returns a potentially non-square matrix formed by applying the covariance  
 * function to pairs of input vectors. To apply the covariance 
 * function to a sigle matrix see cov_exp_quad_same(). This 
 * function instead calculates the covariances between two separate matrices 
 * (hence, "cross"). This version of the function does not require the nugget 
 * parameter as this function implicitly assumes that X1 and X2 are "different"
 * data, in which case the nugget does not apply. If producing a true covariance 
 * matrix with diagonals corresponding variances, then 'cov_exp_quad_same()'
 * should be used instead. 
 *
 * @param X1 Matrix of input vectors. Each input is a row of the matrix.
 * @param X2 Matrix of input vectors. Each input is a row of the matrix. Must 
 *           have same number of columns as X1, but number of rows can differ. 
 * @param rho Vector of length scale parameters. Must have length equal to 
 *            number of columns of X1. Dimension k is scaled by 1/rho[k]^2.
 *            Note that this is different from 'cov_exp_quad()', which scales by 
 *            1/rho[k]^2. This difference is due to the fact that rstan currently
 *            doesn't support the elementwise power function. 
 * @param alpha Marginal standard deviation of the GP. 
 * 
 * @return The matrix formed by applying to covariance function 
 *         to each pair of inputs. Has dimension (rows(X1), rows(X2))
 *         and the (i, j) element of the returned matrix is the 
 *         covariance between X1[i] and X2[j]. 
 */
matrix cov_exp_quad_cross(matrix X1, matrix X2, vector rho, real alpha) {
  int N1 = rows(X1); 
  int N2 = rows(X2); 
  int k = num_elements(rho); 
  vector[k] weights = 1.0 ./ square(rho);
  
  matrix[N1, N2] K; 
  
  for (i in 1:N1) {
    for (j in 1:N2) {
      K[i, j] = alpha^2 * exp(-0.5 * squared_distance_weighted(X1[i], X2[j], weights));
    }
  }

  return(K); 
}

real calc_gp_predictive_mean(vector k_xX, vector K_inv_y, real mu) {
  
  return(mu + dot_product(k_xX, K_inv_y));
  
}

real calc_gp_predictive_var(matrix L, vector k_xX, matrix x_mat, int N, vector rho, real alpha, real sigma) {
  
  vector[N] v = L \ k_xX;
  real k_x = cov_exp_quad_same(x_mat, rho, alpha, sigma)[1, 1] - dot_self(v); 
  
  return(k_x);
  
}

/**
* Evaluates the approximate log-likelihood at a given parameter value. This is an approximation
* to the Gaussian likelihood, where the sum of squares term (the sufficient statistic)
* has been replaced by a fit Gaussian Process (GP) emulator, and then the GP is integrated
* out in order to account for the approximation uncertainty. Additive constants in the 
* log-likelihood are dropped.
*-
* @param L Cholesky factor of the covariance matrix of the design points. 
* @param K_inv_y The vector inv( K(X, X) ) * y, where X is the design matrix and y the 
*        observed response vector. Using the Cholesky factor, this can be calculated as 
*        L' \ (L \ y). 
* @param N The number of design points.
* @param n The number of observed data points (i.e. the dimension of the MVN likelihood).  
* @param k The dimension of the parameter/input space; i.e. the number of calibration parameters.
* @param rho Vector of lengthscale parameters for the kernel. 
* @param alpha Marginal standard deviation parameter for the kernel. 
* @param sigma Nugget (standard deviation) parameter for the kernel. 
* @param mu The GP (constant) mean. 
* @param tau The precision parameter in the Gaussian likelihood. 
* @param x Vector, the parameter value at which to evaluate the log-likelihood. 
*
* @return The Gaussian log-likelihood (up to an additive constant), calculated using 
*         the GP approximation to the sufficient statistic, and marginalized over the 
*         GP. 
*/
real gp_approx(matrix L, vector K_inv_y, matrix X, int N, int n, int k, vector rho, real alpha, real sigma, real mu, real tau, vector x) {
  
  // Temporary: should remove this once I upgrade to cmdstanr and can use argument overloading
  matrix[1, k] x_mat = to_matrix(x, 1, k, 1); 
  
  // Predictive mean
  vector[N] k_xX = to_vector(cov_exp_quad_cross(x_mat, X, rho, alpha));
  real mu_x = mu + dot_product(k_xX, K_inv_y);
  
  // Predictive variance
  vector[N] v = L \ k_xX;
  real k_x = cov_exp_quad_same(x_mat, rho, alpha, sigma)[1, 1] - dot_self(v); 
  
  
  return(0.5 * n * log(tau) - 0.5 * tau * mu_x + 0.125 * square(tau) * k_x); 
  
}


real gp_mean_gaussian_llik(vector x, matrix X, vector K_inv_y, int N, int n, int k, vector rho, real alpha, real mu, real tau) {
  // Temporary: should remove this once I upgrade to cmdstanr and can use argument overloading
  matrix[1, k] x_mat = to_matrix(x, 1, k, 1);
  
  // Predictive mean
  vector[N] k_xX = to_vector(cov_exp_quad_cross(x_mat, X, rho, alpha));
  real gp_pred_mean = calc_gp_predictive_mean(k_xX, K_inv_y, mu); 
  
  return(0.5 * n * log(tau) - 0.5 * tau * gp_pred_mean);
}

vector seq_fun(real start, real end, int N_by) {
 
  real h;
  vector[N_by] out;
  h = (end-start) / (N_by-1);
  for (i in 1:N_by) {
    out[i] = start + (i-1)*h; 
  }
  
  return(out);
  
 }

/**
 * A numerical approximation of the Laplace transform of the log-normal 
 * density function L(s), which is related to the log-normal Moment 
 * Generating Function (MGF) M via L(s) = M(-s). This approximation first cuts
 * of the integrand (which is supported on the entire real line) at calculated
 * cut points, and then applies the trapezoidal rule on the now finite support.
 *
 * @param s Value at which to evaluate the Laplace transform. 
 */
real lognormal_mgf_numerical_approx(real s, real mu, real sigma, int num_eval, real tol, real M) {
  real eps = 0.5 * tol;   
  real cut_lower = mu - sigma * sqrt2() * sqrt(log(M/eps));
  real cut_upper = log(1/s) + log(log(M/eps));
  real dx = (cut_upper - cut_lower) / num_eval; 
  vector[num_eval] x = seq_fun(cut_lower, cut_upper, num_eval); 
  vector[num_eval] fx = exp(-s * exp(x) - 0.5 * square(x - mu) / square(sigma)); 
  // vector[num_eval] integrand_x = exp(-1.0 * (s * exp(x) + 0.5 * (1 / square(sigma)) * square(x - mu))); 
  // return(dx / (sigma * sqrt2() * sqrt(pi())) * (0.5 * (fx[1] + fx[num_eval]) + sum(fx[2:(num_eval - 1)]))); 
  
  return(0.5 * dx * (1 / sqrt(2.0 * pi())) * (1 / sigma) * (fx[1] + fx[num_eval] + 2.0*sum(fx[2:(num_eval-1)]))); 
  
}


real gp_log_approx(matrix L, vector K_inv_y, matrix X, int N, int n, int k, vector rho, 
                   real alpha, real sigma, real mu, real tau, vector x, int num_eval, real tol, real M) {
  
  // Temporary: should remove this once I upgrade to cmdstanr and can use argument overloading
  matrix[1, k] x_mat = to_matrix(x, 1, k, 1); 
  
  // Predictive mean
  vector[N] k_xX = to_vector(cov_exp_quad_cross(x_mat, X, rho, alpha));
  real mu_x = mu + dot_product(k_xX, K_inv_y);
  
  // Predictive variance
  vector[N] v = L \ k_xX;
  real k_x = cov_exp_quad_same(x_mat, rho, alpha, sigma)[1, 1] - dot_self(v); 
  
  // Approximate expectation over GP
  real gp_exp = lognormal_mgf_numerical_approx(0.5 * tau, mu_x, sqrt(k_x), num_eval, tol, M); 
  
  return(0.5 * n * log(tau) + log(gp_exp)); 
  
}


