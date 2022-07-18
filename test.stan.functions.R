# test.stan.functions.R
# Tests the user-defined Stan functions in 'gussian_process_functions.stan'. 
#
# Andrew Roberts
# Working Directory: /projectnb2/dietzelab/arober/test_code

library(rstan)

func_squared_distance_weighted <- 
'
functions {
  real squared_distance_weighted(row_vector x, row_vector y, vector w) {
  real dist = 0.0; 
  int N = num_elements(x);
  
  for(i in 1:N) {
    dist += w[i] * (x[i] - y[i])^2; 
  }
  
  return(dist); 
}
}
'

func_cov_exp_quad_gen_1 <- 
'
functions {
  matrix cov_exp_quad_gen(matrix X, vector rho, real alpha, real sigma) {
    int N = rows(X); 
    int k = num_elements(rho); 
    vector[k] weights = rho .^ -2.0; 
      
    cov_matrix[N] K; 
      
    for (i in 1:(N - 1)) {
      K[i, i] = alpha^2 + sigma^2;
      for (j in (i + 1):N) {
        K[i, j] = alpha^2 * exp(-0.5 * squared_distance_weighted(X[i], X[j], weights));
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = alpha^2 + sigma^2;
    
    return(K); 
  }
}
'




expose_stan_functions(stanc(model_code = func_squared_distance_weighted))
expose_stan_functions(stanc(model_code = func_cov_exp_quad_gen_1))



#
# Test squared_distance_weighted
# 

# Test 1
x <- c(1, 2, 3)
y <- x
all.equal(squared_distance_weighted(x, y, rep(1.0, 3)), sum((x - y)^2))

# Test 2
x <- rnorm(100, 10, 1)
y <- rnorm(100)
all.equal(squared_distance_weighted(x, y, rep(1.0, 100)), sum((x - y)^2))

# Test 3
x <- rnorm(100, 10, 1)
y <- rnorm(100)
w <- runif(100)
all.equal(squared_distance_weighted(x, y, w), sum(w * (x - y)^2))



