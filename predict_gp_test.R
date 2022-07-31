# predict_gp_test.R
# Runs Stan model 'predict_gp_test.stan', which simply calculates GP 
# predictive means and variances.
# This is just a test for the user-defined Stan functions defined in 
# 'gaussian_process_functions.stan'.
#
# Andrew Roberts
# Working Directory: /projectnb2/dietzelab/arober/test_code

library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#
# Define GP functions to test Stan output against
# 

# Squared Exponential (i.e. Gaussian) kernel
K <- function(X1, X2 = NA, rho, alpha) {
  if(is.na(X2[1][1])) X2 <- X1
  C <- matrix(NA, nrow = nrow(X1), ncol = nrow(X2))
  
  for(i in 1:nrow(X1)) {
    for(j in 1:nrow(X2)) {
      C[i, j] <- sum((1 / rho) * (X1[i, ] - X2[j, ])^2)
    }
  }
  
  return(alpha^2 * exp(-0.5 * C))
}

predict_mean <- function(X_pred, X_obs, y_obs, rho, alpha, sigma) {
  eps <- sqrt(.Machine$double.eps) # "Jitter" to ensure positive definite matrices
  cross_cov <- K(X_pred, X_obs, rho, alpha)
  data_cov <- K(X_obs, X_obs, rho, alpha) + diag(rep(eps + sigma^2, nrow(X_obs)))
  return(cross_cov %*% solve(data_cov) %*% y_obs)
}

predict_var<- function(X_pred, X_obs, y_obs, rho, alpha, sigma) {
  eps <- sqrt(.Machine$double.eps) 
  prior_var <- K(X_pred, X_pred, rho, alpha) + diag(rep(eps + sigma^2, nrow(X_pred)))
  cross_cov <- K(X_pred, X_obs, rho, alpha)
  data_cov <- K(X_obs, X_obs, rho, alpha) + diag(rep(eps + sigma^2, nrow(X_obs)))
  return(prior_var - cross_cov %*% solve(data_cov) %*% t(cross_cov))
}


#
# Compile Stan Code
# 

stan.model.path <- 'predict_gp_test.stan'
model <- stan_model(stan.model.path)


#
# Simple 1D Test
#

# Kernel Parameters
rho <- array(1.0, dim = 1)
alpha <- 1.0
sigma <- 0.0

# Observed data
N <- 8
X <- matrix(seq(0, 2*pi, length=N), ncol=1)
y <- sin(X)

# Test data
m <- 100
X_test <- matrix(seq(-0.5, 2*pi + 0.5, length=m), ncol=1)

# Predict values using R code to be compared to Stan output
means <- as.vector(predict_mean(X_test, X, y, rho, alpha, sigma))
cov_pred <- predict_var(X_test, X, y, rho, alpha, sigma)
vars <- diag(cov_pred)

# Stan Model 
stan.list <- list(k = 1, 
                  N = N, 
                  m = m, 
                  X = X, 
                  y = as.vector(y), 
                  X_pred = X_test, 
                  gp_rho = rho, 
                  gp_alpha = alpha, 
                  gp_sigma = sigma)

stan.results <- sampling(model, data = stan.list, iter = 1, chains = 1, seed=596858228, algorithm="Fixed_param")
stan.output <- extract(stan.results)
means.stan <- as.vector(stan.output$mean_pred)
vars.stan <- as.vector(stan.output$var_pred)
K.stan <- stan.output$K_out[1,,]
K.cross.stan <- stan.output$K_cross_out[1,,]

# Test covariance functions
K.r <- K(X, X, rho, alpha) + diag(rep(sqrt(.Machine$double.eps), nrow(X)))
K.cross.r <- K(X, X_test, rho, alpha)

print(paste0("K matrix from R and Stan are equal: ", all.equal(K.r, K.stan)))
print(paste0("K cross matrix from R and Stan are equal: ", all.equal(K.cross.r, K.cross.stan)))

# Test mean and variance predictions
print(paste0("Predictive means from R and Stan are equal: ", all.equal(means, means.stan)))
print(paste0("Predictive variances from R and Stan are equal: ", all.equal(vars, vars.stan)))

# Plots of predicted means
par(mfrow=c(1,2))
matplot(X_test, means, type = 'l', col='gray', lty=1)
matpoints(X, y, pch = 20, cex = 2)

matplot(X_test, means.stan, type = 'l', col='gray', lty=1)
matpoints(X, y, pch = 20, cex = 2)


#
# 2D Test with different kernel parameters and nugget
#

# Kernel Parameters
rho <- c(2.0, 0.5)
alpha <- 10.0
sigma <- 1.0

# Observed data
nx <- 20
x <- seq(0, 2*pi, length = nx)
X_df <- expand.grid(x, x)
X <- as.matrix(X.df)
y <- sin(X[,1] + 4*X[,2])

# Test data
m <- 30
x_test <- seq(-0.5, 2*pi + 0.5, length=m)
X_test_df <- expand.grid(x_test, x_test)
X_test <- as.matrix(X_test_df)

# Predict values using R code to be compared to Stan output
means <- as.vector(predict_mean(X_test, X, y, rho, alpha, sigma))
cov_pred <- predict_var(X_test, X, y, rho, alpha, sigma)
vars <- diag(cov_pred)

# Stan Model 
stan.list <- list(k = 2, 
                  N = nrow(X), 
                  m = nrow(X_test), 
                  X = X, 
                  y = y, 
                  X_pred = X_test, 
                  gp_rho = rho, 
                  gp_alpha = alpha, 
                  gp_sigma = sigma)

stan.results <- sampling(model, data = stan.list, iter = 1, chains = 1, seed=596858228, algorithm="Fixed_param")
stan.output <- extract(stan.results)
means.stan <- as.vector(stan.output$mean_pred)
vars.stan <- as.vector(stan.output$var_pred)
K.stan <- stan.output$K_out[1,,]
K.cross.stan <- stan.output$K_cross_out[1,,]

# Test covariance functions
K.r <- K(X, X, rho, alpha) + diag(rep(sqrt(.Machine$double.eps) + sigma^2, nrow(X)))
K.cross.r <- K(X, X_test, rho, alpha)

print(paste0("K matrix from R and Stan are equal: ", all.equal(K.r, K.stan)))
print(paste0("K cross matrix from R and Stan are equal: ", all.equal(K.cross.r, K.cross.stan)))

# Test mean and variance predictions
print(paste0("Predictive means from R and Stan are equal: ", all.equal(means, means.stan)))
print(paste0("Predictive variances from R and Stan are equal: ", all.equal(vars, vars.stan)))


