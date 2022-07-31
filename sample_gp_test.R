# sample_gp_test.R
# Runs Stan model 'sample_gp_test.stan', which simply samples from a GP prior. 
# This is just a test for the user-defined Stan functions defined in 
# 'gaussian_process_functions.stan'. Specifically, this script tests 
# cov_exp_quad_same() and cov_exp_quad_cross(), depending on which 
# function is specified in sample_gp_test.stan (either can produce the same 
# output in the case of sampling from the prior as only square matrices are 
# required). However, note that the nugget tests can NOT be conducted with 
# cov_exp_quad_cross(), as this function does not have a nugget argument. 
#
# Andrew Roberts
# Working Directory: /projectnb2/dietzelab/arober/test_code

library(rstan)

#
# Setup
#
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan.model.path <- 'sample_gp_test.stan'
model <- stan_model(stan.model.path)


#
# Test 1: simple 1D example
#

N <- 100
X <- matrix(seq(-4, 4, length.out = N), ncol = 1)

# Settings to pass to Stan
stan_list_1 <- list(N = N, 
                   k = 1, 
                   X = X, 
                   gp_rho = array(1.0, dim = 1), 
                   gp_alpha = 1.0, 
                   gp_sigma = 0.0, 
                   gp_mean = 0.0, 
                   gp_sigma_test = 0.0)

# Compile and fit Stan Model
fit1 <- sampling(model, data = stan_list_1, warmup = 0, iter = 25, chains = 1, 
                 seed = 494838, refresh = 4000, algorithm = "Fixed_param")
samples1 <- extract(fit1)
par(mfrow=c(1,2))
matplot(X, t(samples1$y), type = 'l', col='red', lty=1, xlab = 'x', ylab = 'y')
matplot(X, t(samples1$y_baseline), type = 'l', col='blue', lty=1, xlab = 'x', ylab = 'y_baseline')

#
# Test 2: changing marginal standard deviation
#

N <- 100
X <- matrix(seq(-4, 4, length.out = N), ncol = 1)

# Settings to pass to Stan
stan_list_2 <- list(N = N, 
                    k = 1, 
                    X = X, 
                    gp_rho = array(1.0, dim = 1), 
                    gp_alpha = 5.0, 
                    gp_sigma = 0.0, 
                    gp_mean = 0.0, 
                    gp_sigma_test = 0.0)

# Compile and fit Stan Model
fit2 <- sampling(model, data = stan_list_2, warmup = 0, iter = 25, chains = 1, 
                 seed = 494838, refresh = 4000, algorithm = "Fixed_param")
samples2 <- extract(fit2)
par(mfrow=c(1,2))
matplot(X, t(samples2$y), type = 'l', col='red', lty=1, xlab = 'x', ylab = 'y')
matplot(X, t(samples2$y_baseline), type = 'l', col='blue', lty=1, xlab = 'x', ylab = 'y_baseline')


#
# Test 3: Shorten lengthscale; note that Stan's cov_exp_quad() does not allow unique
#         lengthscale parameters for each dimension, but in this 1D case it doesn't
#         make a difference. 

N <- 100
X <- matrix(seq(-4, 4, length.out = N), ncol = 1)

# Settings to pass to Stan
stan_list_3 <- list(N = N, 
                    k = 1, 
                    X = X, 
                    gp_rho = array(0.1, dim = 1), 
                    gp_alpha = 5.0, 
                    gp_sigma = 0.0, 
                    gp_mean = 0.0, 
                    gp_sigma_test = 0.0)

# Compile and fit Stan Model
fit3 <- sampling(model, data = stan_list_3, warmup = 0, iter = 25, chains = 1, 
                 seed = 494838, refresh = 4000, algorithm = "Fixed_param")
samples3 <- extract(fit3)
par(mfrow=c(1,2))
matplot(X, t(samples3$y), type = 'l', col='red', lty=1, xlab = 'x', ylab = 'y')
matplot(X, t(samples3$y_baseline), type = 'l', col='blue', lty=1, xlab = 'x', ylab = 'y_baseline')



# Test 4: Decrease marginal standard deviation and Increase length scale 

N <- 100
X <- matrix(seq(-4, 4, length.out = N), ncol = 1)

# Settings to pass to Stan
stan_list_4 <- list(N = N, 
                    k = 1, 
                    X = X, 
                    gp_rho = array(100.0, dim = 1), 
                    gp_alpha = 0.5, 
                    gp_sigma = 0.0, 
                    gp_mean = 0.0, 
                    gp_sigma_test = 0.0)

# Compile and fit Stan Model
fit4 <- sampling(model, data = stan_list_4, warmup = 0, iter = 25, chains = 1, 
                 seed = 494838, refresh = 4000, algorithm = "Fixed_param")
samples4 <- extract(fit4)
par(mfrow=c(1,2))
matplot(X, t(samples4$y), type = 'l', col='red', lty=1, xlab = 'x', ylab = 'y')
matplot(X, t(samples4$y_baseline), type = 'l', col='blue', lty=1, xlab = 'x', ylab = 'y_baseline')


#
# Test 5: Add non-zero nugget
#

N <- 100
X <- matrix(seq(-4, 4, length.out = N), ncol = 1)
nugget <- 1.0

# Settings to pass to Stan
stan_list_5 <- list(N = N, 
                    k = 1, 
                    X = X, 
                    gp_rho = array(1.0, dim = 1), 
                    gp_alpha = 1.0, 
                    gp_sigma = nugget, 
                    gp_mean = 0.0, 
                    gp_sigma_test = nugget)

# Compile and fit Stan Model
fit5 <- sampling(model, data = stan_list_5, warmup = 0, iter = 25, chains = 1, 
                 seed = 494838, refresh = 4000, algorithm = "Fixed_param")
samples5 <- extract(fit5)
par(mfrow=c(1,2))
matplot(X, t(samples5$y), type = 'l', col='red', lty=1, xlab = 'x', ylab = 'y')
matplot(X, t(samples5$y_baseline), type = 'l', col='blue', lty=1, xlab = 'x', ylab = 'y_baseline')


#
# Test 6: Basic 2D example
#

nx <- 20
x <- seq(0, 2, length = nx)
X <- as.matrix(expand.grid(x, x))

# Settings to pass to Stan
stan_list_6 <- list(N = nrow(X), 
                    k = ncol(X), 
                    X = X, 
                    gp_rho = c(1.0, 1.0), 
                    gp_alpha = 1.0, 
                    gp_sigma = 0.0, 
                    gp_mean = 0.0, 
                    gp_sigma_test = 0.0)

# Compile and fit Stan Model
fit6 <- sampling(model, data = stan_list_6, warmup = 0, iter = 1000, chains = 1, 
                 seed = 494838, refresh = 4000, algorithm = "Fixed_param")
samples6 <- extract(fit6)

j <- 10 # Change this counter to view different samples one at a time
par(mfrow=c(1,2))
persp(x, x, matrix(samples6$y[j, ], ncol=nx), theta=-30, phi=30, xlab="x1", 
      ylab="x2", zlab="y")
persp(x, x, matrix(samples6$y_baseline[j, ], ncol=nx), theta=-30, phi=30, xlab="x1", 
      ylab="x2", zlab="y")

# Summary stats
print(paste0("Min: y = ", min(samples6$y), "; y_baseline = ", min(samples6$y_baseline)))
print(paste0("Max: y = ", max(samples6$y), "; y_baseline = ", max(samples6$y_baseline)))
print(paste0("Mean: y = ", mean(samples6$y), "; y_baseline = ", mean(samples6$y_baseline)))


#
# Test 7: Scaling different dimensions by different lengthscales; note: can't 
#         test against cov_exp_quad() for this one as it doesn't support this
#         feature. 
#

nx <- 20
x <- seq(0, 2, length = nx)
X.df <- expand.grid(x, x)
X <- as.matrix(X.df)

# Settings to pass to Stan
stan_list_7 <- list(N = nrow(X), 
                    k = ncol(X), 
                    X = X, 
                    gp_rho = c(0.1, 10.0), 
                    gp_alpha = 1.0, 
                    gp_sigma = 0.0, 
                    gp_mean = 0.0, 
                    gp_sigma_test = 0.0)

# Compile and fit Stan Model
fit7 <- sampling(model, data = stan_list_7, warmup = 0, iter = 100, chains = 1, 
                 seed = 494838, refresh = 4000, algorithm = "Fixed_param")
samples7 <- extract(fit7)

j <- 5 # Change this counter to view different samples one at a time
persp(x, x, matrix(samples7$y[j, ], ncol=nx), theta=-30, phi=30, xlab="x1", 
      ylab="x2", zlab="y")

# Plotting projections onto each dimension
selector.x1.0 <- X.df$Var1 == 0
selector.x2.0 <- X.df$Var2 == 0
par(mfrow=c(1,2))
matplot(x, t(samples7$y)[selector.x1.0,], type = 'l', col='red', lty=1, xlab = 'x1 = 0, x2')
matplot(x, t(samples7$y)[selector.x2.0,], type = 'l', col='blue', lty=1, xlab = 'x1, x2 = 0')

# Comparing against cov_exp_quad() in the first direction (since cov_exp_quad() used lengthscale
# 0.1 in both dimensions)
par(mfrow=c(1,2))
matplot(x, t(samples7$y)[selector.x2.0,], type = 'l', col='red', lty=1, 
        xlab = 'x1 = 0, x2', ylab = 'y')
matplot(x, t(samples7$y_baseline)[selector.x2.0,], type = 'l', col='blue', lty=1, 
        xlab = 'x1, x2 = 0', ylab = 'y_baseline')


#
# Test 8: 1D test, varying GP mean argument. 
#

N <- 100
X <- matrix(seq(-4, 4, length.out = N), ncol = 1)
gp.mean <- 100.0

# Settings to pass to Stan
stan_list_8 <- list(N = N, 
                    k = 1, 
                    X = X, 
                    gp_rho = array(1.0, dim = 1), 
                    gp_alpha = 1.0, 
                    gp_sigma = 0.0, 
                    gp_mean = gp.mean, 
                    gp_sigma_test = 0.0)

# Compile and fit Stan Model
fit8 <- sampling(model, data = stan_list_8, warmup = 0, iter = 100, chains = 1, 
                 seed = 494838, refresh = 4000, algorithm = "Fixed_param")
samples8 <- extract(fit8)
par(mfrow=c(1,2))
matplot(X, t(samples8$y), type = 'l', col='red', lty=1, xlab = 'x', ylab = 'y')
matplot(X, t(samples8$y_baseline), type = 'l', col='blue', lty=1, xlab = 'x', ylab = 'y_baseline')





