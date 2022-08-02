# parameter_calibration_1d_example.R
# Run MCMC for parameter calibration on toy example with simple scalar function.
# Note that the Stan function f() defined in parameter_calibration_1d_example.stan
# should be equivalent to the function f() defined in this file. This file 
# simulates a dataset and then performs parameter calibration using the 
# sample likelihood and priors used to simulate the dataset. It therefore 
# assumes a perfect model. This file reads 'parameter_calibration_1d_example.stan'
# and 'parameter_calibration_gp_1d_example.stan' to conduct brute force and GP
# approx MCMC, respectively, 
#
# Andrew Roberts
# Working Directory: /projectnb2/dietzelab/arober/pecan_personal

library(rstan)
library(mlegp)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("stan.helper.functions.R")

#
# Settings
#

# The function that is acting as the computer model. Should be equivalent to 
# the function of the same name in the Stan file. 
f <- function(u) {
  return(u)
}

# Hyperparameters (shape and rate) for Gamma hyperprior on tau (precision parameter). Using values 
# taken from an actual PEcAn run. 
tau.shape <- 16
tau.rate <- 1.0

# Number of observations to simulate
n <- 1000

# Gaussian prior parameters on calibration parameter
u.mean <- 0.5
u.sigma <- 0.25

# Directories
base.dir <- '.'
out.dir <- file.path(base.dir, 'output')

#
# Simulate Data, fixing u and tau at their prior means
#
u <- u.mean
tau <- tau.shape / tau.rate
y <- rnorm(n, f(u), 1/sqrt(tau))


#
# Brute Force parameter calibration: try to recover parameters by running MCMC in Stan
#

# Compile Stan code
stan.model.path <- file.path(base.dir, 'parameter_calibration_1d_example.stan')
model <- stan_model(stan.model.path)

# Data to pass to Stan
stan.list <- list(n = n, 
                  y = y, 
                  a = tau.shape, 
                  b = tau.rate, 
                  u_mean = u.mean, 
                  u_sigma = u.sigma)

# MCMC
fit <- sampling(model, stan.list, iter = 50000, chains = 4)
summary(fit)
samples.brute.force <- extract(fit)

#
# Fit GP regression
#

# Design matrix, evaluate model at design points
# TODO: Try using more extreme u values
N <- 10
X <- matrix(seq(qnorm(.01, u.mean, u.sigma), qnorm(.99, u.mean, u.sigma), length = N), ncol = 1)
y.model <- f(X)

# Define the sufficient statistic (SS). For simplicity here, I'm not considering any other explanatory variables
# x that are being conditioned on. If there were, we would have to make sure we were lining up the observed data
# and the process model data correctly conditional on x. So in this case, y_model[1] is just a constant predictor
# of the output variable given the first calibration parameter (design point) X[1,1]. 
SS <- rep(0, N)
for(i in seq(1, N)) {
  SS[i] <- sum((y.model[i] - y)^2)
}

# Scale SS to unit interval
SS.max <- max(SS)
SS.min <- min(SS)
SS <- 1 - (SS.max - SS) / (SS.max - SS.min)
SS <- matrix(SS, ncol = 1)
plot(as.vector(X), as.vector(SS), xlab = 'Design Point', ylab = 'Sufficient Statistic')

# Fit GP
gp.fit <- mlegp(X, SS, nugget.known = 0, constantMean = 1)


#
# Parameter Calibration with Gaussian Process Approximation, propagating uncertainty by 
# integrating out GP.
#

# Compile Stan code
stan.gp.model.path <- file.path(base.dir, 'parameter_calibration_gp_1d_example.stan')
model.gp <- stan_model(stan.gp.model.path)

# Data to pass to Stan
stan.list.gp.params <- create.gp.params.list(gp.fit, "mlegp")
stan.list.other <- list(N = N, 
                        n = n,
                        y = as.vector(SS),
                        X = X, 
                        a = tau.shape, 
                        b = tau.rate, 
                        u_mean = u.mean, 
                        u_sigma = u.sigma)
stan.list.gp <- as.list(c(stan.list.gp.params, stan.list.other))

# MCMC
fit.gp <- sampling(model.gp, stan.list.gp, iter = 50000, chains = 4)
summary(fit.gp)
samples.gp <- extract(fit.gp)
sampler.params <- get_sampler_params(fit.gp, inc_warmup = TRUE)
chain1 <- sampler.params[[1]]

tau.samples <- extract(fit.gp, inc_warmup = TRUE)[["tau"]]



#
# Parameter Calibration with Gaussian Process Approximation, evaluating at mean function.
#
















