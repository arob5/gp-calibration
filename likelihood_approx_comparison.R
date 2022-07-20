# likelihood_approx_comparison.R
# Test code that compares the true Gaussian likelihood to the approximate likelihood
# with GP mean replacing the sufficient statistic, and to approximate likelihood
# where the GP has been integrated out to account for interpolation uncertainty. 
#
# Andrew Roberts
# Working Directory: /projectnb2/dietzelab/arober/pecan_personal

library(mvtnorm)
library(rstan)
library(mlegp)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#
# Settings
#

# The function that is acting as the computer model. Should be equivalent to 
# the function of the same name in the Stan file. 
f <- function(u) {
  return( 10.0 + 4.0*u - u^3 )
}

# Directories
base.dir <- '.'
out.dir <- file.path(base.dir, '..', 'test_output')

# Log unnormalized isotropic multivariate normal density
dmvnorm.log.unnorm <- function(y, u, tau) {
  n <- length(y)
  mu.vec <- rep(f(u), n)
  (n/2)*log(tau) - (tau/2)*sum((y - mu.vec)^2)
}


#
# Define support of calibration parameter u and consider set of fixed values of 
# precision parameter tau for Gaussian likelihood.
#

# Support of calibration parameter
u.min <- -2
u.max <- 2

# Precision parameter values. Note that the model output will be scaled to be 
# in the unit interval. 
sigma.vals <- c(.1, .2)
tau.vals <- 1/sigma.vals^2


#
# Define "true" value of calibration parameter and simulate observed data.
#

# True calibration parameter value
u <- 0

# Simulate observed data
n <- 1000
y.obs <- matrix(NA, nrow = n, ncol = length(sigma.vals))

for(j in seq_along(sigma.vals)) {
  y.obs[,j] <- rnorm(n, f(u), sigma.vals[j])
}


#
# Plot likelihood
#

n.u <- 1000
u.vals <- seq(u.min, u.max, length = n.u)
likelihood.exact <- matrix(NA, nrow = n.u, ncol = length(sigma.vals))

for(i in seq(1, n.u)) {
  for(j in seq_along(tau.vals)) {
    likelihood.exact[i, j] <- dmvnorm.log.unnorm(y.obs[,j], u.vals[i], tau.vals[j])
  }
}

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(matrix(u.vals, ncol = 1), likelihood.exact, type = 'l', 
        lty=1, xlab="u", ylab="Unnormalized Log-Likelihood")
points(u, 0, col="blue", pch=16)
legend("topright", inset=c(-0.3,-0.3), legend=tau.vals, 
       col=seq_along(tau.vals), 
       fill=seq_along(tau.vals), title="tau")



#
# Fit GP regression
#

# Design matrix, evaluate model at design points
N <- 100
X <- matrix(seq(qnorm(.025, u.mean, u.sigma), qnorm(.975, u.mean, u.sigma), length = N), ncol = 1)
y_model <- f(X)

# Define the sufficient statistic (SS). For simplicity here, I'm not considering any other explanatory variables
# x that are being conditioned on. If there were, we would have to make sure we were lining up the observed data
# and the process model data correctly conditional on x. So in this case, y_model[1] is just a constant predictor
# of the output variable given the first calibration parameter (design point) X[1,1]. 
SS <- rep(0, N)
for(i in seq(1, N)) {
  SS[i] <- sum((y_model[i] - y)^2)
}

# Scale SS to unit interval
SS.max <- max(SS)
SS.min <- min(SS)
SS <- 1 - (SS.max - SS) / (SS.max - SS.min)
SS <- matrix(SS, ncol = 1)
plot(as.vector(X), as.vector(SS), xlab = 'Design Point', ylab = 'Sufficient Statistic')

# Fit GP
gp_mle <- mlegp(X, SS, nugget.known = 0, constantMean = 1)



