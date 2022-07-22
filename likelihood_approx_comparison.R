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
  return( (7.0 + 4.0*u + u^2)/147 )
}

f.string <- "(10.0 + 4.0*u + u^3) / 147"

# Directories
base.dir <- '.'
base.out.dir <- file.path(base.dir, '..', 'test_output')

# Log unnormalized isotropic multivariate normal density
dmvnorm.log.unnorm <- function(y, u, tau) {
  n <- length(y)
  mu.vec <- rep(f(u), n)
  (n/2)*log(tau) - (tau/2)*sum((y - mu.vec)^2)
}

# Log unnormalized isotropic multivariate normal density, as function of 
# sufficient statistic
dmvnorm.log.unnorm.SS <- function(SS, tau, n) {
  (n/2)*log(tau) - (tau/2)*SS
}

# Create sub-directory in output directory
tag <- 'design_gaps'
subdir.name <- paste0('llik_approx_comparison_', tag)
out.dir <- file.path(base.out.dir, subdir.name)
dir.create(out.dir)

# Write text file with function f()
file.con <- file(file.path(out.dir, "f.txt"))
writeLines(f.string, file.con)
close(file.con)


#
# Define support of calibration parameter u and consider set of fixed values of 
# precision parameter tau for Gaussian likelihood.
#

# True calibration parameter value
u <- 4

# Support of calibration parameter
u.min <- -2
u.max <- 10
u.extrapolation.width <- 0

n.u <- 1000
u.vals <- seq(u.min - u.extrapolation.width, u.max + u.extrapolation.width, length = n.u)

# Precision parameter values. Note that the model output will be scaled to be 
# in the unit interval. 
sigma.vals <- c(.1, .2)
tau.vals <- 1/sigma.vals^2


#
# Simulate observed data.
#

# Simulate observed data
n <- 1000
y.obs <- matrix(NA, nrow = n, ncol = length(sigma.vals))

for(j in seq_along(sigma.vals)) {
  y.obs[,j] <- rnorm(n, f(u), sigma.vals[j])
}


#
# Plot likelihood
#

likelihood.exact <- matrix(NA, nrow = n.u, ncol = length(sigma.vals))

for(i in seq(1, n.u)) {
  for(j in seq_along(tau.vals)) {
    likelihood.exact[i, j] <- dmvnorm.log.unnorm(y.obs[,j], u.vals[i], tau.vals[j])
  }
}

png(file.path(out.dir, 'exact_llik.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(matrix(u.vals, ncol = 1), likelihood.exact, type = 'l', 
        lty=1, xlab="u", ylab="Unnormalized Log-Likelihood", main = "Exact Log-Likelihood")
legend("right", inset=c(-0.2,-0.3), legend=tau.vals, 
       col=seq_along(tau.vals), 
       fill=seq_along(tau.vals), title="tau")
dev.off()

#
# Fit GP regression
#

# Design matrix, evaluate model at design points
N <- 10
design.points <- c(seq(-2, 0, length=4), seq(0.5, 2, length=7), 
                   seq(2.5, 5, length=3), seq(5.5, 8, length=3), seq(8.5, 10, length=10))
N <- length(design.points)
X <- matrix(design.points, ncol = 1)
y_model <- f(X)

# Define the sufficient statistic (SS). For simplicity here, I'm not considering any other explanatory variables
# x that are being conditioned on. If there were, we would have to make sure we were lining up the observed data
# and the process model data correctly conditional on x. So in this case, y_model[1] is just a constant predictor
# of the output variable given the first calibration parameter (design point) X[1,1]. 
SS <- matrix(NA, nrow = N, ncol = length(tau.vals))
for(i in seq(1, N)) {
  for(j in seq_along(tau.vals)) {
    SS[i, j] <- sum((y_model[i] - y.obs[,j])^2)
  }
}

# Scale to unit interval
# max.SS <- apply(SS, 2, max)
# min.SS <- apply(SS, 2, min)
# SS.scale <- matrix(NA, nrow = N, ncol = length(tau.vals))
# for(j in seq_along(tau.vals)) {
#   SS.scale[, j] <- (SS[, j] - min.SS[j]) / (max.SS[j] - min.SS[j])
# }


matplot(X, SS, xlab = 'Design Point', pch = 16, ylab = 'Sufficient Statistic')

# Fit GP
gp.fits <- vector(mode = 'list', length = length(tau.vals))
for(j in seq_along(tau.vals)) {
  gp.fits[[j]] <- mlegp(X, SS[,j], nugget.known = 0, constantMean = 1)
}

# True likelihood at design points
llik.design <- matrix(NA, nrow = N, ncol = length(tau.vals))
for(i in seq(1, N)) {
  for(j in seq_along(tau.vals)) {
    llik.design[i,j] <- dmvnorm.log.unnorm(y.obs[,j], X[i,1], tau.vals[j])
  }
}

# GP predictive mean at test points
gp.pred.means <- matrix(NA, nrow = n.u, ncol = length(tau.vals))
for(j in seq_along(tau.vals)) {
  gp.pred.means[,j] <- predict.gp(gp.fits[[j]], matrix(u.vals, ncol=1))
}

matplot(matrix(u.vals, ncol=1), gp.pred.means, xlab = 'u', type = 'l', 
        lty = 1, ylab = 'GP Predictive Mean', main = 'GP Predictive Mean of SS')

likelihood.gp.mean <- matrix(NA, nrow = n.u, ncol = length(sigma.vals))
for(i in seq(1, n.u)) {
  for(j in seq_along(tau.vals)) {
    likelihood.gp.mean[i, j] <- dmvnorm.log.unnorm.SS(gp.pred.means[i, j], tau.vals[j], n)
  }
}

png(file.path(out.dir, 'gp_mean_llik.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(matrix(u.vals, ncol = 1), likelihood.gp.mean, type = 'l', 
        lty=1, xlab="u", ylab="Unnormalized Log-Likelihood", 
        main = "GP Mean Approx Log-Likelihood")
# points(u, 0, col="blue", pch=16)
matpoints(X, llik.design, pch=16)
legend("right", inset=c(-0.2,-0.3), legend=tau.vals, 
       col=seq_along(tau.vals), 
       fill=seq_along(tau.vals), title="tau")
dev.off()

#
# Now incorporate uncertainty via "GP integrated out" log-likelihood
#

# Compile Stan code
stan.model.path <- file.path(base.dir, 'likelihood_approx_comparison.stan')
model <- stan_model(stan.model.path)

# Calculate GP integrated out likelihood in Stan (uq = uncertainty quantification)
likelihood.gp.uq <- matrix(NA, nrow = n.u, ncol = length(sigma.vals))
for(j in seq_along(tau.vals)) {
  stan.list <- list(N = N,
                    n = n, 
                    m = n.u, 
                    tau = tau.vals[j],
                    X = X,
                    y = SS[,j],
                    u_vals = matrix(u.vals, ncol=1),
                    gp_rho = array(gp.fits[[j]]$beta, dim = 1),
                    gp_alpha = sqrt(gp.fits[[j]]$sig2),
                    gp_sigma = sqrt(gp.fits[[j]]$nugget), 
                    gp_mean = gp.fits[[j]]$Bhat
                   )
  
  stan.fit <- sampling(model, data = stan.list, warmup = 0, iter = 1, chains = 1, 
                       seed = 494838, refresh = 4000, algorithm = "Fixed_param")
  stan.output <- extract(stan.fit)
  likelihood.gp.uq[, j] <- as.vector(stan.output$u_vals_llik)
}  

png(file.path(out.dir, 'gp_integrated_out_llik.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(matrix(u.vals, ncol = 1), likelihood.gp.uq, type = 'l', 
        lty=1, xlab="u", ylab="Unnormalized Log-Likelihood", 
        main = "GP Integrated Out Approx Log-Likelihood")
# points(u, 0, col="blue", pch=16)
matpoints(X, llik.design, pch=16)
legend("right", inset=c(-0.2,-0.3), legend=tau.vals, 
       col=seq_along(tau.vals), 
       fill=seq_along(tau.vals), title="tau")
dev.off()









