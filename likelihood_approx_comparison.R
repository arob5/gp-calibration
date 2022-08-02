# likelihood_approx_comparison.R
# Test code that compares the true Gaussian likelihood to the approximate likelihood
# with GP mean replacing the sufficient statistic, and to approximate likelihood
# where the GP has been integrated out to account for interpolation uncertainty. 
#
# Andrew Roberts
# Working Directory: /projectnb2/dietzelab/arober/pecan_personal

#
# TODO:
#   - Re-write code to only run with one tau value at a time, and clean up plots/create new legends
#   - Think about GP uncertainty; why does GP seem more confident than I would expect?
#   - Run some actual parameter calibration examples.


library(mvtnorm)
library(rstan)
library(mlegp)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("stan.helper.functions.R")

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

# The function that is acting as the computer model. Should be equivalent to 
# the function of the same name in the Stan file. 
f <- function(u) {
  return(u)
}

f.string <- "u"

# Directories
base.dir <- '.'
# base.out.dir <- file.path(base.dir, '..', 'test_output')
base.out.dir <- file.path(base.dir, 'output')

# Create sub-directory in output directory
tag <- '2'
subdir.name <- paste0('llik_approx_comparison_', tag)
out.dir <- file.path(base.out.dir, subdir.name)
dir.create(out.dir)

# Write text file with function f()
file.con <- file(file.path(out.dir, "f.txt"))
writeLines(f.string, file.con)
close(file.con)


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

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

# Log unnormalized isotropic multivariate normal density, integrated over GP
# approximation of sufficient statistic
dmvnorm.gp.approx <- function(tau, n, gp.mean, gp.var) {
  (n/2)*log(tau) - (tau/2)*gp.mean  + (tau^2/8)*gp.var
}

# Squared Exponential (i.e. Gaussian) kernel (covariance)
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

# Squared Exponential (i.e. Gaussian) kernel (correlation)
C <- function(X1, X2 = NA, rho, alpha) {
  if(is.na(X2[1][1])) X2 <- X1
  cor.mat <- matrix(NA, nrow = nrow(X1), ncol = nrow(X2))
  
  for(i in 1:nrow(X1)) {
    for(j in 1:nrow(X2)) {
      cor.mat[i, j] <- sum((1 / rho) * (X1[i, ] - X2[j, ])^2)
    }
  }
  
  return(exp(-0.5 * cor.mat))
}

predict_mean <- function(X_pred, X_obs, y_obs, rho, alpha, sigma, mu) {
  if(sigma == 0) {
    eps <- sqrt(.Machine$double.eps) # "Jitter" to ensure positive definite matrices
  } else {
    eps <- sigma^2
  }
  
  cross_cov <- K(X_pred, X_obs, rho, alpha)
  data_cov <- K(X_obs, X_obs, rho, alpha) + diag(rep(eps, nrow(X_obs)))
  return(mu + cross_cov %*% solve(data_cov) %*% (y_obs - mu))
}

predict_var<- function(X_pred, X_obs, y_obs, rho, alpha, sigma) {
  if(sigma == 0) {
    eps <- sqrt(.Machine$double.eps) # "Jitter" to ensure positive definite matrices
  } else {
    eps <- sigma^2
  }
  
  prior_var <- K(X_pred, X_pred, rho, alpha) + diag(rep(eps, nrow(X_pred)), nrow = nrow(X_pred))
  cross_corr <- C(X_pred, X_obs, rho, alpha)
  data_cov <- K(X_obs, X_obs, rho, alpha) + diag(rep(eps, nrow(X_obs)))
  
  var.mat <- prior_var - alpha^4 * (cross_corr %*% solve(data_cov) %*% t(cross_corr))
  vars <- diag(var.mat)
  vars[vars < 0] <- 0
  diag(var.mat) <- vars
  
  return(var.mat)
}


# -----------------------------------------------------------------------------
# Simulate data
#   Fixing "true" values of calibration parameter u and precision parameter 
#   tau, then simulating observed data using these values.
# -----------------------------------------------------------------------------

# True calibration parameter value
u <- .5

# Support of calibration parameter
u.min <- 0
u.max <- 1
u.extrapolation.width <- .01

n.u <- 1000
u.vals <- seq(u.min - u.extrapolation.width, u.max + u.extrapolation.width, length = n.u)
X.pred <- matrix(u.vals, ncol=1)

# Precision parameter values.
sigma <- .3
tau <- 1/sigma^2

# Simulate observed data
n <- 1000
y.obs <- matrix(NA, nrow = n, ncol = 1)
y.obs[,1] <- rnorm(n, f(u), sigma)


# Design matrix, evaluate model at design points
N <- 10
design.points <- c(0, .05, .07, .08, .1, .2, .25, .3, 1)
N <- length(design.points)
X <- matrix(design.points, ncol = 1)
y_model <- f(X)


# -----------------------------------------------------------------------------
# Calculate true values of log-likelihood at design points and test points.
# -----------------------------------------------------------------------------

# True log-likelihood at design points
llik.design <- matrix(NA, nrow = N, ncol = 1)
for(i in seq(1, N)) {
  llik.design[i,1] <- dmvnorm.log.unnorm(y.obs[,1], X[i,1], tau)
}

# True log-likelihood at test points
llik.pred <- matrix(NA, nrow = n.u, ncol = 1)
for(i in seq(1, n.u)) {
  llik.pred[i, 1] <- dmvnorm.log.unnorm(y.obs[,1], u.vals[i], tau)
}

png(file.path(out.dir, 'exact_llik.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(matrix(u.vals, ncol = 1), llik.pred, type = 'l', 
        lty=1, xlab="u", ylab="Unnormalized Log-Likelihood", 
        main = "Exact Log-Likelihood", col="red")
dev.off()

# -----------------------------------------------------------------------------
# Fit GP regression
# -----------------------------------------------------------------------------

# Define the sufficient statistic (SS). For simplicity here, I'm not considering any other explanatory variables
# x that are being conditioned on. If there were, we would have to make sure we were lining up the observed data
# and the process model data correctly conditional on x. So in this case, y_model[1] is just a constant predictor
# of the output variable given the first calibration parameter (design point) X[1,1]. 
SS <- matrix(NA, nrow = N, ncol = 1)
for(i in seq(1, N)) {
  SS[i, 1] <- sum((y_model[i] - y.obs[,1])^2)
}

# Similarly calculate true SS at test points for reference in subsequent plots
SS.pred <- matrix(NA, nrow = n.u, ncol = 1)
for(i in seq(1, n.u)) {
  SS.pred[i, 1] <- sum((f(u.vals[i]) - y.obs[,1])^2)
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
gp.fit <- mlegp(X, SS, nugget.known = 0, constantMean = 1)

# -----------------------------------------------------------------------------
# GP predictive mean plots: 
#   - Predictive mean of sufficient statistic with confidence interval
#   - Log-likelihood evaluated at sufficient statistic with confidence interval
# -----------------------------------------------------------------------------

# GP predictive mean and standard error at test points
interval.pct <- .95
p <- 1 - (1 - interval.pct)/2

gp.pred <- predict.gp(gp.fit, X.pred, se.fit = TRUE)
gp.pred.means <- gp.pred$fit
gp.pred.se <- sqrt((gp.pred$se.fit)^2 + gp.fit$nugget)
gp.pred.upper <- qnorm(p, gp.pred.means, gp.pred.se)
gp.pred.lower <- qnorm(p, gp.pred.means, gp.pred.se, lower.tail = FALSE)

png(file.path(out.dir, 'gp_pred_mean_SS.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(X.pred, gp.pred.means, xlab = 'u', type = 'l', lty = 1, 
     ylab = paste0('GP Predictive Mean and ', 100*interval.pct, '% CI'),
     main = 'GP Predictive Mean of SS', 
     ylim = c(min(gp.pred.lower), max(gp.pred.upper)),
     col = 'blue')
points(X, SS, pch = 16)
matlines(X.pred, gp.pred.upper, lty=5, col="gray")
matlines(X.pred, gp.pred.lower, lty=5, col="gray")
dev.off()


matplot(matrix(u.vals, ncol=1), gp.pred.means, xlab = 'u', type = 'l', 
        lty = 1, ylab = paste0('GP Predictive Mean and ', 100*interval.pct, '% interval'),
        main = 'GP Predictive Mean of SS', 
        ylim = c(min(gp.pred.lower), max(gp.pred.upper)), 
        col = "red")
matpoints(X, SS, pch = 16)
matlines(matrix(u.vals, ncol=1), gp.pred.upper, lty=5, col="gray")
matlines(matrix(u.vals, ncol=1), gp.pred.lower, lty=5, col="gray")
matlines(matrix(u.vals, ncol=1), SS.pred, col="blue")


# Likelihood using GP mean approximation
likelihood.gp.mean <- matrix(NA, nrow = n.u, ncol = length(sigma.vals))
likelihood.gp.upper <- matrix(NA, nrow = n.u, ncol = length(tau.vals))
likelihood.gp.lower <- matrix(NA, nrow = n.u, ncol = length(tau.vals))

for(i in seq(1, n.u)) {
  for(j in seq_along(tau.vals)) {
    likelihood.gp.mean[i, j] <- dmvnorm.log.unnorm.SS(gp.pred.means[i, j], tau.vals[j], n)
    likelihood.gp.upper[i, j] <- dmvnorm.log.unnorm.SS(gp.pred.upper[i, j], tau.vals[j], n)
    likelihood.gp.lower[i, j] <- dmvnorm.log.unnorm.SS(gp.pred.lower[i, j], tau.vals[j], n)
  }
}

png(file.path(out.dir, 'gp_mean_llik.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(matrix(u.vals, ncol = 1), likelihood.gp.mean, type = 'l', 
        lty=1, xlab="u", 
        ylab=paste0("Unnormalized Log-Likelihood and ", 100*interval.pct, "% interval"), 
        main = "GP Mean Approx Log-Likelihood")
legend("right", inset=c(-0.2,-0.3), legend=tau.vals, 
       col=seq_along(tau.vals), 
       fill=seq_along(tau.vals), title="tau")
matlines(matrix(u.vals, ncol=1), likelihood.gp.upper, lty=3)
matlines(matrix(u.vals, ncol=1), likelihood.gp.lower, lty=3)
matlines(matrix(u.vals, ncol=1), llik.pred, lty=5)
matpoints(X, llik.design, pch=16)
dev.off()


# -----------------------------------------------------------------------------
# Incorporate uncertainty via "GP integrated out" log-likelihood
# -----------------------------------------------------------------------------

# Compile Stan code
stan.model.path <- file.path(base.dir, 'likelihood_approx_comparison.stan')
model <- stan_model(stan.model.path)

# Calculate GP integrated out likelihood in Stan (uq = uncertainty quantification)
likelihood.gp.uq <- matrix(NA, nrow = n.u, ncol = length(sigma.vals))
stan.output <- vector(mode = 'list', length = length(tau.vals))
for(j in seq_along(tau.vals)) {
  gp.stan.params <- create.gp.params.list(gp.fits[[j]], "mlegp")
  stan.list <- as.list(c(gp.stan.params, 
                         list(N = N,
                              n = n, 
                              N_pred = n.u, 
                              tau = tau.vals[j],
                              X = X,
                              y = SS[,j],
                              u_vals = matrix(u.vals, ncol=1)
                             )
                         )
                        )
  
  stan.fit <- sampling(model, data = stan.list, warmup = 0, iter = 1, chains = 1, 
                       seed = 494838, refresh = 4000, algorithm = "Fixed_param")
  stan.output[[j]] <- extract(stan.fit)
  likelihood.gp.uq[, j] <- as.vector(stan.output[[j]]$u_vals_llik)
}  

png(file.path(out.dir, 'gp_integrated_out_llik.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(matrix(u.vals, ncol = 1), likelihood.gp.uq, type = 'l', 
        lty=1, xlab="u", ylab="Unnormalized Log-Likelihood", 
        main = "GP Integrated Out Approx Log-Likelihood")
matpoints(X, llik.design, pch=16)
legend("right", inset=c(-0.2,-0.3), legend=tau.vals, 
       col=seq_along(tau.vals), 
       fill=seq_along(tau.vals), title="tau")
matlines(matrix(u.vals, ncol=1), llik.pred, lty=5)
dev.off()


# -----------------------------------------------------------------------------
# Tests/Validation: 
# -----------------------------------------------------------------------------

#
# Re-produce predictive means plot using Stan functionality instead of R
#

interval.pct <- .8
p <- 1 - (1 - interval.pct)/2

gp.pred.means.stan <- matrix(NA, nrow = n.u, ncol = length(tau.vals))
gp.pred.se.stan <- matrix(NA, nrow = n.u, ncol = length(tau.vals))
gp.pred.upper.stan <- matrix(NA, nrow = n.u, ncol = length(tau.vals))
gp.pred.lower.stan <- matrix(NA, nrow = n.u, ncol = length(tau.vals))
for(j in seq_along(tau.vals)) {
  gp.pred.means.stan[,j] <- stan.output[[j]]$mean_test
  gp.pred.se.stan[,j] <- sqrt(stan.output[[j]]$var_test)
  gp.pred.upper.stan[,j] <- qnorm(p, gp.pred.means.stan[,j], gp.pred.se.stan[,j])
  gp.pred.lower.stan[,j] <- qnorm(p, gp.pred.means.stan[,j], gp.pred.se.stan[,j], lower.tail = FALSE)
}

png(file.path(out.dir, 'gp_pred_mean_SS_stan.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(matrix(u.vals, ncol=1), gp.pred.means.stan, xlab = 'u', type = 'l', 
        lty = 1, ylab = paste0('GP Predictive Mean and ', 100*interval.pct, '% interval'),
        main = 'GP Predictive Mean of SS', 
        ylim = c(min(gp.pred.lower.stan), max(gp.pred.upper.stan)))
legend("right", inset=c(-0.2,-0.3), legend=tau.vals, 
       col=seq_along(tau.vals), 
       fill=seq_along(tau.vals), title="tau")
matpoints(X, SS, pch = 16)
matlines(matrix(u.vals, ncol=1), gp.pred.upper.stan, lty=3)
matlines(matrix(u.vals, ncol=1), gp.pred.lower.stan, lty=3)
matlines(matrix(u.vals, ncol=1), SS.pred, lty=5)
dev.off()


#
# Re-produce GP integrated out plot using R instead of Stan
#

likelihood.gp.uq.test <- matrix(NA, nrow = n.u, ncol = length(sigma.vals))

for(j in seq_along(tau.vals)) {
  gp.pred <- predict.gp(gp.fits[[j]], matrix(u.vals, ncol=1), se.fit = TRUE)
  gp.pred.mean <- as.vector(gp.pred$fit)
  gp.pred.var <- as.vector(gp.pred$se.fit)^2 + gp.fits[[j]]$nugget
  for(i in seq(1, n.u)) {
    likelihood.gp.uq.test[i, j] <- dmvnorm.gp.approx(tau.vals[j], n, gp.pred.mean[i], gp.pred.var[i])
  }
}

png(file.path(out.dir, 'gp_integrated_out_llik_r_test.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(matrix(u.vals, ncol = 1), likelihood.gp.uq.test, type = 'l', 
        lty=1, xlab="u", ylab="Unnormalized Log-Likelihood", 
        main = "GP Integrated Out Approx Log-Likelihood")
# points(u, 0, col="blue", pch=16)
matpoints(X, llik.design, pch=16)
legend("right", inset=c(-0.2,-0.3), legend=tau.vals, 
       col=seq_along(tau.vals), 
       fill=seq_along(tau.vals), title="tau")
matlines(matrix(u.vals, ncol=1), llik.pred, lty=5)
dev.off()






