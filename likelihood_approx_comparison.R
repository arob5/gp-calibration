# likelihood_approx_comparison.R
# Test code that compares the true Gaussian likelihood to the approximate likelihood
# with GP mean replacing the sufficient statistic, and to approximate likelihood
# where the GP has been integrated out to account for interpolation uncertainty. 
#
# Andrew Roberts
# Working Directory: /projectnb2/dietzelab/arober/pecan_personal
#
# TODO:
#   - Try running example that fixed GP params instead of fitting them
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

# Random seed (for generating random data)
seed <- 10
set.seed(seed)

# The function that is acting as the computer model. Should be equivalent to 
# the function of the same name in the Stan file. 
f <- function(u) {
  return(u)
}

f.string <- "u"

# If TRUE, sets default GP params instead of fitting the params via mlegp
use.default.gp.params <- TRUE
gp.param.defaults <- list(gp_rho = array(0.5), 
                          gp_alpha = 300, 
                          gp_sigma = sqrt(.Machine$double.eps), 
                          gp_mean = 0.0)

# Directories
base.dir <- '.'
# base.out.dir <- file.path(base.dir, '..', 'test_output')
base.out.dir <- file.path(base.dir, 'output')

# Create sub-directory in output directory
tag <- 'pres_7'
subdir.name <- paste0('llik_approx_comparison_', tag)
out.dir <- file.path(base.out.dir, subdir.name)
dir.create(out.dir)


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
      C[i, j] <- sum((1 / rho^2) * (X1[i, ] - X2[j, ])^2)
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
      cor.mat[i, j] <- sum((1 / rho^2) * (X1[i, ] - X2[j, ])^2)
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

predict_var <- function(X_pred, X_obs, y_obs, rho, alpha, sigma) {
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

predict_mean_chol <- function(X_pred, X_obs, y_obs, rho, alpha, sigma, mu) {
  
  # Nugget
  if(sigma == 0) {
    eps <- sqrt(.Machine$double.eps)
  } else {
    eps <- sigma^2
  }
  
  L <- t(chol(K(X_obs, X_obs, rho, alpha) + diag(rep(eps, nrow(X_obs)))))
  K.inv.y <- solve(t(L), solve(L, y_obs - mu))
  cross_cov <- K(X_pred, X_obs, rho, alpha)
  
  return(mu + cross_cov %*% K.inv.y)
  
}


predict_var_chol <- function(X_pred, X_obs, y_obs, rho, alpha, sigma) {

  # Nugget
  if(sigma == 0) {
    eps <- sqrt(.Machine$double.eps)
  } else {
    eps <- sigma^2
  }
  
  L <- t(chol(K(X_obs, X_obs, rho, alpha) + diag(rep(eps, nrow(X_obs)))))
  pred.vars <- vector("numeric", nrow(X_pred))
  for(i in seq_along(pred.vars)) {
    k_x <- K(X_pred[i,,drop=FALSE], X_pred[i,,drop=FALSE], rho, alpha) + eps
    k_Xx <- K(X_obs, X_pred[i,,drop=FALSE], rho, alpha)
    v <- solve(L, k_Xx)
    pred.vars[i] <- k_x - sum(v^2)
  }
  
  return(pred.vars)
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
u.extrapolation.width <- .1

n.u <- 1000
u.vals <- seq(u.min - u.extrapolation.width, u.max + u.extrapolation.width, length = n.u)
X.pred <- matrix(u.vals, ncol=1)

# Precision parameter values.
sigma <- 0.3
tau <- 1/sigma^2

# Simulate observed data
n <- 1000
y.obs <- matrix(NA, nrow = n, ncol = 1)
y.obs[,1] <- rnorm(n, f(u), sigma)


# Design matrix, evaluate model at design points
N <- 4
design.points <- seq(u.min, u.max, length.out = N)
# design.points <- c(0, .05, .07, .08, .1, .2, .25, .3, 1)
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

if(!use.default.gp.params) {
  # Fit GP
  gp.fit <- mlegp(X, SS, nugget.known = 0, constantMean = 1)
  
  # Predictive means and standard errors
  gp.pred.mlegp <- predict.gp(gp.fit, X.pred, se.fit = TRUE)
  gp.means.mlegp <- gp.pred.mlegp$fit
  gp.se.mlegp <- sqrt((gp.pred.mlegp$se.fit)^2 + gp.fit$nugget)
}

# -----------------------------------------------------------------------------
# Run Stan code
# -----------------------------------------------------------------------------

# Compile Stan code
stan.model.path <- file.path(base.dir, 'likelihood_approx_comparison.stan')
model <- stan_model(stan.model.path)

# Run Stan code
if(use.default.gp.params) {
  stan.params.gp <- gp.param.defaults
} else {
  stan.params.gp <- create.gp.params.list(gp.fit, "mlegp")
}

stan.params.other <- list(N = N, n = n, N_pred = n.u, tau = tau, X = X, y = SS[,1], u_vals = X.pred)
stan.params <- as.list(c(stan.params.gp, stan.params.other))
stan.fit <- sampling(model, data = stan.params, warmup = 0, iter = 1, chains = 1, 
                     seed = 494838, refresh = 4000, algorithm = "Fixed_param")

# Save Stan data
stan.output <- extract(stan.fit)
llik.gp.uq <- as.vector(stan.output$u_vals_llik) # GP integrated out likelihood in Stan (uq = uncertainty quantification)
gp.means.stan <- as.vector(stan.output$mean_test)
gp.se.stan <- sqrt(as.vector(stan.output$var_test))

# Tests
gp.means.chol <- predict_mean_chol(X.pred, X, SS, stan.params.gp$gp_rho, 
                                   stan.params.gp$gp_alpha, stan.params.gp$gp_sigma, stan.params.gp$gp_mean)
gp.se.chol <- sqrt(predict_var_chol(X.pred, X, SS, stan.params.gp$gp_rho, stan.params.gp$gp_alpha, stan.params.gp$gp_sigma))
print(paste0("Max abs difference in pred means (R, chol) vs pred means (Stan, chol): ", max(abs(gp.means.chol - gp.means.stan))))
print(paste0("Max abs difference in pred SEs (R, chol) vs pred SEs (Stan, chol): ", max(abs(gp.se.chol - gp.se.stan))))

if(!use.default.gp.params) {
  print(paste0("Max abs difference in pred means (R, chol) vs pred means (mlegp): ", max(abs(gp.means.chol - gp.means.mlegp))))
  print(paste0("Max abs difference in pred SEs (R, chol) vs pred SEs (mlegp): ", max(abs(gp.se.chol - gp.se.mlegp))))
}

# -----------------------------------------------------------------------------
# GP predictive mean plots: 
#   - Predictive mean of sufficient statistic with confidence interval
#   - Log-likelihood evaluated at sufficient statistic with confidence interval
# -----------------------------------------------------------------------------

# GP predictive mean and standard error at test points
interval.pct <- .95
p <- 1 - (1 - interval.pct)/2

gp.pred.upper <- qnorm(p, gp.means.stan, gp.se.stan)
gp.pred.lower <- qnorm(p, gp.means.stan, gp.se.stan, lower.tail = FALSE)

png(file.path(out.dir, 'gp_pred_mean_SS.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(X.pred, gp.means.stan, xlab = 'u', type = 'l', lty = 1, 
     ylab = paste0('GP Predictive Mean and ', 100*interval.pct, '% CI'),
     main = 'GP Predictive Mean of SS', 
     ylim = c(min(gp.pred.lower), max(gp.pred.upper)),
     col = 'blue')
points(X, SS, pch = 16, col="black")
lines(X.pred, gp.pred.upper, lty=1, col="gray")
lines(X.pred, gp.pred.lower, lty=1, col="gray")
lines(X.pred, SS.pred, lty=1, col="red")
legend("right", inset=c(-0.2,-0.3), 
       legend = c("GP pred mean", "Design points", paste0(100*interval.pct, "% CI"), "True SS"), 
       col = c("blue", "black", "gray", "red"), lty = c(1, NA, 1, 1), pch = c(NA, 16, NA, NA))
dev.off()


# Log-likelihood using GP mean approximation
llik.gp.mean <- dmvnorm.log.unnorm.SS(gp.means.stan, tau, n)
llik.gp.upper <- dmvnorm.log.unnorm.SS(gp.pred.upper, tau, n)
llik.gp.lower <- dmvnorm.log.unnorm.SS(gp.pred.lower, tau, n)

png(file.path(out.dir, 'gp_mean_llik.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(X.pred, llik.gp.mean, type="l", lty=1, xlab="u", 
     ylab=paste0("Unnormalized Log-Likelihood and ", 100*interval.pct, "% CI"),
     main = "GP Mean Approx Log-Likelihood", col="blue")
points(X, llik.design, pch=16, col="black")
lines(X.pred, llik.gp.upper, lty=1, col="gray")
lines(X.pred, llik.gp.lower, lty=1, col="gray")
lines(X.pred, llik.pred, lty=1, col="red")
legend("right", inset=c(-0.2,-0.3), 
       legend = c("GP pred mean", "Design points", paste0(100*interval.pct, "% CI"), "True llik"), 
       col = c("blue", "black", "gray", "red"), lty = c(1, NA, 1, 1), pch = c(NA, 16, NA, NA))
dev.off()


# -----------------------------------------------------------------------------
# Incorporate uncertainty via "GP integrated out" log-likelihood
# -----------------------------------------------------------------------------

png(file.path(out.dir, 'gp_integrated_out_llik.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(X.pred, llik.gp.uq, type = "l", lty = 1, xlab = "u", ylab="Unnormalized Log-Likelihood",
     main = "GP Integrated Out Approx Log-Likelihood", col = "blue")
points(X, llik.design, pch=16, col="black")
lines(X.pred, llik.pred, lty=1, col="red")
legend("right", inset=c(-0.2,-0.3), 
       legend = c("GP integrated out llik", "Design points", "True llik"), 
       col = c("blue", "black", "red"), lty = c(1, NA, 1), pch = c(NA, 16, NA))
dev.off()

# Plot showing the standard errors
png(file.path(out.dir, 'gp_pred_std_err.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(X.pred, gp.se.stan, type = "l", lty = 1, xlab = "u", ylab = "GP Predictive Std Err", 
     main = "GP Predictive Std Err at Test Points", col = "blue")
# lines(X.pred, gp.se.mlegp, lty = 1, col = "red")
# legend("right", inset=c(-0.2,-0.3), 
#        legend = c("Cholesky (Stan)", "Inverse (mlegp)"), 
#        col = c("blue", "red"), lty = c(1, 1))
dev.off()

# Plot showing the "penalty term": e^(tau^2/8 * k(u,u))
gp.penalty.term <- (tau^2 / 8) * gp.se.stan^2
png(file.path(out.dir, 'gp_penalty_term.png'), width=600, height=350)
plot(X.pred, gp.penalty.term, type = "l", lty = 1, xlab = "u", ylab = "tau^2/8 * k(u,u)", 
     main = "Log-Likelihood Pentalty Term at Test Points", col = "blue")
dev.off()

# -----------------------------------------------------------------------------
# Print information to file
# -----------------------------------------------------------------------------

file.con <- file(file.path(out.dir, "run_info.txt"))
file.text <- c(paste0("Random seed: ", seed), paste0("Model: f(u) = ", f.string), paste0("True u value: ", u),  
               paste0("tau: ", tau), paste0("n: ", n), paste0("Design points: ", paste0(design.points, collapse = ", ")), 
               paste0("use.default.gp.params: ", use.default.gp.params),
               "GP params:", paste0(names(stan.params.gp), collapse = ", "), 
               paste0(stan.params.gp, collapse = ", "))
writeLines(file.text, file.con)
close(file.con)


# -----------------------------------------------------------------------------
# Tests/Validation: 
# -----------------------------------------------------------------------------

#
# Re-produce predictive means plot using R means instead of Stan
#

gp.pred.upper.test <- qnorm(p, gp.means.chol, gp.se.chol)
gp.pred.lower.test <- qnorm(p, gp.means.chol, gp.se.chol, lower.tail = FALSE)

png(file.path(out.dir, 'gp_pred_mean_SS_r_test.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(X.pred, gp.means.chol, xlab = 'u', type = 'l', 
        lty = 1, ylab = paste0('GP Predictive Mean and ', 100*interval.pct, '% CI'),
        main = 'GP Predictive Mean of SS', 
        ylim = c(min(gp.pred.lower.test), max(gp.pred.upper.test)), col = "blue")
points(X, SS, pch = 16, col = "black")
lines(X.pred, gp.pred.upper.test, lty=1, col="gray")
lines(X.pred, gp.pred.lower.test, lty=1, col="gray")
lines(X.pred, SS.pred, lty=1, col="red")
legend("right", inset=c(-0.2,-0.3), 
       legend = c("GP pred mean", "Design points", paste0(100*interval.pct, "% CI"), "True SS"), 
       col = c("blue", "black", "gray", "red"), lty = c(1, NA, 1, 1), pch = c(NA, 16, NA, NA))
dev.off()

#
# Re-produce predictive means likelihood plot using R means instead of Stan
#

llik.gp.mean.test<- dmvnorm.log.unnorm.SS(gp.means.chol, tau, n)
llik.gp.upper.test <- dmvnorm.log.unnorm.SS(gp.pred.upper.test, tau, n)
llik.gp.lower.test <- dmvnorm.log.unnorm.SS(gp.pred.lower.test, tau, n)

png(file.path(out.dir, 'gp_mean_llik_r_test.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(X.pred, llik.gp.mean.test, type="l", lty=1, xlab="u", 
     ylab=paste0("Unnormalized Log-Likelihood and ", 100*interval.pct, "% CI"),
     main = "GP Mean Approx Log-Likelihood", col="blue")
points(X, llik.design, pch=16, col="black")
lines(X.pred, llik.gp.upper.test, lty=1, col="gray")
lines(X.pred, llik.gp.lower.test, lty=1, col="gray")
lines(X.pred, llik.pred, lty=1, col="red")
legend("right", inset=c(-0.2,-0.3), 
       legend = c("GP pred mean", "Design points", paste0(100*interval.pct, "% CI"), "True llik"), 
       col = c("blue", "black", "gray", "red"), lty = c(1, NA, 1, 1), pch = c(NA, 16, NA, NA))
dev.off()


#
# Re-produce GP integrated out plot using R instead of Stan
#

llik.gp.uq.test <- dmvnorm.gp.approx(tau, n, gp.means.chol, gp.se.chol^2)

png(file.path(out.dir, 'gp_integrated_out_llik_r_test.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(X.pred, llik.gp.uq.test, type = 'l', 
     lty=1, xlab="u", ylab="Unnormalized Log-Likelihood", 
     main = "GP Integrated Out Approx Log-Likelihood", col = "blue")
points(X, llik.design, pch=16, col = "black")
lines(X.pred, llik.pred, lty=1, col = "red")
legend("right", inset=c(-0.2,-0.3), 
       legend = c("GP integrated out llik", "Design points", "True llik"), 
       col = c("blue", "black", "red"), lty = c(1, NA, 1), pch = c(NA, 16, NA))
dev.off()






