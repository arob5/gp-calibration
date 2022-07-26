# likelihood_approx_comparison.R
# Test code that compares the true Gaussian likelihood to the approximate likelihood
# with GP mean replacing the sufficient statistic, and to approximate likelihood
# where the GP has been integrated out to account for interpolation uncertainty. 
#
# Andrew Roberts
# Working Directory: /projectnb2/dietzelab/arober/pecan_personal

#
# TODO: code up the covariance function in R (copy from other file) and try to 
#       re-produce the MLE GP standard errors. 
#


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


# Create sub-directory in output directory
tag <- 'extrap'
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
u.extrapolation.width <- 5

n.u <- 1000
u.vals <- seq(u.min - u.extrapolation.width, u.max + u.extrapolation.width, length = n.u)

# Precision parameter values.
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

#
# Tests 
# TODO: Move these to a separate testing file

# # 1.) Test that kernel gives same results
# gp <- gp.fits[[1]]
# rho <- 1 / (2 * gp$beta)
# alpha <- sqrt(gp$sig2)
# sigma <- sqrt(gp$nugget)
# mu <- gp$Bhat
# 
# K.r <- K(X, rho = rho, alpha = alpha) + diag(rep(sigma^2, nrow(X)))
# K.mlegp <- calcVarMatrix(X, gp$beta, gp$a, gp$nugget, gp$sig2, 0, gp$numObs)
# print(paste0("Cov matrices K(X) equal: ", all.equal(K.r, K.mlegp)))
# 
# # 2.) Test that mean prediction gives same results
# gp.mean.test <- predict_mean(X, X, SS[, 1],
#                              rho = rho, 
#                              alpha = alpha, 
#                              sigma = sigma, 
#                              mu = mu)
# 
# gp.pred <- predict(gp, se.fit = TRUE)
# gp.mean <- gp.pred$fit
# max(abs(as.vector(gp.mean.test) - as.vector(gp.mean)))
# 
# # 3.) Test that standard error prediction gives same results
# c <- calcCorOneObs(X, gp$beta, gp$a, X[1,]) * alpha^2
# c.test <- K(as.matrix(X[1,]), X, rho, alpha)
# print(paste0("Covariance between design and one test point equal: ", all.equal(c, c.test)))
# 
# K.inv.test <- solve(K(X, X, rho, alpha) + diag(rep(sigma^2, nrow(X))))
# print(paste0("Inverse cov matrices equal: ", all.equal(K.inv.test, gp.fits[[1]]$invVarMatrix)))
# 
# v <- calcPredictionError(gp, X[1,], nugget = gp$nugget)
# v.test <- sqrt(predict_var(X[1,,drop=FALSE], X, SS[,1], rho, alpha, sigma))
# print(paste0("Predictive standard errors equal (single test point): ", all.equal(v, v.test)))
# 
# gp.se.test <- sqrt(diag(predict_var(X, X, SS[, 1], rho, alpha, sigma)))
# gp.se <- as.vector(gp.pred$se.fit)
# max(abs(gp.se.test - gp.se))




# True likelihood at design points
llik.design <- matrix(NA, nrow = N, ncol = length(tau.vals))
for(i in seq(1, N)) {
  for(j in seq_along(tau.vals)) {
    llik.design[i,j] <- dmvnorm.log.unnorm(y.obs[,j], X[i,1], tau.vals[j])
  }
}

# GP predictive mean and standard error at test points
interval.pct <- .8
p <- 1 - (1 - interval.pct)/2

gp.pred.means <- matrix(NA, nrow = n.u, ncol = length(tau.vals))
gp.pred.se <- matrix(NA, nrow = n.u, ncol = length(tau.vals))
gp.pred.upper <- matrix(NA, nrow = n.u, ncol = length(tau.vals))
gp.pred.lower <- matrix(NA, nrow = n.u, ncol = length(tau.vals))
for(j in seq_along(tau.vals)) {
  gp.pred <- predict.gp(gp.fits[[j]], matrix(u.vals, ncol=1), se.fit = TRUE)
  gp.pred.means[,j] <- gp.pred$fit
  gp.pred.se[,j] <- sqrt((gp.pred$se.fit)^2 + gp.fits[[j]]$nugget)
  gp.pred.upper[,j] <- qnorm(p, gp.pred.means[,j], gp.pred.se[,j])
  gp.pred.lower[,j] <- qnorm(p, gp.pred.means[,j], gp.pred.se[,j], lower.tail = FALSE)
}

png(file.path(out.dir, 'gp_pred_mean_SS.png'), width=600, height=350)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(matrix(u.vals, ncol=1), gp.pred.means, xlab = 'u', type = 'l', 
        lty = 1, ylab = paste0('GP Predictive Mean and ', 100*interval.pct, '% interval'),
        main = 'GP Predictive Mean of SS', 
        ylim = c(min(gp.pred.lower), max(gp.pred.upper)))
matpoints(X, SS, pch = 16)
matlines(matrix(u.vals, ncol=1), gp.pred.upper)
matlines(matrix(u.vals, ncol=1), gp.pred.lower)
dev.off()

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
stan.output <- vector(mode = 'list', length = length(tau.vals))
for(j in seq_along(tau.vals)) {
  stan.list <- list(N = N,
                    n = n, 
                    m = n.u, 
                    tau = tau.vals[j],
                    X = X,
                    y = SS[,j],
                    u_vals = matrix(u.vals, ncol=1),
                    gp_rho = array(1 / (2*gp.fits[[j]]$beta), dim = 1),
                    gp_alpha = sqrt(gp.fits[[j]]$sig2),
                    gp_sigma = sqrt(gp.fits[[j]]$nugget), 
                    gp_mean = gp.fits[[j]]$Bhat
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
# points(u, 0, col="blue", pch=16)
matpoints(X, llik.design, pch=16)
legend("right", inset=c(-0.2,-0.3), legend=tau.vals, 
       col=seq_along(tau.vals), 
       fill=seq_along(tau.vals), title="tau")
dev.off()


#
# Validating GP integrated out Stan calculations by producing the same plot using 
# R code. 
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
dev.off()

#
# Tests
#

# GP Interpolates (predictive mean at design points)
for(j in seq_along(tau.vals)) {
  print(paste0("GP Interpolates: ", all.equal(stan.output[[j]]$mean_design[1,], SS[,j])))
}

# Predictive variance at design points
head(sqrt(stan.output[[1]]$var_design[1,]))


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
matpoints(X, SS, pch = 16)
matlines(matrix(u.vals, ncol=1), gp.pred.upper.stan)
matlines(matrix(u.vals, ncol=1), gp.pred.lower.stan)
dev.off()


