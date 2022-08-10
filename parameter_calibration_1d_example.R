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

# TODO: 
# - Define bounds on u and tau
# - 0 < sigma^2 < var(data)
# - Bounds on u given by bounds on design points (don't extrapolate)
# Consider scaling response and design prior to fitting GP


library(bayesplot)
library(ggplot2)
library(rstan)
library(mlegp)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("stan.helper.functions.R")
source("helper.functions.R")

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

# Desired standard deviation in the Gaussian model
sigma <- .3

# Hyperparameters (shape and rate) for Gamma hyperprior on tau (precision parameter). 
tau.shape <- 1 / sigma^2
tau.rate <- 1.0

# Number of observations to simulate
n <- 1000

# Gaussian prior parameters on calibration parameter
u.mean <- 0.5
u.sigma <- 0.25

# If TRUE, sets default GP params instead of fitting the params via mlegp
use.default.gp.params <- FALSE
gp.param.defaults <- list(gp_rho = array(0.5), 
                          gp_alpha = 300, 
                          gp_sigma = sqrt(.Machine$double.eps), 
                          gp_mean = 0.0)

# Design matrix, evaluate model at design points for GP regression
# TODO: Try using more extreme u values
N <- 10
X <- matrix(seq(qnorm(.01, u.mean, u.sigma), qnorm(.99, u.mean, u.sigma), length = N), ncol = 1)
N <- nrow(X)
y.model <- f(X)

# Test points at which to evaluate GP predictions to investigate likelihood approximation
N.pred <- 1000
u.pred <- seq(qnorm(.01, u.mean, u.sigma), qnorm(.99, u.mean, u.sigma), length = N.pred)
X.pred <- matrix(u.pred, ncol=1)

# Directories
base.dir <- '.'
# base.out.dir <- file.path(base.dir, '..', 'test_output')
base.out.dir <- file.path(base.dir, 'output')

# Create sub-directory in output directory
tag <- 'pres_3'
subdir.name <- paste0('param_cal_1d_test_', tag)
out.dir <- file.path(base.out.dir, subdir.name)
dir.create(out.dir)


# -----------------------------------------------------------------------------
# Simulate Data, fixing u and tau at their prior means
# -----------------------------------------------------------------------------
u <- u.mean
tau <- tau.shape / tau.rate
y <- rnorm(n, f(u), 1/sqrt(tau))
save.gaussian.llik.plot(y, X.pred, out.dir, "exact_llik.png", f)


# -----------------------------------------------------------------------------
# Brute Force parameter calibration: 
#   Try to recover parameters by running MCMC in Stan
# -----------------------------------------------------------------------------

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
fit.brute.force <- sampling(model, stan.list, iter = 50000, chains = 4, seed = seed)
summary(fit.brute.force)
posterior.brute.force <- as.array(fit.brute.force)
samples.brute.force <- extract(fit.brute.force)

# Posterior uncertainty intervals
color_scheme_set("red")
int.prob <- 0.5
int.outer.prob <- 0.9
mcmc_intervals(posterior.brute.force, 
               pars = c("tau", "u"), 
               prob = int.prob, 
               prob_outer = int.outer.prob, 
               point_est = "median") + 
  ggtitle("Posterior Intervals: Brute force param calibration") + 
  ylab("Parameter") + 
  xlab(paste0("Inner ", 100*int.prob, "%; Outer ", 100*int.outer.prob, "%"))
ggsave(file.path(out.dir, "intervals.brute.force.png"), bg = "white")

# Posterior histogram
mcmc_hist(posterior.brute.force, pars = c("tau", "u")) + 
  ggtitle("Posterior Histogram: Brute force param calibration")
ggsave(file.path(out.dir, "hist.brute.force.png"), bg = "white")

# Posterior kernel density estimates
mcmc_dens(posterior.brute.force, pars = c("tau", "u")) + 
  ggtitle("Posterior Kernel Density Estimates: Brute force param calibration")
ggsave(file.path(out.dir, "dens.brute.force.png"), bg = "white")

# Trace plots
color_scheme_set("mix-blue-red")
mcmc_trace(posterior.brute.force, pars = c("tau", "u"),
           facet_args = list(ncol = 1, strip.position = "left")) + 
  ggtitle("Trace Plots: Brute force param calibration")
ggsave(file.path(out.dir, "trace.brute.force.png"), bg = "white")


# -----------------------------------------------------------------------------
# Fit GP regression
# -----------------------------------------------------------------------------

# Define the sufficient statistic (SS). For simplicity here, I'm not considering any other explanatory variables
# x that are being conditioned on. If there were, we would have to make sure we were lining up the observed data
# and the process model data correctly conditional on x. So in this case, y_model[1] is just a constant predictor
# of the output variable given the first calibration parameter (design point) X[1,1]. 
SS <- rep(0, N)
for(i in seq(1, N)) {
  SS[i] <- sum((y.model[i] - y)^2)
}

# Similarly calculate true SS at test points for reference in subsequent plots
SS.pred <- matrix(NA, nrow = N.pred, ncol = 1)
for(i in seq(1, N.pred)) {
  SS.pred[i, 1] <- sum((f(X.pred[i,1]) - y)^2)
}

# Scale SS to unit interval
# SS.max <- max(SS)
# SS.min <- min(SS)
# SS <- 1 - (SS.max - SS) / (SS.max - SS.min)
# SS <- matrix(SS, ncol = 1)
# plot(as.vector(X), as.vector(SS), xlab = 'Design Point', ylab = 'Sufficient Statistic')

# Fit GP
gp.fit <- mlegp(X, SS, nugget.known = 0, constantMean = 1)

# -----------------------------------------------------------------------------
# Obtain GP predictions at test points for visualizing likelihood approximation
# -----------------------------------------------------------------------------

if(use.default.gp.params) {
  stan.params.gp <- gp.param.defaults
} else {
  stan.params.gp <- create.gp.params.list(gp.fit, "mlegp")
}

# Compile Stan code
stan.llik.approx.path <- file.path(base.dir, 'likelihood_approx_comparison.stan')
stan.llik.approx.model <- stan_model(stan.llik.approx.path)

# Run Stan code
stan.params.other.llik.approx <- list(N = N, n = n, N_pred = N.pred, tau = tau, X = X, y = SS, u_vals = X.pred)
stan.params.llik.approx <- as.list(c(stan.params.gp, stan.params.other.llik.approx))
stan.llik.approx.fit <- sampling(stan.llik.approx.model, data = stan.params.llik.approx, warmup = 0, 
                                 iter = 1, chains = 1, seed = seed, refresh = 4000, algorithm = "Fixed_param")

# Save Stan results
stan.llik.approx.output <- extract(stan.llik.approx.fit)
llik.gp.uq <- as.vector(stan.llik.approx.output$u_vals_llik) # GP integrated out likelihood in Stan (uq = uncertainty quantification)
gp.means.stan <- as.vector(stan.llik.approx.output$mean_test)
gp.se.stan <- sqrt(as.vector(stan.llik.approx.output$var_test))

# Save plots visualizing GP approximation
interval.pct <- .95
save.gp.pred.mean.plot(interval.pct, tau, y, X, X.pred, SS, SS.pred, 
                       gp.means.stan, gp.se.stan, llik.gp.uq, out.dir, "gp_pred_mean", f)

# -----------------------------------------------------------------------------
# Parameter Calibration with Gaussian Process Approximation,  
# propagating uncertainty by integrating out GP.
# -----------------------------------------------------------------------------

# Compile Stan code
stan.gp.model.path <- file.path(base.dir, 'parameter_calibration_gp_1d_example.stan')
model.gp <- stan_model(stan.gp.model.path)

# Data to pass to Stan
stan.list.other <- list(N = N, 
                        n = n,
                        y = as.vector(SS),
                        X = X, 
                        a = tau.shape, 
                        b = tau.rate, 
                        u_mean = u.mean, 
                        u_sigma = u.sigma, 
                        u_lower = min(X), 
                        u_upper = max(X))
stan.list.gp <- as.list(c(stan.params.gp, stan.list.other))

# MCMC
fit.gp <- sampling(model.gp, stan.list.gp, iter = 50000, chains = 4, seed = seed)
summary(fit.gp)
samples.gp <- extract(fit.gp)
posterior.gp<- as.array(fit.gp)

# Posterior uncertainty intervals
color_scheme_set("red")
mcmc_intervals(posterior.gp, 
               pars = c("tau", "u[1]"), 
               prob = int.prob, 
               prob_outer = int.outer.prob, 
               point_est = "median") + 
  ggtitle("Posterior Intervals: GP approx param calibration") + 
  ylab("Parameter") + 
  xlab(paste0("Inner ", 100*int.prob, "%; Outer ", 100*int.outer.prob, "%"))
ggsave(file.path(out.dir, "intervals.gp.png"), bg = "white")

# Posterior histogram
mcmc_hist(posterior.gp, pars = c("tau", "u[1]")) + 
  ggtitle("Posterior Histogram: GP approx param calibration")
ggsave(file.path(out.dir, "hist.gp.png"), bg = "white")

# Posterior kernel density estimates
mcmc_dens(posterior.gp, pars = c("tau", "u[1]")) + 
  ggtitle("Posterior Kernel Density Estimates: GP approx param calibration")
ggsave(file.path(out.dir, "dens.gp.png"), bg = "white")

# Trace plots
color_scheme_set("mix-blue-red")
mcmc_trace(posterior.gp, pars = c("tau", "u[1]"),
           facet_args = list(ncol = 1, strip.position = "left")) + 
  ggtitle("Trace Plots: GP approx param calibration")
ggsave(file.path(out.dir, "trace.gp.png"), bg = "white")


# -----------------------------------------------------------------------------
# Parameter Calibration with Gaussian Process Approximation, 
# evaluating at predictive mean function.
# -----------------------------------------------------------------------------

# Compile Stan code
stan.gp.mean.model.path <- file.path(base.dir, 'parameter_calibration_gp_mean_1d_example.stan')
model.gp.mean <- stan_model(stan.gp.mean.model.path)

# MCMC
fit.gp.mean <- sampling(model.gp.mean, stan.list.gp, iter = 50000, chains = 4, seed = seed)
summary(fit.gp.mean)
posterior.gp.mean <- as.array(fit.gp.mean)

# Posterior uncertainty intervals
color_scheme_set("red")
mcmc_intervals(posterior.gp.mean, 
               pars = c("tau", "u[1]"), 
               prob = int.prob, 
               prob_outer = int.outer.prob, 
               point_est = "median") + 
  ggtitle("Posterior Intervals: GP mean approx param calibration") + 
  ylab("Parameter") + 
  xlab(paste0("Inner ", 100*int.prob, "%; Outer ", 100*int.outer.prob, "%"))
ggsave(file.path(out.dir, "intervals.gp.mean.png"), bg = "white")

# Posterior histogram
mcmc_hist(posterior.gp.mean, pars = c("tau", "u[1]")) + 
  ggtitle("Posterior Histogram: GP mean approx param calibration")
ggsave(file.path(out.dir, "hist.gp.mean.png"), bg = "white")

# Posterior kernel density estimates
mcmc_dens(posterior.gp.mean, pars = c("tau", "u[1]")) + 
  ggtitle("Posterior Kernel Density Estimates: GP mean approx param calibration")
ggsave(file.path(out.dir, "dens.gp.mean.png"), bg = "white")

# Trace plots
color_scheme_set("mix-blue-red")
mcmc_trace(posterior.gp.mean, pars = c("tau", "u[1]"),
           facet_args = list(ncol = 1, strip.position = "left")) + 
  ggtitle("Trace Plots: GP mean approx param calibration")
ggsave(file.path(out.dir, "trace.gp.mean.png"), bg = "white")




# -----------------------------------------------------------------------------
# Print information to file
# -----------------------------------------------------------------------------

file.con <- file(file.path(out.dir, "run_info.txt"))
file.text <- c(paste0("Random seed: ", seed), paste0("Model: f(u) = ", f.string), paste0("True u value: ", u),  
               paste0("tau: ", tau), paste0("n: ", n), paste0("Design points: ", paste0(as.vector(X), collapse = ", ")),
               paste0("tau.shape: ", tau.shape), paste0("tau.rate: ", tau.rate), paste0("u.mean: ", u.mean), 
               paste0("u.sigma: ", u.sigma), 
               paste0("use.default.gp.params: ", use.default.gp.params),
               "GP params:", paste0(names(stan.params.gp), collapse = ", "), 
               paste0(stan.params.gp, collapse = ", "))
writeLines(file.text, file.con)
close(file.con)












