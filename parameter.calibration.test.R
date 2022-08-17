# parameter.calibration.test.R
# Main interface for running test code comparing various algorithms for parameter 
# calibration.
#
# Andrew Roberts

library(bayesplot)
library(ggplot2)
library(rstan)
library(mlegp)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("helper.functions.R")
source("stan.helper.functions.R")
source("gaussian.process.functions.R")
source("mcmc.GP.test.R")


# -----------------------------------------------------------------------------
# Settings
#   - All paths should be relative to base.dir$output.dir. Paths will be 
#     automatically expanded to their full paths in preprocess.settings()
#   - Some fields can be set to NULL and will be automatically filled by 
#     preprocess.settings()
# -----------------------------------------------------------------------------

settings <- list(
  
  # General 
  seed = 5, 
  k = 1, # Dimension of u
  base.dir = getwd(),
  output.dir = "output",
  run.id = Sys.time(), 
  run.description = "Test",
  
  # General MCMC
  n.mcmc.chains = 4, 
  interval.prob = 0.5, 
  interval.prob.outer = 0.9,
  interval.point.est = "median",
  
  # Algorithms
  mcmc.brute.force.stan = TRUE, 
  mcmc.gp.stan = FALSE, 
  mcmc.gp.mean.stan = FALSE, 
  mcmc.pecan = FALSE,
  mcmc.brute.force.stan.path = "parameter_calibration_1d_example.stan",
  mcmc.gp.stan.path = "parameter_calibration_gp_1d_example.stan", 
  mcmc.gp.mean.path = "parameter_calibration_gp_mean_1d_example.stan",
  
  # Likelihood (used to generate synthetic dataset)
  n = 1000,
  N.pred = 1000, # Used for producing plots and as GP prediction test points
  lik.type = "gaussian",
  model = function(u) {u},
  u.true = 0.5,
  sigma.true = 0.3, 
  tau.true = NULL, 
  
  # Priors
  tau.gamma.shape = NULL,
  tau.gamma.rate = 1.0,
  u.gaussian.mean = NULL,
  u.gaussian.sd = 0.25,
  
  # Gaussian Process: used for algorithms gp.stan, gp.mean.stan, and pecan)
  X = NULL, # Manually input design matrix; will override below settings
  N = 5, 
  gp.library = "mlegp", 
  
  # Brute Force algorithm settings
  n.itr.mcmc.brute.force = 50000
  
)

settings <- preprocess.settings(settings)


# -----------------------------------------------------------------------------
# Create output directory for current run
# -----------------------------------------------------------------------------

run.dir <- file.path(settings$output.dir, paste0("param_cal_", settings$run.id))
dir.create(run.dir)
saveRDS(settings, file = file.path(run.dir, "settings.RData"))


# -----------------------------------------------------------------------------
# Generate Synthetic Dataset
# -----------------------------------------------------------------------------

set.seed(settings$seed)
pars <- c("tau", "u")

# Generate observed data
f <- settings$model
y.obs <- rnorm(settings$n, f(settings$u.true), settings$sigma.true)

# Save plot of true likelihood
u.pred <- seq(qnorm(.01, settings$u.true, settings$sigma.true), 
              qnorm(.99, settings$u.true, settings$sigma.true), length = settings$N.pred)
X.pred <- matrix(u.pred, ncol=1)
save.gaussian.llik.plot(y.obs, X.pred, run.dir, "exact_llik.png", f, settings$tau.true)


# -----------------------------------------------------------------------------
# Fit Gaussian Process Regression
# -----------------------------------------------------------------------------

if(any(as.logical(settings[c("mcmc.gp.stan", "mcmc.gp.mean.stan", "mcmc.pecan")]))) {
  X <- settings$X
  N <- settings$N
  y.model <- apply(X, 1, f)
}


# -----------------------------------------------------------------------------
# Brute Force parameter calibration: 
#   Try to recover parameters by running MCMC in Stan
# -----------------------------------------------------------------------------

# Compile Stan code
model.brute.force <- stan_model(settings$mcmc.brute.force.stan.path)

# Data to pass to Stan
stan.list.brute.force <- list(n = settings$n, 
                              y = y.obs, 
                              a = settings$tau.gamma.shape, 
                              b = settings$tau.gamma.rate, 
                              u_mean = settings$u.gaussian.mean, 
                              u_sigma = settings$u.gaussian.sd)

# MCMC
fit.brute.force <- sampling(model.brute.force, stan.list.brute.force, iter = settings$n.itr.mcmc.brute.force, 
                            chains = settings$n.mcmc.chains, seed = settings$seed)
saveRDS(summary(fit.brute.force), file = file.path(run.dir, "summary.brute.force.RData"))

save.posterior.intervals.plot(as.array(fit.brute.force), pars, 
                              file.path(run.dir, "post.int.brute.force.png"), 
                              int.prob = interval.prob, int.outer.prob = interval.prob.outer,
                              point.est = interval.point.est, append.title = "Brute force parameter calibration")




