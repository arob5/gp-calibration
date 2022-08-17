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
# -----------------------------------------------------------------------------

settings <- list(
  
  # General 
  seed = 5, 
  k = 1, # Dimension of u
  base.dir = getwd(),
  output.dir = file.path(base.dir, "output"),
  run.id = Sys.time(), 
  run.description = "",
  
  # Algorithms
  mcmc.brute.force.stan = FALSE, 
  mcmc.gp.stan = FALSE, 
  mcmc.gp.mean.stan = FALSE, 
  mcmc.pecan = FALSE,
  
  # Likelihood (used to generate synthetic dataset)
  n = 1000,
  N.pred = 1000, # Used for producing plots and as GP prediction test points
  lik.type = "gaussian",
  model = function(u) {u},
  u.true = 0.5,
  sigma.true = 0.3, 
  tau.true = 1 / sigma.true^2, 
  
  # Priors
  tau.gamma.shape = tau.true,
  tau.gamma.rate = 1.0,
  u.guassian.mean = u.true,
  u.gaussian.sd = 0.25,
  
  # Gaussian Process: used for algorithms gp.stan, gp.mean.stan, and pecan)
  X = NULL, # Manually input design matrix; will override below settings
  N = 5, 
  library = "mlegp", 
  
  # pecan algorithm settings
  
)


# -----------------------------------------------------------------------------
# Create output directory for current run
# -----------------------------------------------------------------------------

run.dir <- paste0("param_cal_", settings$run.id)
dir.create(file.path(settings$output.dir, run.dir))
save(settings, file = file.path(run.dir, "settings.RData"))


# -----------------------------------------------------------------------------
# Generate Synthetic Dataset
# -----------------------------------------------------------------------------

set.seed(settings$seed)

# Generate observed data
f <- settings$model
y.obs <- rnorm(settings$n, f(settings$u.true), settings$sigma.true)

# Save plot of true likelihood
u.pred <- seq(qnorm(.01, settings$u.true, settings$sigma.true), 
              qnorm(.99, settings$u.true, settings$sigma.true), length = settings$N.pred)
X.pred <- matrix(u.pred, ncol=1)
save.gaussian.llik.plot(y.obs, X.pred, run.dir, "exact_llik.png", f, settings$tau.true)











