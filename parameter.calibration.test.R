# parameter.calibration.test.R
# Main interface for running test code comparing various algorithms for parameter 
# calibration.
#
# Andrew Roberts

# TODO: 
#   - Add functionality to time runs

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

# PEcAn MCMC settings:
# - n.itr: 350000
# - adapt frequency: 250
# - Min scale factor: 0.1
# - AR target: 0.3


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
  mcmc.pecan = TRUE,
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
  u.rng = c(-Inf, Inf),
  
  # Gaussian Process: used for algorithms gp.stan, gp.mean.stan, and pecan)
  X = NULL, # Manually input design matrix; will override below settings
  N = 5, 
  gp.library = "mlegp", 
  gp.plot.interval.pct = .95, 
  
  # Brute Force algorithm settings
  n.itr.mcmc.brute.force = 50000,
  
  # pecan algorithm settings
  n.itr.mcmc.pecan = 10000,
  SS.joint.sample = FALSE, 
  resample.tau = TRUE, 
  adapt.frequency = 250, 
  adapt.min.scale = 0.1, 
  accept.rate.target = 0.3
  
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
# Gaussian Process Regression
# -----------------------------------------------------------------------------

if(any(as.logical(settings[c("mcmc.gp.stan", "mcmc.gp.mean.stan", "mcmc.pecan")]))) {
  # Design points
  X <- settings$X
  N <- settings$N
  y.model <- apply(X, 1, f)
  
  # Sufficient statistic that will be emulated
  SS <- rep(0, N)
  for(i in seq(1, N)) {
    SS[i] <- sum((y.model[i] - y.obs)^2)
  }
  
  # Similarly calculate true SS at test points for reference in subsequent plots
  SS.pred <- rep(0, settings$N.pred)
  for(i in seq(1, settings$N.pred)) {
    SS.pred[i] <- sum((f(X.pred[i]) - y.obs)^2)
  }
  
  # Fit GP regression
  if(settings$gp.library == "mlegp") {
    gp.fit <- mlegp(X, SS, nugget.known = 0, constantMean = 1)
  }
  
  # Map kernel parameters to Stan parameterization
  gp.stan.params <- create.gp.params.list(gp.fit, settings$gp.library)
  gp.obj <- create.gp.obj(gp.fit, settings$gp.library, X, SS)
  saveRDS(gp.obj, file = file.path(run.dir, "gp.obj.RData"))
  
  # Save plots demonstrating GP fit
  if(settings$k == 1) {
    save.gp.pred.mean.plot(gp.obj, X, X.pred, SS, SS.pred, f, settings$gp.plot.interval.pct, run.dir)
  }
  
}


# -----------------------------------------------------------------------------
# Brute Force parameter calibration: 
#   HMC using exact likelihood (no GP approximation)
# -----------------------------------------------------------------------------

if(settings$mcmc.brute.force.stan) {
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
  save.posterior.plots(as.array(fit.brute.force), pars, run.dir, settings$interval.prob, settings$interval.prob.outer,
                       settings$interval.point.est, ".brute.force", "Brute Force")
}


# -----------------------------------------------------------------------------
# PEcAn parameter calibration: 
#   Adaptive Metropolis-within-Gibbs with re-sampling of sufficient statistic
# -----------------------------------------------------------------------------

if(settings$mcmc.pecan) {
  
  # For each chain, randomly choose one of the design points as the initial parameter for MCMC. 
  u.init.indices <- sample(seq_len(settings$N), settings$n.mcmc.chains)
  
  # Set initial proposal variances and initial parameter values
  proposal.vars <- vector(mode = "list", length = settings$n.mcmc.chains)
  u.init <- vector(mode = "list", length = settings$n.mcmc.chains)
  for(i in seq_len(settings$n.mcmc.chains)) {
    proposal.vars[[i]] <- 0.1 * diff(qnorm(c(0.05, 0.95), settings$u.gaussian.mean, settings$u.gaussian.sd))
    u.init[[i]] <- X[u.init.indices[i],,drop = FALSE]
  }
  
  mcmc.pecan.results <- mcmc.GP.test(gp.obj = gp.obj, 
                                     n = settings$n, 
                                     n.itr = settings$n.itr.mcmc.pecan, 
                                     u.rng = settings$u.rng, 
                                     SS.joint.sample = settings$SS.joint.sample, 
                                     resample.tau = settings$resample.tau, 
                                     tau.gamma.shape = settings$tau.gamma.shape, 
                                     tau.gamma.rate = settings$tau.gamma.rate, 
                                     u.prior.mean = settings$u.gaussian.mean, 
                                     u.prior.sd = settings$u.gaussian.sd, 
                                     u0 = u.init[[1]], # TODO: This should be a vector corresponding to different chains
                                     proposal.vars = proposal.vars[[1]],
                                     adapt.frequency = settings$adapt.frequency, 
                                     adapt.min.scale = settings$adapt.min.scale, 
                                     accept.rate.target = settings$accept.rate.target)
  
}











