# parameter.calibration.test.R
# Main interface for running test code comparing various algorithms for parameter 
# calibration.
#
# Andrew Roberts

library(parallel)
library(bayesplot)
library(ggplot2)
library(rstan)
library(mlegp)
library(hetGP)
library(pracma)
library(lhs)

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
  k = 2, # Dimension of u
  base.dir = getwd(),
  output.dir = "output",
  run.id = paste0("2d test", Sys.time()), 
  run.description = "",
  
  # General MCMC
  n.mcmc.chains = 4, 
  warmup.frac = 0.5, 
  interval.prob = 0.5, 
  interval.prob.outer = 0.9,
  interval.point.est = "median",
  
  # Algorithms
  mcmc.brute.force.stan = FALSE, 
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
  model = function(u) {u[1] + u[2]},
  u.true = list(0.5, 1.0),
  # sigma.true = 0.3, 
  tau.true = 1 / 0.3^2, 
  
  # Priors
  tau.gamma.shape = NULL,
  tau.gamma.rate = 1.0,
  u.prior.coef.var = list(0.5, 1.0),
  u.gaussian.mean = NULL,
  u.gaussian.sd = NULL,
  u.rng = NULL,
  
  # Gaussian Process: used for algorithms gp.stan, gp.mean.stan, and pecan)
  X = NULL, # Manually input design matrix; will override below settings
  N = 20, 
  gp.library = "hetGP", 
  log.normal.process = TRUE,
  gp.plot.interval.pct = .95, 
  joint.sample.train.test = FALSE, 
  
  # Brute Force algorithm settings
  n.itr.mcmc.brute.force = 50000,
  
  # GP Stan algorithm settings
  n.itr.mcmc.gp.stan = 1000,
  mgf_num_eval = 1000, 
  mgf_tol = 1e-9, 
  mgf_M = 1,
  
  # pecan algorithm settings
  n.itr.mcmc.pecan = 50000,
  SS.joint.sample = TRUE, 
  resample.tau = FALSE, 
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
pars <- c(paste0("u[", seq(1, settings$k), "]"), "tau")

# Generate observed data
f <- settings$model
y.obs <- generate.observed.data(settings$n, f, settings$u.true, settings$tau.true, settings$lik.type)


# Save plot of true likelihood
if(settings$k == 1) {
  save.gaussian.llik.plot(y.obs, settings$X.pred, run.dir, "exact_llik.png", f, settings$tau.true)
}

# -----------------------------------------------------------------------------
# Gaussian Process Regression
# -----------------------------------------------------------------------------

if(any(as.logical(settings[c("mcmc.gp.stan", "mcmc.gp.mean.stan", "mcmc.pecan")]))) {
  # Design points
  X <- settings$X
  N <- settings$N
  y.model <- apply(X, 1, f)
  
  # Sufficient statistic that will be emulated
  SS <- calc.SS(X, f, y.obs)
  SS.pred <- calc.SS(settings$X.pred, f, y.obs)

  # Fit GP regression
  if(settings$gp.library == "mlegp") {
    gp.fit <- mlegp(X, SS, nugget.known = 0, constantMean = 1)
  } else if(settings$gp.library == "hetGP") {
    gp.fit <- mleHomGP(X, SS, covtype = "Gaussian", known = list(g = .Machine$double.eps))
  } else {
    stop("Invalid GP library: ", settings$gp.library)
  }
  
  # Map kernel parameters to Stan parameterization
  gp.stan.params <- create.gp.params.list(gp.fit, settings$gp.library)
  gp.obj <- create.gp.obj(gp.fit, settings$gp.library, X, SS)
  saveRDS(gp.obj, file = file.path(run.dir, "gp.obj.RData"))
  
  # Save plots demonstrating GP fit
  if(settings$k == 1) {
    save.gp.1d.pred.mean.plot(gp.obj, X, settings$X.pred, SS, SS.pred, f, 
                              settings$gp.plot.interval.pct, file.path(run.dir, "gp_pred_SS.png"))
  }
  
  save.gp.pred.mean.plot(gp.obj, X, settings$X.pred, SS, SS.pred, lognormal.adjustment = FALSE, 
                         log.log.plot = TRUE, file.path(run.dir, "gp_pred_SS_scatter.png")) 

  
  # If response variable in GP regression is log sufficient statistic
  if(settings$mcmc.gp.stan || settings$log.normal.process) {
    log.SS <- log(SS)
    log.SS.pred <- log(SS.pred)
    
    if(settings$gp.library == "mlegp") {
      gp.log.fit <- mlegp(X, log.SS, nugget.known = 0, constantMean = 1)
    } else if(settings$gp.library == "hetGP") {
      gp.log.fit <- mleHomGP(X, log.SS, covtype = "Gaussian", known = list(g = .Machine$double.eps))
    }
    
    gp.log.stan.params <- create.gp.params.list(gp.log.fit, settings$gp.library)
    gp.log.obj <- create.gp.obj(gp.log.fit, settings$gp.library, X, log.SS)
    saveRDS(gp.log.obj, file = file.path(run.dir, "gp.log.obj.RData"))
    
    if(settings$k == 1) {
      save.gp.pred.mean.plot(gp.log.obj, X, X.pred, log.SS, log.SS.pred, f, 
                             settings$gp.plot.interval.pct, file.path(run.dir, "gp_pred_log_SS.png"))
      save.gp.pred.mean.plot(gp.log.obj, X, X.pred, log.SS, log.SS.pred, f, 
                             settings$gp.plot.interval.pct, file.path(run.dir, "gp_pred_exp_log_SS.png"), exp.pred = TRUE)
      save.gp.pred.llik.plot(gp.log.obj, settings$tau.true, settings$n, X, X.pred, log.SS, log.SS.pred,
                             file.path(run.dir, "gp_pred_log_SS_llik.png"))
    }
    
    # TODO: Make sure plot dimensions/scale are the same for this plot and non-log plot
    # Is the LNP really this much of a worse model in this case or is this a bug in the code?
    save.gp.pred.mean.plot(gp.log.obj, X, settings$X.pred, SS, SS.pred, lognormal.adjustment = TRUE, 
                           log.log.plot = TRUE, file.path(run.dir, "lnp_pred_SS_scatter.png")) 

    
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
# GP Stan parameter calibration: 
#   HMC using integrated out GP log approximation
# -----------------------------------------------------------------------------

if(settings$mcmc.gp.stan) {
  pars.gp.stan <- c("u[1]", "tau")
  
  # Compile Stan code
  model.gp.stan <- stan_model(settings$mcmc.gp.stan.path)
  
  # Data to pass to Stan
  stan.list.gp.stan <- list(N = N,
                            n = settings$n, 
                            y = log.SS, 
                            X = X,
                            a = settings$tau.gamma.shape, 
                            b = settings$tau.gamma.rate, 
                            u_mean = settings$u.gaussian.mean, 
                            u_sigma = settings$u.gaussian.sd, 
                            u_lower = min(X), 
                            u_upper = max(X), 
                            gp_rho = gp.log.obj$gp_rho, 
                            gp_alpha = gp.log.obj$gp_alpha, 
                            gp_sigma = gp.log.obj$gp_sigma, 
                            gp_mean = gp.log.obj$gp_mean, 
                            mgf_num_eval = settings$mgf_num_eval, 
                            mgf_tol = settings$mgf_tol, 
                            mgf_M = settings$mgf_M)
  
  # MCMC
  fit.gp.stan <- sampling(model.gp.stan, stan.list.gp.stan, iter = settings$n.itr.mcmc.gp.stan, 
                          chains = settings$n.mcmc.chains, seed = settings$seed)
  saveRDS(summary(fit.gp.stan), file = file.path(run.dir, "summary.gp.stan.RData"))
  save.posterior.plots(as.array(fit.gp.stan), pars.gp.stan, run.dir, settings$interval.prob, settings$interval.prob.outer,
                       settings$interval.point.est, ".gp.stan", "GP Stan")
  
  
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
  
  # Modeling the sufficient statistic either as Gaussian process or log-normal process.
  if(settings$log.normal.process) {
    gp.obj.pecan <- gp.log.obj
  } else {
    gp.obj.pecan <- gp.obj
  }
  
  # Set range of u to prevent extrapolation
  u.range.pecan <- range(X)
  
  # Set up parallel computation
  cl <- makeCluster(settings$n.mcmc.chains)
  clusterExport(cl, ls())
  clusterEvalQ(cl, {
    library(bayesplot)
    library(ggplot2)
    library(rstan)
    library(mlegp)
    library(TruncatedNormal)
    library(mvtnorm)
    library(pracma)
  })
  
  mcmc.pecan.results <- parLapply(cl, seq(1, settings$n.mcmc.chains), 
                                  function(chain) mcmc.GP.test(gp.obj = gp.obj.pecan, 
                                                               n = settings$n, 
                                                               n.itr = settings$n.itr.mcmc.pecan, 
                                                               u.rng = u.range.pecan, 
                                                               SS.joint.sample = settings$SS.joint.sample, 
                                                               resample.tau = settings$resample.tau, 
                                                               tau.gamma.shape = settings$tau.gamma.shape, 
                                                               tau.gamma.rate = settings$tau.gamma.rate, 
                                                               u.prior.mean = settings$u.gaussian.mean, 
                                                               u.prior.sd = settings$u.gaussian.sd, 
                                                               u0 = u.init[[chain]],
                                                               proposal.vars = proposal.vars[[chain]],
                                                               adapt.frequency = settings$adapt.frequency, 
                                                               adapt.min.scale = settings$adapt.min.scale, 
                                                               accept.rate.target = settings$accept.rate.target, 
                                                               log.normal.process = settings$log.normal.process))

  stopCluster(cl)
  samples.pecan <- par.mcmc.results.to.arr(mcmc.pecan.results, c(pars, "SS"), settings$n.itr.mcmc.pecan, settings$warmup.frac)
  saveRDS(mcmc.summary(samples.pecan, pars), file = file.path(run.dir, "summary.pecan.RData"))
  save.posterior.plots(samples.pecan, pars, run.dir, settings$interval.prob, settings$interval.prob.outer,
                       settings$interval.point.est, ".pecan", "PEcAn")
  save.SS.tau.samples.plot(samples.pecan, file.path(run.dir, "SS_tau_samples_pecan.png"))
  
}











