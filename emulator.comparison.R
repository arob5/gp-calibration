# emulator.comparison.R
# Main interface for running test code comparing various GP emulators for modeling the 
# sum of squared errors (or log of this quantity); i.e. the sufficient statistic 
# of the Gaussian likelihood. 
#
# Andrew Roberts

library(ggplot2)
library(mlegp)
library(hetGP)
library(lhs)

source("helper.functions.R")
source("stan.helper.functions.R")
source("gaussian.process.functions.R")
source("mcmc.GP.test.R")

# TODO:
# Likelihood plots (lik vs pred lik) with horizontal error bars on pred like mimicking plot() method from hetGP


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
  k = 2, # Dimension of input space
  base.dir = getwd(),
  output.dir = "output",
  run.id = paste0("", Sys.time()), 
  run.description = "",
  
  # Likelihood (used to generate synthetic dataset)
  n = 1000,
  lik.type = "gaussian",
  model = function(u) {u[1] + u[2]},
  u.true = list(0.5, 1.0),
  sigma.true = 0.3, 
  tau.true = NULL, 
  
  # Priors
  tau.gamma.shape = NULL,
  tau.gamma.rate = 1.0,
  u.prior.coef.var = list(0.5, 1.0),
  u.gaussian.mean = NULL,
  u.gaussian.sd = NULL,
  u.rng = NULL,
  
  # Gaussian Process Regression 
  X = NULL, # Manually input design matrix; will override below settings
  N = 10, 
  N.pred = 1000,
  gp.library = c("hetGP", "mlegp"), 
  log.normal.process = TRUE,
  gp.plot.interval.pct = .95,
  normalize.y = FALSE,
  scale.X = FALSE,
  num.itr.cv = 1
)

settings <- preprocess.settings(settings)


# -----------------------------------------------------------------------------
# Create output directory for current run
# -----------------------------------------------------------------------------

run.dir <- file.path(settings$output.dir, paste0("emulator_comparison_", settings$run.id))
dir.create(run.dir)
saveRDS(settings, file = file.path(run.dir, "settings.RData"))


# -----------------------------------------------------------------------------
# Generate Synthetic Dataset
# -----------------------------------------------------------------------------

set.seed(settings$seed)

# Generate observed data
f <- settings$model
y.obs <- rnorm(settings$n, f(settings$u.true), settings$sigma.true)

# Save plot of true likelihood
if(settings$k == 1) {
  save.gaussian.llik.plot(y.obs, settings$X.pred, run.dir, "exact_llik.png", f, settings$tau.true)
}


# -----------------------------------------------------------------------------
# Main Cross Validation Loop
# -----------------------------------------------------------------------------

# Store cross validation output to be used in computation of metrics, etc.
cv.obj <- vector(mode = "list", length = settings$num.itr.cv)

for(cv in seq_len(settings$num.itr.cv)) {

  # Generate training and test sets
  X.list <- LHS.train.test(settings$N, settings$N.pred, settings$k, settings$u.gaussian.mean, settings$u.gaussian.sd)
  X <- X.list$X
  X.test <- X.list$X.test

  # Calculate sufficient statistic at training and test locations
  SS <- calc.SS(X, f, y.obs)
  SS.test <- calc.SS(X.test, f, y.obs)

  # Calculate true likelihood at training and test locations
  llik <- dmvnorm.log.SS(SS, settings$tau.true, settings$n)
  llik.test <- dmvnorm.log.SS(SS.test, settings$tau.true, settings$n)
  
  # Add info to CV object
  cv.obj[[cv]][c("X", "X.test", "SS", "SS.test", "llik", "llik.test")] <- list(X, X.test, SS, SS.test, llik, llik.test)
  
  # Fits GPs
  gp.fits <- fit.GPs(settings$gp.library, X, SS, log.SS = settings$log.normal.process)
  
  # Map kernel parameters to Stan parameterization
  gp.obj.list <- lapply(seq_along(gp.fits), 
                        function(i) create.gp.obj(gp.fits[[i]], settings$gp.library[i], X, SS, log.y = settings$log.normal.process))
  
  # Predict at test locations
  gp.pred.list <- lapply(seq_along(gp.fits), 
                         function(i) predict_gp(X.test, gp.obj.list[[i]], pred.cov = FALSE,
                                                lognormal.adjustment = settings$log.normal.process))
                                                                    
  # Prediction metrics
  gp.rmse.vec <- sapply(seq_along(gp.fits), function(i) sqrt(sum((gp.pred.list[[i]] - SS.test)^2))) / settings$N.pred
  gp.mae.vec <- sapply(seq_along(gp.fits), function(i) sum(abs(gp.pred.list[[i]] - SS.test))) / settings$N.pred
  
  # Add GP fit and prediction data to CV object
  # TODO: convert log-GP mean/var to original space?
  cv.obj[[cv]][["gp.list"]] <- vector(mode = "list", length = length(gp.fits))
  names(cv.obj[[cv]][["gp.list"]]) <- settings$gp.library
  for(gp in seq_along(gp.fits)) {
    gp.lib <- names(cv.obj[[cv]][["gp.list"]])[[gp]]
    cv.obj[[cv]][[gp.lib]][["gp.obj"]] <- gp.fits[[gp]]
    cv.obj[[cv]][[gp.lib]][["pred.mean"]] <- gp.pred.list[[gp]]$mean
    cv.obj[[cv]][[gp.lib]][["pred.var"]] <- gp.pred.list[[gp]]$var
    
  }
  
}
  
  









  
