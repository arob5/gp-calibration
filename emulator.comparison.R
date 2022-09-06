# emulator.comparison.R
# Main interface for running test code comparing various GP emulators for modeling the 
# sum of squared errors (or log of this quantity); i.e. the sufficient statistic 
# of the Gaussian likelihood. 
#
# Andrew Roberts

library(data.table)
library(colorspace)
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
  k = 1, # Dimension of input space
  base.dir = getwd(),
  output.dir = "output",
  run.id = paste0("test_1d_", Sys.time()), 
  run.description = "1d test.",
  
  # Likelihood (used to generate synthetic dataset)
  n = 1000,
  lik.type = "gaussian",
  model = function(u) {u},
  u.true = list(0.5),
  sigma.true = 0.3, 
  tau.true = NULL, 
  
  # Priors
  tau.gamma.shape = NULL,
  tau.gamma.rate = 1.0,
  u.prior.coef.var = 1.0,
  u.gaussian.mean = NULL,
  u.gaussian.sd = NULL,
  u.rng = NULL,
  
  # Gaussian Process Regression 
  X = NULL, # Manually input design matrix; will override below settings
  N = 5, 
  N.pred = 1000,
  gp.library = c("hetGP", "mlegp"), 
  log.normal.process = TRUE,
  gp.plot.interval.pct = .95,
  normalize.y = FALSE,
  scale.X = FALSE,
  num.itr.cv = 5
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
  u.plot.rng <- qnorm(c(.01, .99), settings$u.gaussian.mean, settings$u.gaussian.sd)
  save.gaussian.llik.plot(y.obs, matrix(seq(u.plot.rng[1], u.plot.rng[2], length.out = 1000), ncol = 1), 
                          run.dir, "exact_llik.png", f, settings$tau.true)
}


# -----------------------------------------------------------------------------
# Main Cross Validation Loop
# -----------------------------------------------------------------------------

# Store cross validation output to be used in computation of metrics, etc.
cv.obj.names <- c(c("X", "X.test", "SS", "SS.test", "llik", "llik.test"), settings$gp.library)
gp.list.names <- c("gp.obj", "pred.mean", "pred.var", "rmse", "mae", "llik.pred.test", "lik.l1.diff")
cv.obj <- vector(mode = "list", length = length(cv.obj.names))
for(gp.lib in settings$gp.library) {
  cv.obj[[gp.lib]] <- vector(mode = "list", length = length(gp.list.names))
  names(cv.obj[[gp.lib]]) <- gp.list.names
}
cv.obj[["num.itr.cv"]] <- settings$num.itr.cv

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
  cv.obj[["X"]][[cv]] <- X
  cv.obj[["X.test"]][[cv]] <- X.test
  cv.obj[["SS"]][[cv]] <- SS
  cv.obj[["SS.test"]][[cv]] <- SS.test
  cv.obj[["llik"]][[cv]] <- llik
  cv.obj[["llik.test"]][[cv]] <- llik.test
  
  # Fits GPs
  gp.fits <- fit.GPs(settings$gp.library, X, SS, log.SS = settings$log.normal.process)
  
  # Map kernel parameters to Stan parameterization
  gp.obj.list <- lapply(seq_along(gp.fits), 
                        function(i) create.gp.obj(gp.fits[[i]], settings$gp.library[i], X, SS, log.y = settings$log.normal.process))
  
  # Predict at test locations
  gp.pred.list <- lapply(seq_along(gp.fits), 
                         function(i) predict_gp(X.test, gp.obj.list[[i]], pred.cov = FALSE,
                                                lognormal.adjustment = settings$log.normal.process))
  
  # Evaludating likelihood at predictions
  llik.pred.test <- lapply(seq_along(gp.fits), 
                           function(i) dmvnorm.log.SS(gp.pred.list[[i]]$mean, settings$tau.true, settings$n))
                                                                    
  # Prediction metrics
  gp.rmse.vec <- sapply(seq_along(gp.fits), function(i) sqrt(sum((gp.pred.list[[i]]$mean - SS.test)^2))) / settings$N.pred
  gp.mae.vec <- sapply(seq_along(gp.fits), function(i) sum(abs(gp.pred.list[[i]]$mean - SS.test))) / settings$N.pred
  gp.lik.l1.diff <- sapply(seq_along(gp.fits), function(i) sum(abs(exp(llik.pred.test[[i]]) - exp(llik.test)))) / settings$N.pred
  
  # Add GP fit and prediction data to CV object
  for(gp in seq_along(settings$gp.library)) {
    gp.lib <- settings$gp.library[[gp]]
    cv.obj[[gp.lib]][["gp.obj"]][[cv]] <- gp.obj.list[[gp]]
    cv.obj[[gp.lib]][["pred.mean"]][[cv]] <- gp.pred.list[[gp]]$mean
    cv.obj[[gp.lib]][["pred.var"]][[cv]] <- gp.pred.list[[gp]]$var
    cv.obj[[gp.lib]][["rmse"]][[cv]] <- gp.rmse.vec[gp]
    cv.obj[[gp.lib]][["mae"]][[cv]] <- gp.mae.vec[gp]
    cv.obj[[gp.lib]][["llik.pred.test"]][[cv]] <- llik.pred.test[[gp]]
    cv.obj[[gp.lib]][["lik.l1.diff"]][[cv]] <- gp.lik.l1.diff[gp]
  }
  
}
  
  









  
