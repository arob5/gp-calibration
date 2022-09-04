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
# Fit GPs
# -----------------------------------------------------------------------------

# Design points
X <- settings$X
N <- settings$N
y.model <- apply(X, 1, f)
  
# Sufficient statistic that will be emulated
SS <- rep(0, N)
for(i in seq(1, N)) {
  SS[i] <- sum((y.model[i] - y.obs)^2)
}
  
# Similarly calculate true SS at test points
SS.pred <- rep(0, settings$N.pred)
for(i in seq(1, settings$N.pred)) {
  SS.pred[i] <- sum((f(settings$X.pred[i]) - y.obs)^2)
}

# If modeling log(SS)
if(settings$log.normal.process) {
  SS <- log(SS)
  SS.pred <- log(SS.pred)
}
  
# Fit GP regression
gp.fits <- vector(mode = "list", length = length(settings$gp.library))
gp.fits.index <- 1
if("mlegp" %in% settings$gp.library) {
  gp.fits[[gp.fits.index]] <- mlegp(X, SS, nugget.known = 0, constantMean = 1)
  gp.fits.index <- gp.fits.index + 1
}

if("hetGP" %in% settings$gp.library) {
  gp.fits[[gp.fits.index]] <- mleHomGP(X, SS, covtype = "Gaussian", known = list(g = .Machine$double.eps))
  gp.fits.index <- gp.fits.index + 1
}
  
if(gp.fits.index != length(gp.fits) + 1) {
  stop("Invalid GP library detected.")
}

# Map kernel parameters to Stan parameterization
gp.stan.params.list <- lapply(seq_along(gp.fits), function(i) create.gp.params.list(gp.fits[[i]], settings$gp.library[i]))
gp.obj.list <- lapply(seq_along(gp.fits), function(i) create.gp.obj(gp.fits[[i]], settings$gp.library[i], X, SS))
for(i in seq_along(gp.obj.list)) {
  saveRDS(gp.obj.list[[i]], file = file.path(run.dir, paste0("gp.obj.", settings$gp.library[i], ".RData")))
}

# Save plots demonstrating GP fit
# if(settings$k == 1) {
#   save.gp.pred.mean.plot(gp.obj, X, settings$X.pred, SS, SS.pred, f, 
#                          settings$gp.plot.interval.pct, file.path(run.dir, "gp_pred_SS.png"))
# }
  

# -----------------------------------------------------------------------------
# Cross Validation
# -----------------------------------------------------------------------------

X.list <- vector(mode = "list", length = settings$num.itr.cv)
X.pred.list <- vector(mode = "list", length = settings$num.itr.cv)

for(cv in seq_len(settings$num.itr.cv)) {

  # Generate train and test sets
  X.combined <- randomLHS(settings$N + settings$N.pred, settings$k)
  X.list[[cv]] <- X.combined[1:settings$N,,drop=FALSE]
  X.pred.list[[cv]] <- X.combined[-(1:settings$N),,drop=FALSE]

  # Apply inverse CDS transform using prior distributions.
  for(j in seq_len(settings$k)) {
    X.list[[cv]][,j] <- qnorm(X.list[[cv]][,j], settings$u.gaussian.mean[[j]], settings$u.gaussian.sd[[j]])
    X.pred.list[[cv]][,j] <- qnorm(X.pred.list[[cv]][,j], settings$u.gaussian.mean[[j]], settings$u.gaussian.sd[[j]])
  }
  
  # Run full model at training and test locations
  y.model <- apply(X.list[[cv]], 1, f)
  y.test <- apply(X.pred.list[[cv]], 1, f)
  
  # Calculate sufficient statistic at training and test locations
  SS <- rep(0, N)
  SS.test <- rep(0, settings$N.pred)
  
  for(i in seq(1, N)) {
    SS[i] <- sum((y.model[i] - y.obs)^2)
  }
  for(i in seq(1, settings$N.pred)) {
    SS.test[i] <- sum((y.test[i] - y.obs)^2)
  }
  
  # If modeling log(SS)
  if(settings$log.normal.process) {
    SS <- log(SS)
    SS.test <- log(SS.test)
  }
  
  # Fit GP regressions
  gp.fits <- vector(mode = "list", length = length(settings$gp.library))
  gp.fits.index <- 1
  if("mlegp" %in% settings$gp.library) {
    gp.fits[[gp.fits.index]] <- mlegp(X.list[[cv]], SS, nugget.known = 0, constantMean = 1)
    gp.fits.index <- gp.fits.index + 1
  }
  
  if("hetGP" %in% settings$gp.library) {
    gp.fits[[gp.fits.index]] <- mleHomGP(X.list[[cv]], SS, covtype = "Gaussian", known = list(g = .Machine$double.eps))
    gp.fits.index <- gp.fits.index + 1
  }
  
  if(gp.fits.index != length(gp.fits) + 1) {
    stop("Invalid GP library detected.")
  }
  
  # Map kernel parameters to Stan parameterization
  gp.stan.params.list <- lapply(seq_along(gp.fits), function(i) create.gp.params.list(gp.fits[[i]], settings$gp.library[i]))
  gp.obj.list <- lapply(seq_along(gp.fits), function(i) create.gp.obj(gp.fits[[i]], settings$gp.library[i], X.list[[cv]], SS))
  
  # Predict at test locations
  gp.pred.list <- lapply(seq_along(gp.fits), function(i) predict_gp(X.pred.list[[cv]], gp.obj.list[[i]], pred.cov = FALSE))
  if(settings$log.normal.process) {
    gp.pred.list <- lapply(gp.pred.list, exp)
  }
  
  # Prediction metrics
  gp.rmse.vec <- sapply(seq_along(gp.fits), sqrt(sum((gp.pred.list[[i]] - SS.test)^2))) / settings$N.pred
  gp.mae.vec <- sapply(seq_along(gp.fits), sum(abs(gp.pred.list[[i]] - SS.test))) / settings$N.pred
  # TODO: likelihood difference metrics
}
  
  









  
