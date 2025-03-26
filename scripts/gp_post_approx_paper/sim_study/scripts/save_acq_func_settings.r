#
# save_acq_func_settings.r
# Saves acquisition settings list, as used in `emulator_seq_design.r`
#
# Andrew Roberts
#

library(data.table)

experiment_tag <- "vsem"
round <- 2L

# Specify where settings will be saved.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")
experiment_dir <- file.path(base_dir, "output", "gp_inv_prob", experiment_tag)
out_dir <- file.path(experiment_dir, paste0("round", round), "design")

# Number of grid points for all integrated uncertainty criteria (grid points
# used to approximate integral over design space). This is not actually stored
# in the settings, but can be used to check that all runs use the same 
# number of grid points (if that is desired).
n_grid <- 500L

# Settings for candidate points for all acquisition criteria (the set of inputs
# that is optimized over in the discrete optimization).
candidate_settings <- list(n_prior=500L, n_mcmc=500L)

# These are additional arguments that will be passed directly to the acquisition
# functions when they are called.
acq_args <- list(adjustment="rectified")

# Create table of settings, one row per a particular acquisition function 
# setup. The flag `integrated` indicates whether or not evaluation of the 
# acquisition function requires approximating an integral over the input space.
# The flag `goal_oriented` indicates whether the acquisition function explicitly
# targets the goal of posterior approximation; practically, this typically 
# indicates whether the acquisition is simply defined wrt the underlying GP,
# or whether it takes into account the likelihood structure. The `integrated`
# and `goal_oriented` flags are not actually settings, but rather labels used
# in organizing/categorizing the acquisition criteria.
#
# The number of integration grid points is partitioned between points sampled
# from the prior and points sampled from an MCMC run, `n_prior` and `n_mcmc`,
# respectively.

l <- list()

# IVAR, pure sequential, prior weights.
l[[1]] <- list(name="llik_IVAR_grid_gp", response_heuristic=NA,
               integrated=TRUE, goal_oriented=FALSE, 
               int_grid_settings=list(n_prior=500L, n_mcmc=0L))

# IVAR, pure sequential, posterior weights.
l[[2]] <- list(name="llik_IVAR_grid_gp", response_heuristic=NA,
               integrated=TRUE, goal_oriented=TRUE, 
               int_grid_settings=list(n_prior=0L, n_mcmc=500L))

# IVAR, pure sequential, prior-posterior mixture weights.
l[[3]] <- list(name="llik_IVAR_grid_gp", response_heuristic=NA,
               integrated=TRUE, goal_oriented=TRUE, 
               int_grid_settings=list(n_prior=100L, n_mcmc=400L))

# IVAR, cl-min batch heuristic, prior weights.
l[[4]] <- list(name="llik_IVAR_grid_gp", response_heuristic="cl-min",
               integrated=TRUE, goal_oriented=FALSE, 
               int_grid_settings=list(n_prior=500L, n_mcmc=0L))

# IVAR, cl-max batch heuristic, prior weights.
l[[5]] <- list(name="llik_IVAR_grid_gp", response_heuristic="cl-max",
               integrated=TRUE, goal_oriented=FALSE, 
               int_grid_settings=list(n_prior=500L, n_mcmc=0L))

# IVAR, cl-min batch heuristic, posterior weights.
l[[6]] <- list(name="llik_IVAR_grid_gp", response_heuristic="cl-min",
               integrated=TRUE, goal_oriented=TRUE, 
               int_grid_settings=list(n_prior=0L, n_mcmc=500L))

# IVAR, cl-max batch heuristic, posterior weights.
l[[7]] <- list(name="llik_IVAR_grid_gp", response_heuristic="cl-max",
               integrated=TRUE, goal_oriented=TRUE, 
               int_grid_settings=list(n_prior=0L, n_mcmc=500L))

# IVAR, cl-min batch heuristic, prior-posterior mixture weights.
l[[8]] <- list(name="llik_IVAR_grid_gp", response_heuristic="cl-min",
               integrated=TRUE, goal_oriented=TRUE, 
               int_grid_settings=list(n_prior=100L, n_mcmc=400L))

# IVAR, cl-max batch heuristic, prior-posterior mixture weights.
l[[9]] <- list(name="llik_IVAR_grid_gp", response_heuristic="cl-max",
               integrated=TRUE, goal_oriented=TRUE, 
               int_grid_settings=list(n_prior=100L, n_mcmc=400L))

# IEVAR, pure sequential, prior weights.
l[[10]] <- list(name="llik_IEVAR_grid", response_heuristic=NA,
               integrated=TRUE, goal_oriented=TRUE, 
               int_grid_settings=list(n_prior=500L, n_mcmc=0L))

# max-var, pure sequential.
l[[11]] <- list(name="llik_neg_var_gp", response_heuristic=NA,
                integrated=FALSE, goal_oriented=FALSE)

# max-var-lik, pure sequential.
l[[12]] <- list(name="llik_neg_var_lik", response_heuristic=NA,
                integrated=FALSE, goal_oriented=TRUE)

#
# Add settings applied to all acquisitions.
#

for(i in seq_along(l)) {
  l[[i]]$candidate_settings <- candidate_settings
  l[[i]]$acq_args <- acq_args
}

# Save to file.
dir.create(out_dir, recursive=TRUE)
saveRDS(l, file.path(out_dir, "acq_settings_list.rds"))




