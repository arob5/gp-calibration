#
# save_acq_func_settings.r
# Saves settings for different acquisition function setups to use in 
# `seq_design_no_mcmc.r`.
#
# Andrew Roberts
#

library(data.table)

experiment_tag <- "vsem"

# Specify where settings will be saved.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")
experiment_dir <- file.path(base_dir, "output", "gp_inv_prob", experiment_tag)
out_dir <- file.path(experiment_dir, "side_analysis", "seq_design_no_mcmc")

# Number of grid points for all integrated uncertainty criteria.
n_grid <- 500L


# NOTE: I believe the below results where using the truncated Gaussian predictive dist.
# NOTE: (these numbers are probably lower bounds, since as the GP acquires more
#        design points prediction will be more expensive).
# - With kergp, llik_IEVAR_grid, 1000 candidates, 500 grid points a single 
#   iteration took about 300 sec. Implying close to 17 hours for 200 acquisitions.
# - Also tried this setup with llik_IENT_grid_gp, which took about 174 sec.
#  - Tried same setup with `llik_IEVAR_grid` but using hetGP instead. Took only
#    about 44.483 seconds. Much faster.

# Create table of settings, one row per a particular acquisition function 
# setup. The flag `integrated` indicates whether or not evaluation of the 
# acquisition function requires approximating an integral over the input space.
# The flag `goal_oriented` indicates whether the acquisition function explicitly
# targets the goal of posterior approximation; practically, this typically 
# indicates whether the acquisition is simply defined wrt the underlying GP,
# or whether it takes into account the likelihood structure.

l <- list()

l[[1]] <- list(name="llik_IEVAR_grid", response_heuristic=NA,
               integrated=TRUE, goal_oriented=TRUE, n_grid=n_grid)
l[[2]] <- list(name="llik_IENT_grid_gp", response_heuristic=NA,
               integrated=TRUE, goal_oriented=FALSE, n_grid=n_grid)
l[[3]] <- list(name="llik_neg_var_lik", response_heuristic=NA,
               integrated=FALSE, goal_oriented=TRUE)
l[[4]] <- list(name="llik_neg_entropy_gp", response_heuristic=NA,
               integrated=FALSE, goal_oriented=FALSE)
l[[5]] <- list(name="llik_neg_var_gp", response_heuristic=NA,
               integrated=FALSE, goal_oriented=FALSE)

dt <- data.table::rbindlist(l, use.names=TRUE, fill=TRUE)
dt[, acq_settings_id := 1:nrow(dt)]
fwrite(dt, file.path(out_dir, "acq_settings.csv"))




