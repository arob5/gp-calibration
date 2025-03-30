#
# seq_design_no_mcmc.r
#
# Evaluating the surrogate sequential design procedures without using any MCMC
# procedures. Practically, this means that we can test the longer term behavior
# of these methods in a way that would be computationally costly with MCMC.
# This also means that candidate and weighting point sets will typically be 
# with respect to the prior or a uniform measure, rather than the approximate
# posterior.
#
# Andrew Roberts
# 

library(ggplot2)
library(data.table)
library(assertthat)
library(docopt)
library(tictoc)

# Temporarily putting some settings here. These are used to subset the 
# acquisition settings data.table so that only the selected settings will 
# be run.
acq_ids <- c(1,2)

# -----------------------------------------------------------------------------
# docopt string for parsing command line arguments.  
# -----------------------------------------------------------------------------

"Usage:
  test_docopt.r [options]
  test_docopt.r (-h | --help)

Options:
  -h --help                                 Show this screen.
  --experiment_tag=<experiment_tag>         Used to locate the base output directory.
  --design_tag=<design_tag>                 Design tag.
  --em_tag=<em_tag>                         Emulator tag.
  --em_id=<em_id>                           Emulator ID.
  --n_candidates=<n_candidates>             Number of candidate points for grid-based optimization.
  --n_batch=<n_batch>                       Number of new design points to acquire.
" -> doc

# Read command line arguments.
cmd_args <- docopt(doc)
experiment_tag <- cmd_args$experiment_tag
design_tag <- cmd_args$design_tag
em_tag <- cmd_args$em_tag
em_id <- cmd_args$em_id
n_candidates <- as.integer(cmd_args$n_candidates)
n_batch <- as.integer(cmd_args$n_batch)

print(paste0("experiment_tag: ", experiment_tag))
print(paste0("design_tag: ", design_tag))
print(paste0("em_tag: ", em_tag))
print(paste0("em_id: ", em_id))
print(paste0("n_candidates: ", n_candidates))
print(paste0("n_batch: ", n_batch))

# ------------------------------------------------------------------------------
# Settings 
# ------------------------------------------------------------------------------

# Base directory: all paths are relative to this directory.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")

# Directory to source code.
src_dir <- file.path(base_dir, "src")

# Define directories and ensure required paths exist.
experiment_dir <- file.path(base_dir, "output", "gp_inv_prob", experiment_tag)
em_dir <- file.path(experiment_dir, "round1", "em", em_tag)
base_design_dir <- file.path(experiment_dir, "round1", "design")
inv_prob_dir <- file.path(experiment_dir, "inv_prob_setup")
base_out_dir <- file.path(experiment_dir, "side_analysis", "seq_design_no_mcmc")
acq_settings_path <- file.path(base_out_dir, "acq_settings.csv")


print(paste0("base_out_dir: ", base_out_dir))

# Source required scripts.
source(file.path(src_dir, "general_helper_functions.r"))
source(file.path(src_dir, "inv_prob_test_functions.r"))
source(file.path(src_dir, "statistical_helper_functions.r"))
source(file.path(src_dir, "plotting_helper_functions.r"))
source(file.path(src_dir, "seq_design.r"))
source(file.path(src_dir, "gp_helper_functions.r"))
source(file.path(src_dir, "gpWrapper.r"))
source(file.path(src_dir, "llikEmulator.r"))
source(file.path(src_dir, "gp_mcmc_functions.r"))
source(file.path(src_dir, "seq_design_gp.r"))
source(file.path(src_dir, "seq_design_for_post_approx.r"))

# Load R project. 
renv::load(base_dir)
print(".libPaths()")
print(.libPaths())
renv::status()

# ------------------------------------------------------------------------------
# Load inverse problem setup and exact MCMC samples.
# ------------------------------------------------------------------------------

inv_prob <- readRDS(file.path(inv_prob_dir, "inv_prob_list.rds"))
test_info_prior <- readRDS(file.path(inv_prob_dir, "test_info_prior.rds"))
test_info_post <- readRDS(file.path(inv_prob_dir, "test_info_post.rds"))
llik_func <- inv_prob$llik_obj$get_llik_func()

# ------------------------------------------------------------------------------
# Load log-likelihood emulator.
# ------------------------------------------------------------------------------

# Fetch design ID used to train the emulator.
em_id_curr <- em_id
em_id_map <- fread(file.path(em_dir, "id_map.csv"))
design_id <- em_id_map[em_id==em_id_curr, design_id]
design_tag <- em_id_map[em_id==em_id_curr, design_tag]

out_dir <- file.path(base_out_dir, em_tag, em_id, design_tag, design_id)
print(paste0("out_dir: ", out_dir))
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# Load emulator and design.
em_llik <- readRDS(file.path(em_dir, em_id, "em_llik.rds"))
design_info <- readRDS(file.path(base_design_dir, design_tag, paste0(design_id, ".rds")))

if(design_info$design_method != "LHS") {
  stop("Currently this script is set up to augment an existing LHS design.")
}

# ------------------------------------------------------------------------------
# Read acquisition function settings.
# ------------------------------------------------------------------------------

acq_settings <- fread(acq_settings_path)
acq_settings <- acq_settings[acq_settings_id %in% acq_ids,, drop=FALSE]

print("Running acquisition settings:")
print(acq_settings)

# ------------------------------------------------------------------------------
# Set up metrics to track within sequential design loop.
# ------------------------------------------------------------------------------

tracking_settings <- list(interval=10L, func_list=list())

tracking_settings$func_list$pw_prior <- function(model) {
  pred <- model$emulator_model$calc_pred_multi_func(list(crps="crps", mse="mse"),
                                                    type="pw", X_new=test_info_prior$input, 
                                                    Y_new=matrix(test_info_prior$llik, ncol=1L))
  pred[, .(mean=mean(y1)), by=func]
}

tracking_settings$func_list$pw_post <- function(model) {
  pred <- model$emulator_model$calc_pred_multi_func(list(crps="crps", mse="mse"),
                                                         type="pw", X_new=test_info_post$input, 
                                                         Y_new=matrix(test_info_post$llik, ncol=1L))
  pred[, .(mean=mean(y1)), by=func]
}

tracking_settings$func_list$agg_prior <- function(model) {
  model$emulator_model$calc_pred_func("log_score", type="agg", 
                                      X_new=test_info_prior$input, 
                                      Y_new=matrix(test_info_prior$llik, ncol=1L),
                                      return_cov=TRUE)
}

tracking_settings$func_list$agg_post<- function(model) {
  model$emulator_model$calc_pred_func("log_score", type="agg", 
                                      X_new=test_info_post$input, 
                                      Y_new=matrix(test_info_post$llik, ncol=1L),
                                      return_cov=TRUE)
}


# ------------------------------------------------------------------------------
# Candidate and Grid Points
# ------------------------------------------------------------------------------

# Finite set of points that will be optimized over. For each em ID/design ID
# we use the same candidate grid for all acquisition functions.
candidate_grid_path <- file.path(out_dir, "candidate_grid.rds")
  
if(file.exists(candidate_grid_path)) {
  print(paste0("Reading candidate points from: ", candidate_grid_path))
  candidate_grid <- readRDS(candidate_grid_path)
} else {
  candidate_grid <- update_LHS_sample(design_info$input, n_batch=n_candidates,
                                      prior_dist_info=inv_prob$par_prior)
  candidate_grid <- candidate_grid[(design_info$N_design+1L):nrow(candidate_grid),,drop=FALSE]
  print(paste0("Writing candidate points to: ", candidate_grid_path))
  saveRDS(candidate_grid, candidate_grid_path)
}

# Same for the grid points used for numerical integration over the input space.
# This is only required for some acquisition functions.
grid_points_path <- file.path(out_dir, "integration_grid_points.rds")

if(file.exists(grid_points_path)) {
  print(paste0("Reading grid points from: ", grid_points_path))
  grid_points <- readRDS(grid_points_path)
} else {
  n_grid <- acq_settings[!is.na(n_grid), n_grid]
  if(length(n_grid) == 0L) { # Acquisition functions do not require grid points.
    grid_points <- NULL
  } else {
    if(!all(n_grid[1] == n_grid)) {
      stop("For now, we require `n_grid` to be the same for all acq funcs.")
    }
    
    grid_points <- get_batch_design("LHS", N_batch=n_grid[1], 
                                    prior_params=inv_prob$par_prior)
    print(paste0("Writing grid points to: ", grid_points_path))
    saveRDS(grid_points, grid_points_path)
  }
}

# ------------------------------------------------------------------------------
# Sequential Design
# ------------------------------------------------------------------------------

# Arguments for `run_seq_design`.
args <- list(model=em_llik, n_batch=n_batch, opt_method="grid", 
             true_func=llik_func, reoptimize_hyperpar=FALSE,
             candidate_grid=candidate_grid, tracking_settings=tracking_settings,
             grid_points=grid_points)

for(i in 1:nrow(acq_settings)) {
  acq_settings_id <- acq_settings[i, acq_settings_id] 
    
  print(paste0("-----> ID: ", acq_settings_id, " ; Name: ", acq_settings[i, name]))
  
  # Append acquisition function settings to function arguments list.
  args$acq_func_name <- acq_settings[i, name]
  args$response_heuristic <- acq_settings[i, response_heuristic]
  if(is.na(args$response_heuristic)) args$response_heuristic <- NULL
  
  # Run sequential design loop.
  tic()
  acq_results <- do.call(run_seq_design, args)
  tictoc_info <- toc()
  
  # Store additional information to save.
  acq_results$runtime <- tictoc_info$toc - tictoc_info$tic
  acq_results$acq_settings <- acq_settings[i,]
  
  # Save results to file. Files are named by their acquisition ID.
  filename <- paste0("acq_results_id_", acq_settings_id, ".rds")
  print(paste0("Saving: ", filename))
  saveRDS(acq_results, file.path(out_dir, filename))
}


