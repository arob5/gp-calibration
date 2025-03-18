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

# ------------------------------------------------------------------------------
# Settings 
# ------------------------------------------------------------------------------

# Base directory: all paths are relative to this directory.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")

# Directory to source code.
src_dir <- file.path(base_dir, "src")

# Settings.
experiment_tag <- "vsem"
design_tag <- "LHS_200"
em_tag <- "llik_quad_mean"

# Define directories and ensure required paths exist.
experiment_dir <- file.path(base_dir, "output", "gp_inv_prob", experiment_tag)
em_dir <- file.path(experiment_dir, "round1", "em", em_tag)
base_design_dir <- file.path(experiment_dir, "round1", "design")
inv_prob_dir <- file.path(experiment_dir, "inv_prob_setup")
out_dir <- file.path(experiment_dir, "side_analysis")

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

# ------------------------------------------------------------------------------
# Sequential Design 
# ------------------------------------------------------------------------------

# Number of candidate points to optimize over.
n_candidates <- 500L

# Number of points to acquire.
n_batch <- 10L

# An an initial test, just considering a single replicate.
em_id_curr <- "1008787650"

# Fetch design ID used to train the emulator.
em_id_map <- fread(file.path(em_dir, "id_map.csv"))
design_id <- em_id_map[em_id==em_id_curr, design_id]
design_tag <- em_id_map[em_id==em_id_curr, design_tag]

# Load emulator and design.
em_llik <- readRDS(file.path(em_dir, em_id_curr, "em_llik.rds"))
design_info <- readRDS(file.path(base_design_dir, design_tag, paste0(design_id, ".rds")))


candidate_grid <- get_batch_design("LHS", N_batch=n_candidates, 
                                   prior_params=inv_prob$par_prior)
acq_func_name <- "llik_neg_var_gp"
llik_func <- inv_prob$llik_obj$get_llik_func()

# TODO: need to add an "evaluation_function" option that computes some evaluation
# criterion/criteria every iteration, or every x iterations.

acq_results <- run_seq_design(em_llik, acq_func_name, n_batch, opt_method="grid",
                              response_heuristic=NULL, true_func=llik_func, 
                              reoptimize_hyperpar=FALSE, 
                              candidate_grid=candidate_grid)





