#
# inv_prob_setup.r
# Part 1 in the inverse problem simulation study workflow. This script:
# (1) loads an inverse problem from those predefined in 
#     `inv_prob_test_functions.r`.
# (2) runs (exact) MCMC to obtain baseline exact posterior samples to be used 
#     for validation later on. These samples will be combined with prior samples
#     and saved together.
# (3) generates two sets of validation points for later use in evaluating 
#     emulators. One set is sampled from the prior and the other from the 
#     exact posterior.
# (4) the information from steps (1) through (3) is all saved to file, for use 
#     in the subsequent parts of the simulation study.
#
# Unlike the scripts in the later parts of the simulation study, this script 
# is not set up to be executed remotely on the cluster, given that its runtime 
# is typically reasonable. This can of course be changed if needed. Note 
# also that this script sets up a single, fixed inverse problem and various 
# associated information. There is no replication done here, as is the case 
# in the subsequent script.
#
# Andrew Roberts
#

library(ggplot2)
library(data.table)
library(assertthat)

# ------------------------------------------------------------------------------
# Settings 
# ------------------------------------------------------------------------------

# Seed for random number generator.
seed <- 36743235
set.seed(seed)

# Used to create the main directory that will store all outputs from the 
# experiment (sim study).
experiment_tag <- "vsem"

# Number of points in each of the two validation sets.
n_test_prior <- 500L
n_test_post <- 500L

# Sampling method to use in constructing the prior validation set.
design_method_test <- "LHS"

# Number of prior samples to include in the samples table that will be saved.
n_samp_prior <- 100000L

# Specifications for exact MCMC.
mcmc_settings <- list(test_label="exact", mcmc_func_name="mcmc_bt_wrapper", 
                      sampler="DEzs", n_itr=200000L, try_parallel=FALSE,  
                      n_chain=4L, defer_ic=TRUE)

# Starting iteration defining the end of the burn-in/warm-up for exact MCMC.
burn_in_start <- 100000L

# ------------------------------------------------------------------------------
# Setup 
# ------------------------------------------------------------------------------

# Filepath definitions.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")
src_dir <- file.path(base_dir, "src")
experiment_dir <- file.path(base_dir, "output", "gp_inv_prob", experiment_tag)
out_dir <- file.path(experiment_dir, "inv_prob_setup")

# Source required files.
source(file.path(src_dir, "general_helper_functions.r"))
source(file.path(src_dir, "inv_prob_test_functions.r"))
source(file.path(src_dir, "statistical_helper_functions.r"))
source(file.path(src_dir, "plotting_helper_functions.r"))
source(file.path(src_dir, "seq_design.r"))
source(file.path(src_dir, "gp_helper_functions.r"))
source(file.path(src_dir, "gpWrapper.r"))
source(file.path(src_dir, "llikEmulator.r"))
source(file.path(src_dir, "mcmc_helper_functions.r"))
source(file.path(src_dir, "gp_mcmc_functions.r"))

# Create experiment directories.
dir.create(experiment_dir, recursive=TRUE)
dir.create(out_dir)

# Save settings.
setup_settings <- list(seed=seed, experiment_tag=experiment_tag, 
                       n_test_prior=n_test_prior, n_test_post=n_test_post,
                       design_method_test=design_method_test, 
                       n_samp_prior=n_samp_prior, mcmc_settings=mcmc_settings,
                       burn_in_start=burn_in_start)
saveRDS(setup_settings, file=file.path(out_dir, "setup_settings.RData"))

# ------------------------------------------------------------------------------
# Inverse problem setup 
# ------------------------------------------------------------------------------

inv_prob <- get_vsem_test_1(default_conditional=FALSE, 
                            default_normalize=TRUE)
saveRDS(inv_prob, file=file.path(out_dir, "inv_prob_list.RData"))

# ------------------------------------------------------------------------------
# Exact MCMC 
# ------------------------------------------------------------------------------

# MCMC sampling using exact likelihood.
mcmc_list <- run_mcmc(inv_prob$llik_obj, inv_prob$par_prior, mcmc_settings)

# Table storing MCMC samples.
samp_dt <- mcmc_list$samp
mcmc_metadata <- mcmc_list$output_list

# Attach prior samples.
prior_samp <- sample_prior(inv_prob$par_prior, n=n_samp_prior)
samp_dt <- append_samples_mat(samp_dt, prior_samp, param_type="par", 
                              test_label="prior")

# Save to file.
saveRDS(samp_dt, file=file.path(out_dir, "samp_exact.RData"))
saveRDS(mcmc_metadata, file=file.path(out_dir, "samp_exact_metadata.RData"))

# ------------------------------------------------------------------------------
# Construct validation test points.
# ------------------------------------------------------------------------------

# Validation inputs sampled from prior.
test_info_prior <- get_init_design_list(inv_prob, design_method_test, n_test_prior)
saveRDS(test_info_prior, file=file.path(out_dir, "test_info_prior.RData"))

# Validation inputs sub-sampled from true posterior.
exact_mcmc_samp <- select_mcmc_samp_mat(samp_dt, test_label="exact",
                                        param_type="par", itr_start=burn_in_start)

test_info_post <- get_init_design_list(inv_prob, "subsample",
                                       N_design=n_test_post, 
                                       design_candidates=exact_mcmc_samp)
saveRDS(test_info_post, file=file.path(out_dir, "test_info_post.RData"))
