#
# analyze_approx_mcmc.r
# This script is intended to be run after `run_approx_mcmc.r` and 
# `postprocess_approx_mcmc.r` to provide summaries of the approximate MCMC 
# output and compare the various approximate samples to the exact MCMC samples.
# This script assumes that the MCMC output saved to file has already been 
# post-processed and mixing issues have been addressed.
# 
# Andrew Roberts
#

library(ggplot2)
library(data.table)
library(assertthat)
library(docopt)

experiment_tag <- "vsem"
round <- 1L

# ------------------------------------------------------------------------------
# Settings 
# ------------------------------------------------------------------------------

# Base directory: all paths are relative to this directory.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")

# Directory to source code.
src_dir <- file.path(base_dir, "src")

# Define directories and ensure required paths exist.
experiment_dir <- file.path(base_dir, "output", "gp_inv_prob", experiment_tag)


mcmc_settings_path <- file.path(experiment_dir, "mcmc_approx_settings.rds")
inv_prob_path <- file.path(experiment_dir, "inv_prob_setup", "inv_prob_list.rds")

# Source required scripts.
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
source(file.path(src_dir, "sim_study_functions.r"))

# Load R project. 
# renv::load(base_dir)
# print(".libPaths()")
# print(.libPaths())
# renv::status()

# ------------------------------------------------------------------------------
# Read and compile MCMC output.  
# ------------------------------------------------------------------------------

# Load MCMC tags.
# TODO

mcmc_tag <- "marginal-rect"
em_tag <- "llik_quad_mean"
design_tag <- "LHS_200"


info <- get_samp_dt_reps_agg(experiment_dir, round, mcmc_tag, em_tag, design_tag, 
                             only_valid=TRUE, format_long=FALSE)
samp_dt_reps_agg <- info$samp

# ------------------------------------------------------------------------------
# Plot results.
# ------------------------------------------------------------------------------

hist(samp_dt_reps_agg[param_name=="KEXT", mean], 10)






