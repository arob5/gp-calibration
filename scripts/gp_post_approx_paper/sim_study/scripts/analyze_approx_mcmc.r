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

# -----------------------------------------------------------------------------
# docopt string for parsing command line arguments.  
# -----------------------------------------------------------------------------

"Usage:
  test_docopt.r [options]
  test_docopt.r (-h | --help)

Options:
  -h --help                                 Show this screen.
  --experiment_tag=<experiment_tag>         Used to locate the base output directory.
  --run_id=<run_id>                         Used to define output directory.
  --em_dir=<out_dir>                        llikEmulator base directory.
  --em_id=<em_id>                           llikEmulator sub-directory.
" -> doc

# ------------------------------------------------------------------------------
# Settings 
# ------------------------------------------------------------------------------

# Base directory: all paths are relative to this directory.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")

# Directory to source code.
src_dir <- file.path(base_dir, "src")

# Set variables controlling filepaths.
experiment_tag <- "vsem"
run_id <- "test_run6"
em_dir <- "init_emulator/LHS_250"
em_id <- "1018157756"

print(paste0("experiment_tag: ", experiment_tag))
print(paste0("run_id: ", run_id))
print(paste0("em_dir: ", em_dir))
print(paste0("em_id: ", em_id))

# Define directories and ensure required paths exist.
experiment_dir <- file.path(base_dir, "output", "gp_inv_prob", experiment_tag)
em_id_dir <- file.path(experiment_dir, em_dir, em_id)
out_dir <- file.path(experiment_dir, run_id, em_dir, em_id)

print(paste0("experiment_dir: ", experiment_dir))
print(paste0("em_id_dir: ", em_id_dir))
print(paste0("out_dir: ", out_dir))

mcmc_settings_path <- file.path(experiment_dir, "mcmc_approx_settings.rds")
inv_prob_path <- file.path(experiment_dir, "inv_prob_setup", "inv_prob_list.rds")

print("-----> Checking required files exist:")
assert_that(file.exists(mcmc_settings_path))
assert_that(file.exists(inv_prob_path))

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

# Load R project. 
renv::load(base_dir)
print(".libPaths()")
print(.libPaths())
renv::status()

# Chain weights.
chain_weights <- calc_chain_weights(samp_list$info)



