#
# seq_design.r
#
# This script runs one round of sequential acquisition of design points. It 
# therefore should only be run for rounds >= 2 (i.e., not the initial design
# round). The user specifies:
#  (i.) the experiment ID.
#  (i.) the round.
#  (iii.) An MCMC tag and ID.
#  (iv.) Sequential design settings.
# 
# The round is an integer; e.g., `2` for the first round of sequential 
# acquisition after the initial design round (round 1). In the round 2 example.
# all outputs from this run will be stored in 
# `<experiment_dir>/round2/design/<acq_id>/<mcmc_tag>/<mcmc_id>`. The MCMC 
# tag/ID and acquisition ID are discussed below.
#
# The MCMC tag/ID correspond to an MCMC run from the previous round. Recall that
# MCMC IDs are unique within tags, so both the tag and ID are required to specify
# a unique run from the previous round. This will be used to identify the 
# current emulator that will be updated. It will also be used to read MCMC
# output from the current posterior approximation, which can be used to specify
# candidate points for discrete optimization, and weights used in some
# acquisition criteria that integrate over the design space.
#
# Perhaps the most important sequential design setting is the number of points
# to acquire, which is typically called `n_batch` in the code. Currently,
# this is passed as a command line argument to this script. While it is not
# explicitly enforced, the code is set up so that all sequential design runs 
# within a single round of sequential design will use the same `n_batch` value.
# The other major portion of the settings is the acquisition function settings.
# This is specified by passing and acquisition function ID "acq_id" via a 
# command line argument. The acquisition function ID corresponds to the index
# of an element of the list, 
# `<experiment_dir>/<round_dir>/design/acq_settings_list.rds`, which must
# be saved in advance of running this script.
#
# Andrew Roberts
#

# TODO:
# - For IEVAR criteria, check the "plugin" argument. I might have to set this to TRUE always.
#      -> Conclusion: only relevant for forward model emulator, so don't worry about it for now.
# - Implement batch heuristics.
# - Implement log-normal adjustments:
#      -> Ensure that IVAR is actually using the adjustment.
#      -> For IEVAR, include adjustment heuristic by applying the adjustment post-hoc.
# - Make sure I am passing the "rectified" adjustment to all of the algorithms.
# - Add prior to IEVAR so it targets log posterior density, not just llik. Currently
#   doesn't matter since I'm using uniform priors.

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
  --round=<round>                           Integer specifying the round.
  --mcmc_tag=<mcmc_tag>                     MCMC tag from previous round.
  --mcmc_id=<mcmc_id>                       MCMC ID from previous round.
  --n_batch=<n_batch>                       Number of points to acquire.
  --acq_id=<acq_id>                         Acquisition function ID.
" -> doc

# ------------------------------------------------------------------------------
# Settings 
# ------------------------------------------------------------------------------

# Base directory: all paths are relative to this directory.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")

# Directory to source code.
src_dir <- file.path(base_dir, "src")

# Read command line arguments.
cmd_args <- docopt(doc)
experiment_tag <- cmd_args$experiment_tag
round <- as.integer(cmd_args$round)
mcmc_tag <- cmd_args$mcmc_tag
mcmc_id <- cmd_args$mcmc_id
n_batch <- cmd_args$n_batch
acq_id <- cmd_args$acq_id

print(paste0("experiment_tag: ", experiment_tag))
print(paste0("round: ", round))
print(paste0("mcmc_tag: ", mcmc_tag))
print(paste0("mcmc_id: ", mcmc_id))
print(paste0("n_batch: ", n_batch))
print(paste0("acq_id: ", acq_id))

# Define directories.
prev_round <- round - 1L
experiment_dir <- file.path(base_dir, "output", "gp_inv_prob", experiment_tag)
round_dir <- file.path(experiment_dir, paste0("round", round))
design_dir <- file.path(round_dir, "design")
acq_dir <- file.path(design_dir, paste0("acq_", acq_id))
mcmc_input_dir <- file.path(experiment_dir, paste0("round", prev_round),
                            "mcmc", mcmc_tag, mcmc_id)
out_dir <- file.path(acq_dir, mcmc_tag, mcmc_id) # Outputs fromn this file saved here.

# File path to acquisition settings.
acq_settings_path <- file.path(design_dir, "acq_settings_list.rds")

print(paste0("MCMC input dir (from previous round): ", mcmc_input_dir))
print(paste0("Output directory: ", out_dir))
print(paste0("Acq settings path: ", acq_settings_path))

if(!file.exists(acq_settings_path)) {
  stop("Acquisition settings path does not exist.")
}

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
source(file.path(src_dir, "seq_design.r"))
source(file.path(src_dir, "seq_design_gp.r"))
source(file.path(src_dir, "seq_design_for_post_approx.r"))

# Create output directory.
dir.create(out_dir, recursive=TRUE)








