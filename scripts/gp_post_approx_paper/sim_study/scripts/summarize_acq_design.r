# summarize_acq_design.r
#
# Currently, this is very rough code for reading sequential design results,
# conditioning the new emulators, and summarizing emulator error metrics.
# This file should almost be certainly split into multiple steps (e.g., the 
# emulators should probably be updated at the end of the sequential design
# step).
#
# Andrew Roberts

library(data.table)
library(ggplot2)

experiment_tag <- "vsem"
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")

# Directory to source code.
src_dir <- file.path(base_dir, "src")

# Base experiment directory.
experiment_dir <- file.path(base_dir, "output", "gp_inv_prob", experiment_tag)

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
source(file.path(src_dir, "sim_study_functions.r"))


round <- 2L
acq_id <- 11L
mcmc_tag_prev <- "mcwmh-joint-rect"
mcmc_id_prev <- "1008787650"





load_acq_data <- function(experiment_dir, round, acq_id, mcmc_tag_prev,
                          mcmc_id_prev) {
  
  results_path <- file.path(experiment_dir, paste0("round", round), 
                            "design", paste0("acq_", acq_id), mcmc_tag_prev,
                            mcmc_id_prev, "acq_results.rds")
  acq_results <- readRDS(results_path)
  return(acq_results)
}

load_acq_data_reps <- function(experiment_dir, round, acq_id, mcmc_tag_prev,
                               design_tag_prev, em_tag_prev) {
  .NotYetImplemented()
}



