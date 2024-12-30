#
# postprocess_approx_mcmc.r
# This script is intended to be run after `run_approx_mcmc.r` to identify 
# potential issues with MCMC runs. See `run_approx_mcmc.r` for details on 
# specifying filepaths, etc., as this script uses the same setup. This script 
# should be run before `analyze_approx_mcmc.r`, as this latter script assumes 
# individual chains are well-mixed, burn-in has already been dropped, issues 
# corrected, etc.
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

# Set variables controlling filepaths.
experiment_tag <- "vsem"
run_id <- "mcmc_round1"
em_dir <- "init_emulator/LHS_250"

print(paste0("experiment_tag: ", experiment_tag))
print(paste0("run_id: ", run_id))
print(paste0("em_dir: ", em_dir))

# Define directories and ensure required paths exist.
experiment_dir <- file.path(base_dir, "output", "gp_inv_prob", experiment_tag)
base_out_dir <- file.path(experiment_dir, run_id, em_dir)

print(paste0("experiment_dir: ", experiment_dir))
print(paste0("base_out_dir: ", base_out_dir))

mcmc_settings_path <- file.path(experiment_dir, "mcmc_approx_settings.rds")
inv_prob_path <- file.path(experiment_dir, "inv_prob_setup", "inv_prob_list.rds")

print("-----> Checking required files exist:")
assert_that(file.exists(mcmc_settings_path))
assert_that(file.exists(inv_prob_path))

mcmc_settings_list <- readRDS(mcmc_settings_path)

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
# renv::load(base_dir)
# print(".libPaths()")
# print(.libPaths())
# renv::status()

# ------------------------------------------------------------------------------
# Define helper function to process a single MCMC run.
# ------------------------------------------------------------------------------

process_test_label <- function(lbl, em_id, rhat_threshold=1.05) {
  # A chain is defined as "valid" if the within-chain split R-hat values for 
  # all parameters do not exceed `rhat_threshold`. `n_itr` is the total number  
  # of iterations with non-NA values across all valid chains. An MCMC run 
  # (test_label) is marked as valid if it has at least one valid chain and 
  # at least one valid iteration.

  # Check to see if MCMC output exists.
  out_dir <- file.path(base_out_dir, em_id)
  samp_list_path <- file.path(out_dir, paste0("mcmc_samp_", lbl, ".rds"))
  if(!file.exists(samp_list_path)) {
    mcmc_summary <- data.table(em_id=em_id, test_label=lbl, n_valid_chains=0L, 
                               max_rhat=NA, n_valid_itr=0L, status="missing")
    return(list(summary=mcmc_summary, rhat=NULL))
  }
  
  # Read MCMC list.
  samp_list <- readRDS(samp_list_path)

  # Check to see if error occurred during run.
  err_occurred <- !is.null(samp_list$output_list[[1]]$condition)
  if(err_occurred) {
    mcmc_summary <- data.table(em_id=em_id, test_label=lbl, n_valid_chains=0L,
                               max_rhat=NA, n_valid_itr=0L, status="err")
    return(list(summary=mcmc_summary, rhat=NULL))
  }

  # MCMC samples for the current test label.
  samp_dt <- select_mcmc_samp(samp_list$samp, test_labels=lbl)
  
  # Within-chain Rhat statistics.
  id_cols <- c("test_label", "chain_idx", "param_type", "param_name")
  rhat_dt <- calc_R_hat(samp_dt, within_chain=TRUE)$R_hat_vals
  rhat_dt <- rhat_dt[, .SD, .SDcols=c(id_cols, "R_hat")]
  rhat_max_dt <- rhat_dt[, .(rhat=max(R_hat, na.rm=TRUE)), by=.(test_label, chain_idx)]
  rhat_max_dt[, valid := (rhat <= rhat_threshold)]
  
  # Compute MCMC run summary.
  n_valid_chains <- rhat_max_dt[, sum(valid, na.rm=TRUE)]
  max_rhat <- rhat_max_dt[, max(rhat)]
  valid_chains <- rhat_max_dt[valid==TRUE, chain_idx]
  if(length(valid_chains) > 0) {
    n_valid_itr <- samp_dt[(chain_idx %in% valid_chains) & !is.na(sample), 
                            .N, by=.(param_type, param_name)][, unique(N)]
  } else {
    n_valid_itr <- 0L
  }
  
  mcmc_summary <- data.table(em_id = em_id, 
                             test_label = lbl,
                             n_valid_chains = n_valid_chains,
                             max_rhat = max_rhat, 
                             n_valid_itr = n_valid_itr,
                             status = "valid")
  mcmc_summary[(n_valid_chains < 1) | (n_valid_itr < 1), status := "invalid"]
  rhat_dt[, em_id := em_id]
  
  return(list(summary=mcmc_summary, rhat=rhat_dt))
}


# ------------------------------------------------------------------------------
# Compute summaries for each MCMC run within `base_out_dir`.
# ------------------------------------------------------------------------------

em_id_dirs <- setdiff(list.files(base_out_dir), 
                      c("mcmc_summary.csv", "rhat_chain_summary.csv"))
lbls <- names(mcmc_settings_list)
mcmc_summary <- data.table(em_id=character(), test_label=character(), 
                           n_valid_chains=integer(), max_rhat=numeric(), 
                           n_valid_itr=integer(), status=character())
rhat_chain_summary <- data.table(em_id=character(), test_label=character(), 
                                 chain_idx=integer(), param_type=character(), 
                                 param_name=character(), R_hat=numeric())

for(em_id in em_id_dirs) {
  print(paste0("-----> Processing: ", em_id))
  
  for(lbl in lbls) {
    cat("\t", lbl, "\n")
    results <- process_test_label(lbl, em_id, rhat_threshold=1.05)
    mcmc_summary <- rbindlist(list(mcmc_summary, results$summary), 
                              use.names=TRUE, fill=TRUE)
    
    if(!is.null(results$rhat)) {
      rhat_chain_summary <- rbindlist(list(rhat_chain_summary, results$rhat),
                                      use.names=TRUE)
    }
  }
}

# ------------------------------------------------------------------------------
# Save summaries to file.
# ------------------------------------------------------------------------------

fwrite(mcmc_summary, file.path(base_out_dir, "mcmc_summary.csv"))
fwrite(rhat_chain_summary, file.path(base_out_dir, "rhat_chain_summary.csv"))









