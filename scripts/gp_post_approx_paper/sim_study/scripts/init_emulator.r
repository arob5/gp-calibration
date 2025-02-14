#
# init_emulator.r
# Part 2 in the inverse problem simulation study workflow. This script:
# (1) loads inverse problem data saved to file in `inv_prob_setup.r`. 
# (2) Generates a (random) initial design to train the emulator.
# (3) Fits a Gaussian process (GP) emulator, and uses the fit emulator to 
#     define a log-likelihood surrogate object (llikEmulator).
# (4) Generate GP and llikEmulator predictions at validation points, and compute
#     various metrics summarizing emulator fit.
# (5) the information from steps (1) through (4) is all saved to file, for use 
#     in the subsequent parts of the simulation study.
#
# This script is set up to be run remotely on the cluster. The intention is to 
# run this script (with different random seeds) on a large number of remote jobs
# in order to create replicate initial designs. Each initial design replicate 
# and its associated emulator is saved to a separate subdirectory. For now, a 
# single GP specification is used here, but this could be changed in the future 
# to allow comparison across different GP models.
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
  --run_id=<run_id>                         Directory name to save outputs.
  --n_design=<n_design>                     Number of design points.
  --design_method=<design_method>           Algorithm to generate design points.
  --n_rep=<n_rep>                           Number of replicate designs.
" -> doc

# ------------------------------------------------------------------------------
# Settings 
# ------------------------------------------------------------------------------

# Used to identify the location of the experiment, which has already been 
# initiated by `inv_prob_setup.r`.
experiment_tag <- "vsem"

print("--------------------Running `init_emulator.r` --------------------")
print(paste0("Experiment tag: ", experiment_tag))

# ------------------------------------------------------------------------------
# Setup 
# ------------------------------------------------------------------------------

# Read command line arguments.
cmd_args <- docopt(doc)
n_design <- as.integer(cmd_args$n_design)
design_method <- cmd_args$design_method
n_rep <- as.integer(cmd_args$n_rep)
run_id <- cmd_args$run_id
print(paste0("run_id: ", run_id))

# Filepath definitions.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")
src_dir <- file.path(base_dir, "src")
experiment_dir <- file.path(base_dir, "output", "gp_inv_prob", experiment_tag)
setup_dir <- file.path(experiment_dir, "inv_prob_setup")
base_out_dir <- file.path(experiment_dir, "init_emulator", run_id)

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

# Load R project. 
renv::load(base_dir)
print(".libPaths()")
print(.libPaths())
renv::status()

# Create output directories.
dir.create(base_out_dir, recursive=TRUE)

# Load inverse problem setup information.
inv_prob <- readRDS(file.path(setup_dir, "inv_prob_list.rds"))
test_info_prior <- readRDS(file.path(setup_dir, "test_info_prior.rds"))
test_info_post <- readRDS(file.path(setup_dir, "test_info_post.rds"))

# Log-likelihood bounds.
llik_bounds <- inv_prob$llik_obj$get_llik_bounds()

# Settings list.
settings <- list(run_id=run_id, experiment_tag=experiment_tag, 
                 n_design=cmd_args$n_design, 
                 design_method=cmd_args$design_method)

# ------------------------------------------------------------------------------
# Define function that creates a single replicate initial design and 
# accompanying emulator.
# ------------------------------------------------------------------------------

construct_init_em <- function(seed) {
  print("-------------- Initiating new design --------------")
  print("-----> Creating output directory")
  set.seed(seed)
  out_dir <- file.path(base_out_dir, seed)
  dir.create(out_dir, recursive=TRUE)
  settings_seed <- settings
  settings_seed$seed <- seed
  saveRDS(settings_seed, file=file.path(out_dir, "init_em_settings.rds"))
  print(paste0("Seed: ", seed))
  print(paste0("Directory: ", out_dir))

  # Generate design.
  print("-----> Generating design")

  # Generate training design. 
  design_info <- get_init_design_list(inv_prob, design_method, n_design)
  saveRDS(design_info, file=file.path(out_dir, "design_info.rds"))
  
  print("-----> Fit emulator")
  
  # Fit GP for log-likelihood.
  gp_obj <- gpWrapperKerGP(design_info$input, matrix(design_info$llik, ncol=1), 
                           scale_input=TRUE, normalize_output=TRUE)
  gp_obj$set_gp_prior("Gaussian", "quadratic", include_noise=FALSE)
  gp_obj$fit(multistart=10, trace=TRUE)
  print(gp_obj$summarize())

  # Instantiate and save log-likelihood emulator object.
  em_llik <- llikEmulatorGP("em_llik", gp_obj, default_conditional=FALSE, 
                            default_normalize=TRUE, lik_par=inv_prob$sig2_model, 
                            llik_bounds=llik_bounds)
  saveRDS(em_llik, file=file.path(out_dir, "em_llik.rds"))
  
  print("-----> Emulator predictions at test points")
  em_pred_list_prior <- em_llik$predict_emulator(test_info_prior$input, 
                                                 return_cov=TRUE)
  em_pred_list_post <- em_llik$predict_emulator(test_info_post$input,
                                                return_cov=TRUE)
  
  saveRDS(em_pred_list_prior, file=file.path(out_dir, "em_pred_list_prior.rds"))
  saveRDS(em_pred_list_post, file=file.path(out_dir, "em_pred_list_post.rds"))
  
  print("-----> llikEmulator error measures")
  err_types <- c("mse", "wmse", "mae", "wmae")
  em_llik_errs_prior <- em_llik$calc_lik_approx_pw_err(llik_true=test_info_prior$llik, 
                                                       input=test_info_prior$input,
                                                       err_type=err_types, 
                                                       em_pred_list=em_pred_list_prior,
                                                       return_type="data.table")
  
  em_llik_errs_post <- em_llik$calc_lik_approx_pw_err(llik_true=test_info_post$llik, 
                                                      input=test_info_post$input,
                                                      err_type=err_types, 
                                                      em_pred_list=em_pred_list_post,
                                                      return_type="data.table")
  
  saveRDS(em_llik_errs_prior, file=file.path(out_dir, "em_llik_errs_prior.rds"))
  saveRDS(em_llik_errs_post, file=file.path(out_dir, "em_llik_errs_post.rds"))
  
  print("-----> GP emulator error measures")
  gp_pw_err_types <- c("mae", "mse", "crps", "log_score")
  gp_agg_err_types <- c("mah", "log_score")
  
  gp_pw_errs_prior <- em_llik$emulator_model$calc_pred_multi_func(gp_pw_err_types, 
                                                                  type="pw", 
                                                                  pred_list=em_pred_list_prior, 
                                                                  Y_new=matrix(test_info_prior$llik, ncol=1))
  gp_pw_errs_post <- em_llik$emulator_model$calc_pred_multi_func(gp_pw_err_types, 
                                                                 type="pw", 
                                                                 pred_list=em_pred_list_post, 
                                                                 Y_new=matrix(test_info_post$llik, ncol=1))
  gp_agg_errs_prior <- em_llik$emulator_model$calc_pred_multi_func(gp_agg_err_types, 
                                                                   type="agg", 
                                                                   pred_list=em_pred_list_prior, 
                                                                   Y_new=matrix(test_info_prior$llik, ncol=1))
  gp_agg_errs_post <- em_llik$emulator_model$calc_pred_multi_func(gp_agg_err_types, 
                                                                  type="agg", 
                                                                  pred_list=em_pred_list_post, 
                                                                  Y_new=matrix(test_info_post$llik, ncol=1))
  
  saveRDS(gp_pw_errs_prior, file=file.path(out_dir, "gp_pw_errs_prior.rds"))
  saveRDS(gp_pw_errs_post, file=file.path(out_dir, "gp_pw_errs_post.rds"))
  saveRDS(gp_agg_errs_prior, file=file.path(out_dir, "gp_agg_errs_prior.rds"))
  saveRDS(gp_agg_errs_post, file=file.path(out_dir, "gp_agg_errs_post.rds"))
}

# ------------------------------------------------------------------------------
# Create replicate designs
# ------------------------------------------------------------------------------

print("-------------- Generating replicate designs --------------")

# Seeds for each replicate design.
seeds <- sample.int(n=.Machine$integer.max, size=cmd_args$n_rep)

# Generate replicate designs.
for(seed in seeds) construct_init_em(seed)



