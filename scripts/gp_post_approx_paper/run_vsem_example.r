#
# run_vsem_example.r
# Bayesian inversion for Very Simple Ecosystem Model (VSEM). This file 
# is written to be able to run on a computing cluster. 
#
# Andrew Roberts
#

# TODO:
# - Perhaps the best way to do this is to have a `run_vsem_example_setup.r` 
#   script that is responsible for saving the inverse problem list and 
#   generating validation data that will be used for testing all emulators 
#   (one set of validation data based on the prior and one set based on the 
#   exact MCMC output).
# - Then a `run_vsem_example_emulator.r` script that reads in the inverse
#   problem setup, creates an initial design, and then fits an emulator. This 
#   script can be called in a distributed fashion using different design seeds
#   to create an ensemble of initial designs and their corresponding emulators.
# - Then a `run_vsem_example_mcmc.r` script that reads in a specific MCMC alg
#   and design/emulator and runs MCMC. An MCMC settings object can be saved 
#   to file beforehand so this file can read settings from that file.
#   Results can be saved for each alg/design combination.
# - In addition, I can write a file that runs sequential design for a large 
#   number of iterations (replicated over each intitial design) using methods 
#   that do not require running MCMC after each iteration.

# TODO:
# - save samp_dt for exact mcmc separately 
# - update the function that runs a bunch of different mcmc algs; wrapper 
#   function should maybe be called something like "run_mcmc_llik_em()", which 
#   translates an mcmc tag into the algorithm to be called. Then can write a 
#   wrapper around this called `run_mcmc_llik_em_approx()` which basically just
#   loops over `run_mcmc_llik_em()`. Save `samp_dt` separately for each to 
#   avoid one massive data.table. 









# TODO:
#    - Need to check llikEmulatorGP and forward emulator and ensure that the 
#      normalization/conditional parameters are properly defined and resulting 
#      in the correct posterior approximations.
#    - Write a function that acts on `par_prior` and returns bounds for each parameter. 
#      Will be -Inf/Inf for unbounded. 
#    - Talk with Mike about reasonable priors to use in this experiment. 
#      e.g., VSEM default parameter range for LAR is enormous, and would
#      probably benefit from something like a log-normal/half-normal/Gamma 
#      prior rather than a uniform one. 
#    - Set up the settings here so that the VSEM parametrer names can be passed
#      in from the bash file; right now `dim_par` and `dim_output` do nothing. 
#    - Add hyperparameter constraints on other kergp models (not just Gaussian+Quadratic). 
#    - Add test points along coordinate axes.
#    - Add basis function forward model emulator.
#    - Figure out weird sampling behavior for approximate posteriors. 

# -----------------------------------------------------------------------------
# docopt string for parsing command line arguments.  
# -----------------------------------------------------------------------------

"Usage:
  test_docopt.r <run_id> <output_dir> [options]
  test_docopt.r (-h | --help)

Options:
  -h --help                                 Show this screen.
  --N_design=<N_design>                     Number of design points.
  --design_method=<design_method>           Algorithm to generate design points. 
  --N_design_test=<N_design_test>           Number of validation points for emulators. 
  --design_method_test=<design_method_test> Algorithm to generate validation points. 
  --mcmc_tags=<mcmc_tag,mcmc_tag>           GP-approx MCMC algorithms. 
  --N_mcmc=<N_mcmc>                         Number of MCMC iterations. 
" -> doc

# -----------------------------------------------------------------------------
# MCMC settings  
# -----------------------------------------------------------------------------

# Helper function will extract "tag" and "mcmc_func_name"; everything else 
# will be treated as an argument to `run_mcmc()` or the MCMC function.

# TODO: where to add cov_prop_init?
# TODO: add marginal acceptance prob.

common_settings <- list(n_itr=100000L, try_parallel=FALSE)
common_bt_settings <- list(mcmc_func_name="mcmc_bt_wrapper", sampler="DEzs", 
                           n_chain=4L, defer_ic=TRUE)

mcmc_settings_list <- list(
  list(mcmc_tag="quantile7", common_settings, common_bt_settings,
       approx_type="quantile", alpha=0.7),
  list(mcmc_tag="quantile8", common_settings, common_bt_settings,
       approx_type="quantile", alpha=0.8),
  list(mcmc_tag="quantile9", common_settings, common_bt_settings,
       approx_type="quantile", alpha=0.9),
  list(mcmc_tag="mean", common_settings, common_bt_settings,
       approx_type="mean"),
  list(mcmc_tag="marginal", common_settings, common_bt_settings,
       approx_type="marginal"),
  list(mcmc_tag="mcwmh-joint", common_settings, mcmc_func_name="mcmc_noisy_llik",
       use_joint=TRUE),
  list(mcmc_tag="mcwmh-ind", common_settings, mcmc_func_name="mcmc_noisy_llik",
       use_joint=FALSE)
)


# -----------------------------------------------------------------------------
# Setup 
# -----------------------------------------------------------------------------

print("---------------------------------- Setup ----------------------------------")

interactive_mode <- interactive()
print(paste0("Interactive mode: ", interactive_mode))

# Base directory (i.e., project directory). Required for loading the R project (renv). 
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")

# Load R project. 
renv::load(base_dir)
print(".libPaths()")
print(.libPaths())
renv::status()

library(lhs)
library(ggplot2)
library(viridis)
library(parallel)
library(gridExtra)
library(data.table)
library(docopt)

if(interactive()) {
  arguments <- list(par_cal_names=c("KEXT","LAR"), N_design=20,
                    N_design_test=200, design_method="LHS",
                    design_method_test="LHS", N_mcmc=100000L,
                    mcmc_tags=c("gp-mean", "gp-marg", "mcwmh-joint", "mcwmh-ind"),
                    run_id="test", output_dir="")
} else {
  # Read command line arguments. 
  arguments <- docopt(doc)
}
 
run_id <- arguments$run_id
output_dir <- arguments$output_dir
required_settings <- c("N_design", "N_design_test", "design_method",
                       "design_method_test", "N_mcmc", "mcmc_tags")
settings <- arguments[required_settings]

# Source and output paths.  
src_dir <- file.path(base_dir, "src")

# Source scripts. 
source(file.path(src_dir, "gp_helper_functions.r"))
source(file.path(src_dir, "gpWrapper.r"))
source(file.path(src_dir, "inv_prob_test_functions.r"))
source(file.path(src_dir, "general_helper_functions.r"))
source(file.path(src_dir, "statistical_helper_functions.r"))
source(file.path(src_dir, "plotting_helper_functions.r"))
source(file.path(src_dir, "llikEmulator.r"))
source(file.path(src_dir, "gp_emulator_functions.r"))
source(file.path(src_dir, "mcmc_helper_functions.r"))
source(file.path(src_dir, "sim_study_functions.r"))
source(file.path(src_dir, "mcmc_calibration_functions.r"))
source(file.path(src_dir, "seq_design.r"))
source(file.path(src_dir, "seq_design_gp.r"))
source(file.path(src_dir, "seq_design_for_post_approx.r"))
source(file.path(src_dir, "gp_mcmc_functions.r"))

# Seeds for random number generator for different portions of analysis. 
design_seed <- 26423

# Design settings for emulator construction. 
if(is.null(settings$N_design)) {
  settings$N_design <- 10L*settings$dim_par
} else {settings$N_design <- as.integer(settings$N_design)}

if(is.null(settings$design_method)) settings$design_method <- "LHS"

# Design settings for emulator validation. 
if(is.null(settings$N_design_test)) {
  settings$N_design_test <- min(100L*settings$dim_par, 500L)
} else {settings$N_design_test <- as.integer(settings$N_design_test)}

if(is.null(settings$design_method_test)) settings$design_method_test <- "LHS"

# Global likelihood normalization settings. 
default_conditional <- FALSE
default_normalize <- TRUE

# GP-accelerated MCMC algorithms to run. 
if(is.null(settings$mcmc_tags)) {
  settings$mcmc_tags <- c("gp-mean", "gp-marg", "gp-quantile", 
                          "mcwmh-joint", "mcwmh-ind")
} else {
  settings$mcmc_tags <- unlist(strsplit(settings$mcmc_tags, ","))
}

# Settings for MCMC runs.  
if(is.null(settings$N_mcmc)) {settings$N_mcmc <- 50000L
} else {settings$N_mcmc <- as.integer(settings$N_mcmc)}

# Create global variables with the setting values. 
missing_settings <- names(settings)[sapply(settings, is.null)]
if(length(missing_settings) > 0) stop("Missing settings: ", missing_settings)
N_design <- settings$N_design
design_method <- settings$design_method
N_design_test <- settings$N_design_test
design_method_test <- settings$design_method_test
mcmc_tags <- settings$mcmc_tags
N_mcmc <- settings$N_mcmc

# If `output_dir` already exists, append timestep and create new directory.
# Otherwise create the directory. 
if(dir.exists(output_dir)) {
  stop("`output_dir` already exists.")
} else {
  dir.create(output_dir, recursive=TRUE)
}
save(settings, file=file.path(output_dir, "settings.RData"))

# Print settings. 
print("--------------------- User specified settings -------------------------")

print("--> General settings")
print(paste0("run ID: ", run_id))
print(paste0("output directory: ", output_dir))
print(paste0("Design sampling seed (design_seed): ", design_seed))
print(paste0("Likelihood densities conditional: ", default_conditional))
print(paste0("Likelihood densities normalized: ", default_normalize))

print("--> Inverse problem")

print("--> Design")
print(paste0("Number design points: ", N_design))
print(paste0("Design method: ", design_method))
print(paste0("Number test points: ", N_design_test))
print(paste0("Test design method: ", design_method_test))

print("--> MCMC settings")
print(paste0("Number iterations: ", N_mcmc))
print(paste0("GP-accelerated MCMC algs: ", paste(mcmc_tags, collapse=", ")))


# -----------------------------------------------------------------------------
# Bayesian inverse problem setup. 
# -----------------------------------------------------------------------------

print("--------------------- Bayesian Inverse Problem Setup -------------------------")
print("Creating VSEM inverse problem.")
inv_prob <- get_vsem_test_1(default_conditional=default_conditional, 
                            default_normalize=default_normalize)
llik_obj_exact <- inv_prob$llik_obj
par_prior <- inv_prob$par_prior
save(inv_prob, file=file.path(output_dir, "inv_prob_list.RData"))

# -----------------------------------------------------------------------------
# Initial design and test points.  
# -----------------------------------------------------------------------------

print("--------------------- Constructing design -------------------------")

set.seed(design_seed)

# Generate training design. 
design_info <- get_init_design_list(inv_prob, design_method, N_design)
save(design_info, file=file.path(output_dir, "design_info.RData"))

# Generate validation design.
test_info <- get_init_design_list(inv_prob, design_method_test, N_design_test)
save(test_info, file=file.path(output_dir, "test_info.RData"))

# -----------------------------------------------------------------------------
# Exact MCMC. 
# -----------------------------------------------------------------------------

print("--------------------- Running exact MCMC ---------------------")

# Initial condition and covariance proposal for MCMC.
# TODO: could also use joint Gaussian approximation that Meng is working on to 
# set these quantities here. 
mcmc_par_init <- design_info$input[which.max(design_info$llik + design_info$lprior),]
cov_prop_init <- cov(design_info$input) # TODO: replace with IS estimate
print(paste0("MCMC initial value:", paste(mcmc_par_init, collapse=", ")))
print("Initial proposal covariance:")
print(cov_prop_init)

# MCMC sampling using exact likelihood.
mcmc_exact_list <- run_mcmc("mcmc_bt_wrapper", llik_obj_exact, defer_ic=TRUE,
                            sampler="DEzs", n_chain=4L, par_prior=par_prior, 
                            n_itr=N_mcmc, test_label="exact", try_parallel=FALSE)

# Table storing MCMC samples.
samp_dt_exact <- mcmc_exact_list$samp
mcmc_list_exact <- list(exact=mcmc_exact_list$output_list)

# Save exact MCMC samples.
fwrite(samp_dt_exact, file.path(output_dir, "mcmc_samp_exact.csv"))
save(mcmc_list_exact, file=file.path(output_dir, "mcmc_list_exact.RData"))

# -----------------------------------------------------------------------------
# Fitting Emulators.   
# -----------------------------------------------------------------------------

print("--------------------- Fitting Emulators -------------------------")

llik_em_list <- list()

# # Log-likelihood emulator, constant mean. 
# print("-> Fitting em_llik_const")
# em_llik_gp <- gpWrapperHet(design_info$input, matrix(design_info$llik, ncol=1), 
#                            normalize_output=TRUE, scale_input=TRUE)
# em_llik_gp$fit("Gaussian", "constant", estimate_nugget=FALSE)
# llik_em_list[["em_llik_const"]] <- llikEmulatorGP("em_llik_const", em_llik_gp, default_conditional=default_conditional, 
#                                                   default_normalize=default_normalize, lik_par=inv_prob$sig2_model, 
#                                                   use_fixed_lik_par=TRUE)
# 
# 
# # Log-likelihood emulator, quadratic mean.
# print("-> Fitting em_llik_quad")
# em_llik_gp2 <- gpWrapperKerGP(design_info$input, matrix(design_info$llik, ncol=1), 
#                               normalize_output=TRUE, scale_input=TRUE)
# em_llik_gp2$fit("Gaussian", "quadratic", estimate_nugget=FALSE, 
#                 optimFun="nloptr::nloptr", trace=TRUE, multistart=10)
# llik_em_list[["em_llik_quad"]] <- llikEmulatorGP("em_llik_quad", em_llik_gp2, default_conditional=default_conditional, 
#                                                  default_normalize=default_normalize, lik_par=inv_prob$sig2_model, 
#                                                  use_fixed_lik_par=TRUE)

# Log-likelihood emulator, constant mean, Gaussian plus quadratic kernel. 
print("-> Fitting em_llik_const_GaussQuad")
em_llik_gp3 <- gpWrapperKerGP(design_info$input, matrix(design_info$llik, ncol=1), 
                              normalize_output=TRUE, scale_input=TRUE)
em_llik_gp3$fit("Gaussian_plus_Quadratic", "constant", estimate_nugget=FALSE, 
                optimFun="nloptr::nloptr", trace=TRUE, multistart=10)
llik_em_list[["em_llik_const_GaussQuad"]] <- llikEmulatorGP("em_llik_const_GaussQuad", em_llik_gp3, 
                                                            default_conditional=default_conditional, 
                                                            default_normalize=default_normalize, 
                                                            lik_par=inv_prob$sig2_model, use_fixed_lik_par=TRUE)

# Save log-likelihood emulators. 
save(llik_em_list, file=file.path(output_dir, "llik_em_list.RData"))

# -----------------------------------------------------------------------------
# Log-Likelihood predictions at test points.   
# -----------------------------------------------------------------------------

# Emulator predictions. 
print("--------------------- Emulator Predictions at Test Points -------------------------")
emulator_pred_list <- list()
for(em_name in names(llik_em_list)) {
  print(paste0("Predicting with emulator: ", em_name))
  emulator_pred_list[[em_name]] <- llik_em_list[[em_name]]$predict_emulator(test_info$input)
}
save(emulator_pred_list, file=file.path(output_dir, "emulator_pred_list.RData"))

# -----------------------------------------------------------------------------
# Run MCMC Samplers.  
# -----------------------------------------------------------------------------

print("--------------------- GP-Accelerated MCMC -------------------------")

llik_em <- llik_em_list[[1]]
mcmc_mcwmh_list <- run_mcmc("mcmc_noisy_llik", llik_em, par_prior=par_prior,
                            n_chain=4L, n_itr=N_mcmc, test_label="mcwmh-joint",
                            try_parallel=FALSE, use_joint=TRUE, 
                            cov_prop=cov_prop_init)

samp_dt_approx <- combine_samp_dt(samp_dt, mcmc_mcwmh_list$samp)
mcmc_list_approx <- mcmc_mcwmh_list$output_list

fwrite(samp_dt_approx, file.path(output_dir, "mcmc_samp_approx.csv"))
save(mcmc_list_approx, file=file.path(output_dir, "mcmc_list_approx.RData"))

# run_approx_mcmc <- TRUE
# if(run_approx_mcmc) {
#   for(em_name in names(llik_em_list)) {
#     print(paste0("----- ", em_name, " -----"))
#     mcmc_info_list <- run_approx_mcmc_comparison(inv_prob_list=inv_prob, 
#                                                  llik_em_obj=llik_em_list[[em_name]], 
#                                                  mcmc_tags=mcmc_tags, save_dir=output_dir, 
#                                                  samp_dt=samp_dt, mcmc_list=mcmc_list, 
#                                                  test_label_suffix=em_name, N_itr=N_mcmc, 
#                                                  par_init=mcmc_par_init, cov_prop=cov_prop_init,
#                                                  overwrite=TRUE)
#     
#     samp_dt <- mcmc_info_list$samp
#     mcmc_list <- mcmc_info_list$mcmc_list
#   }
# }

print("-------------------------- End of Script ------------------------------")




