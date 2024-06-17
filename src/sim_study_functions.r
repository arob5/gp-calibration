# ---------------------------------------------------------------------------------------------------
# sim_study_functions.r
# Simulation study for comparing different approaches to (batch) sequential design with GP emulator. 
# Dependencies: mcmc_calibration_functions.r, gp_emulator_functions.r, 
#               sequential_design_optimiation.r
#
# Andrew Roberts
# ---------------------------------------------------------------------------------------------------

library(lhs)
library(hetGP)
library(mlegp)
library(RJSONIO)
library(ggplot2)
library(viridis)
library(gridExtra)
library(data.table)
library(BayesianTools)

# -----------------------------------------------------------------------------
# Generic functions for running simulation experiments. 
# -----------------------------------------------------------------------------

init_experiment <- function(experiment_id, experiment_type, parent_dir, global_seed, set_seed=FALSE, ...) {
  # The `...` should be used to pass in experiment type-specific settings. 
  
  print(paste0("Creating experiment: ", experiment_id))
  
  experiment_config <- list(base_settings=list(), base_paths=list())
  experiment_config$base_settings$global_seed <- global_seed
  
  # Set the experiment type.
  experiment_tag <- get_experiment_type_tag(experiment_type)
  experiment_config$base_settings$experiment_type <- experiment_type
  experiment_config$base_settings$experiment_type_tag <- experiment_tag
  
  # Create directory. If the run ID already exists, throw error. 
  experiment_dir <- file.path(parent_dir, experiment_id)
  if(file.exists(experiment_dir)) stop("Experiment directory ", experiment_dir, " already exists.")
  experiment_config$base_settings$experiment_id <- experiment_id
  experiment_config$base_paths$experiment_dir <- experiment_dir
  
  # Get list of experiment type-specific settings. 
  experiment_config$experiment_type_settings <- init_experiment_type(experiment_config, ...)
  
  # Create main experiment directory and subdirectories. 
  dir.create(experiment_dir)
  print(paste0("Created experiment directory: ", experiment_dir))
  settings_dir <- file.path(experiment_dir, "settings")
  output_dir <- file.path(experiment_dir, "output")
  code_dir <- file.path(experiment_dir, "code")
  dir.create(settings_dir)
  dir.create(output_dir)
  dir.create(code_dir)
  experiment_config$base_paths$settings_dir <- settings_dir
  experiment_config$base_paths$output_dir <- output_dir
  experiment_config$base_paths$code_dir <- code_dir
  
  # Save experiment config as JSON file. 
  experiment_config_json <- toJSON(experiment_config)
  experiment_config_path <- file.path(settings_dir, "experiment_config.json") 
  write(experiment_config_json, experiment_config_path)
  print(paste0("Saved experiment config settings: ", experiment_config_path))
  experiment_config$base_paths$experiment_config_path <- experiment_config_path
  
  # Optionally set global seed. 
  if(set_seed) {
    set.seed(experiment_config$global_seed)
    print(paste0("Global seed set to ", experiment_config$base_settings$global_seed))
  }
  
  return(experiment_config)
  
}


load_experiment <- function(experiment_id, parent_dir, config_filename="experiment_config.json", 
                            set_seed=FALSE, ...) {
                  
  # Check that the experiment exists.
  experiment_dir <- file.path(parent_dir, experiment_id)
  if(!file.exists(experiment_dir)) stop("Experiment ", experiment_dir, " not found.")
  
  # Check that experiment config exists. If so load and convert to R list. 
  experiment_config_path <- file.path(experiment_dir, "settings", config_filename)
  if(!file.exists(experiment_config_path)) {
    stop("Config for experiment ", experiment_id, " not found at path ", experiment_config_path)
  }
  
  experiment_config <- fromJSON(experiment_config_path, simplify=FALSE)
  
  # Validation checks for base settings must pass or error is thrown. 
  validate_experiment_base_settings(experiment_config)
  status <- list(base_settings_validated=TRUE)
  
  # Validation checks for core experiment type settings must similarly pass. 
  validate_experiment_type_core_settings(experiment_config)
  status$core_experiment_type_settings_validated <- TRUE
  
  # Specialized load function for the specific experiment type.  
  type_return_list <- load_experiment_type(experiment_config, status)
  status <- type_return_list$status
  
  # Optionally set global seed. 
  if(set_seed) {
    set.seed(experiment_config$global_seed)
    print(paste0("Global seed set to ", experiment_config$base_settings$global_seed))
  }
  
  return(list(config=experiment_config, status=status, obj=type_return_list$obj))
  
}

init_experiment_type <- function(config, ...) {
  tag <- config$base_settings$experiment_type_tag 
  
  if(tag != "none") {
    experiment_type_settings <- get(paste0("init_experiment_type_", tag))(config, ...)
    return(experiment_type_settings)
  }
  
  return(NULL)
}


load_experiment_type <- function(config, status) {
  # There is a minimal set of base settings required for all experiments. 
  # If any of these are missing an error is thrown. Beyond that, experiment 
  # status is defined and tracked differently for each experiment type. 
  
  # Validation checks for base settings must pass before moving on to specifics
  # for the experiment type. 
  assert_that(status$base_settings_validated, 
              msg="Base settings must pass validation before running `load_experiment_type()`.")

  # Load specific experiment type. 
  tag <- config$base_settings$experiment_type_tag 
  if(tag != "none") {
    type_load_list <- get(paste0("load_experiment_type_", tag))(config)
    status[[tag]] <- type_load_list$status
    experiment_type_objects <- type_load_list$obj
    return(list(status=status, obj=experiment_type_objects))
  }
  
  return(list(status=status, obj=NULL))
  
}


validate_experiment_base_settings <- function(experiment_config) {
  
  required_base_settings <- c("experiment_id", "experiment_type", "global_seed")
  missing_setting_names <- setdiff(required_base_settings, names(experiment_config$base_settings))
  
  # Ensure required settings are found as names in the list. 
  if(length(missing_setting_names) > 0) {
    stop("`experiment config missing required settings: ", missing_setting_names)
  }
  
  # Ensure the corresponding list elements are not missing. 
  setting_is_missing <- sapply(required_base_settings, 
                               function(setting_name) is.na(experiment_config$base_settings[[setting_name]]) ||
                                                      is.null(experiment_config$base_settings[[setting_name]]))
  if(any(setting_is_missing)) {
    stop("`experiment` config has NA or NULL settings: ", 
         paste(required_base_settings[setting_is_missing], sep=", "))
  }
  
  # Ensure the experiment type and the tag are a matching pair. 
  tag <- get_experiment_type_tag(experiment_config$base_settings$experiment_type)
  if(tag != experiment_config$base_settings$experiment_type_tag) {
    stop("`experiment_config$experiment_type_tag` == ",  experiment_config$experiment_type_tag, 
         " does not match expected tag ", tag)
  }
  
}


validate_experiment_type_core_settings <- function(config) {
  
  tag <- config$base_settings$experiment_type_tag 
  if(tag != "none") get(paste0("validate_experiment_type_core_settings_", tag))(config$experiment_type_settings)
}


get_experiment_type_tag <- function(experiment_type) {
  # Defines valid experiment types. The special type "none" means that the 
  # experiment only uses the base functions. Most experiment types have a 
  # specialized suite of functions, which are identified by the suffix 
  # "_<experiment_type_tag>". 
  
  # This defines the valid experiment types and their associated tags. 
  # The names and values give the type full name and short name (i.e., "tag), 
  # respectively. 
  experiment_type_mapping <- c(none="none", llik_emulator_seq_design="lesd")
  
  if(!(experiment_type %in% names(experiment_type_mapping))) {
    stop("Experiment type ", experiment_type, " not recognized.")
  }
  
  return(unname(experiment_type_mapping[experiment_type]))
  
}

 
# -----------------------------------------------------------------------------
# Experiment Type: llik_emulator_seq_design (lesd, for short). 
#
# Simulation experiment functions for testing log-likelihood emulators
# for accelerated MCMC in Bayesian inverse problems, and sequential design 
# methods for improving the emulators. This experiment type fixes a specific 
# Bayesian inverse problem (likelihood, data, forward model, priors), so 
# the experiments concern varying the emulation methods, sequential design 
# algorithms, and/or MCMC algorithms. 
#
# This experiment type features multiple "rounds" corresponding to rounds 
# of sequential design. Round 0 is the "initial design", round 1 is the 
# first round of (batch) sequential design, and so on. Each round features 
# 3 "stages": setup, design, and sample. Setup involves defining config 
# parameters to specify which methods to test, how many replication to 
# run, etc. The design stage runs algorithms which acquire a batch of 
# design points in parameter space, which typically requires outputs from 
# previous rounds as input. The sample stage runs MCMC algorithms using
# the log-likelihood emulators which have been updated in the previous 
# design stage.
# -----------------------------------------------------------------------------

# Required settings which are fixed across rounds: 
# define_Bayesian_inverse_problem.r: 
#   Script which first initializes the experiment and then creates the 
#   llikEmulatorExactGauss object (storing the likelihood and forward model), 
#   and the data.frames defining the priors (have option to fix Sig_eps, in 
#   which case a prior is not required).
#
#   The global seed will be set anytime an experiment is loaded (maybe this 
#   isn't a good idea). But I should definitely make use of local seeds which 
#   are set before steps involving randomization. For example, lesd should 
#   require the config to contain a data seed which is set at the top of 
#   llikEmulatorExactGauss, ensuring that the same data Y is generated 
#   anytime this file is run. Will have to think about how to handle 
#   the seeds for the design and sampling stages as well; e.g. best to 
#   have a fixed seed for each of these that is the same across all rounds, 
#   or have round specific seeds for each? 

init_experiment_type_lesd <- function(config, inverse_problem_seed, fix_lik_par, ...) {
  # `inverse_problem_seed` is the integer random seed which is set in prior to calling 
  # `define_Bayesian_inverse_problem()` which ensures the synthetic data generation, etc 
  #  generates the same exact inverse problem each time the experiment is loaded. 
  #  `fix_lik_par` is TRUE/FALSE, determines whether the likelihood parameters 
  #  (e.g. observation covariance) should be fixed or learned. 
  
  settings <- list()
  assert_that((inverse_problem_seed %% 1) == 0, msg="`inverse_problem_seed` must be integer") # Ensure this is an integer
  settings$inverse_problem_seed <- inverse_problem_seed
  
  assert_that(is.logical(fix_lik_par))
  settings$fix_lik_par <- fix_lik_par
  
  return(settings)
  
}


validate_experiment_type_core_settings_lesd <- function(exp_type_config) {
  
  required_core_settings <- c("inverse_problem_seed", "fix_lik_par")
  missing_setting_names <- setdiff(required_core_settings, names(exp_type_config))
  
  # Ensure required settings are found as names in the list. 
  if(length(missing_setting_names) > 0) {
    stop("`experiment type config missing required settings: ", missing_setting_names)
  }
  
  # Ensure the corresponding list elements are not missing. 
  setting_is_missing <- sapply(required_core_settings, 
                               function(setting_name) is.na(exp_type_config[[setting_name]]) ||
                                                      is.null(exp_type_config[[setting_name]]))
  if(any(setting_is_missing)) {
    stop("experiment type config has NA or NULL settings: ", 
         paste(required_core_settings[setting_is_missing], sep=", "))
  }
  
  assert_that((exp_type_config$inverse_problem_seed %% 1) == 0, 
              msg="`inverse_problem_seed` must be integer") # Ensure this is an integer
  assert_that(is.logical(exp_type_config$fix_lik_par))
  
}


load_experiment_type_lesd <- function(config) {
  
  status <- list()
  
  #
  # Setup phase: checking the files that will be fixed throughout the entire experiment 
  # and must be correctly specified before proceeding with the experiment. 
  #
  
  status$setup <- list()
  
  # Function defining Bayesian inverse problem setup. 
  status$setup$inverse_problem <- list(file_exists=FALSE, file_validated=FALSE)
  inv_prob_script_path <- file.path(config$base_paths["code_dir"], "define_Bayesian_inverse_problem.r")
  inv_prob_file_exists <- file.exists(inv_prob_script_path)
  status$setup$inverse_problem$file_exists <- inv_prob_file_exists
  if(!inv_prob_file_exists) return(list(status=status)) 
  
  # Source file, then run function to load the objects defining the Bayesian inverse problem. 
  source(inv_prob_script_path)
  func_name <- "define_Bayesian_inverse_problem"
  inv_prob_func_exists <- exists(func_name, envir=.GlobalEnv)
  if(!inv_prob_func_exists) return(list(status=status))
  inv_prob_list <- get(func_name)(config)
  status$setup$inverse_problem$file_validated <- ensure_valid_inv_prob_list_lesd(inv_prob_list, 
                                                                                 config$experiment_type_settings$fix_lik_par)

  # File defining valid sequential design algorithms. 
  # File defining valid sampling algorithms. 
  # File defining valid likelihood emulation methods. 

  return(list(status=status, obj=inv_prob_list))
  
}


ensure_valid_inv_prob_list_lesd <- function(inv_prob_list, fix_lik_par) {
  # TODO: should update this to provide more info on what specifically is not 
  # satisfying the requirements. There are also currently no constraints on 
  # `lik_par_prior_params`
  
  required_elements <- c("llik_exact", "par_prior_params", "ground_truth")
  
  # Ensure it is a list with correct names. 
  if(!is.list(inv_prob_list)) return(FALSE)
  if(length(setdiff(required_elements, names(inv_prob_list))) > 0) return(FALSE)
  if(!fix_lik_par && !("lik_par_prior_params" %in% names(inv_prob_list))) return(FALSE)
  
  # Ensure `llik_exact` is a reference class which inherits from class `llikEmulator` and 
  # has attribute `exact_llik == TRUE`. 
  if(!inherits(inv_prob_list$llik_exact, "llikEmulator")) return(FALSE)
  if(!inv_prob_list$llik_exact$exact_llik) return(FALSE)
  
  # Ensure prior are provided for each calibration parameter. 
  if(!isTRUE(setequal(inv_prob_list$llik_exact$input_names, rownames(inv_prob_list$par_prior_params)))) return(FALSE)
  
  # Ensure ground truth is a list with the true parameters. 
  if(!is.list(inv_prob_list$ground_truth)) return(FALSE)
  if(length(setdiff(c("par_true", "lik_par_true"), names(inv_prob_list$ground_truth))) > 0) return(FALSE)
  
  
  return(TRUE)
}


# -----------------------------------------------------------------------------
# Convenience functions for testing different emualtor/MCMC methods to 
# approximately solve a Bayesian inverse problem. 
# -----------------------------------------------------------------------------

run_approx_mcmc_comparison <- function(inv_prob_list, llik_em_obj, mcmc_tags, save_dir=NULL, 
                                       samp_dt=NULL, mcmc_list=NULL, test_label_suffix=NULL, ...) {
  # Runs different sampling algorithms to try to approximate the posterior of the 
  # parameters in the inverse problem defined by the `inv_prob_list` object.
  # Saves a csv file containing the samples. 
  # TODO: Need to generalize to allow sampling likelihood parameters. 
  #
  # Args:
  #    inv_prob_list: list, with elements "fwd", "par_prior", "llik_obj".
  #    llik_em_obj: Object of class that inherits from `llikEmulator`. This 
  #                 log-likelihood emulator is passed to all of the approximate 
  #                 GP-accelerated MCMC algorithms. 
  #    mcmc_tags: character, vector of MCMC valid test labels that determine 
  #               which MCMC algorithms will be run. Valid options include 
  #               "gp-mean", "gp-marg", "mcwmh-joint", "mcwmh-ind", "pseudo-marg", 
  #               "acc-prob-marg". 
  #    samp_dt: data.table, with colnames "test_label", "param_type", "param_name",   
  #             "itr", and "sample". If non-NULL, then this function will append
  #              to this data.table. Otherwise a new one will be
  #             created. 
  #    mcmc_list: list, an optional list storing MCMC information returned by 
  #               the MCMC functions other than the samples themselves (e.g., 
  #               proposal covariance). If provided, will be appended to. Otherwise
  #               a new one will be created. 
  #    test_label_suffix: character(1), an optional string to append to the end  
  #                       of each MCMC tag when defining the test labels. A 
  #                       hyphen will automatically be added between the MCMC
  #                       tag and test label. Otherwise the test labels will be set 
  #                       to the MCMC tags.
  #    ...: additional arguments will be passed to the MCMC functions; e.g., "par_init", 
  #         "N_itr", proposal adaptation parameters, etc. 
  #
  # Returns: 
  #    Invisibly returns a list with named elements: 
  #        - "samp": data.table containing the MCMC samples.
  #        - "mcmc_list": list containing MCMC output other than samples.  
  #    Side effect: If `save_dir` specifies a valid directory, then 
  #                 saves the data.table as a csv file and the list as a RData file.  
  
  # Create new data.table for samples if an existing one is not provided. 
  if(is.null(samp_dt)) {
    samp_dt <- data.table(test_label=character(), param_type=character(), 
                          param_name=character(), itr=integer(), sample=numeric())
  }
  if(is.null(mcmc_list)) {
    mcmc_list <- list()
  }
  
  # Run MCMC algorithms. 
  valid_mcmc_tags <- c("gp-mean", "gp-marg", "mcwmh-joint", "mcwmh-ind", "pseudo-marg", "acc-prob-marg")
  invalid_tags <- setdiff(mcmc_tags, valid_mcmc_tags)
  if(length(invalid_tags) > 0) {
    message("Unsupported MCMC tags will not be run: ", invalid_tags)
  }
  
  # Add underscore before the suffix. 
  if(!is.null(test_label_suffix)) test_label_suffix <- paste0("-", test_label_suffix)
  
  # Output directory: Create if it doesn't exist. If files already exist, append timestamp 
  # to filenames to avoid overwriting. 
  samp_filename <- "mcmc_samp"
  list_filename <- "mcmc_list"
  if(!is.null(save_dir)) {
    if(!dir.exists(save_dir)) dir.create(save_dir)
    timestamp <- as.character(Sys.time())
    if(file.exists(file.path(save_dir, paste0(samp_filename, ".csv")))) samp_filename <- paste(samp_filename, timestamp, sep="_")
    if(file.exists(file.path(save_dir, paste0(list_filename, ".csv")))) list_filename <- paste(list_filename, timestamp, sep="_")
  }
  samp_filename <- paste0(samp_filename, ".csv")
  list_filename <- paste0(list_filename, ".RData")
  
  # Run MCMC algorithms. 
  if("gp-mean" %in% mcmc_tags) {
    mcmc_output <- mcmc_gp_unn_post_dens_approx(llik_emulator=llik_em_obj, par_prior_params=inv_prob_list$par_prior, 
                                                approx_type="mean", ...)
    lbl <- paste0("gp-mean", test_label_suffix)
    samp_dt <- append_mcmc_output(samp_dt, mcmc_output$samp, test_label=lbl)
    mcmc_list[[lbl]] <- mcmc_output[setdiff(colnames(mcmc_output), "samp")]
  }
  
  if("gp-marg" %in% mcmc_tags) {
    mcmc_output <- mcmc_gp_unn_post_dens_approx(llik_emulator=llik_em_obj, par_prior_params=inv_prob_list$par_prior, 
                                                approx_type="marginal", ...)
    lbl <- paste0("gp-marg", test_label_suffix)
    samp_dt <- append_mcmc_output(samp_dt, mcmc_output$samp, test_label=lbl)
    mcmc_list[[lbl]] <- mcmc_output[setdiff(colnames(mcmc_output), "samp")]
  }
  
  if("mcwmh-joint" %in% mcmc_tags) {
    mcmc_output <- mcmc_gp_noisy(llik_emulator=llik_em_obj, par_prior_params=inv_prob$par_prior, 
                                 mode="MCMH", use_gp_cov=TRUE, ...)
                                 
    lbl <- paste0("mcwmh-joint", test_label_suffix)
    samp_dt <- append_mcmc_output(samp_dt, mcmc_output$samp, test_label=lbl)
    mcmc_list[[lbl]] <- mcmc_output[setdiff(colnames(mcmc_output), "samp")]
  }

  if("mcwmh-ind" %in% mcmc_tags) {
    mcmc_output <- mcmc_gp_noisy(llik_emulator=llik_em_obj, par_prior_params=inv_prob$par_prior, 
                                 mode="MCMH", use_gp_cov=FALSE, ...)
    
    lbl <- paste0("mcwmh-ind", test_label_suffix)
    samp_dt <- append_mcmc_output(samp_dt, mcmc_output$samp, test_label=lbl)
    mcmc_list[[lbl]] <- mcmc_output[setdiff(colnames(mcmc_output), "samp")]
  }
  
  if("pseudo-marg" %in% mcmc_tags) {
    mcmc_output <- mcmc_gp_noisy(llik_emulator=llik_em_obj, par_prior_params=inv_prob$par_prior, 
                                 mode="pseudo-marg", ...)
    lbl <- paste0("pseudo-marg", test_label_suffix)
    samp_dt <- append_mcmc_output(samp_dt, mcmc_output$samp, test_label=lbl)
    mcmc_list[[lbl]] <- mcmc_output[setdiff(colnames(mcmc_output), "samp")]
  }
  
  if("pseudo-marg" %in% mcmc_tags) {
    mcmc_output <- mcmc_gp_noisy(llik_emulator=llik_em_obj, par_prior_params=inv_prob$par_prior, 
                                 mode="pseudo-marg", ...)
    lbl <- paste0("pseudo-marg", test_label_suffix)
    samp_dt <- append_mcmc_output(samp_dt, mcmc_output$samp, test_label=lbl)
    mcmc_list[[lbl]] <- mcmc_output[setdiff(colnames(mcmc_output), "samp")]
  }
  
  
  
}





















# -----------------------------------------------------------------------------
# Setup:
#    - Specify forward model (i.e., computer model or simulator). 
#    - Set true values of calibration parameters. 
#    - Set true values of likelihood parameters. 
#    - Set prior distributions for calibration and likelihood parameters. 
# -----------------------------------------------------------------------------

# # Computer model and synthetic data generation. 
# computer_model_data <- generate_vsem_test_case(4)
# print(computer_model_data$ref_pars[computer_model_data$pars_cal_sel,])
# 
# # Priors on calibration parameters. 
# theta_prior_params <- computer_model_data$ref_pars[computer_model_data$pars_cal_sel,]
# theta_prior_params[, "dist"] <- c("Uniform", "Uniform")
# theta_prior_params[,"param1"] <- c(1.3, 0.4) 
# theta_prior_params[,"param2"] <- c(1.7, 0.6)
# theta_prior_params <- theta_prior_params[, c("dist", "param1", "param2")]
# print(paste0(rep("-", 30), " Calibration parameters priors ", rep("-", 30)))
# print(theta_prior_params)
# 
# # Priors on likelihood parameters. 
# sig2_eps_prior_info <- get_IG_priors_numerical_test(sig2_true = diag(computer_model_data$Sig_eps), 
#                                                     bias_frac = c(0.1, -0.15), bins = 50,
#                                                     coef_var = c(0.3, 0.5), return_prior_plots = TRUE, 
#                                                     output_variables = computer_model_data$output_vars)
# sig2_eps_prior_params <- sig2_eps_prior_info$prior
# print(paste0(rep("-", 30), " Likelihood parameters priors ", rep("-", 30)))
# print(sig2_eps_prior_params)


# -----------------------------------------------------------------------------
# GP-Accelerated MCMC Helper Functions
# -----------------------------------------------------------------------------

get_1d_linear_Gaussian_approx_post_density <- function(data_seed, design_seed, G, sig2_eps, Sig0_theta, N_design,  
                                                       emulator_settings, N_obs=100, design_method="LHS") {
  
  # List to store outputs. 
  plot_list <- list()
  
  # -----------------------------
  # Linear Gaussian Model Setup. 
  # -----------------------------
  
  linear_Gaussian_info <- generate_linear_Gaussian_test_data(data_seed, N_obs=N_obs, 
                                                             D=1, Sig_theta=Sig0_theta, 
                                                             G=G, sig2_eps=sig2_eps)
  computer_model_data <- linear_Gaussian_info$computer_model_data
  theta_prior_params <- linear_Gaussian_info$theta_prior_params
  linear_Gaussian_info$true_posterior$SSR <- get_computer_model_SSR(computer_model_data, 
                                                                    theta_vals=linear_Gaussian_info$true_posterior$mean, 
                                                                    na.rm=TRUE)
  
  # Plot simulated data/ground truth from linear Gaussian model. 
  sim_data <- as.data.frame(computer_model_data$data_obs)
  sim_data$ref <- computer_model_data$data_ref
  sim_data$time <- 1:N_obs
  plot_list$sim_data <- ggplot(sim_data) + geom_point(mapping=aes(x=time, y=y)) + 
                        geom_line(mapping=aes(x=time, y=ref), color="red") + 
                        xlab("t") + ylab("output") + ggtitle("Ground Truth and Observed Data")
                   
  # -----------------------------------
  # Design points and plot grid points  
  # -----------------------------------
  
  # Define input bounds at 1st and 99th percentiles of Gaussian prior. Used for grid bounds 
  # for plotting, as well as to truncate the prior for the GP-accelerated MCMC algorithms. 
  input_bounds <- matrix(c(min=qnorm(.01, theta_prior_params$param1, theta_prior_params$param2), 
                           max=qnorm(.99, theta_prior_params$param1, theta_prior_params$param2)), ncol=1)
  rownames(input_bounds) <- c("min", "max")
  colnames(input_bounds) <- computer_model_data$pars_cal_names
  
  # The GP-accelerated MCMC algorithms utilize the truncated prior. 
  theta_prior_params_trunc <- truncate_prior_theta(theta_prior_params, input_bounds)
  
  # Initial Design. 
  design_settings <- data.frame(N_design=N_design, design_method=design_method, design_seed=design_seed)
  init_design_info <- get_input_output_design(N_points = design_settings$N_design,
                                              design_method = design_settings$design_method, 
                                              scale_inputs = TRUE,
                                              param_ranges = input_bounds,  
                                              computer_model_data = computer_model_data, 
                                              theta_prior_params = theta_prior_params)
  init_design_info$lpost <- calc_lpost_theta_product_lik(computer_model_data = computer_model_data, 
                                                         theta_vals = init_design_info$inputs, 
                                                         vars_obs = sig2_eps, 
                                                         SSR = init_design_info$outputs,
                                                         na.rm=TRUE, theta_prior_params=theta_prior_params, 
                                                         return_list=FALSE)
  
  # Grid of points spread across prior (for plotting).  
  prior_grid_info <- get_input_output_design(N_points = 1000,
                                             design_method = "grid", 
                                             scale_inputs = TRUE,
                                             param_ranges = init_design_info$input_bounds,
                                             computer_model_data = computer_model_data, 
                                             theta_prior_params = theta_prior_params, 
                                             design_seed = design_settings$design_seed)
  prior_grid_info$lpost <- calc_lpost_theta_product_lik(computer_model_data = computer_model_data, 
                                                        theta_vals = prior_grid_info$inputs, 
                                                        vars_obs = diag(computer_model_data$Sig_eps), 
                                                        SSR = prior_grid_info$outputs,
                                                        na.rm = TRUE, theta_prior_params=theta_prior_params, 
                                                        return_list = FALSE)
  
  # Plot true posterior and design points. 
  prior_grid_df <- data.frame(input=prior_grid_info$inputs[,1], lpost=prior_grid_info$lpost)
  design_df <- data.frame(input=init_design_info$inputs[,1], lpost=init_design_info$lpost)
  plot_list$true_post <- ggplot() + geom_line(aes(x=input, y=lpost), prior_grid_df, color="red") + 
                         geom_point(aes(x=input, y=lpost), design_df) + 
                         xlab("u") + ylab("Unnormalized log posterior") + 
                         ggtitle("Prior Design and True lpost values")

  
  # -----------------------------------
  # Fit Emulator   
  # -----------------------------------
  
  # Fit emulators on initial design. 
  gp_fits <- fit_independent_GPs(X_train=init_design_info$inputs_scaled, Y_train=init_design_info$outputs_normalized, 
                                 gp_lib=emulator_settings$gp_lib, gp_kernel=emulator_settings$kernel)$fits
  emulator_info_list <- list(gp_fits=gp_fits, input_bounds=init_design_info$input_bounds, 
                             output_stats=init_design_info$output_stats, settings=emulator_settings)
  
  # Induced log joint density (i.e. log unnormalized posterior density) emulator. 
  lpost_emulator <- get_lpost_emulator_obj(emulator_info_list=emulator_info_list, design_info_list=init_design_info, 
                                           computer_model_data=computer_model_data, sig2_eps=sig2_eps, 
                                           theta_prior_params=theta_prior_params_trunc, center_output=TRUE, scale_output=TRUE)
  
  # Plot initial emulator fit.
  pred_grid <- predict_lpost_emulator(inputs_new_scaled=prior_grid_info$inputs_scaled, lpost_emulator=lpost_emulator, 
                                      unscale=TRUE, uncenter=TRUE)
  plot_list$gp_fit <- plot_gp_fit_1d(prior_grid_info$inputs, prior_grid_info$lpost, init_design_info$inputs,
                                     init_design_info$lpost, pred_grid$mean, pred_grid$var, 
                                     vertical_line=linear_Gaussian_info$true_posterior$mean,
                                     xlab="u", ylab="Unnormalized log posterior", main_title="GP lpost Predictions", 
                                     CI_prob=0.95)
  
  # ------------------------------------
  # Approximate (log) posterior density    
  # ------------------------------------
  
  # Joint (u, phi) posterior density. 
  phi_grid <- seq(min(prior_grid_info$outputs), max(prior_grid_info$outputs), length.out=50)
  u_phi_grid <- expand.grid(prior_grid_info$inputs, phi_grid)
  colnames(u_phi_grid) <- c("u", "phi")
  
  lprior_u <- calc_lprior_theta(matrix(u_phi_grid$u, ncol=1), theta_prior_params_trunc)
  gp_pred_list <- predict_independent_GPs(X_pred=prior_grid_info$inputs_scaled, gp_obj_list=emulator_info_list$gp_fits, 
                                          gp_lib=emulator_settings$gp_lib, 
                                          denormalize_predictions=TRUE, output_stats=emulator_info_list$output_stats)
  gp_pred_means <- sapply(gp_pred_list, function(x) x$mean)
  gp_pred_vars <- sapply(gp_pred_list, function(x) x$var_comb)
  phi_pred_mean_grid <- expand.grid(gp_pred_means, phi_grid)
  phi_pred_var_grid <- expand.grid(gp_pred_vars, phi_grid)
  
  u_phi_grid$gp_pred_mean <- phi_pred_mean_grid$Var1
  u_phi_grid$gp_pred_var <- phi_pred_mean_grid$Var2
  u_phi_grid$a <- u_phi_grid$gp_pred_mean - 0.5 * (u_phi_grid$gp_pred_var/sig2_eps)
  
  u_phi_grid$lpost_u_phi <- lprior_u - 0.5*log(u_phi_grid$gp_pred_var) - 0.5*u_phi_grid$gp_pred_mean^2/u_phi_grid$gp_pred_var + 
                            0.5*u_phi_grid$a^2/u_phi_grid$gp_pred_var - 0.5*(u_phi_grid$phi-u_phi_grid$a)^2/u_phi_grid$gp_pred_var
  
  # Marginal u posterior density (phi marginalized out). 
  a <- gp_pred_means - 0.5 * gp_pred_vars/sig2_eps
  lpost_marg_u <- calc_lprior_theta(prior_grid_info$inputs, theta_prior_params_trunc) - 0.5*(gp_pred_means^2 - a^2)/gp_pred_vars
  
  # Plot the densities. 
  plot_list$u_phi_lpost <- get_2d_heatmap_plot(X=u_phi_grid[,c("u","phi")], y=u_phi_grid$lpost_u_phi, raster=TRUE, 
                                               param_names=c("u", "phi"), bigger_is_better=TRUE, main_title="(u,phi) posterior log density", 
                                               legend="Log Post (u,phi) Density", 
                                               point_coords=c(linear_Gaussian_info$true_posterior$mean, linear_Gaussian_info$true_posterior$SSR))
  
  u_marg_post_df <- data.frame(input=prior_grid_info$inputs[,1], lpost=lpost_marg_u[,1])
  plot_list$u_marg_lpost <- ggplot(u_marg_post_df, aes(x=input, y=lpost)) + geom_line() + 
                            xlab("u") + ylab("log E_{phi}[p(u|Sigma,Y)]") + ggtitle("Log Post u Marginal Density")
  
  # Collect objects to return. 
  obj_list <-   list(computer_model_data=computer_model_data, linear_Gaussian_info=linear_Gaussian_info, 
                     prior_params=theta_prior_params, prior_params_trunc=theta_prior_params_trunc, 
                     init_design_info=init_design_info, prior_grid_info=prior_grid_info, 
                     lpost_emulator=lpost_emulator, gp_pred_grid_list=gp_pred_list)
  
  return(list(obj=obj_list, plots=plot_list))

}


run_gp_mcmc_tests <- function(run_settings_list, computer_model_data=NULL, theta_prior_params=NULL,   
                              emulator_info_list=NULL, theta_init=NULL, N_chain=4, N_itr=2000, burn_ins=0.5*N_itr,  
                              learn_sig_eps=FALSE, return_cov_prop_scale=TRUE, return_SSR_samp=TRUE, 
                              cov_prop_init=NULL, ...) {
  # Currently this assumed fixed sig2_eps. 
  #
  # Args:
  #    computer_model_data: list, the computer model data list. 
  #    theta_prior_params: list, containing prior info. 
  #    emulator_info_list: list, the emulator info list.  
  #    TODO: run_settings_list: list
  #    algs: character, vector of names of MCMC algorithms to run. Valid options 
  #          are "ind_gp_gibbs", "ind_gp_marg", "ind_gp_joint", "ind_gp_trajectory". 
  #    ...: other named arguments are passed to the MCMC functions (e.g. arguments
  #         to control proposal covariance adaptation, etc.). 
  #
  
  # Parameters whose samples will be returned. 
  param_types <- "theta"
  if(learn_sig_eps) param_types <- c(param_types, "sig_eps")
  if(return_cov_prop_scale) param_types <- c(param_types, "cov_prop_scale")
  if(return_SSR_samp) param_types <- c(param_types, "SSR")
  
  # Set burn-ins. 
  test_labels <- sapply(run_settings_list, function(run) run$test_label)
  if(is.null(burn_ins)) burn_ins <- rep(1, length(test_labels))
  else if((length(burn_ins)==1) && length(test_labels)>1) burn_ins <- rep(burn_ins, length(test_labels))
  burn_ins <- setNames(burn_ins, test_labels)


  # Run MCMC algorithms.
  for(j in seq_along(test_labels)) {
    # Get MCMC function.
    test_label <- test_labels[j]
    alg_name <- run_settings_list[[j]]$alg
    mcmc_func <- get(paste0("mcmc_calibrate_", alg_name))
    
    # Get arguments to MCMC function. 
    mcmc_args <- get_mcmc_func_args_list(run_settings_list[[j]], computer_model_data, 
                                         theta_prior_params, emulator_info_list, theta_init, N_itr,
                                         learn_sig_eps, sig2_eps_init, cov_prop_init)
    mcmc_output <- do.call(mcmc_func, mcmc_args)

    # Format MCMC output. 
    col_sel <- intersect(param_types, names(mcmc_output))
    mcmc_samp_dt_alg <- format_mcmc_output(samp_list=mcmc_output[col_sel], test_label=test_label)
    mcmc_samp_dt_alg <- select_mcmc_samp(mcmc_samp_dt_alg, burn_in_start=burn_ins[test_label])
    
    if(j == 1) mcmc_samp_dt <- copy(mcmc_samp_dt_alg)
    else mcmc_samp_dt <- rbindlist(list(mcmc_samp_dt, mcmc_samp_dt_alg), use.names=TRUE)
                   
  }
  
  return(list(mcmc_samp_dt=mcmc_samp_dt, burn_ins=burn_ins))
  
}


get_mcmc_func_args_list <- function(run_settings_list, computer_model_data=NULL, theta_prior_params=NULL, 
                                    emulator_info_list=NULL, theta_init=NULL, N_itr=NULL, learn_sig_eps=NULL, 
                                    sig2_eps_init=NULL, cov_prop_init=NULL) {
  # A helper function to `run_gp_mcmc_tests()`, which returns a named list of MCMC function arguments 
  # for a specific GP-MCMC run. The run-specific MCMC settings in `run_settings_list` take 
  # precedent, but if a required argument is missing in this list, it then populated by one of the 
  # other arguments of this function, e.g. `computer_model_data`. Required arguments for all 
  # GP-MCMC functions are "computer_model_data", "theta_prior_params", "emulator_info_list". If 
  # "learn_sig_eps" is FALSE, then "sig2_eps_init" is also required. Optional arguments that are 
  # not provided here will just use the defaults for the GP-MCMC function in question. Certain 
  # GP-MCMC functions have their own particular required arguments are well. Additional elements 
  # of the list `run_settings_list` are allowed, and will be ignored by the GP-MCMC functions. 
  #
  # Args:
  #    run_settings_list: list, containing run-specific settings and arguments for a MCMC run. 
  #    Remaining arguments are used to populate settings in `run_settings_list` which are NULL. 
  #    The main use case here envisions many different runs, each with their own specific settings. 
  #    Settings that are fixed across all runs can then be provided in these additional arguments
  #    so they don't need to be set for each specific run. 
  #
  # Returns:
  #    list, a modified version of `run_settings_list`, where missing elements have been
  #    populated by the other arguments passed to this function. 
  
  # Set global arguments which have not been specified in the algorithm-specific settings. 
  if(is.null(run_settings_list$computer_model_data)) run_settings_list$computer_model_data <- computer_model_data
  if(is.null(run_settings_list$theta_prior_params)) run_settings_list$theta_prior_params <- theta_prior_params
  if(is.null(run_settings_list$emulator_info_list)) run_settings_list$emulator_info_list <- emulator_info_list
  if(is.null(run_settings_list$theta_init)) run_settings_list$theta_init <- theta_init
  if(is.null(run_settings_list$N_itr)) run_settings_list$N_itr <- N_itr
  if(is.null(run_settings_list$learn_sig_eps)) run_settings_list$learn_sig_eps <- learn_sig_eps
  if(is.null(run_settings_list$sig2_eps_init)) run_settings_list$sig2_eps_init <- sig2_eps_init
  if(is.null(run_settings_list$cov_prop_init)) run_settings_list$cov_prop_init <- cov_prop_init
  
  
  # Ensure required arguments are present. 
  if(is.null(run_settings_list$computer_model_data)) stop("`run_settings_list` is missing `computer_model_data`")
  if(is.null(run_settings_list$theta_prior_params)) stop("`run_settings_list` is missing `theta_prior_params`")
  if(is.null(run_settings_list$emulator_info_list)) stop("`run_settings_list` is missing `emulator_info_list`")
  if(!isTRUE(run_settings_list$learn_sig_eps) && is.null(run_settings_list$sig2_eps_init))
    stop("`run_settings_list` requires element `sig2_eps_init` when `learn_sig_eps` is not TRUE.")
  
  return(run_settings_list)
  
}


# -----------------------------------------------------------------------------
# Saving/loading tests results.  
# -----------------------------------------------------------------------------

create_mcmc_run <- function(settings, base_run_dir="output", set_global_variables=FALSE, set_seed=TRUE) {
  # TODO: add argument checking to ensure all required arguments are present. 
  
  # Extract Git hash. 
  settings$git_hash <- system("git rev-parse HEAD", intern=TRUE)
  
  # Validate settings. 
  validate_mcmc_run_settings(settings)
  
  # Create directory. If the run ID already exists, throw error. 
  run_id <- settings$run_id
  if(is.null(run_id)) stop("`run_id` missing from settings.")
  run_id_dir <- file.path(base_run_dir, run_id)
  if(file.exists(run_id_dir)) stop("Run ID path ", run_id_dir, " already exists.")
  dir.create(run_id_dir)
  
  # Create sub-directories for plots and analysis. 
  plot_dir <- file.path(run_id_dir, "plots")
  analysis_dir <- file.path(run_id_dir, "analysis")
  dir.create(plot_dir)
  dir.create(analysis_dir)
  
  # Add paths to settings. 
  settings$run_id_dir <- run_id_dir
  settings$analysis_dir <- analysis_dir
  settings$plot_dir <- plot_dir
  
  # Save top-level settings as JSON file. 
  settings_json <- toJSON(settings)
  write(settings_json, file.path(run_id_dir, "settings.json"))
  
  # Optionally set global seed. 
  if(set_seed) {
    set.seed(settings$global_seed)
    print(paste0("Global seed set to ", settings$global_seed))
  }
  
  # Optionally store settings as global variables. 
  if(set_global_variables) {
    convert_list_to_global_vars(settings)
    print("`settings` elements stored as global variables.")
  }
  
  return(invisible(settings))
  
}


load_mcmc_run_data <- function(run_id, base_run_dir="output", set_global_variables=FALSE, set_seed=TRUE) {
  # TODO: add argument checking to ensure all required arguments are present. 
  
  # Check that the run exists.
  run_id_dir <- file.path(base_run_dir, run_id)
  if(!file.exists(run_id_dir)) stop("Run ID ", run_id_dir, " not found.")
  
  # Check that top-level settings exists. If so load and convert to R list. 
  settings_path <- file.path(run_id_dir, "settings.json")
  if(!file.exists(settings_path)) stop("Settings for run ", run_id, " not found at path ", settings_path)
  settings <- fromJSON(settings_path)
  
  # Check required settings are present. 
  validate_mcmc_run_settings(settings)
  
  # Optionally set global seed. 
  if(set_seed) {
    set.seed(settings$global_seed)
    cat("Global seed set to", settings$global_seed)
  }
  
  # Optionally store settings as global variables. 
  if(set_global_variables) {
    convert_list_to_global_vars(settings)
    print("`settings` elements stored as global variables.")
  }
  
}


convert_list_to_global_vars <- function(settings) {
  
  # Argument checks.
  if(class(settings) != "list") stop("`settings` must be a list.")
  var_names <- names(settings)
  if(is.null(var_names) || any(var_names == "")) {
    stop("`settings` missing names for some or all elements.")
  }
  
  # Ensure no variables already exist in global environment. 
  for(i in seq_along(settings)) {
    var_name <- var_names[i]
    if(exists(var_name, where=.GlobalEnv)) 
      stop("Variable ", var_name, " already exists in `.GlobalEnv`. No global variables were created.")
  }
  
  # Create global variables. 
  for(i in seq_along(settings)) assign(var_names[i], settings[[i]], envir=.GlobalEnv)
  
}


validate_mcmc_run_settings <- function(settings) {
  
  # Must be named list. 
  if(class(settings) != "list") stop("`settings` must be a list.")
  var_names <- names(settings)
  if(is.null(var_names) || any(var_names == "")) {
    stop("`settings` missing names for some or all elements.")
  }
  
  # Must contain set of required settings. 
  required_settings <- c("run_id", "global_seed", "data_seed", "design_seed", "N_obs", 
                         "sig2_eps", "N_design", "design_method", "N_mcmc", "git_hash")
  missing_settings <- setdiff(required_settings, var_names)
  extra_settings <- setdiff(var_names, required_settings)
  
  if(length(extra_settings) > 0) message("Non-required settings found: ", paste(extra_settings, collapse=", "))
  if(length(missing_settings) > 0) stop("Missing required settings: ", paste(missing_settings, collapse=", "))
  
}


# -----------------------------------------------------------------------------
# Other helper functions. 
# -----------------------------------------------------------------------------

get_IG_priors_numerical_test <- function(sig2_true, bias_frac=0, coef_var=0.1, return_prior_plots=FALSE, 
                                         output_variables=NULL, bins=50) {
  # A convenience function to generate inverse gamma priors on the variance parameters in the product Gaussian 
  # likelihood. The priors are constructed by setting the mean of the inverse gamma priors to 
  # sig2_true * (1 + bias_frac), so that a `bias_frac` of 0 implies the prior will be centered on the true value. 
  # The spread of the prior is determined by passing the coefficient of variation `coef_var`. Together, the
  # mean and coefficient of variation determine the shape and scale parameters of the inverse gamma distributions. 
  # Optionally, this function also returns histograms of samples from the prior distributions. 
  #
  # Args:
  #    sig2_true: numeric(p), vector of length p containing the true variance parameters for each of the p outputs. 
  #    bias_frac: numeric, described above. Can be a vector of length p, or length 1, in which case the same 
  #               value will be used for each of the p outputs. 
  #    coef_var: numeric, the coefficient of variation. Can be a vector of length p, or length 1, in which case the same 
  #              value will be used for each of the p outputs.
  #    return_prior_plots: logical(1), if TRUE returns list of prior histogram plots in addition to the prior specifications. 
  #    output_variables: character(p), vector of output variable names. Used to label the plots. 
  #    bins: integer(1), number of bins in the histogram plots. 
  #
  # Returns:
  #    If `return_prior_plots` is FALSE, returns list containing the prior information. Otherwise, returns a list 
  #    where the prior info list is the first element, and the second element is a list of prior histogram plots. 
  
  # Parameters for inverse gamma priors. 
  prior_means <- (1 + bias_frac) * sig2_true
  IG_shape_params <- 2 + (1 / coef_var^2)
  IG_scale_params <- sig2_true * (IG_shape_params - 1)
  
  if(any(IG_shape_params < 2)) stop("Inverse Gamma shape parameter(s) less than 2, meaning variance is infinite.")
  
  sig_eps_prior_params <- list(dist = "IG", 
                               IG_shape = IG_shape_params, 
                               IG_scale = IG_scale_params)
  
  if(!return_prior_plots) return(sig_eps_prior_params)
  
  # Sample from prior and produce histograms. 
  N_samp <- 10000
  p <- length(sig2_true)
  prior_sig_eps_samples <- matrix(NA, nrow = N_samp, ncol = p)
  if(is.null(output_variables)) output_variables <- paste0("output", 1:p)
  colnames(prior_sig_eps_samples) <- output_variables
  
  for(i in seq_len(N_samp)) {
    prior_sig_eps_samples[i,] <- sample_prior_Sig_eps(sig_eps_prior_params)
  }
  
  plts <- vector(mode = "list", length = p)
  for(j in seq_len(p)) {
    x_col <- output_variables[j]
    plts[[j]] <- ggplot(as.data.frame(prior_sig_eps_samples), aes(x = .data[[x_col]])) + 
      geom_histogram(bins = bins) + 
      ggtitle(paste0("Prior Samples: sig2, ", x_col)) + 
      ylab(x_col) + 
      geom_vline(xintercept = sig2_true[j], color = "red")
  }
  
  return(list(prior=sig_eps_prior_params, plots=plts))
  
}




