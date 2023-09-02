

mcmc_info <- list(N_mcmc_exact = 50000, N_mcmc_approx = 50000, learn_sig_eps = TRUE, burn_in_exact = 10000, burn_in_approx = 10000)
init_design_settings <- list(N_design = 10, design_method = "LHS", N_design_reps = 1, batch_sizes = c(5, 10, 15))
emulator_settings <- data.frame(gp_lib = c("hetGP"), 
                                kernel = "Gaussian", 
                                transformation_method = c("truncated"),
                                emulator_target = "SSR",
                                scale_X = TRUE, 
                                normalize_y = TRUE)

# Also need to pass in information to determine data used to compute validation metrics. e.g. allow option to 
# sub-sample exact posterior samples. 
validation_settings <- list(sample_metric_settings  = list(metrics = c("mean", "cov"), param_types = c("theta")), 
                            emulator_metric_settings = list(metrics = c("crps", "rmse"), target = "lpost"))

acquisition_settings <- data.frame(acq_method = "IVAR_lpost", 
                                   candidate_method = "", 
                                   integrate_method = "approx_post", 
                                   batch_method = c("exact", "kriging_believer", "constant_liar_optimist", "constant_liar_pessimist"))

# Each acquisition setting should be combination of: acquisition function, candidate point method, integrate point method, 
# batch acquisition strategy. As baseline should have batch acquisition strategy "exact", which runs model after each 
# iteration. Should allow multiple batch sizes so that we can investigate the effect of batch size when considering errors. 
# For all of these tests, integrate point method will use support point sub-sample of approx posterior points. For baseline, 
# also have integrate points based on prior. For low-dimensional examples, can use prior grid as candidate points to produce 
# nice plots, but for more expensive runs can use a prior-approx posterior mixture (again using support points). Should also 
# probably just set `input_bounds` in advance and keep it the same for all designs. Also need consistent way of determining 
# point estimates to use for lpost_emulator: I think the current init estimates approach for round0 is reasonable. Then for 
# other rounds use mean of approximate MCMC output. Should probably weight emulator-based metrics by 1.) prior grid, 
# 2.) approx posterior samples, 3.) exact posterior samples.  
#
#
# If running on lots of random initial designs, it might be better to batch out over designs on the cluster, while using multicore
# on the acquisition settings. 


# For now only allow one emulator and iterate over acquisition settings, but should generalize to multiple emulators eventually
# (e.g. to be able to compare loss emulation to basis function approach). 
seq_design_batch_sim_test <- function(computer_model_data, prior_info, mcmc_info, emulator_settings, init_design_settings, 
                                      acquisition_settings, validation_settings, random_seed) {
  
  # ---------------------------------------------------------------------------------------------------------
  # Should save a file containing all input information including random seed so that output is reproducible. 
  # ---------------------------------------------------------------------------------------------------------
  # TODO: Will batching out jobs mess up random seed? If so could have separate random seed generated for each job, and 
  # saved to file. 
  
  set.seed(random_seed)
  
  
  # ---------------------------------------------------------------------------------------------------------
  # Run Exact MCMC and format samples. 
  # ---------------------------------------------------------------------------------------------------------
  
  # Run MCMC. 
  time_start <- proc.time()
  mcmc_exact_list <- mcmc_calibrate_product_lik(computer_model_data = computer_model_data, 
                                                theta_prior_params = theta_prior_params, 
                                                learn_sig_eps = TRUE,
                                                sig_eps_prior_params = sig_eps_prior_params,
                                                N_mcmc = N_mcmc_exact)
  mcmc_exact_runtime <- (proc.time() - time_start)[["elapsed"]]
  print(paste0("Exact MCMC runtime: ", mcmc_exact_runtime, " seconds."))
  
  # Format MCMC samples. 
  mcmc_samp_dt <- format_mcmc_output(samp_list = mcmc_exact_list[c("theta", "sig_eps")], test_label = "exact")
  burn_ins <- c(exact = mcmc_info$burn_in_exact)
  
  # Store formatted samples of each parameter type with burn-in dropped. 
  samp_exact_theta <- select_mcmc_samp(mcmc_samp_dt, burn_in_start = burn_ins, test_labels = "exact", param_types = "theta")[, .(param_name, sample)]
  samp_exact_theta <- as.matrix(unstack(samp_exact_theta, sample ~ param_name))[, computer_model_data$pars_cal_names]
  samp_exact_sig_eps <- select_mcmc_samp(mcmc_samp_dt, burn_in_start = burn_ins, test_labels = "exact", param_types = "sig_eps")[, .(param_name, sample)]
  samp_exact_sig_eps<- as.matrix(unstack(samp_exact_sig_eps, sample ~ param_name))[, computer_model_data$pars_cal_names]
  
  # ---------------------------------------------------------------------------------------------------------
  # Experiments with random initial designs. 
  # ---------------------------------------------------------------------------------------------------------
  
  # TODO: replace for loop with batching out jobs on the cluster. 
  for(i in seq_len(init_design_settings$N_design_reps)) {
    run_batch_design_test(computer_model_data, prior_info, mcmc_info, emulator_settings, init_design_settings, 
                          acquisition_settings, validation_settings)
  }

}


run_batch_design_test <- function(computer_model_data, prior_info, mcmc_info, emulator_settings, init_design_settings, 
                                  acquisition_settings, validation_settings, validation_data, input_bounds, 
                                  save_method = "file", dir_save_name = NULL) {
                                  
  # ---------------------------------------------------------------------------------------------------------
  # Generate initial design and compute initial parameter estimates based on design. 
  # ---------------------------------------------------------------------------------------------------------
  
  # TODO: pass input_bounds in here. 
  init_design_info <- get_input_output_design(N_points = init_design_settings$N_design,
                                              design_method = init_design_settings$design_method, 
                                              scale_inputs = TRUE,
                                              computer_model_data = computer_model_data, 
                                              theta_prior_params = prior_info$theta_prior_params)
  init_design_info$init_estimates <- get_init_param_estimates(init_design_info, computer_model_data, 
                                                              prior_info$sig_eps_prior_params)

  # ---------------------------------------------------------------------------------------------------------
  # Fit Emulators on Initial Design 
  # ---------------------------------------------------------------------------------------------------------
  
  # Fit emulators on initial design. 
  gp_fits <- fit_independent_GPs(X_train = init_design_info$inputs_scaled, Y_train = init_design_info$outputs_normalized, 
                                 gp_lib = emulator_settings$gp_lib, gp_kernel = emulator_settings$kernel)$fits
  emulator_info_list <- list(gp_fits = gp_fits, input_bounds = init_design_info$input_bounds, 
                             output_stats = init_design_info$output_stats, settings = emulator_settings)
  
  # Induced log joint density (i.e. log unnormalized posterior density) emulator. 
  lpost_emulator <- get_lpost_emulator_obj(emulator_info_list = emulator_info_list, design_info_list = init_design_info, 
                                           computer_model_data = computer_model_data, sig2_eps = init_design_info$init_estimates$best_sig2_eps, 
                                           theta_prior_params = theta_prior_params)
  

}






