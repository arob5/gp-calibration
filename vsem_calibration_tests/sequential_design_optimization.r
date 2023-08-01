#
# sequential_design_optimization.r
# Functions related to sequential design for GPs and Bayesian Optimization. 
#
# Andrew Roberts
#

run_sequential_design_optimization <- function(acquisition_settings, init_design_settings, emulator_settings, computer_model_data, 
                                               sig2_eps, theta_prior_params, theta_grid_ref = NULL, optimize_sig_eps = FALSE, 
                                               sig_eps_prior_params = NULL) {
  # TODO: generalize so this works with log-normal process. 
  # If `optimize_sig_eps` is TRUE, then treat `sig2_eps` as the initial condition. 
  # TODO: make sure input/output scaling is done correctly. 
  # TODO: Need to clarify what `lpost_max` is tracking; I believe it should be conditional on the most recent 
  #       value of sig2_eps. 
  # TODO: track objective values, both in fixed sig_eps and optimized sig_eps setting. 
  
  # TODO: optimize sig2_eps

  # Create initial design and fit GP. 
  design_info <- get_input_output_design(N_points = init_design_settings$N_design, 
                                         design_method = init_design_settings$design_method, 
                                         computer_model_data = computer_model_data, 
                                         theta_prior_params = theta_prior_params, 
                                         transformation_method = emulator_settings$transformation_method)

  gp_fits <- fit_independent_GPs(X_train = design_info$inputs_scaled, 
                                 Y_train = design_info$outputs_normalized, 
                                 gp_lib = emulator_settings$gp_lib, 
                                 gp_kernel = emulator_settings$kernel)$fits
  
  emulator_info_list <- list(gp_fits = gp_fits, 
                             input_bounds = design_info$input_bounds, 
                             output_stats = design_info$output_stats, 
                             settings = emulator_settings)
  
  # Current observed objective values (i.e. log posterior values). 
  design_objective_vals <- calc_lpost_theta_product_lik(theta_vals = design_info$inputs, 
                                                        computer_model_data = computer_model_data,  
                                                        SSR = design_info$outputs, 
                                                        vars_obs = sig2_eps, 
                                                        na.rm = TRUE, 
                                                        return_list = FALSE, 
                                                        theta_prior_params = theta_prior_params)
  design_best_idx_curr <- which.max(design_objective_vals)
  design_best_idx <- c(design_best_idx_curr)
  lpost_max <- design_objective_vals[design_best_idx] # TODO: probably want to track lpost_max over time as well. 
  new_design_idx_curr <- 1
  
  # Sequential Design/Bayesian Optimization loop. 
  for(i in seq_len(acquisition_settings$N_opt_iter)) {
    
    # Optimize sig_eps, if specified. 
    if(optimize_sig_eps) {
      sig2_eps <- optimize_sig_eps_cond_post(SSR_theta = design_info$outputs[design_best_idx_curr], 
                                             sig_eps_prior_params = sig_eps_prior_params,
                                             n_obs = computer_model_data$n_obs)
    }
    
    # Obtain new batch of design points (without running the forward model).  
    inputs_new <- batch_acquisition_opt_one_step(emulator_info_list = emulator_info_list, 
                                                 acquisition_settings = acquisition_settings, 
                                                 design_input_curr = design_info$inputs, 
                                                 design_objective_curr = design_objective_vals, 
                                                 design_best_idx = design_best_idx_curr,  
                                                 computer_model_data = computer_model_data, 
                                                 sig2_eps = sig2_eps, 
                                                 theta_prior_params = theta_prior_params, 
                                                 theta_grid_ref = theta_grid_ref)

    # Run forward model at input points in batch.  
    # TODO: this should be parallelized.
    SSR_new <- get_computer_model_SSR(computer_model_data = computer_model_data, theta_vals = inputs_new, na.rm = TRUE)
    lpost_new <- calc_lpost_theta_product_lik(computer_model_data = computer_model_data, 
                                              theta_vals = inputs_new, 
                                              SSR = SSR_new,
                                              vars_obs = sig2_eps, 
                                              na.rm = TRUE, 
                                              theta_prior_params = theta_prior_params, 
                                              return_list = FALSE)
    
    # Update design and keep track of current MAP estimate. If conditional on a new sig_eps, then `design_best_idx_curr`
    # is always updated. If sig_eps is fixed, then may not be updated. 
    design_best_idx_new <- which.max(lpost_new)
    if(optimize_sig_eps || (lpost_new[design_best_idx_new] > lpost_max)) {
      design_best_idx_curr <- length(design_info$outputs) + design_best_idx_new
      lpost_max <- lpost_new[design_best_idx_new]
    }
    design_best_idx <- c(design_best_idx, design_best_idx_curr)
    design_info$inputs <- rbind(design_info$inputs, inputs_new)
    design_info$outputs <- rbind(design_info$outputs, SSR_new)

    # Index of the first new design point in batch added to the design.  
    new_design_idx_curr <- new_design_idx_curr + acquisition_settings$batch_size
    
    # Update GP (including hyperparameter estimates). 
    # TODO: modify `update_independent_GPs()` so that it can optionally modify hyperparameter estimates or not. 
    emulator_info_list$gp_fits <- update_independent_GPs(gp_fits = emulator_info_list$gp_fits, 
                                                         gp_lib = emulator_settings$gp_lib, 
                                                         X_new = inputs_new, 
                                                         Y_new = SSR_new, 
                                                         input_bounds = emulator_info_list$input_bounds, 
                                                         output_stats = emulator_info_list$output_stats) 
  }
  
  return(list(emulator_info_list = emulator_info_list, design_best_idx))
  
}


optimize_sig_eps_cond_post <- function(SSR_theta, sig_eps_prior_params, n_obs) {
  # Computes the closed form optimization of the conditional posterior p(Sig_eps | theta, Y) for the product likelihood 
  # with inverse gamma priors on the variance parameters. 
  #
  # Args:
  #    SSR_theta: numeric(), vector of length equal to the number of outputs P, containing the squared L2 errors for each output computed using 
  #               the calibration input `theta` being conditioned upon in the conditional posterior. The value `theta` is not required to be 
  #               explicitly passed here, only SSR(theta). 
  #    sig_eps_prior_params: list, must contain elements named "IG_shape" and "IG_scale" which 
  #                          each correspond to P-length vectors storing the parameters for the independent Inverse Gamma priors on each
  #                          variance parameter.
  #    n_obs: numeric(), vector of length equal to the number of outputs P. The number of observations for each output. 
  
  (0.5 * SSR_theta + sig_eps_prior_params$IG_scale) / (0.5 * n_obs + sig_eps_prior_params$IG_shape + 1)
  
}
                           

batch_acquisition_opt_one_step <- function(lpost_emulator, acquisition_settings) {
  # Acquires a batch of new input points via a greedy, heuristic approach. Currently only supports the kriging 
  # believer heuristic, but should be generalized to other heuristics as well (e.g. constant liar). Returns the 
  # batch of (scaled) input points, but does not run the forward model at these new inputs or update the GPs. 
  #
  # Args:
  #    emulator_info_list: list, the emulator information list. 
  #    acquisition_settings: list, the acquisition settings list. 
  #    computer_model_data: list, the computer model data list. 
  #    theta_prior_params: data.frame, defining prior distributions on the calibration parameters. Required by some  
  #                        acquisition functions. 
  #
  # Returns:
  #    matrix of dimension (batch size, parameter dimension) containing the batch of scaled input points. 
  
  # Set constant liar value to use as responses in pseudo-updates, if relevant. 
  if(grepl("constant_liar", acquisition_settings$batch_method)) {
    lpost_emulator$constant_liar_value <- get_constant_liar_value(lpost_emulator, acquisition_settings$batch_method)
  }
  
  # Acquire batch of input points (without running forward model). 
  for(b in 1:acquisition_settings$batch_size) {
    
    # Optimize sequential acquisition function.  
    input_new_scaled <- acquisition_opt_one_step(lpost_emulator = lpost_emulator, acquisition_settings = acquisition_settings)

    # Update lpost emulator using batch heuristic method. Does not affect underlying GPs, including their hyperparameters. 
    lpost_emulator <- pseudo_update_lpost_emulator(lpost_emulator, inputs_new_scaled = input_new_scaled,
                                                   pseudo_update_method = acquisition_settings$batch_method)
  }
  
  return(lpost_emulator)
  
}


acquisition_opt_one_step <- function(lpost_emulator, acquisition_settings, max_objective_curr = NULL) {
  # Performs a single one-point acquisition by optimizing the specified acquisition function. Returns 
  # the newly acquired (scaled) design point. 
  #
  # Args:
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`.  
  #    acquisition_settings: list, the acquisition settings list. 
  #    max_objective_curr: the current observed maximum objective value. Passed to some acquisition functions 
  #                        (e.g. expected improvement) as a threshold value. Could generalize this later to 
  #                        allow the threshold setting to be specified in `acquisition_settings` and then 
  #                        compute the threshold rather than assuming it is the maximum objective. 
  #
  # Returns:
  #    matrix of dimension 1xd, the (scaled) input returned by optimizing the acquisition. 
  
  # Select next point by optimizing acquisition function.
  if(acquisition_settings$opt_method == "grid") {
    theta_new <- optimize_acquisition_grid(acquisition_type = acquisition_settings$acquisition_type, 
                                           theta_grid_candidate = acquisition_settings$theta_grid_candidate, 
                                           lpost_emulator = lpost_emulator,
                                           N_MC_samples = N_MC_samples, 
                                           theta_grid_integrate = acquisition_settings$theta_grid_integrate, 
                                           N_subsample_candidate = acquisition_settings$N_subsample_candidate, 
                                           N_subsample_integrate = acquisition_settings$N_subsample_integrate, 
                                           threshold_lpost = max_objective_curr)
  } else {
    stop("Invalid acquisition optimization method: ", opt_method)
  }
  
  return(theta_new)
  
}


acquisition_opt_one_step_old <- function(emulator_info_list, acquisition_settings, computer_model_data, 
                                         max_objective_curr = NULL, sig2_eps = NULL, theta_prior_params = NULL) {
  # Performs a single one-point acquisition by optimizing the specified acquisition function. Returns 
  # the newly acquired (scaled) design point. 
  #
  # Args:
  #    emulator_info_list: list, the emulator information list. 
  #    acquisition_settings: list, the acquisition settings list. 
  #    computer_model_data: list, the computer model data list. 
  #    max_objective_curr: the current observed maximum objective value. Passed to some acquisition functions 
  #                        (e.g. expected improvement) as a threshold value. Could generalize this later to 
  #                        allow the threshold setting to be specified in `acquisition_settings` and then 
  #                        compute the threshold rather than assuming it is the maximum objective. 
  #    sig2_eps: numeric, vector of length P = number of output variables, containing the likelihood observation variances.
  #              Required by most of the acquisition functions. 
  #    theta_prior_params: data.frame, defining prior distributions on the calibration parameters. Required by some  
  #                        acquisition functions. 
  #
  # Returns:
  #    matrix of dimension 1xd, the (scaled) input returned by optimizing the acquisition. 
  
  # Select next point by optimizing acquisition function.
  if(acquisition_settings$opt_method == "grid") {
    theta_new <- optimize_acquisition_grid(acquisition_type = acquisition_settings$acquisition_type, 
                                           theta_grid_candidate = acquisition_settings$theta_grid_candidate, 
                                           emulator_info_list = emulator_info_list,
                                           computer_model_data = computer_model_data,
                                           theta_prior_params = theta_prior_params, 
                                           sig2_eps = sig2_eps, 
                                           N_MC_samples = N_MC_samples, 
                                           theta_grid_integrate = acquisition_settings$theta_grid_integrate, 
                                           N_subsample_candidate = acquisition_settings$N_subsample_candidate, 
                                           N_subsample_integrate = acquisition_settings$N_subsample_integrate, 
                                           threshold_lpost = max_objective_curr)
  } else {
    stop("Invalid acquisition optimization method: ", opt_method)
  }
                              
  return(theta_new)

}


optimize_acquisition_grid <- function(acquisition_type, theta_grid_candidate, lpost_emulator, N_MC_samples = NULL, theta_grid_integrate = NULL,
                                      N_subsample_candidate = NULL, N_subsample_integrate = NULL, threshold_lpost = NULL) { 
  # Returns a new (scaled) design point given by optimizing the acquisition function over a finite set of candidate (i.e. grid) points. 
  # All input points (candidate or integration points) are assumed to already be properly scaled. 
  #
  # Args:
  #    acquisition_type: character, used to select the acquisition function. Acquisition function naming convention is 
  #                      "acquisition_<acquisition_type>". 
  #    theta_grid_candidate: matrix of shape (# candidate points, d=dimension of input space). The points should already be scaled. 
  #                          The acquisition function will be evaluated at each candidate point and the arg max over the points returned. 
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`.
  #    All remaining arguments are fed to the acquisition function; the specific required arguments depend on the particular acquisition function. 
  #
  # Returns:
  #    matrix of dimension 1xd, the (scaled) input in `theta_grid_candidate` with the largest acquisition value. 
  
  # Get acquisition function. 
  acquisition_func <- get(paste0("acquisition_", acquisition_type))
  
  # If specified, obtain sub-sample of candidate and/or integration points. 
  if(!is.null(N_subsample_candidate)) theta_grid_candidate <- theta_grid_candidate[sample(1:nrow(theta_grid_candidate), size = N_subsample_candidate, replace = FALSE),, drop=FALSE]
  if(!is.null(N_subsample_integrate)) theta_grid_integrate <- theta_grid_integrate[sample(1:nrow(theta_grid_integrate), size = N_subsample_integrate, replace = FALSE),, drop=FALSE]
  
  # Evaluate acquisition function on grid of reference inputs. 
  acquisition_vals_grid <- acquisition_func(theta_vals = theta_grid_candidate, 
                                            lpost_emulator = lpost_emulator,
                                            threshold_lpost = threshold_lpost, 
                                            N_MC_samples = N_MC_samples, 
                                            theta_grid_integrate = theta_grid_integrate)
  
  # Return scaled input in reference grid that maximizes the acquisition. 
  return(theta_grid_candidate[which.max(acquisition_vals_grid),,drop = FALSE])
  
}


optimize_acquisition_grid_old <- function(acquisition_type, theta_grid_candidate, emulator_info_list, computer_model_data, 
                                          theta_prior_params = NULL, sig2_eps = NULL, N_MC_samples = NULL, theta_grid_integrate = NULL, 
                                          N_subsample_candidate = NULL, N_subsample_integrate = NULL, threshold_lpost = NULL) { 
  # Returns a new (scaled) design point given by optimizing the acquisition function over a finite set of candidate (i.e. grid) points. 
  # All input points (candidate or integration points) are assumed to already be properly scaled. 
  #
  # Args:
  #    acquisition_type: character, used to select the acquisition function. Acquisition function naming convention is 
  #                      "acquisition_<acquisition_type>". 
  #    theta_grid_candidate: matrix of shape (# candidate points, d=dimension of input space). The points should already be scaled. 
  #                          The acquisition function will be evaluated at each candidate point and the arg max over the points returned. 
  #    emulator_info_list: character(1), the emulator information list. 
  #    computer_model_data: character(1), the computer model data list. 
  #    All remaining arguments are fed to the acquisition function; the specific required arguments depend on the particular acquisition function. 
  #
  # Returns:
  #    matrix of dimension 1xd, the (scaled) input in `theta_grid_candidate` with the largest acquisition value. 
  
  # Get acquisition function. 
  acquisition_func <- get(paste0("acquisition_", acquisition_type))
  
  # If specified, obtain sub-sample of candidate and/or integration points. 
  if(!is.null(N_subsample_candidate)) theta_grid_candidate <- theta_grid_candidate[sample(1:nrow(theta_grid_candidate), size = N_subsample_candidate, replace = FALSE),, drop=FALSE]
  if(!is.null(N_subsample_integrate)) theta_grid_integrate <- theta_grid_integrate[sample(1:nrow(theta_grid_integrate), size = N_subsample_integrate, replace = FALSE),, drop=FALSE]
  
  # Evaluate acquisition function on grid of reference inputs. 
  acquisition_vals_grid <- acquisition_func(theta_vals = theta_grid_candidate, 
                                            emulator_info_list = emulator_info_list,
                                            computer_model_data = computer_model_data, 
                                            theta_prior_params = theta_prior_params, 
                                            sig2_eps = sig2_eps, 
                                            threshold_lpost = threshold_lpost, 
                                            N_MC_samples = N_MC_samples, 
                                            theta_grid_integrate = theta_grid_integrate)
  
  # Return scaled input in reference grid that maximizes the acquisition. 
  return(theta_grid_candidate[which.max(acquisition_vals_grid),,drop = FALSE])

}


pseudo_update_lpost_emulator <- function(lpost_emulator, input_new_scaled, pseudo_update_method, input_new_unscaled = NULL) {
  
  # Determine response value to use in pseudo-update. 
  if(pseudo_update_method == "kriging_believer") {
    output_lpost_new <- NULL
  } else if(grepl("constant_liar", pseudo_update_method)) {
    if(is.null(lpost_emulator$constant_liar_value)) stop("`constant_liar_value` missing for pseudo-update.")
    output_lpost_new <- lpost_emulator$constant_liar_value
  } else {
    stop("Invalid value for `pseudo_update_method`: ", pseudo_update_method)
  } 
  
  lpost_emulator <- update_lpost_emulator(lpost_emulator, inputs_new_scaled = input_new_scaled, outputs_lpost_new = output_lpost_new)
  
  return(lpost_emulator)
  
}


get_constant_liar_value <- function(lpost_emulator, constant_liar_method) {
  
  if(pseudo_update_method == "constant_liar_pessimist") {
    return(min(lpost_emulator$outputs_lpost))
  } else if(pseudo_update_method == "constant_liar_optimist") {
    return(max(lpost_emulator$outputs_lpost))
  } else if(pseudo_upodate_method == "constant_liar_mean") {
    return(mean(lpost_emulator$outputs_lpost))
  } else {
    stop("Invalid constant liar method: ", constant_liar_method)
  }
  
}
                             

# ---------------------------------------------------------------------------------
# Acquisition Functions:
#    - For both optimization and design. 
#    - All acquisition functions assume that arguments related to the input/ 
#      calibration space (e.g. `theta_vals`, `theta_grid_integrate`) are scaled.
#    - All acquisition functions are designed to be maximized. 
#    - Acquisition function naming convention is "acquisition_<acquisition_type>"
# ---------------------------------------------------------------------------------

acquisition_EI_lpost <- function(theta_vals, lpost_emulator, threshold_lpost, ...) {
  # Closed-form implementation of the expected improvement (EI) acquisition function applied to an lpost emulator 
  # which is a Gaussian Process. 
  #
  # Args:
  #    theta_vals: matrix of dimension M x D, each row an input location at which to compute the value of the acquisition function. 
  #                The inputs are assumed to already be properly scaled. 
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`. 
  #    threshold_lpost: the threshold objective value to use in the EI calculation, typically the current observed maximum of the 
  #                     unnormalized log posterior. 
  #
  # Returns:
  #    numeric, vector of length equal to length of `theta_vals`; the EI acquisition function evaluations at inputs `theta_vals`. 
  
  # Predictive mean and variance of lpost emulator evaluated at inputs `theta_vals`. 
  lpost_pred_list <- predict_lpost_emulator(inputs_new_scaled = theta_vals, lpost_emulator = lpost_emulator, include_nugget = TRUE)
  mu <- lpost_pred_list$mean
  sig <- sqrt(lpost_pred_list$var)
  
  improvement <- mu - threshold_lpost
  EI <- improvement * pnorm(improvement / sig) + sig * dnorm(improvement / sig)
  
  return(EI)
  
}


acquisition_EI_lpost_old <- function(theta_vals, emulator_info_list, computer_model_data, 
                                     theta_prior_params, sig2_eps, threshold_lpost, gp_pred_list = NULL, ...) {
  
  # Predictive mean and variance of log posterior approximation evaluated at inputs `theta_vals`. 
  lpost_pred_list <- predict_lpost_GP_approx(theta_vals_scaled = theta_vals, emulator_info_list = emulator_info_list,
                                             sig2_eps = sig2_eps, theta_prior_params = theta_prior_params, 
                                             N_obs = computer_model_data$n_obs, include_nugget = TRUE, 
                                             gp_pred_list = gp_pred_list, return_vals = c("mean", "var"))
  
  mu <- lpost_pred_list$mean
  sig <- sqrt(lpost_pred_list$var)
  
  improvement <- mu - threshold_lpost
  EI <- improvement * pnorm(improvement / sig) + sig * dnorm(improvement / sig)
                                            
  return(EI)
  
}


acquisition_EI_lpost_MC <- function(theta_vals, emulator_info_list, computer_model_data, 
                                    theta_prior_params, sig2_eps, threshold_lpost, N_MC_samples = 1000, gp_pred_list = NULL, ...) {
  # A Monte Carlo (MC) approximation to the expected improvement (EI) criterion targeting the approximation 
  # to the unnormalized log posterior density. This is the MC analog of `acquisition_EI_lpost()`; the latter assumes 
  # that the predictive distributions of the GPs targeting the SSR functions are Gaussian, while this function allows 
  # for approximate computation in the case of non-Gaussian predictive distributions. 
  #
  # Args:
  #    theta_vals: matrix of dimension M x D, each row an input location at which to compute the value of the acquisition function. 
  #                The inputs are assumed to already be properly scaled. 
  #    emulator_info_list: list, the emulator information list. 
  #    computer_model_data: list, the computer model information list. 
  #    theta_prior_params: data.frame, defining prior distributions on the calibration parameters. 
  #    sig2_eps: numeric, vector of length P = number of output variables, containing the likelihood observation variances.
  #    threshold_lpost: the threshold objective value to use in the EI calculation, typically the current observed maximum.
  #    N_MC_samples: integer, the number of samples to draw from the log posterior density approximation at each input in order 
  #                  to approximate PI. 
  #    gp_pred_list: A list, as returned by `predict_independent_GPs()`. This allows the predictive means and 
  #                  variances of the underlying GP emulators evaluated at the input locations at the to be passed 
  #                  if they have already been computed. If NULL, they are computed here. 
  #    
  # Returns:
  #    numeric vector of length M, containing the evaluations of the acquisition at the M input points in `theta_vals`.
  
  # Have this return a matrix of dim N_MC_samples x nrow(theta_vals)
  lpost_samp <- sample_GP_lpost_theta(theta_vals_scaled = theta_vals, 
                                      emulator_info_list = emulator_info_list,
                                      computer_model_data = computer_model_data, 
                                      theta_prior_params = theta_prior_params, 
                                      sig2_eps = sig2_eps, 
                                      N_samples = N_MC_samples, 
                                      gp_pred_list = gp_pred_list)
                                    
  # Monte Carlo estimates of acquisition at each point. 
  acq_EI_MC_estimates <- lpost_samp - threshold_lpost
  acq_EI_MC_estimates[acq_EI_MC_estimates < 0] <- 0
  
  return(colMeans(acq_EI_MC_estimates))
  
}


acquisition_PI_lpost_MC <- function(theta_vals, emulator_info_list, computer_model_data, theta_prior_params, sig2_eps, 
                                    threshold_lpost, N_MC_samples = 1000, gp_pred_list = NULL, ...) {
  # A Monte Carlo (MC) approximation to the probability of improvement (PI) criterion targeting the approximation 
  # to the unnormalized log posterior density. This is the MC analog of `acquisition_PI_lpost()`; the latter assumes 
  # that the predictive distributions of the GPs targeting the SSR functions are Gaussian, while this function allows 
  # for approximate computation in the case of non-Gaussian predictive distributions. 
  #
  # Args:
  #    theta_vals: matrix of dimension M x D, each row an input location at which to compute the value of the acquisition function. 
  #                The inputs are assumed to already be properly scaled. 
  #    emulator_info_list: list, the emulator information list. 
  #    computer_model_data: list, the computer model information list. 
  #    theta_prior_params: data.frame, defining prior distributions on the calibration parameters. 
  #    sig2_eps: numeric, vector of length P = number of output variables, containing the likelihood observation variances.
  #    threshold_lpost: the threshold objective value to use in the PI calculation, typically the current observed maximum.
  #    N_MC_samples: integer, the number of samples to draw from the log posterior density approximation at each input in order 
  #                  to approximate PI. 
  #    gp_pred_list: A list, as returned by `predict_independent_GPs()`. This allows the predictive means and 
  #                  variances of the underlying GP emulators evaluated at the input locations at the to be passed 
  #                  if they have already been computed. If NULL, they are computed here. 
  #    
  # Returns:
  #    numeric vector of length M, containing the evaluations of the acquisition at the M input points in `theta_vals`.
  
  # Have this return a matrix of dim N_MC_samples x nrow(theta_vals)
  lpost_samp <- sample_GP_lpost_theta(theta_vals_scaled = theta_vals, 
                                      emulator_info_list = emulator_info_list,
                                      computer_model_data = computer_model_data, 
                                      theta_prior_params = theta_prior_params, 
                                      sig2_eps = sig2_eps, 
                                      N_samples = N_MC_samples, 
                                      gp_pred_list = gp_pred_list)
  
  # Monte Carlo estimates of acquisition at each point. 
  acq_EI_MC_estimates <- lpost_samp - threshold_lpost
  acq_EI_MC_estimates <- colMeans(acq_EI_MC_estimates > 0)

  return(acq_EI_MC_estimates)
  
}


acquisition_EIVAR_lpost <- function(theta_vals, emulator_info_list, sig2_eps, theta_grid_integrate, ...) {
  # Implements the expected integrated variance (EIVAR) criteria that targets the log posterior in the loss emulation setting. 
  # In this case the inner two integrals of EIVAR are available in closed form. The outer integral, the expectation over the input 
  # space is approximated by a finite sum over grid points `theta_grid_integrate`.
  #
  # Args:
  #    theta_vals: matrix of dimension M x D, each row an input location at which to compute the value of the acquisition function. 
  #                The inputs are assumed to already be properly scaled. 
  #    emulator_info_list: list, the emulator information list. 
  #    theta_prior_params: data.frame, defining prior distributions on the calibration parameters. 
  #    sig2_eps: numeric, vector of length P = number of output variables, containing the likelihood observation variances.
  #    theta_grid_integrate: matrix, of dimension M_integrate x D. The set of inputs used to approximate the integral over the 
  #                          input space required in evaluating the EIVAR criterion. 
  # 
  # Returns:
  #    numeric vector of length M, containing the evaluations of the acquisition at the M input points in `theta_vals`. 
  #    Technically returns the negative of EIVAR to align with the acquisition convention that bigger is better. 
  
  # Handle case of single input. 
  if(is.null(nrow(theta_vals))) theta_vals <- matrix(theta_vals, nrow = 1)
  
  # Vector to store EIVAR estimates at inputs `theta_vals`. 
  EIVAR_est <- vector(mode = "numeric", length = nrow(theta_vals))
  
  # Current GP predictive distribution. 
  gp_fits_curr <- emulator_info_list$gp_fits
  
  for(i in 1:nrow(theta_vals)) {
    
    # Update variance by conditioning on theta evaluation value. Should not affect `emulator_info_list` outside of local function scope. 
    emulator_info_list$gp_fits <- update_independent_GPs(gp_fits = gp_fits_curr, gp_lib = emulator_info_list$settings$gp_lib, 
                                                         X_new = theta_vals[i,,drop=FALSE], Y_new = NULL, update_hyperparameters = FALSE)
                                                  
    # Compute unnormalized log posterior approximation predictive variance at each theta grid location. 
    lpost_pred_var_grid <- predict_lpost_GP_approx(theta_vals_scaled = theta_grid_integrate, emulator_info_list = emulator_info_list, 
                                                   sig2_eps = sig2_eps, include_nugget = TRUE, return_vals = "var")$var
    
    # Estimate EIVAR via discrete sum over theta grid locations. 
    EIVAR_est[i] <- -1.0 * mean(lpost_pred_var_grid)
    
  }
  
  return(EIVAR_est)
  
}


acquisition_VAR_lpost <- function(theta_vals, emulator_info_list, sig2_eps, ...) {
  # Implements the acquisition which is simply defined as the variance of the unnormalized log posterior approximation in  
  # the loss emulation setting. 
  #
  # Args:
  #    theta_vals: matrix of dimension M x D, each row an input location at which to compute the value of the acquisition function. 
  #                The inputs are assumed to already be properly scaled. 
  #    emulator_info_list: list, the emulator information list. 
  #    sig2_eps: numeric, vector of length P = number of output variables, containing the likelihood observation variances. 
  # 
  # Returns:
  #    numeric vector of length M, containing the evaluations of the acquisition at the M input points in `theta_vals`. 
  
  # Handle case of single input. 
  if(is.null(nrow(theta_vals))) theta_vals <- matrix(theta_vals, nrow = 1)
  
  # Compute unnormalized log posterior approximation predictive variance at each theta grid location. 
  lpost_pred_var_grid <- predict_lpost_GP_approx(theta_vals_scaled = theta_vals, emulator_info_list = emulator_info_list, 
                                                 sig2_eps = sig2_eps, return_vals = "var", include_nugget = TRUE, include_sig_eps_prior = FALSE)
  
  return(lpost_pred_var_grid$var)
  
}


acquisition_VAR_post <- function(theta_vals, emulator_info_list, computer_model_data, theta_prior_params, sig2_eps, ...) {
  # Implements the acquisition which is simply defined as the variance of the unnormalized posterior density approximation 
  # in the loss emulation setting. Due to the typical large dynamic range in predictive mean/variance values, the log of the 
  # acquisition is returned (i.e. the log predictive variance of the unnormalized posterior density approximation). Let 
  # mu and sig2 denote the mean and variance of the predictive distribution of the log posterior approximation (lpost). Then 
  # the log variance of the posterior approximation is log[exp(sig2) - 1] + 2*mu + sig2. Direct computation of the first term 
  # is problematic when the range of sig2 values is very large. However, when sig2 is large (say, > 100) then 
  # log[exp(sig2) - 1] ~ sig2 is a very good approximation. This approximation is applied at this cutoff for numerical stability.
  # Thus, for large values of sig2, the returned value is 2(sig2 + mu). Note that this is similar to the upper confidence bound (UCB)
  # acquisition for the random field approximation of the LOG posterior; however, the UCB uses the predictive standard deviation instead
  # of the variance. The use of the variance in 2(sig2 + mu) often means that the predictive variance dominates this function in situations
  # when the GP is fairly uncertain. 
  #
  # Args:
  #    theta_vals: matrix of dimension M x D, each row an input location at which to compute the value of the acquisition function. 
  #                The inputs are assumed to already be properly scaled. 
  #    emulator_info_list: list, the emulator information list. 
  #    computer_model_data: list, the computer model information list.
  #    theta_prior_params: data.frame, defining prior distributions on the calibration parameters. 
  #    sig2_eps: numeric, vector of length P = number of output variables, containing the likelihood observation variances. 
  # 
  # Returns:
  #    numeric vector of length M, containing the evaluations of the acquisition at the M input points in `theta_vals`. 
  
  # Handle case of single input. 
  if(is.null(nrow(theta_vals))) theta_vals <- matrix(theta_vals, nrow = 1)
  
  # Compute unnormalized log posterior approximation predictive variance at each theta grid location. 
  lpost_pred_list<- predict_lpost_GP_approx(theta_vals_scaled = theta_vals, 
                                            emulator_info_list = emulator_info_list, 
                                            sig2_eps = sig2_eps, 
                                            theta_prior_params = theta_prior_params, 
                                            N_obs = computer_model_data$n_obs,
                                            include_nugget = TRUE, 
                                            include_sig_eps_prior = FALSE)
  mu <- lpost_pred_list$mean
  sig2 <- lpost_pred_list$var
  

  # Numerically stable computation of the log variance. 
  idx_approx_sel <- (sig2 >= 100)
  term_1 <- vector(mode = "numeric", length = length(sig2))
  term_1[!idx_approx_sel] <- log(exp(sig2[!idx_approx_sel]) - 1)
  term_1[idx_approx_sel] <- sig2[idx_approx_sel] # Apply approximation. 
  
  log_pred_var_post <- term_1 + 2*mu + sig2

  
  return(log_pred_var_post)
  
}






