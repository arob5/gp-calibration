#
# seq_design_gp.r
#
# Functions related to sequential design and optimization for Gaussian processes
# (GPs). Currently, this concerns design and optimization directly for GPs, in 
# addition to exponentiated GPs (i.e., log-normal processes). The methods 
# implemented here fall under the framework of optimizing some objective function 
# to select the design points. We refer to the objective function as an 
# "acquisition function" ("acq" for short), but this function is also called
# a design criterion or refinement function in the literature. 
#
# Depends: These functions are designed to work with objects from the 
#          gpWrapper class. 
#
# Andrew Roberts
#

library(matrixStats)


acquire_batch_input_sequentially <- function(gp, acq_func_name, N_batch, model_response_heuristic, 
                                             opt_method, f_exact=NULL, reoptimize_hyperpar=FALSE, ...) {
  # This function implements an acquisition optimization-based sequential design loop, sequentially 
  # requiring design points (inputs) serially one at a time. If `model_response_heuristic` is "none" then 
  # this means the full forward model `f_exact` must be run after each new input is acquired. Other options  
  # for `model_response_heuristic` implement heuristics to avoid forward model evaluations. 
  # `model_response_heuristic` can be "KB", "CL_pessimist", "CL_optimist", or "none"; the latter runs the 
  # full forward model. `opt_method` currently only allows "grid". This function can act as a typical 
  # purely sequential design look; it can also act as a single step of a batch sequential design 
  # loop, provided that a model response heuristic is specified. 
  #
  # Returns:
  #     Note that `gp_updated` will contain the pseudo model responses if a model response 
  #     heuristic is used. 
  
  # Should verify that `gp` is of class gpWrapper with `dim_Y==1`. 
  # TODO: validate_args_acquire_batch_input_sequentially()

  assert_that(is_gp(gp))
  assert_that(gp$Y_dim==1L)
  
  # TODO: need to add `reoptimize_hyperpar` for gpWrapperKerGP. 
  if((model_response_heuristic != "none") && (reoptimize_hyperpar)) {
    message("`reoptimize_hyperpar` is TRUE but `model_response_heuristic` is not none.")
  }
  
  # Make a copy to avoid modifying the emulator provided in argument. 
  gp_copy <- gp$copy(shallow=FALSE)
  
  # Objects to store the acquired inputs and the associated model (perhaps pseudo) responses. 
  input_batch <- matrix(nrow=N_batch, ncol=gp_copy$X_dim, dimnames=list(NULL,gp_copy$X_names))
  response_batch <- rep(NA_real_, N_batch)
  
  for(i in 1:N_batch) {
    # Acquire new input point. 
    input_new <- optimize_acq_single_input(acq_func_name, gp_copy, opt_method, ...)
    input_batch[i,] <- input_new
    
    # Acquire model response or pseudo model response at acquired input. 
    response_new <- get_acq_model_response(input_new, model_response_heuristic, gp, f_exact, ...)
    response_batch[i] <- response_new
    
    # Update emulator. 
    gp_copy$update(matrix(input_new, nrow=1), matrix(response_new, nrow=1), 
                   update_hyperpar=reoptimize_hyperpar, ...)
  }
  
  return(list(input_batch=input_batch, response_batch=response_batch, gp_updated=gp_copy))
}


optimize_acq_single_input <- function(acq_func_name, gp, opt_method, candidate_grid=NULL, ...) {
  # Optimizes an acquisition function to return the single input minimizing the 
  # acquisition function. 
  
  # Define objective function for the optimization. 
  acq_func <- get(paste0("acq_", acq_func_name))
  
  # Dispatch to the correct optimization algorithm. 
  if(opt_method == "grid") input_new <- minimize_objective_grid(acq_func, candidate_grid, gp=gp, ...)
  else stop("`opt_method` ", opt_method, " not supported.")
  
  return(input_new)
}


evaluate_gp_acq_func_vectorized <- function(acq_func, input_mat, gp=NULL, ...) {
  apply(input_mat, 1, function(input) acq_func(matrix(input, nrow=1), gp=gp, ...))
}


get_acq_model_response <- function(input, model_response_heuristic, 
                                   gp=NULL, f_exact=NULL, ...) {
  
  if(model_response_heuristic == "none") return(f_exact(input))
  else .NotYetImplemented()
}



# -----------------------------------------------------------------------------
# Acquisition Functions
#
# Acquisition functions must be called as `acq_<acq_func_name>(input, ...)`, 
# where `input` is a numeric D-length vector or 1xD matrix representing a 
# parameter value. 
#
# -----------------------------------------------------------------------------

acq_IVAR_grid <- function(input, gp, grid_points, weights=1/nrow(grid_points), ...) {
  # A grid-based (sample sum approximation) of the integrated variance criterion 
  # for GPs (also known as integrated mean squared prediction error). When `input` 
  # is a matrix with more than 1 row, then the acquisition will be computed
  # in batch mode, meaning that it considers conditioning on the entire 
  # batch of inputs. Note that this is different from the function simply 
  # being vectorized across multiple inputs. For the latter, use 
  # `evaluate_acq_func_vectorized()`. 
  
  # TODO: validate_args_acq_IVAR_grid()
  
  N_grid <- nrow(grid_points)
  if(length(weights)==1L) weights <- rep(weights, N_grid)
  gp_copy <- gp$copy(shallow=FALSE)
  
  # Condition the GP on the new batch of inputs `input`. Since the conditional variance 
  # does not depend on the response, the associated batch response is just set to a 
  # vector of zeros. 
  pseudo_response <- matrix(0, nrow=nrow(input), ncol=1)
  gp_copy$update(input, pseudo_response, update_hyperpar=FALSE, ...)
  
  # Evaluate conditional variance at grid points. 
  pred_cond <- gp_copy$predict(grid_points, return_mean=FALSE, return_var=TRUE, ...)
  
  # Return the weighted sum of conditional variances. 
  return(sum(drop(pred_cond$var) * weights))
}


acq_IEVAR_grid <- function(input, gp, grid_points, weights=1/nrow(grid_points), log_scale=TRUE, ...) {
  # This function targets exploration for an exponentiated GP. It can be thought 
  # of as an integrated mean squared prediction error criterion for 
  # log-normal processes. The outer integral over the design space is 
  # approximated with a discrete sum over grid points `grid_points` which are given 
  # weights `weights`. When `input` is a matrix with more than 1 row, then the
  # acquisition will be computed in batch mode, meaning that it considers conditioning
  # on the entire batch of inputs. Note that this is different from the function simply 
  # being vectorized across multiple inputs. For the latter, use 
  # `evaluate_acq_func_vectorized()`
  
  assert_that(gp$Y_dim==1L)
  
  log_evar <- gp$calc_expected_exp_cond_var(grid_points, input, log_scale=TRUE, ...)[,1]
  log_summands <- log_evar + log(weights)
  log_IEVAR <- matrixStats::logSumExp(log_summands)
  
  if(log_scale) return(log_IEVAR)
  return(exp(log_IEVAR))
}


acq_IENT_grid <- function(input, gp, grid_points, weights=1/nrow(grid_points), ...) {
  # A grid-based (sample sum approximation) of the integrated conditional entropy 
  # criterion for GPs. When `input` is a matrix with more than 1 row, 
  # then the acquisition will be computed in batch mode,
  # meaning that it considers conditioning on the entire 
  # batch of inputs. Note that this is different from the function simply 
  # being vectorized across multiple inputs. For the latter, use 
  # `evaluate_acq_func_vectorized()`. 
  
  N_grid <- nrow(grid_points)
  if(length(weights)==1) weights <- rep(weights, N_grid)
  gp_copy <- gp$copy(shallow=FALSE)
  
  # Condition the GP on the new batch of inputs `input`. Since the conditional entropy 
  # does not depend on the response, the associated batch response is just set to a 
  # vector of zeros. 
  pseudo_response <- matrix(0, nrow=nrow(input), ncol=1)
  gp_copy$update(input, pseudo_response, update_hyperpar=FALSE, ...)
  
  # Evaluate conditional entropy at grid points. 
  integrand_vals <- -1 * acq_neg_entropy(grid_points, gp_copy, ...)
  
  # Return the weighted sum of conditional variances. 
  return(sum(drop(integrand_vals) * weights))
  
}


acq_IEENT_grid <- function(input, gp, grid_points, weights=1/nrow(grid_points), log_scale=TRUE, ...) {
  # This function targets exploration for the exponentiated GP. It can be thought 
  # of as an integrated conditional entropy criterion for 
  # log-normal processes. The outer integral over the design space is approximated with a discrete
  # sum over grid points `grid_points` which are given weights `weights`. 
  # When `input` is a matrix with more than 1 row, then the
  # acquisition will be computed in batch mode, meaning that it considers conditioning
  # on the entire batch of inputs. Note that this is different from the function simply 
  # being vectorized across multiple inputs. For the latter, use 
  # `evaluate_acq_func_vectorized()`. 
  
  # TODO: validate_args_acq_IEENT_grid()
  
  N_grid <- nrow(grid_points)
  if(length(weights)==1) weights <- rep(weights, N_grid)
  gp_copy <- gp$copy(shallow=FALSE)
  
  # Emulator predictions at acquisition evaluation locations.
  pred <- gp_copy$predict(input, return_mean=TRUE, return_var=FALSE, ...)
  
  # Update the GP model, treating the predictive mean as the observed 
  # response at the acquisition evaluation locations. 
  gp_copy$update(input, pred$mean, update_hyperpar=FALSE, ...)
  
  # Predict with the conditional ("cond") GP (i.e., the updated GP) at the grid locations. 
  # Convert to the exponentiated scale to obtain log-normal predictive quantities. 
  pred_cond <- gp_copy$predict(grid_points, return_mean=TRUE, return_var=TRUE, ...)
  
  # Return the weighted sum of conditional entropy evaluations.
  summands <- pred_cond$mean + 0.5*log(pred_cond$var) + 0.5*log(2*pi) + 0.5
  return(sum(drop(summands) * weights))

}


acq_neg_var <- function(input, gp, ...) {
  # Simply returns the negative predictive variance at the input point 
  # `input`. The negative is due to the fact that the framework assumes
  # that acquisition functions are always minimized, so to implement 
  # the "maximum variance" acquisition it must be negated here. 
  
  -1 * gp$predict(input, return_mean=FALSE, return_var=TRUE, ...)$var[,1]
}

acq_neg_entropy <- function(input, gp, ...) {
  # The negative entropy of the GP predictive distribution at the input point
  # `input`. 
  
  -0.5 * log(2*pi*gp$predict(input, return_mean=FALSE, return_var=TRUE, ...)$var[,1]) - 0.5
}


acq_neg_exp_entropy <- function(input, gp, ...) {
  # Returns the negative entropy of the log-normal process 
  # (LNP) `exp(gp)` at the input points `input`. This is technically 
  # the negative entropy (computed with respect to the natural log) 
  # up to a multiplicative constant. Multiplying the returned value
  # by C^2, where C = log2(e) gives the exact negative entropy. 
  
  gp_pred <- gp$predict(input, return_mean=TRUE, return_var=TRUE, ...)
  lnp_pred <- convert_Gaussian_to_LN(mean_Gaussian=gp_pred$mean, var_Gaussian=gp_pred$var,
                                     return_mean=FALSE, return_var=TRUE, log_scale=TRUE)
  
  -lnp_pred$log_mean - 0.5 * lnp_pred$log_var - 0.5*log(2*pi) - 0.5
}


acq_neg_exp_var <- function(input, gp, log_scale=TRUE, ...) {
  # Returns the negative predictive variance of the log-normal process 
  # (LNP) `exp(gp)` at the input points `input`.
  # The negative is due to the fact that the framework assumes
  # that acquisition functions are always minimized, so to implement 
  # the "maximum LNP variance" acquisition it must be negated here. 
  
  gp_pred <- gp$predict(input, return_mean=TRUE, return_var=TRUE, ...)
  lnp_pred <- convert_Gaussian_to_LN(mean_Gaussian=gp_pred$mean, var_Gaussian=gp_pred$var,
                                     return_mean=FALSE, return_var=TRUE, log_scale=TRUE)
  
  if(log_scale) return(-lnp_pred$log_var)
  return(-exp(lnp_pred$log_var))
}


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


get_init_param_estimates <- function(design_info, computer_model_data, sig_eps_prior_params, verbose = TRUE) {
  # Produces estimates for calibration and likelihood parameters given a set of design data and priors on 
  # each set of parameters. The returned estimate for the calibration parameter is selected from the finite 
  # set of design points, so the forward model is not run in this function at all.
  # This function essentially performs a mini coordinate ascent over the discrete set of calibration parameter
  # values in the design, alternating between calibration parameter and likelihood parameter updates. 
  # This function, among other uses, provides a cheap way to set the initial values for MCMC sampling. 
  #
  # Args:
  #    design_info: list, design info list containing input design points and corresponding observed outputs. 
  #    computer_model_data: list, computer model data list. 
  #    sig_eps_prior_params: list, containing prior information on likelihood parameters. 
  #    verbose: logical, whether or not to print results after each iteration. 
  #
  # Returns:
  #    list, containing elements "best_input_idx", "best_input", "best_sig2_eps", and "lpost_max". 
  
  # Start by fixing calibration parameter at value within design that minimizes sum of squared errors across all 
  # output variables. 
  best_input_idx_prev <- which.min(rowSums(design_info$outputs))
  # theta_curr_SSR_min <- design_info$inputs[theta_curr_SSR_min_idx,,drop=FALSE]
  if(verbose) print(paste0("Init idx: ", best_input_idx_prev))
  
  for(i in 1:100) {
    
    # Likelihood parameter update: Closed form MLE estimate of likelihood parameters given calibration parameters. 
    sig2_eps_estimate <- optimize_sig_eps_cond_post(design_info$outputs[best_input_idx_prev,], sig_eps_prior_params, computer_model_data$n_obs)
    
    # Calibration parameter update: Compute log joint posterior density over design inputs with MLE estimate of likelihood 
    # parameters and set calibration parameters to the maximizing value. 
    lpost_design <- calc_lpost_theta_product_lik(theta_vals = design_info$inputs, 
                                                 computer_model_data = computer_model_data, 
                                                 SSR = design_info$outputs, 
                                                 vars_obs = sig2_eps_estimate, 
                                                 na.rm = TRUE, 
                                                 return_list = FALSE,
                                                 theta_prior_params = theta_prior_params)
    best_input_idx_curr <- which.max(lpost_design)
    
    if(verbose) print(paste0("Curr idx: ", best_input_idx_curr, " --- sig2: ", paste0(sig2_eps_estimate, collapse = ",")))
    
    if(best_input_idx_curr == best_input_idx_prev) break
    best_input_idx_prev <- best_input_idx_curr
    
  }
  
  best_input <- design_info$inputs[best_input_idx_curr,,drop=FALSE]
  lpost_max <- lpost_design[best_input_idx_curr]
  
  
  return(list(best_input_idx = best_input_idx_curr, best_input = best_input, best_sig2_eps = sig2_eps_estimate, lpost_max = lpost_max))
  
}

                           
batch_acquisition_opt_one_step <- function(lpost_emulator, acquisition_settings, verbose = TRUE) {
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
    input_new_scaled <- acquisition_opt_one_step(lpost_emulator = lpost_emulator, acquisition_settings = acquisition_settings, verbose = verbose)

    # Update lpost emulator using batch heuristic method. Does not affect underlying GPs, including their hyperparameters. 
    lpost_emulator <- pseudo_update_lpost_emulator(lpost_emulator, input_new_scaled = input_new_scaled,
                                                   pseudo_update_method = acquisition_settings$batch_method)
  }
  
  return(lpost_emulator)
  
}


acquisition_opt_one_step <- function(lpost_emulator, acquisition_settings, verbose = TRUE) {
  # Performs a single one-point acquisition by optimizing the specified acquisition function. Returns 
  # the newly acquired (scaled) design point. 
  #
  # Args:
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`.  
  #    acquisition_settings: list, the acquisition settings list. 
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
                                           verbose = verbose)
                                           
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
                                      N_subsample_candidate = NULL, N_subsample_integrate = NULL, verbose = TRUE) { 
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
                                            N_MC_samples = N_MC_samples, 
                                            theta_grid_integrate = theta_grid_integrate, 
                                            verbose = verbose)
  
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
  
  if(constant_liar_method == "constant_liar_pessimist") {
    return(min(lpost_emulator$outputs_lpost))
  } else if(constant_liar_method == "constant_liar_optimist") {
    return(max(lpost_emulator$outputs_lpost))
  } else if(constant_liar_method == "constant_liar_mean") {
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

acquisition_EI_lpost <- function(theta_vals, lpost_emulator, threshold_lpost = NULL, ...) {
  # Closed-form implementation of the expected improvement (EI) acquisition function applied to an lpost emulator 
  # which is a Gaussian Process. 
  #
  # Args:
  #    theta_vals: matrix of dimension M x D, each row an input location at which to compute the value of the acquisition function. 
  #                The inputs are assumed to already be properly scaled. 
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`. 
  #    threshold_lpost: the threshold objective value to use in the EI calculation, typically the current observed maximum of the 
  #                     unnormalized log posterior. If NULL, set to maximum observed value of the unnormalized log posterior, as 
  #                     stored in `lpost_emulator`. 
  #
  # Returns:
  #    numeric, vector of length equal to length of `theta_vals`; the EI acquisition function evaluations at inputs `theta_vals`. 
  
  # Threshold to use in EI calculation. 
  if(is.null(threshold_lpost)) threshold_lpost <- max(lpost_emulator$outputs_lpost)
  
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


# TODO: speed up by pre-compoting k(X_integrate, X_design). This matrix can then easily be modified when conditioning on an additional 
#       candidate point. 
acquisition_EIVAR_lpost <- function(theta_vals, lpost_emulator, theta_grid_integrate, verbose = TRUE, ...) {
  # Implements the expected integrated variance (EIVAR) criteria that targets the log posterior in the loss emulation setting. 
  # In this case the inner two integrals of EIVAR are available in closed form. The outer integral, the expectation over the input 
  # space is approximated by a finite sum over grid points `theta_grid_integrate`.
  #
  # Args:
  #    theta_vals: matrix of dimension M x D, each row an input location at which to compute the value of the acquisition function. 
  #                The inputs are assumed to already be properly scaled. 
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`.
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
  
  for(i in 1:nrow(theta_vals)) {
    
    # Update variance by conditioning on theta evaluation value. Should not affect `lpost_emulator` outside of local function scope.
    lpost_emulator_temp <- update_lpost_emulator(lpost_emulator, inputs_new_scaled = theta_vals[i,,drop=FALSE], outputs_lpost_new = NULL)
    
    # Compute unnormalized log posterior approximation predictive variance at each theta grid location. 
    lpost_pred_var_grid <- predict_lpost_emulator(inputs_new_scaled = theta_grid_integrate, lpost_emulator = lpost_emulator_temp, return_vals = "var", 
                                                  include_nugget = TRUE, verbose = verbose, unscale = FALSE, uncenter = FALSE)$var
    
    # Estimate EIVAR via discrete sum over theta grid locations. 
    EIVAR_est[i] <- -1.0 * mean(lpost_pred_var_grid)
    
  }
  
  return(EIVAR_est)
  
}


acquisition_EIVAR_lpost_old <- function(theta_vals, emulator_info_list, sig2_eps, theta_grid_integrate, ...) {
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


acquisition_IVAR_post <- function(theta_vals, lpost_emulator, theta_grid_integrate, verbose = TRUE, include_nugget = TRUE, ...) {
  # Implements the expected integrated variance (EIVAR) criteria that targets the log posterior in the loss emulation setting. 
  # In this case the inner two integrals of EIVAR are available in closed form. The outer integral, the expectation over the input 
  # space is approximated by a finite sum over grid points `theta_grid_integrate`. In the code, I use the variable naming 
  # convention "int", "can", and "n" to refer to the integration points (`theta_grid_integrate`), candidate points 
  # (`theta_vals`), and current design inputs (`lpost_emulator$inputs$inputs_scaled`), respectively. 
  #
  # Args:
  #    theta_vals: matrix of dimension M x D, each row an input location at which to compute the value of the acquisition function. 
  #                The inputs are assumed to already be properly scaled. 
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`.
  #    theta_grid_integrate: matrix, of dimension M_integrate x D. The set of inputs used to approximate the integral over the 
  #                          input space required in evaluating the EIVAR criterion. 
  # 
  # Returns:
  #    numeric vector of length M, containing the evaluations of the acquisition at the M input points in `theta_vals`. 
  #    Technically returns the negative of EIVAR to align with the acquisition convention that bigger is better. 
  #
  # TODO: 
  #    - should have a way to also pass in unscaled inputs, or to compute them here. This function ends up unscaling 
  #      multiple times, which is wasteful.
  #    - should pass `include_nugget` to all acquisition functions. 
  #    - in the greedy batch procedure, there are computations I can pass in without re-computing; 
  #      e.g. all of the prior quantities and `K_int_n_Kinv`.
  #    - pull out terms that can be pulled out of the integral. 
  #    - maybe compute unscaled inputs once at the beginning 
  
  n <- nrow(lpost_emulator$inputs_lpost$inputs_scaled)
  N_int <- nrow(theta_grid_integrate)
  
  # Handle case of single input. 
  if(is.null(nrow(theta_vals))) theta_vals <- matrix(theta_vals, nrow = 1)
  
  # Vector to store log IVAR estimates at inputs `theta_vals`. 
  log_IVAR_est <- vector(mode = "numeric", length = nrow(theta_vals))
  
  # Compute lpost prior mean and kernel evaluations. 
  lpost_mu0_int <- calc_lpost_mean(lpost_emulator, inputs_scaled = theta_grid_integrate)
  lpost_k0_int_n <- calc_lpost_kernel(lpost_emulator, theta_grid_integrate, lpost_emulator$inputs_lpost$inputs_scaled)
  lpost_k0_int_can <- calc_lpost_kernel(lpost_emulator, theta_grid_integrate, theta_vals)
  
  # Compute predictive means conditional on current design points. This is used to average over the unknown response when 
  # conditioning on new points in `theta_vals`. 
  lpost_curr_pred_can <- predict_lpost_emulator(theta_vals, lpost_emulator, return_vals = c("mean", "var"), verbose = verbose, 
                                                include_nugget = include_nugget, unscale = TRUE, uncenter = FALSE)

  for(i in 1:nrow(theta_vals)) {
    
    # Update variance by conditioning on theta evaluation value. Should not affect `lpost_emulator` outside of local function scope.
    lpost_emulator_temp <- update_lpost_emulator(lpost_emulator, inputs_new_scaled = theta_vals[i,,drop=FALSE], outputs_lpost_new = NULL)
    repeated_design_input <- ifelse(nrow(lpost_emulator_temp$inputs_lpost$inputs_scaled) == n, TRUE, FALSE)

    # Compute unnormalized log posterior approximation predictive mean and variance at each theta grid location. The variance prediction is exact, 
    # while the mean prediction uses the plug-in kriging believer approach.
    lpost_pred_int <- predict_lpost_emulator(inputs_new_scaled = theta_grid_integrate, lpost_emulator = lpost_emulator_temp, return_vals = c("mean", "var"), 
                                             include_nugget = include_nugget, verbose = verbose, prior_mean_vals_new = lpost_mu0_int, 
                                             unscale = TRUE, uncenter = FALSE)
    lpost_pred_mean_KB_int <- lpost_pred_int$mean
    lpost_pred_var_int <- lpost_pred_int$var
    
    if(!repeated_design_input) {
      # Computations used both in kernel term for conditioned GP and averaged unknown response term. 
      B <- (lpost_k0_int_n %*% lpost_emulator_temp$K_inv[1:n, n+1, drop=FALSE]) + (lpost_k0_int_can[,i,drop=FALSE] * lpost_emulator_temp$K_inv[n+1, n+1]) 
        
      # Term capturing the uncertainty in the value of the unknown response. 
      log_uncertainty_term <- lpost_curr_pred_can$var[i] * drop(B)^2
    }
    
    # Numerically stable calculation of log[exp(k) - 1] term. 
    idx_approx_sel <- (lpost_pred_var_int >= 100)
    log_exp_term <- vector(mode = "numeric", length = length(lpost_pred_var_int))
    log_exp_term[!idx_approx_sel] <- log(exp(lpost_pred_var_int[!idx_approx_sel]) - 1)
    log_exp_term[idx_approx_sel] <- lpost_pred_var_int[idx_approx_sel] # Apply approximation. 
    
    # Sum terms to compute log expected variance for current candidate point over all integration points. 
    if(!repeated_design_input) {
      log_IVAR <- lpost_pred_var_int + log_exp_term + 2 * (lpost_pred_mean_KB_int + lpost_curr_pred_can$var[i] * drop(B)^2)
    } else {
      log_IVAR <- lpost_pred_var_int + log_exp_term + 2 * lpost_pred_mean_KB_int
    }
    
    # Approximate integral by summing values over integration points. Use numerically stable computation 
    # of log-sum-exp. Not normalizing sum, as division by number of integration points does not affect 
    # the optimization over candidate points. 
    log_IVAR_est[i] <- -1.0 * (matrixStats::logSumExp(log_IVAR) - log(N_int))
    
  }
  
  return(log_IVAR_est)
  
}


acquisition_IVAR_post_old <- function(theta_vals, lpost_emulator, theta_grid_integrate, verbose = TRUE, include_nugget = TRUE, ...) {
  # Implements the expected integrated variance (EIVAR) criteria that targets the log posterior in the loss emulation setting. 
  # In this case the inner two integrals of EIVAR are available in closed form. The outer integral, the expectation over the input 
  # space is approximated by a finite sum over grid points `theta_grid_integrate`. In the code, I use the variable naming 
  # convention "int", "can", and "n" to refer to the integration points (`theta_grid_integrate`), candidate points 
  # (`theta_vals`), and current design inputs (`lpost_emulator$inputs$inputs_scaled`), respectively. 
  #
  # Args:
  #    theta_vals: matrix of dimension M x D, each row an input location at which to compute the value of the acquisition function. 
  #                The inputs are assumed to already be properly scaled. 
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`.
  #    theta_grid_integrate: matrix, of dimension M_integrate x D. The set of inputs used to approximate the integral over the 
  #                          input space required in evaluating the EIVAR criterion. 
  # 
  # Returns:
  #    numeric vector of length M, containing the evaluations of the acquisition at the M input points in `theta_vals`. 
  #    Technically returns the negative of EIVAR to align with the acquisition convention that bigger is better. 
  #
  # TODO: 
  #    - should have a way to also pass in unscaled inputs, or to compute them here. This function ends up unscaling 
  #      multiple times, which is wasteful.
  #    - should pass `include_nugget` to all acquisition functions. 
  #    - in the greedy batch procedure, there are computations I can pass in without re-computing; 
  #      e.g. all of the prior quantities and `K_int_n_Kinv`.
  #    - pull out terms that can be pulled out of the integral. 
  
  n <- nrow(lpost_emulator$inputs_lpost$inputs_scaled)
  N_int <- nrow(theta_grid_integrate)
  
  # Handle case of single input. 
  if(is.null(nrow(theta_vals))) theta_vals <- matrix(theta_vals, nrow = 1)
  
  # Vector to store IVAR estimates at inputs `theta_vals`. 
  log_IVAR_est <- vector(mode = "numeric", length = nrow(theta_vals))
  
  # Compute lpost prior mean and kernel evaluations. 
  lpost_mu0_int <- calc_lpost_mean(lpost_emulator, inputs_scaled = theta_grid_integrate)
  lpost_mu0_can <- calc_lpost_mean(lpost_emulator, inputs_scaled = theta_vals)
  lpost_k0_int_n <- calc_lpost_kernel(lpost_emulator, theta_grid_integrate, lpost_emulator$inputs_lpost$inputs_scaled)
  lpost_k0_int_can <- calc_lpost_kernel(lpost_emulator, theta_grid_integrate, theta_vals)
  
  # Compute predictive means conditional on current design points. This is used to average over the unknown response when 
  # conditioning on new points in `theta_vals`. 
  lpost_curr_pred_can <- predict_lpost_emulator(theta_vals, lpost_emulator, return_vals = c("mean", "var"), verbose = verbose, 
                                                include_nugget = include_nugget)
                                            
  for(i in 1:nrow(theta_vals)) {
    
    # Update variance by conditioning on theta evaluation value. Should not affect `lpost_emulator` outside of local function scope.
    lpost_emulator_temp <- update_lpost_emulator(lpost_emulator, inputs_new_scaled = theta_vals[i,,drop=FALSE], outputs_lpost_new = NULL)
    
    # Compute unnormalized log posterior approximation predictive variance at each theta grid location. 
    lpost_pred_var_int <- predict_lpost_emulator(inputs_new_scaled = theta_grid_integrate, lpost_emulator = lpost_emulator_temp, return_vals = "var", 
                                                 include_nugget = include_nugget, verbose = verbose)$var
    
    # Computations used both in kernel term for conditioned GP and averaged unknown response term. 
    A <- (lpost_k0_int_n %*% lpost_emulator_temp$K_inv[1:n, 1:n, drop=FALSE]) + 
         (lpost_k0_int_can[,i,drop=FALSE] %*% lpost_emulator_temp$K_inv[n+1, 1:n, drop=FALSE])
    B <- (lpost_k0_int_n %*% lpost_emulator_temp$K_inv[1:n, n+1, drop=FALSE]) + 
         (lpost_k0_int_can[,i,drop=FALSE] %*% lpost_emulator_temp$K_inv[n+1, n+1])
    
    # Numerically stable calculation of log[exp(k) - 1] term. 
    idx_approx_sel <- (lpost_pred_var_int >= 100)
    log_exp_term <- vector(mode = "numeric", length = length(lpost_pred_var_int))
    log_exp_term[!idx_approx_sel] <- log(exp(lpost_pred_var_int[!idx_approx_sel]) - 1)
    log_exp_term[idx_approx_sel] <- lpost_pred_var_int[idx_approx_sel] # Apply approximation. 
    
    # Observed data term. 
    log_obs_term <- drop(A %*% (matrix(lpost_emulator$outputs_lpost, ncol=1) - lpost_emulator$mu0_design))
    
    # Current predictive mean term (conditioning on current n-point design). 
    log_curr_pred_mean_term <- drop(B * (lpost_curr_pred_can$mean[i] - lpost_mu0_can[i]))
    
    # Current predictive variance term (conditioning on current n-point design).
    log_curr_pred_var_term <- lpost_curr_pred_can$var[i] * drop(B)^2

    # Sum terms to compute log IVAR for current candidate point over all integration points. 
    log_IVAR <- lpost_pred_var_int + log_exp_term + 2 * (lpost_mu0_int + log_obs_term + log_curr_pred_mean_term + log_curr_pred_var_term)
    
    # Approximate integral by summing values over integration points. Use numerically stable computation 
    # of log-sum-exp.
    log_IVAR_est[i] <- -1.0 * (matrixStats::logSumExp(log_IVAR) - log(N_int))
    
  }
  
  return(log_IVAR_est)
  
}


acquisition_IVAR_post_MC <- function(theta_vals, lpost_emulator, theta_grid_integrate, N_MC_samples = 1000,
                                     verbose = TRUE, include_nugget = TRUE, ...) {
  # NOTE: this function is intended only for validating `acquisition_IVAR_post()`. It does not fit into the current 
  # acquisition function code standards, but can easily be updated so that it does. Updates required: 1.) remove the 
  # argument `computer_model_data` which is not an argument accepted by acquisition functions and 2.) generate 
  # samples from a function `sample_lpost_emulator()` (not yet implemented) instead of using `sample_GP_lpost_theta`. 
  
  n <- length(lpost_emulator$outputs_lpost)
  N_int <- nrow(theta_grid_integrate)
  
  # Vector to store log IVAR estimates at inputs `theta_vals`. 
  log_IVAR_est <- vector(mode = "numeric", length = nrow(theta_vals))
  
  # Pre-compute prior mean evaluations at integration points. 
  prior_mean_vals_int <- calc_lpost_mean(lpost_emulator, inputs_scaled = theta_grid_integrate)
  
  for(i in 1:nrow(theta_vals)) {
    
    # Condition lpost emulator on current candidate point.  
    lpost_emulator_temp <- update_lpost_emulator(lpost_emulator, inputs_new_scaled = theta_vals[i,,drop=FALSE], outputs_lpost_new = NULL)
    
    # Compute unnormalized log posterior approximation predictive variance at each theta grid location. The predictive mean will also 
    # be required, but the predictive mean depends on the unknown response values. Thus, the predictive mean is approximated via 
    # Monte Carlo below. 
    lpost_pred_var_int <- predict_lpost_emulator(inputs_new_scaled = theta_grid_integrate, lpost_emulator = lpost_emulator_temp, return_vals = "var", 
                                                 include_nugget = include_nugget, verbose = verbose)$var
    
    # Compute the portion of the predictive lpost (log) variance that does not depend on the predictive mean of the lpost emulator. 
    idx_approx_sel <- (lpost_pred_var_int >= 100)
    log_exp_term <- vector(mode = "numeric", length = length(lpost_pred_var_int))
    log_exp_term[!idx_approx_sel] <- log(exp(lpost_pred_var_int[!idx_approx_sel]) - 1)
    log_exp_term[idx_approx_sel] <- lpost_pred_var_int[idx_approx_sel] # Apply approximation. 
    log_var_term1 <- log_exp_term + lpost_pred_var_int
    
    # Generate samples from current GP predictive distribution at current candidate point, since the post emulator is log-normal and 
    # hence depends on the response value. 
    lpost_samp <- drop(sample_lpost_emulator(inputs_new_scaled = theta_vals[i,,drop=FALSE], 
                                             lpost_emulator = lpost_emulator, N_samples = N_MC_samples))
    
    # Approximating the predictive mean term. 
    log_post_pred_var <- vector(mode = "numeric", length = N_MC_samples * nrow(theta_grid_integrate))
    idx <- 1
    
    for(k in seq_len(N_MC_samples)) {

      # Add sampled response value. 
      lpost_emulator_temp$outputs_lpost[n+1] <- lpost_samp[k]
      
      # Compute predictive mean at integration points, conditional on sampled response value. 
      lpost_pred_mean_int <- predict_lpost_emulator(inputs_new_scaled = theta_grid_integrate, lpost_emulator = lpost_emulator_temp, return_vals = "mean", 
                                                    include_nugget = include_nugget, verbose = verbose, prior_mean_vals_new = prior_mean_vals_int)$mean

      # Append predictive means. 
      end_idx <- idx + length(lpost_pred_mean_int) - 1
      log_post_pred_var[idx:end_idx] <- 2 * lpost_pred_mean_int + log_var_term1
      idx <- end_idx+1
    
    }
    
    # Approximate double integral by summing values over integration points (averaging over input space) and Monte Carlo samples 
    # averaging over unknown response. Use numerically stable computation 
    # of log-sum-exp. Not normalizing sum, as division by number of integration points does not affect 
    # the optimization over candidate points. Do normalize the sum approximating the inner integral so that the output
    # of this function can be compared with `acquisition_IVAR_post()`. 
    log_IVAR_est[i] <- -1.0 * (matrixStats::logSumExp(log_post_pred_var) - log(N_MC_samples * N_int))
    
  }

  return(log_IVAR_est)
  
}


acquisition_VAR_lpost <- function(theta_vals, lpost_emulator, ...) {
  # Implements the acquisition which is simply defined as the variance of the unnormalized log posterior approximation in  
  # the loss emulation setting. 
  #
  # Args:
  #    theta_vals: matrix of dimension M x D, each row an input location at which to compute the value of the acquisition function. 
  #                The inputs are assumed to already be properly scaled. 
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`.
  # 
  # Returns:
  #    numeric vector of length M, containing the evaluations of the acquisition at the M input points in `theta_vals`. 
  
  # Handle case of single input. 
  if(is.null(nrow(theta_vals))) theta_vals <- matrix(theta_vals, nrow = 1)
  
  # Compute unnormalized log posterior approximation predictive variance at each theta grid location.
  lpost_pred_var_grid <- predict_lpost_emulator(inputs_new_scaled = theta_vals, lpost_emulator = lpost_emulator, return_vals = "var", 
                                                include_nugget = TRUE)
  
  return(lpost_pred_var_grid$var)
  
}


acquisition_VAR_lpost_old <- function(theta_vals, emulator_info_list, sig2_eps, ...) {
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


acquisition_VAR_post <- function(theta_vals, lpost_emulator, ...) {
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
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`. 
  # 
  # Returns:
  #    numeric vector of length M, containing the evaluations of the acquisition at the M input points in `theta_vals`. 
  
  # Handle case of single input. 
  if(is.null(nrow(theta_vals))) theta_vals <- matrix(theta_vals, nrow = 1)
  
  # Compute unnormalized log posterior approximation predictive variance at each theta grid location.
  lpost_pred_list <- predict_lpost_emulator(inputs_new_scaled = theta_vals, lpost_emulator = lpost_emulator, include_nugget = TRUE)
  mu <- lpost_pred_list$mean
  sig2 <- lpost_pred_list$var
  
  # Numerically stable computation of the log predictive variance. 
  idx_approx_sel <- (sig2 >= 100)
  term_1 <- vector(mode = "numeric", length = length(sig2))
  term_1[!idx_approx_sel] <- log(exp(sig2[!idx_approx_sel]) - 1)
  term_1[idx_approx_sel] <- sig2[idx_approx_sel] # Apply approximation. 
  
  log_pred_var_post <- term_1 + 2*mu + sig2
  
  return(log_pred_var_post)
  
}


acquisition_VAR_post_old <- function(theta_vals, emulator_info_list, computer_model_data, theta_prior_params, sig2_eps, ...) {
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


acquisition_IEVAR_post <- function(theta_vals, lpost_emulator, theta_grid_integrate, verbose = TRUE, 
                                   include_nugget = TRUE, ...) {
  
  log_IEVAR_vals <- vector(mode = "numeric", length = nrow(theta_vals))
  
  for(i in seq_len(nrow(theta_vals))) {
    log_IEVAR_vals[i] <- matrixStats::logSumExp(calc_log_EVAR_post(lpost_emulator, theta_vals[i,,drop=FALSE], theta_grid_integrate))
  }
  
  return(-1.0 * (log_IEVAR_vals - log(nrow(theta_grid_integrate))))
  
}


calc_log_EVAR_post <- function(lpost_emulator, input_candidate, inputs_integrate) {
  
  # TODO:
  #    - Add argument checks. 
  #    - Walk through and pull out steps that can pre-computed for efficiency (will definitely want to pass in 
  #      `prior_mean_vals_integrate` to this function to be passed to `predict_lpost_emulator`). 
  
  # Predictions using GP conditioned on current design.  
  pred_curr_cand <- predict_lpost_emulator(input_candidate, lpost_emulator, return_vals = c("mean", "var", "cov"),
                                           inputs_new_scaled_2 = inputs_integrate, unscale = TRUE, uncenter = FALSE)

  # Update GP using kriging believer approach. 
  lpost_emulator_KB <- update_lpost_emulator(lpost_emulator, inputs_new_scaled = input_candidate, 
                                             outputs_lpost_new = pred_curr_cand$mean, outputs_scaled = FALSE, 
                                             outputs_centered = TRUE)
  
  # Predict with kriging believer GP. 
  pred_KB_int <- predict_lpost_emulator(inputs_integrate, lpost_emulator_KB, return_vals = c("mean", "var"),
                                        unscale = TRUE, uncenter = FALSE)
  
  # Compute log EVAR evaluations at the integration inputs. 
  inflation_factor <- 2 * drop(pred_curr_cand$cov)^2 / pred_curr_cand$var
  log_EVAR_vals <- convert_to_post_emulator_log_moments(pred_KB_int$mean, pred_KB_int$var, return_vals = "log_var")$log_var + inflation_factor
  
  return(log_EVAR_vals)
  
}


acquisition_IEVAR_post_test_mode <- function(theta_vals, lpost_emulator, theta_grid_integrate, verbose = TRUE, 
                                   include_nugget = TRUE, ...) {
  
  log_IEVAR_vals <- vector(mode = "numeric", length = nrow(theta_vals))
  avg_log_KB_vars <- vector(mode = "numeric", length = nrow(theta_vals))
  avg_inflation_term_vals <- vector(mode = "numeric", length = nrow(theta_vals))
  
  for(i in seq_len(nrow(theta_vals))) {
    IEVAR_info <- calc_log_EVAR_post_test_mode(lpost_emulator, theta_vals[i,,drop=FALSE], theta_grid_integrate)
    log_IEVAR_vals[i] <- matrixStats::logSumExp(IEVAR_info$log_EVAR) - log(nrow(theta_grid_integrate))
    avg_log_KB_vars[i] <- matrixStats::logSumExp(IEVAR_info$log_var_term) - log(nrow(theta_grid_integrate))
    avg_inflation_term_vals[i] <- matrixStats::logSumExp(IEVAR_info$inflation_factor) - log(nrow(theta_grid_integrate))
  }
  
  return(list(log_IEVAR_vals = -log_IEVAR_vals, avg_log_KB_vars = -avg_log_KB_vars, avg_inflation_term_vals = -avg_inflation_term_vals))
  
}


calc_log_EVAR_post_test_mode <- function(lpost_emulator, input_candidate, inputs_integrate) {
  
  # TODO:
  #    - Add argument checks. 
  #    - Walk through and pull out steps that can pre-computed for efficiency (will definitely want to pass in 
  #      `prior_mean_vals_integrate` to this function to be passed to `predict_lpost_emulator`). 
  
  # Predictions using GP conditioned on current design.  
  pred_curr_cand <- predict_lpost_emulator(input_candidate, lpost_emulator, return_vals = c("mean", "var", "cov"),
                                           inputs_new_scaled_2 = inputs_integrate, unscale = TRUE, uncenter = FALSE)
  
  # Update GP using kriging believer approach. 
  lpost_emulator_KB <- update_lpost_emulator(lpost_emulator, inputs_new_scaled = input_candidate, 
                                             outputs_lpost_new = pred_curr_cand$mean, outputs_scaled = FALSE, 
                                             outputs_centered = TRUE)
  
  # Predict with kriging believer GP. 
  pred_KB_int <- predict_lpost_emulator(inputs_integrate, lpost_emulator_KB, return_vals = c("mean", "var"),
                                        unscale = TRUE, uncenter = FALSE)
  
  # Compute log EVAR evaluations at the integration inputs. 
  inflation_factor <- 2 * drop(pred_curr_cand$cov)^2 / pred_curr_cand$var
  log_var_term <- convert_to_post_emulator_log_moments(pred_KB_int$mean, pred_KB_int$var, return_vals = "log_var")$log_var
  log_EVAR_vals <- log_var_term + inflation_factor
  
  # Further specify the quantities making up the inflation term. 
  kn_u <- predict_lpost_emulator(inputs_integrate, lpost_emulator, return_vals = "var")$var
  rho_n <- drop(pred_curr_cand$cov) / (sqrt(kn_u * pred_curr_cand$var))
  
  return(list(log_EVAR = log_EVAR_vals, log_var_term = log_var_term, inflation_factor = inflation_factor, 
              kn_u = kn_u, rho_n = rho_n))
  
}


get_IVAR_post_vars <- function(theta_candidate, lpost_emulator, theta_grid_integrate, verbose = TRUE, include_nugget = TRUE, ...) {
  # A function for testing purposes. Instead of integrating over the computed variance values, this function returns the 
  # vector of values. 
  
  n <- nrow(lpost_emulator$inputs_lpost$inputs_scaled)
  N_int <- nrow(theta_grid_integrate)
  
  # Handle case of single input. 
  if(is.null(nrow(theta_candidate))) theta_candidate <- matrix(theta_candidate, nrow = 1)
  
  # Compute lpost prior mean and kernel evaluations. 
  lpost_mu0_int <- calc_lpost_mean(lpost_emulator, inputs_scaled = theta_grid_integrate)
  lpost_k0_int_n <- calc_lpost_kernel(lpost_emulator, theta_grid_integrate, lpost_emulator$inputs_lpost$inputs_scaled)
  lpost_k0_int_can <- calc_lpost_kernel(lpost_emulator, theta_grid_integrate, theta_candidate)
  
  # Compute predictive means conditional on current design points. This is used to average over the unknown response when 
  # conditioning on new points in `theta_vals`. 
  lpost_curr_pred_can <- predict_lpost_emulator(theta_candidate, lpost_emulator, return_vals = c("mean", "var"), verbose = verbose, 
                                                include_nugget = include_nugget, unscale = TRUE, uncenter = FALSE)
  
  # Update variance by conditioning on theta evaluation value. Should not affect `lpost_emulator` outside of local function scope.
  lpost_emulator_temp <- update_lpost_emulator(lpost_emulator, inputs_new_scaled = theta_candidate, outputs_lpost_new = NULL)
  repeated_design_input <- ifelse(nrow(lpost_emulator_temp$inputs_lpost$inputs_scaled) == n, TRUE, FALSE)
    
  # Compute unnormalized log posterior approximation predictive mean and variance at each theta grid location. The variance prediction is exact, 
  # while the mean prediction uses the plug-in kriging believer approach.
  lpost_pred_int <- predict_lpost_emulator(inputs_new_scaled = theta_grid_integrate, lpost_emulator = lpost_emulator_temp, return_vals = c("mean", "var"), 
                                           include_nugget = include_nugget, verbose = verbose, prior_mean_vals_new = lpost_mu0_int, 
                                           unscale = TRUE, uncenter = FALSE)
  lpost_pred_mean_KB_int <- lpost_pred_int$mean
  lpost_pred_var_int <- lpost_pred_int$var
    
  if(!repeated_design_input) {
    # Computations used both in kernel term for conditioned GP and averaged unknown response term. 
    B <- (lpost_k0_int_n %*% lpost_emulator_temp$K_inv[1:n, n+1, drop=FALSE]) + (lpost_k0_int_can * lpost_emulator_temp$K_inv[n+1, n+1]) 
      
    # Term capturing the uncertainty in the value of the unknown response. 
    log_uncertainty_term <- lpost_curr_pred_can$var * drop(B)^2
  }
    
  # Numerically stable calculation of log[exp(k) - 1] term. 
  idx_approx_sel <- (lpost_pred_var_int >= 100)
  log_exp_term <- vector(mode = "numeric", length = length(lpost_pred_var_int))
  log_exp_term[!idx_approx_sel] <- log(exp(lpost_pred_var_int[!idx_approx_sel]) - 1)
  log_exp_term[idx_approx_sel] <- lpost_pred_var_int[idx_approx_sel] # Apply approximation. 
    
  # Sum terms to compute log expected variance for current candidate point over all integration points. 
  if(!repeated_design_input) {
    log_IVAR <- lpost_pred_var_int + log_exp_term + 2 * (lpost_pred_mean_KB_int + log_uncertainty_term)
  } else {
    log_IVAR <- lpost_pred_var_int + log_exp_term + 2 * lpost_pred_mean_KB_int
  }

  return(log_IVAR)
  
}




