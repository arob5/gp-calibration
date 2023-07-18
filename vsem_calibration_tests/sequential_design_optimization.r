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
                           

batch_acquisition_opt_one_step <- function(emulator_info_list, acquisition_settings, design_input_curr, design_objective_curr,
                                           design_best_idx, computer_model_data, sig2_eps, theta_prior_params, theta_grid_ref = NULL) {

  # Deep copy of emulator. 
  # TODO: write this function. 
  gp_fits <- copy_independent_GPs(emulator_info_list$gp_fits)
  
  # Object to store batch of inputs. 
  inputs_new <- matrix(nrow = acquisition_settings$batch_size, ncol = ncol(design_input_curr))
  
  # Acquire batch of input points (without running forward model). 
  for(b in 1:acquisition_settings$batch_size) {

    # Optimize sequential acquisition function.  
    theta_new <- acquisition_opt_one_step(emulator_info_list = emulator_info_list, 
                                          acquisition_settings = acquisition_settings, 
                                          design_input_curr = design_input_curr, 
                                          design_objective_curr = design_objective_curr, 
                                          design_best_idx = design_best_idx, 
                                          computer_model_data = computer_model_data, 
                                          sig2_eps = sig2_eps, 
                                          theta_prior_params = theta_prior_params, 
                                          theta_grid_ref = theta_grid_ref)
    inputs_new[b,] <- theta_new
    
    # Update GP (using e.g. kriging believer, constant liar, etc.) 
    # TODO: write this function
    emulator_info_list$gp_fits <- pseudo_update_independent_GPs(gp_fits = emulator_info_list$gp_fits, 
                                                                gp_lib = emulator_settings$gp_lib, 
                                                                X_new = theta_new, 
                                                                pseudo_update_method = acquisition_settings$batch_heuristic, 
                                                                input_bounds = emulator_info_list$input_bounds, 
                                                                output_stats = emulator_info_list$output_stats)
    
  }
  
  return(inputs_new)
  
}


acquisition_opt_one_step <- function(emulator_info_list, acquisition_settings, design_input_curr, design_objective_curr, design_best_idx,  
                                     computer_model_data, sig2_eps, theta_prior_params, theta_grid_opt = NULL, theta_grid_integrate = NULL) {
  # Note that this is all conditional on fixed likelihood parameters.                               
  

  # Select next point by optimizing acquisition function.
  if(acquisition_settings$opt_method == "grid") {
    theta_new <- optimize_acquisition(acquisition_type = acquisition_settings$acquisition_type, 
                                      opt_method = acquisition_settings$opt_method, 
                                      emulator_info_list = emulator_info_list, 
                                      computer_model_data = computer_model_data, 
                                      theta_prior_params = theta_prior_params, 
                                      sig2_eps = sig2_eps, 
                                      design_input_curr = design_input_curr, 
                                      design_objective_curr = design_objective_curr,
                                      design_best_idx = design_best_idx,
                                      theta_grid_opt = theta_grid_opt, 
                                      theta_grid_integrate = theta_grid_integrate, 
                                      N_MC_samples = acquisition_settings$N_MC_samples)
  } else {
    stop("Invalid acquisition optimization method: ", opt_method)
  }
                              
  return(theta_new)

}


optimize_acquisition_grid <- function(acquisition_type, theta_grid_opt, emulator_info_list, computer_model_data, 
                                      theta_prior_params, sig2_eps, design_input_curr, design_objective_curr, 
                                      design_best_idx, N_MC_samples = NULL, theta_grid_integrate = NULL) {
  
  # Acquisition functions all assume that inputs are already scaled. Scale them here. 
  theta_grid_opt <- scale_input_data(theta_grid_opt, input_bounds = emulator_info_list$input_bounds)
  if(!is.null(theta_grid_integrate)) theta_grid_integrate <- scale_input_data(theta_grid_integrate, input_bounds = emulator_info_list$input_bounds) 
  
  # Get acquisition function. 
  acquisition_func <- get(paste0("acquisition_", acquisition_type))
  
  # Evaluate acquisition function on grid of reference inputs. 
  acquisition_vals_grid <- acquisition_func(theta_vals = theta_grid_opt, 
                                            emulator_info_list = emulator_info_list,
                                            computer_model_data = computer_model_data, 
                                            theta_prior_params = theta_prior_params, 
                                            sig2_eps = sig2_eps, 
                                            design_input_curr = design_input_curr,
                                            design_objective_curr = design_objective_curr, 
                                            design_best_idx = design_best_idx,
                                            N_MC_samples = N_MC_samples, 
                                            theta_grid_integrate = theta_grid_integrate)
  
  # Select input in reference grid that maximizes the acquisition. 
  max_idx <- which.max(acquisition_vals_grid)

  return(theta_grid_ref[max_idx,])

}


optim.EI <- function(f, ninit, end)
{
  ## initialization
  X <- randomLHS(ninit, 2)
  y <- f(X)
  gpi <- newGPsep(X, y, d=0.1, g=1e-6, dK=TRUE)
  da <- darg(list(mle=TRUE, max=0.5), randomLHS(1000, 2))
  mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
  
  ## optimization loop of sequential acquisitions
  maxei <- c()
  for(i in (ninit+1):end) {
    solns <- EI.search(X, y, gpi)
    m <- which.max(solns$val)
    maxei <- c(maxei, solns$val[m])
    xnew <- as.matrix(solns[m,3:4])
    ynew <- f(xnew)
    updateGPsep(gpi, xnew, ynew)
    mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
    X <- rbind(X, xnew)
    y <- c(y, ynew)
  }
  
  ## clean up and return
  deleteGPsep(gpi)
  return(list(X=X, y=y, maxei=maxei))
}


# ---------------------------------------------------------------------------------
# Acquisition Functions:
#    - For both optimization and design. 
#    - All acquisition functions assume that arguments related to the input/ 
#      calibration space (e.g. `theta_vals`, `theta_grid_integrate`) are scaled.
# ---------------------------------------------------------------------------------

acquisition_EI_MC <- function(theta_vals, emulator_info_list, computer_model_data, 
                              theta_prior_params, sig2_eps, design_input_curr, design_objective_curr, 
                              design_best_idx, N_MC_samples = 1000, gp_pred_list = NULL) {
  # Note that lpost_min should be the maximum posterior over both theta and Sigma, not over the conditional posterior. 
  # `theta_vals` should be scaled. 
  
  # Have this return a matrix of dim N_MC_samples x nrow(theta_vals)
  lpost_samp <- sample_GP_lpost_theta(theta_vals_scaled = theta_vals, 
                                    emulator_info_list = emulator_info_list,
                                    computer_model_data = computer_model_data, 
                                    theta_prior_params = theta_prior_params, 
                                    sig2_eps = sig2_eps, 
                                    N_samples = N_MC_samples, 
                                    gp_pred_list = gp_pred_list)
                                    
  # Monte Carlo estimates of acquisition at each point. 
  alpha_EI_MC_estimates <- lpost_samp - design_objective_curr[design_best_idx]
  alpha_EI_MC_estimates[alpha_EI_MC_estimates < 0] <- 0
  
  return(colMeans(alpha_EI_MC_estimates))
  
}


acquisition_PI_MC <- function(theta_vals, emulator_info_list, computer_model_data, 
                              theta_prior_params, sig2_eps, design_input_curr, design_objective_curr,
                              design_best_idx, N_MC_samples = 1000, gp_pred_list = NULL) {
  # Note that lpost_min should be the maximum posterior over both theta and Sigma, not over the conditional posterior. 
  # `theta_vals` should be scaled. 
  
  # Have this return a matrix of dim N_MC_samples x nrow(theta_vals)
  lpost_samp <- sample_GP_lpost_theta(theta_vals_scaled = theta_vals, 
                                      emulator_info_list = emulator_info_list,
                                      computer_model_data = computer_model_data, 
                                      theta_prior_params = theta_prior_params, 
                                      sig2_eps = sig2_eps, 
                                      N_samples = N_MC_samples, 
                                      gp_pred_list = gp_pred_list)
  
  # Monte Carlo estimates of acquisition at each point. 
  alpha_EI_MC_estimates <- lpost_samp - design_objective_curr[design_best_idx]
  alpha_EI_MC_estimates <- colMeans(alpha_EI_MC_estimates > 0)

  return(alpha_EI_MC_estimates)
  
}


# TODO: I think it makes more sense to scale inputs prior to passing to acq functions, so all acq functions assume they are already 
#       scaled. Look into updating this and where the inputs should be scaled. 
acquisition_EIVAR_lpost <- function(theta_vals, emulator_info_list, theta_prior_params, sig2_eps, theta_grid_integrate = theta_vals) {
  # Implements the expected integrated variance (EIVAR) criteria that targets the log posterior in the loss emulation setting. 
  # In this case the inner two integrals of EIVAR are available in closed form. The outer integral, the expectation over the input 
  # space is approximated by a finite sum over grid points `theta_grid_integrate`. By default, the grid values used are the same 
  # values at which the acquisition function is to be evaluated; in the case that the acquisition is only being evaluated at one 
  # or a handful of values, then this default for `theta_grid_integrate` should certainly be overwritten. 
  
  # Handle case of single input. 
  if(is.null(nrow(theta_vals))) theta_vals <- matrix(theta_vals, nrow = 1)
  
  # Vector to store EIVAR estimates at inputs `theta_vals`. 
  EIVAR_est <- vector(mode = "numeric", length = nrow(theta_vals))
  
  for(i in 1:nrow(theta_vals)) {
    
    # Update variance by conditioning on theta evaluation value.
    gp_fits_conditioned <- update_independent_GPs(gp_fits = emulator_info_list$gp_fits, gp_lib = emulator_info_list$settings$gp_lib, 
                                                  X_new = theta_vals[i,,drop=FALSE], Y_new = NULL, update_hyperparamters = FALSE)
                                                  
    # Compute unnormalized log posterior approximation predictive variance at each theta grid location. 
    lpost_pred_var_grid <- predict_lpost_GP_approx(theta_vals_scaled = theta_grid_integrate, emulator_info_list = emulator_info_list, 
                                                   sig2_eps = sig2_eps, include_nugget = TRUE, include_sig_eps_prior = FALSE)
    
    # Estimate EIVAR via discrete sum over theta grid locations. 
    EIVAR_est[i] <- mean(lpost_pred_var_grid)
    
  }
  
  return(EIVAR_est)
  
}










