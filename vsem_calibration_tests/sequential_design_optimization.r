#
# sequential_design_optimization.r
# Functions related to sequential design for GPs and Bayesian Optimization. 
#
# Andrew Roberts
#

Bayes_opt <- function(Bayes_opt_settings, init_design_settings, emulator_settings, computer_model_data, 
                      sig2_eps, theta_prior_params, theta_grid_ref = NULL) {
  # TODO: generalize so this works with log-normnal process. 
  
  # Create initial design and fit GP. 
  design_info <- get_input_output_design(N_points = init_design_settings$N_design, 
                                         design_method = init_design_settings$design_method, 
                                         computer_model_data = computer_model_data, 
                                         theta_prior_params = theta_prior_params, 
                                         transformation_method = emulator_settings$transformation_method)
  design_inputs <- design_info$inputs
  
  gp_fits <- fit_independent_GPs(X_train = design_info$inputs_scaled, 
                                 Y_train = design_info$outputs_normalized, 
                                 gp_lib = emulator_settings$gp_lib, 
                                 gp_kernel = emulator_settings$kernel)$fits
  
  emulator_info_list <- list(gp_fits = gp_fits, 
                             input_bounds = design_info$input_bounds, 
                             output_stats = design_info$output_stats, 
                             settings = emulator_settings)
  
  # Current observed objective values (i.e. log posterior values). 
  design_objective_vals <- calc_lpost_theta_product_lik(theta_vals = design_inputs, 
                                                        computer_model_data = computer_model_data,  
                                                        SSR = design_info$outputs, 
                                                        vars_obs = sig2_eps, 
                                                        na.rm = TRUE, 
                                                        return_list = FALSE, 
                                                        theta_prior_params = theta_prior_params)
  design_best_idx <- which.max(design_objective_vals)
  
  # Bayesian Optimization loop. 
  for(i in seq_len(Bayes_opt_settings$N_opt_iter)) {
    
    design_idx_curr <- init_design_settings$N_design + i
    
    # Obtain new design point and corresponding objective value. 
    opt_results <- bayes_opt_one_step(emulator_info_list = emulator_info_list, 
                                      Bayes_opt_settings = Bayes_opt_settings, 
                                      design_input_curr = design_inputs, 
                                      design_objective_curr = design_objective_vals, 
                                      design_best_idx = design_best_idx[i],  
                                      computer_model_data = computer_model_data, 
                                      sig2_eps = sig2_eps, 
                                      theta_prior_params = theta_prior_params, 
                                      theta_grid_ref = theta_grid_ref)
    
    design_inputs <- opt_results$input
    design_objective_vals <- opt_results$objective
    design_best_idx <- c(design_best_idx, opt_results$idx)

    # Update GP. 
    emulator_info_list$gp_fits <- update_independent_GPs(gp_fits = emulator_info_list$gp_fits, 
                                                         gp_lib = emulator_settings$gp_lib, 
                                                         X_new = opt_results$input[design_idx_curr,,drop=FALSE], 
                                                         Y_new = opt_results$SSR_new, 
                                                         input_bounds = emulator_info_list$input_bounds, 
                                                         output_stats = emulator_info_list$output_stats) 
  }
  
  return(list(design_inputs = design_inputs, objective_vals = design_objective_vals, 
              best_idx = design_best_idx, emulator_info_list = emulator_info_list))
  
}


bayes_opt_one_step <- function(emulator_info_list, Bayes_opt_settings, design_input_curr, design_objective_curr, design_best_idx,  
                               computer_model_data, sig2_eps, theta_prior_params, theta_grid_ref = NULL) {
  # Note that this is all conditional on fixed likelihood parameters.                               
  
  # Select next point by optimizing acquisition function.
  theta_new <- optimize_acquisition(acquisition_type = Bayes_opt_settings$acquisition_type, 
                                    opt_method = Bayes_opt_settings$opt_method, 
                                    emulator_info_list = emulator_info_list, 
                                    computer_model_data = computer_model_data, 
                                    theta_prior_params = theta_prior_params, 
                                    sig2_eps = sig2_eps, 
                                    design_input_curr = design_input_curr, 
                                    design_objective_curr = design_objective_curr,
                                    design_best_idx = design_best_idx,
                                    theta_grid_ref = theta_grid_ref, 
                                    N_MC_samples = Bayes_opt_settings$N_MC_samples)
                                    
  # Run forward model at new point. 
  SSR_new <- get_computer_model_SSR(computer_model_data = computer_model_data, theta_vals = theta_new, na.rm = TRUE)
  lpost_new <- calc_lpost_theta_product_lik(computer_model_data = computer_model_data, 
                                            theta_vals = theta_new, 
                                            SSR = SSR_new,
                                            vars_obs = sig2_eps, 
                                            na.rm = TRUE, 
                                            theta_prior_params = theta_prior_params, 
                                            return_list = FALSE)
  
  # Add new point to design. 
  design_input_curr <- rbind(design_input_curr, matrix(theta_new, nrow=1))
  design_objective_curr <- c(design_objective_curr, lpost_new)
  if(lpost_new > design_objective_curr[design_best_idx]) design_best_idx <- length(design_objective_curr)
  
  return(list(input = design_input_curr, objective = design_objective_curr, 
              idx = design_best_idx, SSR_new = SSR_new))

}


optimize_acquisition <- function(acquisition_type, opt_method, emulator_info_list, computer_model_data, 
                                 theta_prior_params, sig2_eps, design_input_curr, design_objective_curr, 
                                 design_best_idx, theta_grid_ref = NULL, N_MC_samples = NULL) {
  
  if(opt_method == "grid") {
    theta_new <- optimize_acquisition_grid(acquisition_type = acquisition_type, 
                                           emulator_info_list = emulator_info_list, 
                                           computer_model_data = computer_model_data, 
                                           theta_prior_params = theta_prior_params, 
                                           sig2_eps = sig2_eps, 
                                           design_input_curr = design_input_curr,
                                           design_objective_curr = design_objective_curr, 
                                           design_best_idx = design_best_idx,
                                           theta_grid_ref = theta_grid_ref, 
                                           N_MC_samples = N_MC_samples)
  } else {
    stop("Invalid acquisition optimization method: ", opt_method)
  }
  
  return(theta_new)
  
}


optimize_acquisition_grid <- function(acquisition_type, theta_grid_ref, emulator_info_list, computer_model_data, 
                                      theta_prior_params, sig2_eps, design_input_curr, design_objective_curr, 
                                      design_best_idx, N_MC_samples) {
  
  # Get acquisition function. 
  acquisition_func <- get(paste0("acquisition_", acquisition_type))
  
  # Evaluate acquisition function on grid of reference inputs. 
  acquisition_vals_grid <- acquisition_func(theta_grid_ref = theta_grid_ref, 
                                            emulator_info_list = emulator_info_list,
                                            computer_model_data = computer_model_data, 
                                            theta_prior_params = theta_prior_params, 
                                            sig2_eps = sig2_eps, 
                                            design_input_curr = design_input_curr,
                                            design_objective_curr = design_objective_curr, 
                                            design_best_idx = design_best_idx,
                                            N_MC_samples = N_MC_samples)
  
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


# ------------------------------------------------------------------------------
# Acquisition Functions:
#    - For both optimization and design. 
# ------------------------------------------------------------------------------

acquisition_EI_MC <- function(theta_grid_ref, emulator_info_list, computer_model_data, 
                              theta_prior_params, sig2_eps, design_input_curr, design_objective_curr, 
                              design_best_idx, N_MC_samples = 1000, gp_pred_list = NULL) {
  # Note that lpost_min should be the maximum posterior over both theta and Sigma, not over the conditional posterior. 
  # `theta_vals` should be unscaled. 
  
  # Have this return a matrix of dim N_MC_samples x nrow(theta_vals)
  lpost_samp <- samp_GP_lpost_theta(theta_vals = theta_grid_ref, 
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


acquisition_PI_MC <- function(theta_grid_ref, emulator_info_list, computer_model_data, 
                              theta_prior_params, sig2_eps, design_input_curr, design_objective_curr,
                              design_best_idx, N_MC_samples = 1000, gp_pred_list = NULL) {
  # Note that lpost_min should be the maximum posterior over both theta and Sigma, not over the conditional posterior. 
  # `theta_vals` should be unscaled. 
  
  # Have this return a matrix of dim N_MC_samples x nrow(theta_vals)
  lpost_samp <- samp_GP_lpost_theta(theta_vals = theta_grid_ref, 
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












