#
# sequential_design_optimization.r
# Functions related to sequential design for GPs and Bayesian Optimization. 
#
# Andrew Roberts
#

bayes_opt_one_step <- function(emulator_info_list, sig2_eps, ) {
  
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

acquisition_EI_MC <- function(emulator_info_list, theta_vals, lpost_min, computer_model_data, theta_prior_params, sig2_eps, N_MC_samples = 1000) {
  # Note that lpost_min should be the minimum posterior over both theta and Sigma, not over the conditional posterior. 
  
  
  lprior_vals <- apply(theta_vals, 1, function(theta) calc_lprior_theta(theta, theta_prior_params))
  theta_vals <- scale_input_data(theta_vals, input_bounds = emulator_info_list$input_bounds)
  
  gp_pred_list <- predict_independent_GPs(X_pred = theta_vals, 
                                          gp_obj_list = emulator_info_list$gp_fits,  
                                          gp_lib = emulator_info_list$settings$gp_lib, 
                                          denormalize_predictions = TRUE, 
                                          output_stats = emulator_info_list$output_stats)

  SSR_samp <- sample_independent_GPs_pointwise(gp_pred_list = gp_pred_list, 
                                               transformation_methods = emulator_info_list$settings$transformation_method,
                                               include_nugget = TRUE)
  
  lpost_vals <- calc_lpost_theta_product_lik(computer_model_data = computer_model_data, 
                                             SSR = SSR_samp, 
                                             vars_obs = sig2_eps, 
                                             na.rm = TRUE, 
                                             theta_prior_params = theta_prior_params, 
                                             lprior_vals = lprior_vals,
                                             return_list = FALSE)
  
  
  ###########
  
  # Have this return a matrix of dim N_MC_samples x nrow(theta_vals)
  lpost_samp <- samp_GP_lpost_theta(theta_vals = theta_vals, 
                                    emulator_info_list = emulator_info_list,
                                    computer_model_data = computer_model_data, 
                                    theta_prior_params = theta_prior_params, 
                                    sig2_eps = sig2_eps, 
                                    N_samples = N_MC_samples)
  
  
  # Monte Carlo estimates of acquisition at each point. 
  alpha_EI_MC_estimates <- colMeans(pmax(0, lpost_min - lpost_samp))
  
  return(alpha_EI_MC_estimates)
  
}














