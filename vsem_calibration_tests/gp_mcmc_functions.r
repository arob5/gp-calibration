#
# gp_mcmc_functions.r
# MCMC algorithms utilizing a Gaussian process (GP) emulator to accelerate inference. 
#

# -----------------------------------------------------------------------------
# Helper Functions. 
# -----------------------------------------------------------------------------

sample_emulator_cond <- function(input_scaled, emulator_info_list, cond_type, sig2_eps = NULL) {
  # Draws a sample either from p(phi|u) or p(phi|u,Sigma,Y). The former is the GP predictive 
  # distribution for the model-data misfit evaluated at input u. The latter is the conditional 
  # posterior of the model-data misfit value, given input u and variance parameters Sigma. 
  # Currently only supports drawing a single sample, and always includes the nugget variances.
  # Also, sets a floor of zero on the sampled values to prevent negative values. 
  #
  # Args:
  #    input_scaled: numeric(d), the calibration parameter value that serves as the input to the GPs.
  #    emulator_info_list: list, the emulator info list. 
  #    cond_type: character(1), specifies which conditional to sample from. If "prior", 
  #               samples from the GP predictive distribution p(phi|u). If "post", samples
  #               from the GP conditional posterior distribution p(phi|u,Sigma,Y).
  #    sig2_eps: numeric(p), the vector of variance parameters. Only required if cond_type is "post". 
  #
  # Returns:
  #    matrix, of dimension 1xp, with p being the number of output variables. 
  
  # Compute GP predictions. 
  gp_pred_list <- predict_independent_GPs(X_pred = input_scaled, gp_obj_list = emulator_info_list$gp_fits, 
                                          gp_lib = emulator_info_list$settings$gp_lib, include_cov_mat = FALSE, 
                                          denormalize_predictions = TRUE, output_stats = emulator_info_list$output_stats)
  
  # Sample from GP predictive distribution. 
  SSR_samples <- sample_independent_GPs_pointwise(gp_pred_list, transformation_methods = emulator_info_list$settings$transformation_method, 
                                                  include_nugget = TRUE)
  
  # Shift mean if sample from conditional posterior is desired. 
  if(cond_type == "post") {
    gp_pred_vars <- sapply(gp_pred_list, function(x) x$var_comb)
    SSR_samples <- pmax(SSR_samples - matrix(0.5*(1/sig2_eps)*gp_pred_vars, nrow=1), 0)
  }
  
  return(SSR_samples)
  
}


# -----------------------------------------------------------------------------
# Specific MCMC algorithms. 
# -----------------------------------------------------------------------------

# emulator_info: list with "gp_fits", "output_stats", "settings", and "input_bounds"
# TODO: allow joint sampling, incorporating covariance between current and proposal. "output_stats" must be 
# on the correct scale (e.g. it should be on log scale for LNP).
mcmc_calibrate_ind_GP <- function(computer_model_data, theta_prior_params, emulator_info,
                                  theta_init = NULL, sig2_eps_init = NULL, learn_sig_eps = FALSE, sig_eps_prior_params = NULL, 
                                  N_mcmc = 50000, adapt_frequency = 1000, adapt_min_scale = 0.1, accept_rate_target = 0.24, 
                                  proposal_scale_decay = 0.7, Cov_prop_init_diag = 0.1, adapt_cov_method = "AM", 
                                  adapt_scale_method = "MH_ratio", adapt_init_threshold = 3) {
  # `theta_prior_params` should already be truncated, if desired. 
  
  if(learn_sig_eps && sig_eps_prior_params$dist != "IG") {
    stop("`mcmc_calibrate_ind_GP()` requires inverse gamma priors on observation variances.")
  }
  
  # Number observations in time series, number output variables, and dimension of parameter space.
  p <- length(computer_model_data$output_vars)
  d <- length(computer_model_data$pars_cal_names)
  
  # Ensure parameters in `theta_prior_params` are sorted correctly based on ordering in `computer_model_data`. 
  theta_prior_params <- theta_prior_params[computer_model_data$pars_cal_names,]
  
  # Objects to store samples.
  theta_samp <- matrix(nrow = N_mcmc, ncol = d)
  colnames(theta_samp) <- computer_model_data$pars_cal_names
  sig2_eps_samp <- matrix(nrow = N_mcmc, ncol = p)
  colnames(sig2_eps_samp) <- computer_model_data$output_vars
  
  # Set initial conditions. 
  if(is.null(theta_init)) {
    theta_init <- sample_prior_theta(theta_prior_params)
  }
  
  if(learn_sig_eps) {
    if(is.null(sig2_eps_init)) {
      sig2_eps_init <- sample_prior_Sig_eps(sig_eps_prior_params)
    }
  } else {
    if(is.null(sig_eps_init)) stop("Value for `sig_eps_init` must be provided when `learn_sig_eps` is FALSE.")
  }
  
  theta_samp[1,] <- theta_init
  sig2_eps_samp[1,] <- sig2_eps_init
  
  theta_scaled_curr <- scale_input_data(theta_samp[1,,drop=FALSE], input_bounds=emulator_info$input_bounds)
  sig2_eps_curr <- sig2_eps_init
  lprior_theta_curr <- calc_lprior_theta(theta_init, theta_prior_params)
  SSR_curr <- sample_emulator_cond(theta_scaled_curr, emulator_info, cond_type="post", sig2_eps=sig2_eps_curr)
  
  # Proposal covariance.
  Cov_prop <- diag(Cov_prop_init_diag, nrow = d)
  L_prop <- t(chol(Cov_prop))
  log_scale_prop <- 0
  accept_count <- 0
  samp_mean <- theta_init
  
  for(itr in seq(2, N_mcmc)) {
    
    #
    # Metropolis step for theta.
    #
    
    # theta proposals.
    theta_prop <- theta_samp[itr-1,] + (exp(log_scale_prop) * L_prop %*% matrix(rnorm(d), ncol = 1))[,1]
    
    # Immediately reject if proposal is outside of prior bounds (i.e. prior density is 0). If this occurs on the first 
    # iteration we let the normal calculations proceed since we need to initialize `gp_pred_list`. In this case, the 
    # calculations will still return an acceptance probability of 0. After the first iteration, there is no need to 
    # waste computation if we know the acceptance probability will be 0. 
    if((itr > 2) && 
       (any(theta_prop < theta_prior_params[["bound_lower"]], na.rm = TRUE) ||
        any(theta_prop > theta_prior_params[["bound_upper"]], na.rm = TRUE))) {
      
      theta_samp[itr,] <- theta_samp[itr-1,]
      alpha <- 0
      
    } else {
      
      # Accept-Reject step. 
      lpost_theta_curr <- calc_lpost_theta_product_lik(computer_model_data, lprior_vals=lprior_theta_curr, SSR=SSR_curr, 
                                                       vars_obs=sig2_eps_curr, normalize_lik=FALSE, na.rm=TRUE, return_list=FALSE)
      lpost_theta_prop_list <- calc_lpost_theta_product_lik(computer_model_data, theta_vals=matrix(theta_prop, nrow=1), SSR=SSR_curr, 
                                                            vars_obs = sig2_eps_curr, normalize_lik=FALSE, na.rm=TRUE, theta_prior_params=theta_prior_params, 
                                                            return_list=TRUE)
      alpha <- min(1.0, exp(lpost_theta_prop_list$lpost - lpost_theta_curr))
      
      if(runif(1) <= alpha) {
        theta_samp[itr,] <- theta_prop
        theta_scaled_curr <- scale_input_data(matrix(theta_prop,nrow=1), input_bounds=emulator_info$input_bounds)
        lprior_theta_curr <- lpost_theta_prop_list$lprior
        accept_count <- accept_count + 1 
      } else {
        theta_samp[itr,] <- theta_samp[itr-1,]
      }
      
    }
    
    # Adapt proposal covariance matrix and scaling term.
    adapt_list <- adapt_cov_proposal(Cov_prop, log_scale_prop, theta_samp, itr, accept_count, alpha, samp_mean, adapt_frequency, 
                                     adapt_cov_method, adapt_scale_method, accept_rate_target, adapt_min_scale,  
                                     proposal_scale_decay, adapt_init_threshold, L = L_prop)
    Cov_prop <- adapt_list$C
    L_prop <- adapt_list$L
    log_scale_prop <- adapt_list$log_scale
    samp_mean <- adapt_list$samp_mean
    accept_count <- adapt_list$accept_count
    
    
    #
    # Gibbs step for phi. 
    #
    
    SSR_curr <- sample_emulator_cond(theta_scaled_curr, emulator_info, cond_type="post", sig2_eps=sig2_eps_curr)
    
    #
    # Gibbs step for sig2_eps. 
    #
    
    if(learn_sig_eps) {
      sig2_eps_curr <- sample_cond_post_Sig_eps(SSR=SSR_curr, Sig_eps_prior_params=sig_eps_prior_params, n_obs=computer_model_data$n_obs)
      sig2_eps_samp[itr,] <- sig2_eps_curr
    } else {
      sig2_eps_samp[itr,] <- sig2_eps_init
    }
    
  }
  
  return(list(theta = theta_samp, sig2_eps = sig2_eps_samp, Cov_prop = Cov_prop, scale_prop = exp(log_scale_prop)))
  
}



