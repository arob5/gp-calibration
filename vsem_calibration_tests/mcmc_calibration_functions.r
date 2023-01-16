#
# mcmc_calibration_functions.r
# Functions defining likelihoods and priors, and implementions of MCMC for parameter
# calibration of computer model. Used for tests with Very Simple Ecosystem (VSEM)
# model. 
#
# Andrew Roberts
#

library(LaplacesDemon)

llik_Gaussian <- function(par, Sig_eps, par_ref, par_cal_sel, output_vars, PAR, data_obs) {
  # Unnormalized Gaussian log-likelihood. Assumes independence across time, but allows 
  # for correlation across the p outputs by specifying the p x p covariance matrix 
  # Sig_eps. Runs the full forward model to obtain model outputs and calculate likelihood. 
  #
  # Args:
  #    par: numeric vector, calibration parameters at which to evaluate log-likelihood. 
  #    Sig_eps: matrix, dimension pxp, covariance between output variables (assumed 
  #             fixed across time). Rownames and colnames must be set to names of 
  #             the output variables. 
  #    par_ref: data.frame, rownames should correspond to parameters of computer model. 
  #             Must contain column named "best". Parameters that are fixed at nominal 
  #             values (not calibrated) are set to their values given in the "best" column. 
  #    par_cal_sel: integer vector, selects the rows of 'par_ref' that correspond to 
  #                 parameters that will be calibrated. 
  #    output_vars: character vector, used to the select the outputs to be considered in 
  #                 the likelihood; e.g. selects the correct sub-matrix of 'Sig_eps' and the 
  #                 correct columns of 'data_obs'. 
  #    PAR: numeric vector, time series of photosynthetically active radiation used as forcing
  #         term in VSEM. 
  #    data_obs: data.table, dimension n x p (n = length of time series, p = number outputs).
  #              Colnames set to output variable names. 
  #
  # Returns:
  #   Unnormalized log-likelihood across all observations and over the output variables
  #   specified in 'output_vars'. 
  
  # Parameters not calibrated are fixed at default values
  theta <- par_ref$best
  theta[par_cal_sel] <- par
  
  # Run forward model, re-scale NEE.
  pred_model <- as.data.table(VSEM(theta, PAR))[, ..output_vars]
  if("NEE" %in% output_vars) {
    pred_model[, NEE := NEE*1000]
  }
  
  # Evaluate sum (y_i - f_i)^T Sig_eps^{-1} (y_i - f_i) over i.
  Sig_eps <- Sig_eps[output_vars, output_vars]
  L <- t(chol(Sig_eps))
  model_errs <- data_obs[, ..output_vars] - pred_model
  log_quadratic_form <- sum(forwardsolve(L, t(model_errs))^2)
  
  return(-0.5 * log_quadratic_form)
  
}


llik_Gaussian_err <- function(model_errs, Sig_eps, output_vars) {
  # A version of llik_Gaussian() that is parameterized in terms of the n x p model 
  # error matrix Y - f(theta) and the observation covariance Sig_eps. This presumes
  # the forward model has already been run, unlike llik_Gaussian(),  which runs
  # the model in order to calculate the likelihood. 
  #
  # Args:
  #     model_errs: matrix of dimensions n x p, where n = number observations in time 
  #                 series and p = number output variables. 
  #     Sig_eps: matrix of dimensions p x p, the covariance matrix used in the 
  #              multivariate Gaussian likelihood calculation. 
  #
  # Returns:
  #    numeric, the unnormalized log-likelihood. 
  
  Sig_eps <- Sig_eps[output_vars, output_vars]
  L <- t(chol(Sig_eps))
  log_quadratic_form <- sum(forwardsolve(L, t(model_errs))^2)
  
  return(-0.5 * log_quadratic_form)
  
}


run_VSEM <- function(par, par_ref, par_cal_sel, PAR) {
  # Runs the VSEM model using specified parameter settings, returning the outputs 
  # of the model. 
  #
  # Args:
  #    par: numeric vector, values of calibration parameters used to run the model. 
  #    par_ref: data.frame, rownames should correspond to parameters of computer model. 
  #             Must contain column named "best". Parameters that are fixed at nominal 
  #             values (not calibrated) are set to their values given in the "best" column. 
  #    par_cal_sel: integer vector, selects the rows of 'par_ref' that correspond to 
  #                 parameters that will be calibrated. 
  #    PAR: numeric vector, time series of photosynthetically active radiation used as forcing
  #         term in VSEM. 
  #
  # Returns:
  #   matrix of dimensions n x p where n is the length of the time series and p is the number
  #   of output variables. 
  
  # Parameters not calibrated are fixed at default values
  theta <- par_ref$best
  theta[par_cal_sel] <- par
  
  # Run forward model, re-scale NEE.
  pred_model <- as.data.table(VSEM(theta, PAR))[, ..output_vars]
  if("NEE" %in% output_vars) {
    pred_model[, NEE := NEE*1000]
  }
  
  return(pred_model)
  
}


calc_lprior_theta <- function(theta, theta_prior_params) {
  # Evaluates the log prior density on calibration functions at specific values of the settings.  
  #
  # Args:
  #    theta: numeric vector, the value of the calibration parameters at which to evaluate the prior density. 
  #    theta_prior_params: data.frame, with columns "dist", "param1", and "param2". The ith row of the data.frame
  #                        should correspond to the ith entry of 'theta'. Currently, accepted values of "dist" are 
  #                        "Gaussian" (param1 = mean, param2 = std dev) and "Uniform" (param1 = lower, param2 = upper).
  #
  # Returns:
  #    The prior density evaluation log p(theta). Assumes prior independence, so the log-prior is the sum of the log-prior
  #    evaluations for each entry of 'theta'. 
  
  lprior <- 0
  
  theta_prior_params[["val"]] <- theta
  Gaussian_priors <- theta_prior_params[theta_prior_params$dist == "Gaussian",]
  Uniform_priors <- theta_prior_params[theta_prior_params$dist == "Uniform",]
  
  if(nrow(Gaussian_priors) > 0) {
    lprior <- lprior + sum(dnorm(Gaussian_priors$val, Gaussian_priors$param1, Gaussian_priors$param2, log = TRUE))
  }
  
  if(nrow(Uniform_priors) > 0) {
    lprior <- lprior + sum(dunif(Uniform_priors$val, Uniform_priors$param1, Uniform_priors$param2, log = TRUE))
  }
  
  return(lprior)  

}


# TODO: 
#    - Look into need for truncated Gaussian in PEcAn algorithm
#    - Pass in priors instead of hard-coding
mcmc_calibrate <- function(f, llik_func, lprior, par_cal_sel, par_ref, data_obs, output_vars, PAR,
                           theta_init, theta_prior_params, theta_bounds = NA, Sig_eps_init, Sig_eps_prior_params = NA,
                           N_mcmc, learn_Sig_eps) {
  # MCMC with no GP approximations. Can handle correlated outputs. 

  # Number observations in time series, number output variables, and dimension of parameter space
  n <- nrow(data_obs)
  p <- length(output_vars)
  d <- length(theta_sel)

  # Objects to store samples
  par_cal_names <- rownames(par_ref)[par_cal_sel]
  theta_samp <- matrix(nrow = N_mcmc, ncol = par_cal_sel)
  colnames(theta_samp) <- par_cal_names
  Sig_eps_samp <- matrix(nrow = N_mcmc, ncol = 0.5*p*(p+1)) # Each row stores lower triangle of Sig_eps

  # Set initial values
  theta_samp[1,] <- theta_init
  Sig_eps_samp[1,] <- Sig_eps_init
  Sig_eps_curr <- Sig_eps_init
  model_errs_curr <- data_obs - run_VSEM(theta_init, par_ref, par_cal_sel, PAR)
  log_theta_post_curr <- llik_Gaussian_err(model_errs_curr, Sig_eps_curr, output_vars) + calc_lprior_theta(theta_init, theta_prior_params)
  
  # Proposal covariance
  Cov_prop <- diag(1, nrow = d)
  L_prop <- t(chol(Cov_prop))

  for(itr in seq(2, N_mcmc)) {

    #
    # Metropolis step for theta
    #

    # theta proposal
    theta_prop <- theta_samp[itr-1,] + L_prop %*% rnorm(p)

    # Calculate log-likelihood for proposal
    model_errs_prop <- data_obs - run_VSEM(theta_prop, par_ref, par_cal_sel, PAR)
    llik_prop <- llik_Gaussian_err(model_errs_prop, Sig_eps_curr, output_vars)

    # Accept-Reject Step
    log_theta_post_prop <- llik_prop + calc_lprior_theta(theta_prop, theta_prior_params)
    if(runif(1) < exp(log_theta_post_prop - log_theta_post_curr)) {
      theta_samp[itr,] <- theta_prop
      log_theta_post_curr <- log_theta_post_prop
      model_errs_curr <- model_errs_prop
    } else {
      theta_samp[itr,] <- theta_samp[itr-1,]
    }

    #
    # Gibbs step for Sig_eps
    #

    if(learn_Sig_eps) {
      Sig_eps_curr <- sample_Sig_eps()
      Sig_eps_curr <- sample_Sig_eps(model_errs_curr, Sig_eps_prior_params, n)
      Sig_eps_samp[itr,] <- lower.tri(Sig_eps_curr, diag = TRUE)
    }

  }

}


sample_Sig_eps <- function(model_errs_curr, Sig_eps_prior_params, n) {
  
  inv_wishart_scale <- crossprod(model_errs_curr) + Sig_eps_prior_params$scale_matrix
  inv_wishart_dof <- n + Sig_eps_prior_params$dof
  
  
}














