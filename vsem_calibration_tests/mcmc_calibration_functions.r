#
# mcmc_calibration_functions.r
# Functions defining likelihoods and priors, and implementions of MCMC for parameter
# calibration of computer model. Used for tests with Very Simple Ecosystem (VSEM)
# model. 
#
# Andrew Roberts
#

llik_Gaussian <- function(par, Sig_eps, par_ref, par_cal_sel, output_vars, PAR, data_obs) {
  # Unnormalized Gaussian log-likelihood. Assumes independence across time, but allows 
  # for correlation across the p outputs by specifying the p x p covariance matrix 
  # Sig_eps. 
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
  
  # Evaluate sum (y_i - f_i)^T Sig_eps (y_i - f_i) over i.
  Sig_eps <- Sig_eps[output_vars, output_vars]
  L <- t(chol(Sig_eps))
  model_errs <- data_obs[, ..output_vars] - pred_model
  log_quadratic_form <- sum(colSums(forwardsolve(L, t(model_errs))^2))
  
  return(-0.5 * log_quadratic_form)
  
}


