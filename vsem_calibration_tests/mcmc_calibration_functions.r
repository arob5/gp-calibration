#
# mcmc_calibration_functions.r
# Functions defining likelihoods and priors, and implementions of MCMC for parameter
# calibration of computer model. Used for tests with Very Simple Ecosystem (VSEM)
# model. 
#
# Andrew Roberts
#

library(truncnorm)
library(LaplacesDemon)

llik_product_Gaussian <- function(theta_vals = NULL, par_ref = NULL, par_cal_sel = NULL, PAR_data = NULL, data_obs = NULL, output_vars = NULL, 
                                  SSR = NULL, vars_obs, n_obs, normalize = TRUE, na.rm = FALSE, sum_output_lliks = TRUE) {
  # Evaluate Gaussian log-likelihood assuming independence across time and across
  # output variables. In order to calculate the log-likelihood, the forward model must be run and then 
  # the L2 difference between the forward model outputs and the observed data evaluated. This function can 
  # perform all of those steps, or the L2 error can be passed directly using the `SSR` argument to avoid 
  # running the forward model when it is not necessary. For the forward model to be run, the arguments 
  # `theta_vals`, `par_ref`, `par_cal_sel`, `PAR_data`, `data_obs`, and `output_vars` must all be passed. Otherwise, these are 
  # not needed if `SSR` is provided. This function is vectorized in that it can evaluate 
  # the log-likelihood at multiple input parameter values, `theta_vals`. 
  #
  # Args:
  #    theta_vals: matrix, of dimension M x d, where M is the number of calibration parameter values
  #                and d is the number of claibration parameters. These are the parameter inputs at which 
  #                to run the forward to evaluate log-likelihood. 
  #    par_ref: data.frame, rownames should correspond to parameters of computer model. 
  #             Must contain column named "best". Parameters that are fixed at nominal 
  #             values (not calibrated) are set to their values given in the "best" column. 
  #    par_cal_sel: integer vector, selects the rows of 'par_ref' that correspond to 
  #                 parameters that will be calibrated. 
  #    PAR_data: numeric vector, time series of photosynthetically active radiation used as forcing
  #         term in VSEM.
  #    data_obs: data.table, dimension n x p (n = length of time series, p = number outputs).
  #              Colnames set to output variable names. 
  #    output_vars: character vector, used to the select the outputs to be considered in 
  #                 the likelihood.
  #    SSR: matrix, of dimension M x p. The (m, j) entry is the sum of squared errors for the 
  #         jth output variable and mth input parameter. 
  #    vars_obs: numeric(p), vector of observation/noise variances for each output. 
  #    n_obs: integer(p), the number of observations for each output.
  #    normalize: logical(1), whether or not to include the portion of the Gaussian normalization constant
  #               that doesn't depend on the theta or variance parameters. 
  #    na.rm: logical(1), if TRUE ignores missing data when calculating the sum of squared errors. 
  #    sum_output_lliks: logical(1), if TRUE returns matrix of dimension M x p containing the log-likelihood 
  #                      evaluations for each input parameter and each output separately. Otherwise, sums across 
  #                      the rows to calculate the overall likelihood, and thus returns a numeric vector of length M. 
  #
  # Returns:
  #    See `sum_output_lliks` in the Args section above. 
  
  # Run forward model. 
  if(is.null(SSR)) {
    model_outputs_list <- run_VSEM(theta_vals, par_ref, par_cal_sel, PAR_data, output_vars = output_vars) 
    SSR <- calc_SSR(data_obs, model_outputs_list, na.rm = na.rm)
  }
  
  p <- length(vars_obs)
  llik_outputs <- matrix(nrow = nrow(SSR), ncol = ncol(SSR))
  
  for(j in seq(1, p)) {
    llik_outputs[,j] <- -0.5 * n_obs[j] * log(vars_obs[j]) - 0.5 * SSR[,j] / vars_obs[j]
    if(normalize) llik_outputs[,j] <- llik_outputs[,j] - 0.5 * n_obs[j] * log(2*pi)
  }
  
  # Either return likelihood for each output seperately, or return single likelihood across all outputs. 
  if(sum_output_lliks) {
    return(rowSums(llik_outputs))
  }
  
  return(llik_outputs)
  
}


llik_Gaussian_SSR <- function(SSR, vars_obs, n_obs, normalize = TRUE) {
  # Evaluate Gaussian log-likelihood assuming independence across time and across
  # output variables. This is a convenience function that is parameterized in terms
  # of the sum of squared errors for each output variable. This function evaluates 
  # the log-likelihood for each output separately and returns the evaluations in 
  # a vector. `llik_product_Gaussian_SSR()` is essentially the same function 
  # but sums the log-likelihoods, thus returning the log-likelihood over 
  # all outputs. 
  #
  # Args:
  #    SSR: numeric(p), vector of sum of squared errors for each output, where p is 
  #         the number of output variables. 
  #    vars_obs: numeric(p), vector of observation/noise variances for each output. 
  #    n_obs: integer(p), the number of observations for each output.
  #
  # Returns:
  #    numeric(p), the log-likelihood evaluations for the p outputs. 
  
  p <- length(vars_obs)
  llik_outputs <- vector(mode = "numeric", length = p)
  
  for(j in seq(1, p)) {
    llik_outputs[j] <- -0.5 * n_obs[j] * log(vars_obs[j]) - 0.5 * SSR[j] / vars_obs[j]
    if(normalize) llik_outputs[j] <- llik_outputs[j] - 0.5 * n_obs[j] * log(2*pi)
  }
  
  return(llik_outputs)
  
}


llik_product_Gaussian_SSR <- function(SSR, vars_obs, n_obs, normalize = TRUE) {
  # Evaluate Gaussian log-likelihood assuming independence across time and across
  # output variables. This is a convenience function that is parameterized in terms
  # of the sum of squared errors for each output variable. 
  #
  # Args:
  #    SSR: numeric(p), vector of sum of squared errors for each output, where p is 
  #         the number of output variables. 
  #    vars_obs: numeric(p), vector of observation/noise variances for each output. 
  #    n_obs: integer(p), the number of observations for each output. 
  #
  # Returns:
  #    numeric(1), the log-likelihood. 
  
  p <- length(vars_obs)
  
  llik <- -0.5 * sum(n_obs * log(vars_obs)) - 0.5 * sum(SSR / vars_obs)
  if(normalize) {
    llik <- llik - 0.5*log(2*pi) * sum(n_obs)
  }
  
  return(llik)
  
}


calc_SSR <- function(data_obs, model_outputs_list, na.rm = TRUE) {
  # Computes the sum of squared residuals (SSR) between model runs and observed 
  # data on a per-output basis. Can handle multiple model runs (e.g. one per 
  # design point for emulation) or outputs from single model run (e.g. as required
  # during MCMC).
  #
  # Args:
  #    data_obs: data.table of dimension n x p, where n is the length of the time 
  #              series outputs from the model, and p is the number of output variables.
  #              This is Y in my typical notation. 
  #    model_outputs_list: either a matrix of dimension n x p corresponding to the 
  #                        model outputs f(theta) from a single model run. Or a list 
  #                        {f(theta_1), ..., f(theta_N)} of such matrices, collecting
  #                        the outputs from multiple model runs. 
  #    na.rm: logical(1), whether or not to remove NA values from the sum of squares 
  #           calculation, passed to the `colSums()` functions. Default is TRUE. 
  # 
  # Returns:
  #    matrix of dimension N x p, where N is the number of model runs and p is the 
  #    number of output variables. The (i, j) entry of the matrix is the SSR
  #    for the jth output of the ith model run, i.e. ||Y_j - f(j, theta_i)||^2.
  
  # If model only run at single set of calibration parameter values
  if(is.matrix(model_outputs_list)) {
    model_outputs_list <- list(model_outputs_list)
  }
  
  N_param_runs <- length(model_outputs_list)
  N_outputs <- ncol(model_outputs_list[[1]])
  SSR <- matrix(nrow = N_param_runs, ncol = N_outputs)
  
  for(i in seq_along(model_outputs_list)) {
    SSR[i,] <- colSums(data_obs - model_outputs_list[[i]], na.rm = na.rm)^2
  }
  colnames(SSR) <- colnames(data_obs)
  
  return(SSR)
  
}

# TODO: 
#    This function is intended to be the analog of `llik_product_Gaussian` but for the likelihood that models correlation 
#    across output errors. It is a bit tricky to account for missing data, need to debug this as it doesn't appear to be working. 
# llik_Gaussian_correlated_outputs <- function(theta_vals, Sig_eps, par_ref, par_cal_sel, output_vars, PAR_data, data_obs, n_obs, normalize = TRUE) {
#   # Unnormalized Gaussian log-likelihood. Assumes independence across time, but allows 
#   # for correlation across the p outputs by specifying the p x p covariance matrix 
#   # Sig_eps. Runs the full forward model to obtain model outputs and calculate likelihood. 
#   #
#   # Args:
#   #    theta_vals: matrix, of dimension M x d, where M is the number of calibration parameter values
#   #                and d is the number of claibration parameters. These are the parameter inputs at which 
#   #                to run the forward to evaluate log-likelihood. 
#   #    Sig_eps: matrix, dimension pxp, covariance between output variables (assumed 
#   #             fixed across time). Rownames and colnames must be set to names of 
#   #             the output variables. 
#   #    par_ref: data.frame, rownames should correspond to parameters of computer model. 
#   #             Must contain column named "best". Parameters that are fixed at nominal 
#   #             values (not calibrated) are set to their values given in the "best" column. 
#   #    par_cal_sel: integer vector, selects the rows of 'par_ref' that correspond to 
#   #                 parameters that will be calibrated. 
#   #    output_vars: character vector, used to the select the outputs to be considered in 
#   #                 the likelihood; e.g. selects the correct sub-matrix of 'Sig_eps' and the 
#   #                 correct columns of 'data_obs'. 
#   #    PAR: numeric vector, time series of photosynthetically active radiation used as forcing
#   #         term in VSEM. 
#   #    data_obs: data.table, dimension n x p (n = length of time series, p = number outputs).
#   #              Colnames set to output variable names. 
#   #
#   # Returns:
#   #   Unnormalized log-likelihood across all observations and over the output variables
#   #   specified in 'output_vars'. 
#   
#   n_max <- nrow(data_obs)
#   
#   # Run forward model.
#   VSEM_outputs <- run_VSEM(theta_vals, par_ref, par_cal_sel, PAR_data, output_vars = output_vars) 
#   
#   # Compute forward model error with respect to observed data (weighted L2 error). 
#   model_errs_list <- lapply(VSEM_outputs, function(model_outputs) data_obs[, output_vars] - model_outputs)  
#   
#   # Cholesky factor of covariance matrix over outputs
#   Sig_eps <- Sig_eps[output_vars, output_vars]
#   L <- t(chol(Sig_eps))
#   
#   # log det(Sigma) term. This part is the same for each theta. 
#   log_det_term <- sum(n_obs * diag(L))
#   
#   # L2 error term. Different for each theta. 
#   M <- nrow(theta_vals)
#   log_L2_term <- vector(mode = "numeric", length = M)
#   for(i in seq(1, n_max)) {
#     non_missing_data_sel <- !is.na(data_obs[i,])
#     L_i <- L[non_missing_data_sel, non_missing_data_sel]
#     
#     for(m in seq(1, M)) {
#       err_i <- model_errs_list[[m]][i, non_missing_data_sel]
#       log_L2_term[m] <- log_L2_term[m] + sum((forwardsolve(L_i, err_i))^2)
#     }
#   }
#   
#   # Compute log likelihood up to normalizing constant. 
#   llik <- -log_det_term - 0.5 * log_L2_term
#  
#   # Optionally normalize the density. Normalization constant also doesn't depend on theta. 
#   if(normalize) {
#     llik <- llik - 0.5 * log(2*pi) * sum(n_obs)
#   }
#   
#   # TEST
#   llik2 <- sapply(model_errs_list, function(model_errs) llik_Gaussian_err(model_errs, Sig_eps, output_vars = NA, normalize = TRUE))
#   
#   return(list(llik = llik, llik2 = llik2))
#   
# }


llik_Gaussian_err <- function(model_errs, Sig_eps, output_vars = NA, normalize = TRUE) {
  # A version of llik_Gaussian() that is parameterized in terms of the n x p model 
  # error matrix Y - f(theta) and the observation covariance Sig_eps. This presumes
  # the forward model has already been run, unlike llik_Gaussian(),  which runs
  # the model in order to calculate the likelihood. 
  #
  # Args:
  #     model_errs: matrix of dimensions n x p, where n = number observations in time 
  #                 series and p = number output variables. This is the model error matrix
  #                 Y - f(theta). 
  #     Sig_eps: matrix of dimensions p x p, the covariance matrix used in the 
  #              multivariate Gaussian likelihood calculation. 
  #     output_vars: character vector, containing names of output variables to be selected.
  #                  If provided, uses row and column names of `Sig_eps` to select sub-matrix
  #                  associated with the output variables, as well as column names of `model_errs`.
  #                  If NA, uses entire matrix. 
  #    normalize: logical, if TRUE returns the normalized log-density, which requires calculation of 
  #               the determinant term. Otherwise, excludes this term, thus returning an unnormalized
  #               log density (where Sig_eps is treated as a constant). Default is TRUE. 
  #
  # Returns:
  #    numeric, the unnormalized log-likelihood. 
  
  if(!is.na(output_vars)) {
    Sig_eps <- Sig_eps[output_vars, output_vars]
    model_errs <- model_errs[, output_vars]
  }
  L <- t(chol(Sig_eps))
  log_quadratic_form <- sum(forwardsolve(L, t(model_errs))^2)

  if(normalize) {
    return(-0.5 * log_quadratic_form - 0.5 * prod(dim(model_errs)) * log(2*pi) - nrow(model_errs) * sum(log(diag(L))))
  }
  
  return(-0.5 * log_quadratic_form)
  
}


run_VSEM <- function(theta_vals, par_ref, par_cal_sel, PAR_data, output_vars = c("NEE", "Cv", "Cs", "CR")) {
  # Runs the VSEM model using the specified parameter settings, thus returning the outputs of the forward model. 
  # If multiple input parameters are passed, then the model is run at each input and the results are returned in a 
  # list, one for each input. If a single input parameter is passed, this function reduces to `run_VSEM_single_input()`.
  #
  # Args:
  #    theta_vals: matrix of dimension M x d, the values of calibration parameters used to run the model. Each row is 
  #                a setting in the d-dimensional space of calibration parameters. 
  #    par_ref: data.frame, rownames should correspond to parameters of computer model. 
  #             Must contain column named "true_value". Parameters that are fixed at nominal 
  #             values (not calibrated) are set to their values given in the "true_value" column. 
  #    par_cal_sel: integer vector, selects the rows of 'par_ref' that correspond to 
  #                 parameters that will be calibrated. 
  #    PAR_data: numeric vector, time series of photosynthetically active radiation used as forcing
  #              term in VSEM. 
  #    output_vars: character vector, selects columns of output from VSEM model. Default returns
  #                 all four columns. 
  #
  # Returns:
  #    If run at multiple input points (i.e. `theta_vals` has M > 1 rows) then a list of M matrices as returned, where each 
  #    matrix is the return value of `run_VSEM_single_input()` at the respective input value. If M = 1, then a single matrix 
  #    is returned not in a list. 
  
  if(isTRUE(nrow(theta_vals) > 1)) {
    return(lapply(theta_vals, function(theta) run_VSEM_single_input(theta, par_ref, par_cal_sel, PAR_data, output_vars)))
  } else {
    return(run_VSEM_single_input(theta_vals, par_ref, par_cal_sel, PAR_data, output_vars))
  }
  
}


run_VSEM_single_input <- function(par_val, par_ref, par_cal_sel, PAR_data, output_vars = c("NEE", "Cv", "Cs", "CR")) {
  # Runs the VSEM model using specified parameter setting, returning the outputs of the model. 
  #
  # Args:
  #    theta: numeric vector, values of calibration parameters used to run the model. 
  #    par_ref: data.frame, rownames should correspond to parameters of computer model. 
  #             Must contain column named "true_value". Parameters that are fixed at nominal 
  #             values (not calibrated) are set to their values given in the "true_value" column. 
  #    par_cal_sel: integer vector, selects the rows of 'par_ref' that correspond to 
  #                 parameters that will be calibrated. 
  #    PAR_data: numeric vector, time series of photosynthetically active radiation used as forcing
  #         term in VSEM. 
  #    output_vars: character vector, selects columns of output from VSEM model. Default returns
  #                 all four columns. 
  #
  # Returns:
  #   matrix of dimensions n x p where n is the length of the time series and p is the number
  #   of output variables. 

  # Parameters not calibrated are fixed at default values
  theta <- par_ref$true_value
  theta[par_cal_sel] <- par_val
  
  # Run forward model, re-scale NEE.
  VSEM_output <- as.matrix(VSEM(theta, PAR_data))
  if("NEE" %in% output_vars) {
    VSEM_output[, "NEE"] <- VSEM_output[, "NEE"] * 1000
  }
  
  # Compute LAI, if included in output variables. LAI is simply LAR times the above-ground vegetation pool
  # at time t. 
  if("LAI" %in% output_vars) {
    VSEM_output <- cbind(VSEM_output, theta[rownames(par_ref) == "LAR"] * VSEM_output[, "Cv"])
    colnames(VSEM_output)[ncol(VSEM_output)] <- "LAI"
  }
  
  # Select only the output variables 
  VSEM_output <- VSEM_output[, output_vars]
  
  return(VSEM_output)
  
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


calc_lpost_theta_product_lik <- function(lprior_vals = NULL, llik_vals = NULL, theta_vals = NULL, par_ref = NULL, par_cal_sel = NULL, PAR_data = NULL, data_obs = NULL, 
                                         output_vars = NULL, SSR = NULL, vars_obs = NULL, n_obs = NULL, normalize_lik = TRUE, na.rm = FALSE,
                                         theta_prior_params = NULL, return_list = TRUE) {
  # Evaluates the exact log-posterior density, up to the normalizing constant (model evidence). 
  # This functions provides the option to calculate the log prior and log likelihood from scrath, at the given parameter values `theta_vals` 
  # and other necessary arguments for computing these quantities. Or if these quantities have already been calculated, they can be passed 
  # directly. 
  #
  # Args:
  #    See `llik_product_Gaussian()` for args required to calculate log likelihood. 
  #    See `calc_lprior_theta()` for args required to calculate log prior. 
  #    theta_vals: matrix of dimension M x d, where d is the dimension of the input space. The vector of calibration parameter values 
  #                at which to evaluate the posterior. Not required if prior and likelihood are provided as arguments.   
  #    lprior_vals: numeric(M), vector of log prior evaluations. If NULL, these will be computed. 
  #    llik_vals: numeric(M). vector of log likelihood evaluations. If NULL, these will be computed. 
  #    return_list: logical(1), if TRUE returns list which contains log prior, likelihood, and posterior evaluations. Otherwise just 
  #                 returns numeric vector of log posterior evaluations. 
  #
  # Returns:
  #    See `return_list` argument above. 
  
  # Prior
  if(is.null(lprior_vals)) {
    lprior_vals <- sapply(theta_vals, function(theta) calc_lprior_theta(theta, theta_prior_params))
  }
  
  # Likelihood
  if(is.null(llik_vals)) {
    llik_vals <- llik_product_Gaussian(theta_vals = theta_vals, par_ref = par_ref, par_cal_sel = par_cal_sel, PAR_data = PAR_data, data_obs = data_obs, 
                                       output_vars = output_vars, SSR = SSR, vars_obs = vars_obs, n_obs = n_obs, normalize = normalize_lik, 
                                       na.rm = na.rm, sum_output_lliks = TRUE)
  }
  
  if(!return_list) return(lprior_vals + llik_vals)
  
  return(list(lpost = lprior_vals + llik_vals, 
              lprior = lprior_vals, 
              llik = llik_vals))
  
}


# TODO: update this; I believe data_obs should now be treated as a matrix
mcmc_calibrate <- function(par_ref, par_cal_sel, data_obs, output_vars, PAR,
                           theta_init = NA, theta_prior_params, 
                           learn_Sig_eps = FALSE, Sig_eps_init = NA, Sig_eps_prior_params = list(), diag_cov = FALSE, 
                           N_mcmc, adapt_frequency, adapt_min_scale, accept_rate_target, proposal_scale_decay, proposal_scale_init) {
  # MCMC implementation for VSEM carbon model. Accommodates Gaussian likelihood, possibly with correlations 
  # between the different output variables, but assumes independence across time. Samples from posterior over both  
  # calibration parameters (theta) and observation covariance (Sig_eps), or just over theta if `learn_Sig_eps` is FALSE.
  # Allows for arbitrary prior over theta, but assumes Inverse Wishart prior on Sig_eps (if treated as random).
  # MCMC scheme is adaptive random walk Metropolis. This MCMC algorithm does not involve any model emulation/likelihood 
  # approximation. 
  #
  # Args:
  #    par_ref: data.frame, rownames should correspond to parameters of computer model. 
  #             Must contain column named "best". Parameters that are fixed at nominal 
  #             values (not calibrated) are set to their values given in the "best" column. 
  #    par_cal_sel: integer vector, selects the rows of 'par_ref' that correspond to 
  #                 parameters that will be calibrated. 
  #    data_obs: matrix, dimension n x p (n = length of time series, p = number outputs).
  #              Colnames set to output variable names. 
  #    output_vars: character vector, used to the select the outputs to be considered in 
  #                 the likelihood; e.g. selects the correct sub-matrix of 'Sig_eps' and the 
  #                 correct columns of 'data_obs'. 
  #    PAR: numeric vector, time series of photosynthetically active radiation used as forcing
  #         term in VSEM.
  #    theta_init: numeric vector of length p, the initial value of the calibration parameters to use in MCMC. If NA, samples
  #                the initial value from the prior. 
  #    theta_prior_params: data.frame, with columns "dist", "param1", and "param2". The ith row of the data.frame
  #                        should correspond to the ith entry of 'theta'. Currently, accepted values of "dist" are 
  #                        "Gaussian" (param1 = mean, param2 = std dev) and "Uniform" (param1 = lower, param2 = upper).
  #    learn_Sig_eps: logical, if TRUE treats observation covariance matrix as random and MCMC samples from joint 
  #                   posterior over Sig_eps and theta. Otherwise, fixes Sig_eps at value `Sig_eps_init`.
  #    Sig_eps_init: matrix, p x p covariance matrix capturing dependence between output variables. If `learn_Sig_eps` is 
  #                  TRUE then `Sig_eps_init` can either be set to the initial value used for MCMC, or set to NA in which 
  #                  case the initial value will be sampled from the prior. If `learn_Sig_eps` is FALSE, then a non-NA value 
  #                  is required, and treated as the fixed nominal value of Sig_eps_init. 
  #    Sig_eps_prior_params: list, defining prior on Sig_eps. See `sample_prior_Sig_eps()` for details. Only required if `learn_Sig_eps` is TRUE.
  #    diag_cov: logical, if TRUE constrains Sig_eps to be diagonal. If the prior distribution is specified to be product Inverse Gamma then this 
  #              is automatically set to TRUE. Default is FALSE. 
  #    N_mcmc: integer, the number of MCMC iterations. 
  #    adapt_frequency: integer, number of iterations in between each covariance adaptation. 
  #    adapt_min_scale: numeric scalar, used as a floor for the scaling factor in covariance adaptation, see `adapt_cov_proposal()`.
  #    accept_rate_target: numeric scalar, the desired acceptance rate, see `adapt_cov_proposal()`.
  #    proposal_scale_decay: Controls the exponential decay in the adjustment made to the scale of the proposal covariance as a function of the 
  #                          number of iterations. 
  #    proposal_scale_init: numeric, the proposal covariance is initialized to be diagonal, with `proposal_scale_init` along the diagonal. 
  #
  # Returns:
  #    list, with named elements "theta" and "Sig_eps". The former is a matrix of dimension N_mcmc x p with the MCMC samples of 
  #    theta stored in the rows. The latter is of dimension N_mcmc x p(p+1)/2, where each row stores the lower triangle of the 
  #    MCMC samples of Sig_eps, ordered column-wise, using lower.tri(Sig_eps, diag = TRUE). If `learn_Sig_eps` is FALSE, then 
  #    the first row stores the fixed value of Sig_eps, and the remaining rows are NA. 
  
  # Number observations in time series, number output variables, and dimension of parameter space
  n <- nrow(data_obs)
  p <- length(output_vars)
  d <- length(par_cal_sel)

  # Objects to store samples
  par_cal_names <- rownames(par_ref)[par_cal_sel]
  theta_samp <- matrix(nrow = N_mcmc, ncol = par_cal_sel)
  colnames(theta_samp) <- par_cal_names
  if(learn_Sig_eps && isTRUE(Sig_eps_prior_params$dist == "IG")) {
    diag_cov <- TRUE
  }
  if(diag_cov) {
    Sig_eps_samp <- matrix(nrow = N_mcmc, ncol = p)
  } else {
    Sig_eps_samp <- matrix(nrow = N_mcmc, ncol = 0.5*p*(p+1)) # Each row stores lower triangle of Sig_eps
    stop("Need to fix correlated Gaussian likelihood function to work with missing data.")
  }

  # Set initial values
  if(is.na(theta_init)) {
    theta_init <- sample_prior_theta(theta_prior_params)
  }
  if(learn_Sig_eps && is.na(Sig_eps_init)) {
    Sig_eps_init <- sample_prior_Sig_eps(Sig_eps_prior_params)
  } else {
    Sig_eps_init <- diag(1, nrow = p, ncol = p)
  }
  
  theta_samp[1,] <- theta_init
  if(diag_cov) {
    Sig_eps_samp[1,] <- diag(Sig_eps_init)
  } else {
    Sig_eps_samp[1,] <- lower.tri(Sig_eps_init, diag = TRUE)
  }
  
  Sig_eps_curr <- Sig_eps_init
  model_errs_curr <- data_obs[, output_vars] - run_VSEM(theta_init, par_ref, par_cal_sel, PAR, output_vars) # SSR here instead? 
  lprior_theta_curr <- calc_lprior_theta(theta_init, theta_prior_params)
  
  # Proposal covariance
  Cov_prop <- diag(proposal_scale_init, nrow = d)
  log_scale_prop <- 0
  L_prop <- t(chol(Cov_prop))
  accept_count <- 0
  
  for(itr in seq(2, N_mcmc)) {

    #
    # Metropolis step for theta
    #

    # Adapt proposal covariance matrix and scaling term
    if((itr > 3) && ((itr - 1) %% adapt_frequency) == 0) {
      # Cov_prop <- adapt_cov_proposal(Cov_prop, theta_samp[(itr - adapt_frequency):(itr - 1),,drop=FALSE],
      #                                adapt_min_scale, accept_count / adapt_frequency, accept_rate_target)
      if(accept_count == 0) {
        Cov_prop <- adapt_min_scale * Cov_prop
      } else {
        Cov_prop <- stats::cov(theta_samp[(itr - adapt_frequency):(itr - 1),,drop=FALSE])
      }
      L_prop <- t(chol(Cov_prop))
      accept_count <- 0
    }
    
    # theta proposals
    theta_prop <- (theta_samp[itr-1,] + sqrt(exp(log_scale_prop)) * L_prop %*% matrix(rnorm(d), ncol = 1))[1,]

    # Calculate log-likelihoods for current and proposed theta
    model_errs_prop <- data_obs[, output_vars] - run_VSEM(theta_prop, par_ref, par_cal_sel, PAR, output_vars)
    llik_curr <- llik_Gaussian_err(model_errs_curr, Sig_eps_curr) 
    llik_prop <- llik_Gaussian_err(model_errs_prop, Sig_eps_curr)

    # Accept-Reject Step
    lprior_theta_prop <- calc_lprior_theta(theta_prop, theta_prior_params)
    log_theta_post_curr <- llik_curr + lprior_theta_curr
    log_theta_post_prop <- llik_prop + lprior_theta_prop
    alpha <- min(1.0, exp(log_theta_post_prop - log_theta_post_curr))
    if(runif(1) <= alpha) {
      theta_samp[itr,] <- theta_prop
      lprior_theta_curr <- lprior_theta_prop
      model_errs_curr <- model_errs_prop
      accept_count <- accept_count + 1 
    } else {
      theta_samp[itr,] <- theta_samp[itr-1,]
    }

    # Adapt scaling term for proposal
    log_scale_prop <- log_scale_prop + (1 / itr^proposal_scale_decay) * (alpha - accept_rate_target)
    
    #
    # Gibbs step for Sig_eps
    #

    if(learn_Sig_eps) {
      Sig_eps_curr <- sample_cond_post_Sig_eps(model_errs_curr, Sig_eps_prior_params, n)
      if(diag_cov) {
        Sig_eps_samp[itr,] <- diag(Sig_eps_curr)
      } else {
        Sig_eps_samp[itr,] <- lower.tri(Sig_eps_curr, diag = TRUE)
      }
    }

  }
  
  return(list(theta = theta_samp, Sig_eps = Sig_eps_samp))

}


sample_prior_theta <- function(theta_prior_params) {
  # Return sample from prior distribution on the calibration parameters (theta). 
  #
  # Args:
  #    theta_prior_params: data.frame, with columns "dist", "param1", and "param2". The ith row of the data.frame
  #                        should correspond to the ith entry of 'theta'. Currently, accepted values of "dist" are 
  #                        "Gaussian" (param1 = mean, param2 = std dev) and "Uniform" (param1 = lower, param2 = upper).
  #
  # Returns:
  #    numeric vector of length equal to number of rows of `theta_prior_params`, the prior sample. 
  
  theta_samp <- vector(mode = "numeric", length = nrow(theta_prior_params))
  
  for(i in seq_along(theta_samp)) {
    if(theta_prior_params[i, "dist"] == "Gaussian") {
      theta_samp[i] <- rnorm(1, theta_prior_params[i, "param1"], theta_prior_params[i, "param2"])
    } else if(theta_prior_params[i, "dist"] == "Uniform") {
      theta_samp[i] <- runif(1, theta_prior_params[i, "param1"], theta_prior_params[i, "param2"])
    }
  }

  return(theta_samp)  
  
}


sample_prior_Sig_eps <- function(Sig_eps_prior_params) {
  # Returns sample from prior distribution on the p x p observation covariance matrix Sig_eps
  #
  # Args:
  #    Sig_eps_prior_params: list, must contain element named "dist" specifying the prior distribution. Current accepted values are 
  #                          "IW" (Inverse Wishart) or "IG" (independent Inverse Gamma priors). Depending on value of "dist", 
  #                          must also contain either either 1.) names "scale_matrix" and "dof" 
  #                          that are the arguments of the Inverse Wishart distribution, or 2.) names "IG_shape" and "IG_scale" which 
  #                          each correspond to p-length vectors storing the parameters for the independent Inverse Gamma priors on each
  #                          variance parameter. Note that dist "IG" constrains Sig_eps to be diagonal, while "IW" does not. 
  #
  # Returns:
  #    matrix, p x p positive definite matrix sampled from the prior p(Sig_eps).
  
  if(Sig_eps_prior_params$dist == "IW") {
    return(LaplacesDemon::rinvwishart(nu = Sig_eps_prior_params$dof, S = Sig_eps_prior_params$scale_matrix))
  } else if(Sig_eps_prior_params$dist == "IG") {
    p <- length(Sig_eps_prior_params$IG_shape)
    sig2_eps_vars <- vector(mode = "numeric", length = p)
    for(j in seq_len(p)) {
      sig2_eps_vars[j] <- LaplacesDemon::rinvgamma(1, shape = Sig_eps_prior_params$IG_shape[j], scale = Sig_eps_prior_params$IG_scale[j])
    }
    return(diag(sig2_eps_vars))
  }
  
}


sample_cond_post_Sig_eps <- function(model_errs, Sig_eps_prior_params, n) {
  # Return sample of the observation covariance matrix Sig_eps, drawn from the conditional  
  # posterior p(Sig_eps|theta, Y). Under the model assumptions, this conditional posterior 
  # has an inverse Wishart or Inverse Gamma product distribution. This function does not explicitly take theta as an 
  # argument; rather, it assumes the forward model f(theta) has already been run and 
  # the error matrix Y - f(theta) is provided by the argument `model_errs`.
  #
  # Args:
  #     model_errs: matrix of dimensions n x p, where n = number observations in time 
  #                 series and p = number output variables. This is the model error matrix
  #                 Y - f(theta). 
  #    Sig_eps_prior_params: list, must contain element named "dist" specifying the prior distribution. Current accepted values are 
  #                          "IW" (Inverse Wishart) or "IG" (independent Inverse Gamma priors). Depending on value of "dist", 
  #                          must also contain either either 1.) names "scale_matrix" and "dof" 
  #                          that are the arguments of the Inverse Wishart distribution, or 2.) names "IG_shape" and "IG_scale" which 
  #                          each correspond to p-length vectors storing the parameters for the independent Inverse Gamma priors on each
  #                          variance parameter. Note that dist "IG" constrains Sig_eps to be diagonal, while "IW" does not. 
  #    n: The number of observations in Y (i.e. the length of the time series). Should 
  #       equal the number of rows of `model_errs`.
  #
  # Returns:
  #    matrix, p x p positive definite matrix sampled from the prior distribution p(Sig_eps|theta, Y).
  
  if(Sig_eps_prior_params$dist == "IW") {
    inv_wishart_scale <- crossprod(model_errs_curr) + Sig_eps_prior_params$scale_matrix
    inv_wishart_dof <- n + Sig_eps_prior_params$dof
    return(LaplacesDemon::rinvwishart(nu = inv_wisharat_dof, S = inv_wishart_scale))
  } else if(Sig_eps_prior_params$dist == "IG") {
    n <- nrow(model_errs)
    p <- ncol(model_errs)
    output_l2_errs <- colSums(model_errs^2)
    sig2_eps_vars <- vector(mode = "numeric", length = p)
    for(j in seq_len(p)) {
      sig2_eps_vars[j] <- LaplacesDemon::rinvgamma(1, shape = 0.5*n + Sig_eps_prior_params$IG_shape[j], scale = 0.5*output_l2_errs[j] + Sig_eps_prior_params$IG_scale[j])
    }
    return(diag(sig2_eps_vars))
  }
  
}


adapt_cov_proposal <- function(cov_proposal, sample_history, min_scale, accept_rate, accept_rate_target) {
  # Returns an adapted covariance matrix to be used in adaptive MCMC scheme. Computes the new
  # covariance matrix by considering the sample correlation calculated from previous parameter
  # samples, which is scaled by a factor determined by the acceptance rate and target acceptance
  # rate. 
  #
  # Args:
  #    cov_proposal: matrix, p x p positive definite covariance matrix. 
  #    sample_history: matrix, l x p, where l is number of previous parameter samples used 
  #                    in sample correlation calculation. Each row of the matrix is a previous 
  #                    sample. 
  #    min_scale: numeric scalar, used as a floor for the scaling factor. 
  #    accept_rate: numeric scalar, the MCMC accept rate over the l-length parameter history. 
  #    accept_rate_target: numeric scalar, the desired acceptance rate. 
  #
  # Returns:
  #    matrix, p x p covariance matrix. 
  
  if(accept_rate == 0) {
    return(min_scale * cov_proposal)
  } else {
    cor_estimate <- stats::cor(sample_history)
    scale_factor <- max(accept_rate / accept_rate_target, min_scale)
    stdev <- apply(sample_history, 2, stats::sd)
    scale_mat <- scale_factor * diag(stdev, nrow = length(stdev))
    return(scale_mat %*% cor_estimate %*% scale_mat)
  }
}


LHS_train_test <- function(N_train, N_test, prior_params, joint = TRUE, extrapolate = TRUE, order_1d = FALSE, true_theta_samples = NULL) {
  # Generate a set training (design) points and test/validation points for evaluating a Gaussian 
  # process fit using Latin Hypercube Sampling (LHS). Can generate the train and test sets jointly, 
  # meaning they are sampled together in the LHS procedure. Or they can be sampled independently 
  # using two separate LHS schemes. Also allows the option to ensure that the test points are 
  # within the extent of the training points, if only the interpolation accuracy of the GP 
  # is of interest. In this case, the truncation means that the test points will not exactly
  # be distributed as specified in `prior_params`.
  #
  # Args:
  #    N_train: integer, the number of training points to sample. 
  #    N_test: integer, the number of test points to sample. 
  #    prior_params: data.frame containing the prior distribution information of the input 
  #                  parameters, with each row corresponding to a parameter. See `calc_lprior_theta()`
  #                  for the requirements of this data.frame. 
  #    joint: logical, if TRUE jointly samples train/test. Default is TRUE. 
  #    extrapolate: logical, if TRUE no truncation is performed for the test samples. Otherwise 
  #                 truncation is performed to guarantee the samples are within the extent of the 
  #                 training set. Default is TRUE. 
  #    order_1d: logical(1), only relevant if the dimension of the input space (i.e. the number of 
  #              calibration parameters) is one-dimensional. In this case, if `order_1d` is TRUE then 
  #              the design points will be sorted in increasing order. This is convenient when 
  #              plotting 1d GP plots. 
  #    true_theta_samples: matrix, containing samples from the true posterior over the calibration parameters. If this is provided 
  #                        (non-NULL) it overrides the other test data settings. The train data is unaffected. 
  #
  # Returns: 
  #    list, with named elements "X_train" and "X_test", storing the N_train x d and 
  #    N_test x d matrices of samples, respectively. 
  
  # The dimension of the input space
  d <- nrow(prior_params)
  
  # Generate train and test sets
  if(joint) {
    X_combined <- randomLHS(N_train + N_test, d)
    X_train <- X_combined[1:N_train,,drop=FALSE]
    X_test <- X_combined[-(1:N_train),,drop=FALSE]
  } else {
    X_train <- randomLHS(N_train, d)
    X_test <- randomLHS(N_test, d)
  }
  
  # Apply inverse CDS transform using prior distributions.
  for(j in seq_len(d)) {

    if(prior_params[j, "dist"] == "Uniform") {
      
      # Train
      X_train[,j] <- qunif(X_train[,j], prior_params[j, "param1"], prior_params[j, "param2"])
      
      # Test
      if(extrapolate) {
        lower_bound_test <- prior_params[j, "param1"]
        upper_bound_test <- prior_params[j, "param2"]
      } else {
        lower_bound_test <- min(X_train[,j])
        upper_bound_test <- max(X_train[,j])
      }
      X_test[,j] <- qunif(X_test[,j], lower_bound_test, upper_bound_test)
      
    } else if(prior_params[j, "dist"] == "Gaussian") {
      
      # Train
      X_train[,j] <- qnorm(X_train[,j], mean = prior_params[j, "param1"], sd = prior_params[j, "param2"])
      
      # Test
      if(extrapolate) {
        lower_bound_test <- -Inf
        upper_bound_test <- Inf
      } else {
        lower_bound_test <- min(X_train[,j])
        upper_bound_test <- max(X_train[,j])
      }
      X_test[,j] <- qtruncnorm(X_test[,j], a = lower_bound_test, b = upper_bound_test, mean = prior_params[j, "param1"], sd = prior_params[j, "param2"])
    }
    
  }
  
  if(!is.null(true_theta_samples)) {
    X_test <- true_theta_samples
  }
  
  if(order_1d && (d == 1)) {
    X_train <- X_train[order(X_train),,drop=FALSE]
    X_test <- X_test[order(X_test),,drop=FALSE]
  }
  
  return(list(X_train = X_train, X_test = X_test))
  
}


grid_train_test <- function(N_train, N_test, prior_params, extrapolate, tail_prob_excluded = 0.01, true_theta_samples = NULL) {
  # Note that currently the "extrapolate" argument makes no difference, since grid points are placed
  # at the edges of the defined bounds. This can be changed if needed, by adding an "extrapolate_prob" 
  # argument or something like that. 

  # The dimension of the input space.
  d <- nrow(prior_params)
  
  # Number of points marginally for each dimension.
  N_train_per_dim <- N_train^(1/d)
  N_test_per_dim <- N_test^(1/d)
  if((round(N_train_per_dim) != N_train_per_dim) || (round(N_test_per_dim) != N_test_per_dim)) {
    stop("N_train and N_test must have an integer d^th root; d = ", d)
  }
  
  # Grid of inputs
  X_train_marginals <- matrix(nrow = N_train_per_dim, ncol = d)
  X_test_marginals <- matrix(nrow = N_test_per_dim, ncol = d)
  
  for(j in seq_len(d)) {
    
    if(prior_params[j, "dist"] == "Uniform") {
      
      # Train
      X_train_marginals[,j] <- seq(prior_params[j, "param1"], prior_params[j, "param2"], length.out = N_train_per_dim)
      
      # Test
      if(extrapolate) {
        lower_bound_test <- prior_params[j, "param1"]
        upper_bound_test <- prior_params[j, "param2"]
      } else {
        lower_bound_test <- min(X_train_marginals[,j])
        upper_bound_test <- max(X_train_marginals[,j])
      }
      X_test_marginals[,j] <- seq(lower_bound_test, upper_bound_test, length.out = N_test_per_dim)
      
    } else if(prior_params[j, "dist"] == "Gaussian") {
      
      # Train
      left_bound <- qnorm(tail_prob_excluded/2, prior_params[j, "param1"], prior_params[j, "param2"])
      right_bound <- left_bound + 2*abs(left_bound)
      X_train_marginals[,j] <- seq(left_bound, right_bound, length.out = N_train_per_dim)
      
      # Test
      if(extrapolate) {
        lower_bound_test <- left_bound
        upper_bound_test <- right_bound
      } else {
        lower_bound_test <- min(X_train_marginals[,j])
        upper_bound_test <- max(X_train_marginals[,j])
      }
      X_test_marginals[,j] <- seq(lower_bound_test, upper_bound_test, length.out = N_test_per_dim)
    }
    
  }
  
  # Create grids
  X_train <- as.matrix(expand.grid(lapply(seq(1, d), function(j) X_train_marginals[,j])))
  
  if(is.null(true_theta_samples)) {
    X_test <- as.matrix(expand.grid(lapply(seq(1, d), function(j) X_test_marginals[,j])))
  } else {
    X_test <- true_theta_samples
  }
  colnames(X_train) <- rownames(theta_prior_params)
  colnames(X_test) <- rownames(theta_prior_params)
  
  return(list(X_train = X_train, X_test = X_test))
  
}


get_train_test_data <- function(N_train, N_test, prior_params, extrapolate, ref_pars, pars_cal_sel,
                                data_obs, PAR, output_vars, scale_X, normalize_Y, log_SSR, 
                                joint_LHS = FALSE, method = "LHS", order_1d = FALSE, true_theta_samples = NULL) {
  # Generates training (design) dataset and test dataset via Latin Hypercube sampling or 
  # via a grid-based approach. Note that the LHS method takes into account the prior, while 
  # the grid-based approach simply creates a grid of evenly spaced points within the bounds
  # defined by the prior. 
  # Input points are sampled from the space of calibration parameters, while the outputs 
  # are the squared L2 errors between the observed data and VSEM outputs, or the log 
  # of this quantity. Optionally pre-processes the data by scaling the inputs and 
  # normalizing the outputs. 
  #
  # Args:
  #    N_train: Number of design points. 
  #    N_test: Number of test points. 
  #    prior_params: data.frame containing the prior distribution information of the input 
  #                  parameters, with each row corresponding to a parameter. See `calc_lprior_theta()`
  #                  for the requirements of this data.frame. 
  #    extrapolate: logical, if TRUE no truncation is performed for the test samples. Otherwise 
  #                 truncation is performed to guarantee the samples are within the extent of the 
  #                 training set. Default is TRUE. 
  #    ref_pars: data.frame, rownames should correspond to parameters of computer model. 
  #             Must contain column named "best". Parameters that are fixed at nominal 
  #             values (not calibrated) are set to their values given in the "best" column. 
  #    par_cal_sel: integer vector, selects the rows of 'ref_pars' that correspond to 
  #                 parameters that will be calibrated. 
  #    data_obs: matrix, dimension n x p (n = length of time series, p = number outputs).
  #              Colnames set to output variable names. 
  #    output_vars: character vector, used to the select the outputs to be considered in 
  #                 the likelihood; e.g. selects the correct sub-matrix of 'Sig_eps' and the 
  #                 correct columns of 'data_obs'. 
  #    PAR: numeric vector, time series of photosynthetically active radiation used as forcing
  #         term in VSEM.
  #    log_SSR: logical(1), if TRUE includes log-transformed outputs as well, in addition to the
  #             outputs on the original scale. 
  #    joint_LHS: logical, only relevant if `method` is LHS. If TRUE jointly samples train/test in a single LHS 
  #               sample. Othwerwise, samples the train and test data separately. Default is FALSE. 
  #    method: character(1), character string indicating the sampling or deterministic method used to 
  #            generate the train/test points. Valid options are "LHS" (latin hypercube sampling, the default), 
  #            or "grid". 
  #    order_1d: logical(1), only relevant if the dimension of the input space (i.e. the number of 
  #              calibration parameters) is one-dimensional. In this case, if `order_1d` is TRUE then 
  #              the design points will be sorted in increasing order, with the training response 
  #              values ordered accordingly. This is convenient when plotting 1d GP plots. 
  #    true_theta_samples: matrix, containing samples from the true posterior over the calibration parameters. If this is provided 
  #                        (non-NULL) it overrides the other test data settings. The train data is unaffected. 
  #
  # Returns:
  #    list, with names "X_train", "X_test", "X_train_preprocessed", "X_test_preprocessed", 
  #    "Y_train", "Y_train_preprocessed", "Y_test", "input_bounds", and "output_stats".
  #    The latter two encode the transformations applied to the inputs and outputs, respectively 
  #    (see prep_GP_training_data()). 
  
  # Generate train and test input datasets
  if(method == "LHS") {
    X_train_test <- LHS_train_test(N_train, N_test, prior_params, joint_LHS, extrapolate, order_1d, true_theta_samples = true_theta_samples)
  } else if(method == "grid") {
    X_train_test <- grid_train_test(N_train, N_test, prior_params, extrapolate, order_1d, true_theta_samples = true_theta_samples)
  }
  X_train <- X_train_test$X_train
  X_test <- X_train_test$X_test
  
  # Run VSEM on train and test sets to obtain outputs
  model_outputs_list_train <- lapply(X_train, function(theta) run_VSEM(theta, ref_pars, pars_cal_sel, PAR, output_vars))
  model_outputs_list_test <- lapply(X_test, function(theta) run_VSEM(theta, ref_pars, pars_cal_sel, PAR, output_vars))
  Y_train <- calc_SSR(data_obs[, output_vars], model_outputs_list_train, na.rm = TRUE)
  Y_test <- calc_SSR(data_obs[, output_vars], model_outputs_list_test, na.rm = TRUE)

  # Pre-process training data: scale inputs and normalize outputs
  GP_data_preprocessed <- prep_GP_training_data(X_train, Y_train, scale_X, normalize_Y)
  X_train_preprocessed <- GP_data_preprocessed$X
  Y_train_preprocessed <- GP_data_preprocessed$Y
  
  # Pre-process testing input data
  X_test_preprocessed <- scale_input_data(X_test, GP_data_preprocessed$input_bounds)
  
  # Pre-process testing output data
  Y_test_preprocessed <- normalize_output_data(Y = Y_test, output_stats = GP_data_preprocessed$output_stats)
  
  output_list <- list(X_train = X_train, X_test = X_test, X_train_preprocessed = X_train_preprocessed, 
                      X_test_preprocessed = X_test_preprocessed, Y_train = Y_train, 
                      Y_train_preprocessed = Y_train_preprocessed, Y_test = Y_test,
                      Y_test_preprocessed = Y_test_preprocessed,
                      input_bounds = GP_data_preprocessed$input_bounds, 
                      output_stats = GP_data_preprocessed$output_stats)
  
  # Include log-transformed L2 error
  if(log_SSR) {
    output_list[["log_Y_train"]] <- log(output_list$Y_train)
    output_list[["log_Y_test"]] <- log(output_list$Y_test)
    output_list[c("log_Y_train_preprocessed", "log_output_stats")] <- prep_GP_training_data(Y = output_list$log_Y_train, normalize_Y = TRUE)[c("Y", "output_stats")]
    output_list[["log_Y_test_preprocessed"]] <- normalize_output_data(Y = output_list$log_Y_test, output_stats = output_list$log_output_stats)
  }
  
  return(output_list)
  
}


evaluate_GP_emulators <- function(emulator_settings, N_iter, N_design_points, N_test_points, 
                                  theta_prior_params, ref_pars, pars_cal_sel, data_obs, PAR, output_vars, 
                                  joint = TRUE, extrapolate = FALSE) {
  
  # Create list to store results
  nrow_test_results <- N_iter * nrow(emulator_settings)
  rmse_cols <- paste0("rmse_", output_vars)
  rmse_scaled_cols <- paste0("rmse_scaled_", output_vars)
  fit_time_cols <- paste0("fit_time_", output_vars)
  emulator_settings_cols <- c("gp_lib", "target", "kernel", "scale_X", "normalize_y")
  colnames_test_results <- c("iter", emulator_settings_cols, rmse_cols, rmse_scaled_cols, fit_time_cols)
  test_results <- data.frame(matrix(nrow = nrow_test_results, ncol = length(colnames_test_results)))
  colnames(test_results) <- colnames_test_results
  idx <- 1
  
  any_log_SSR <- any(emulator_settings[,"target"] == "log_SSR")
  
  # Loop over design/test datasets
  for(iter in seq_len(N_iter)) {
    # Create design/test sets
    train_test_data <- get_train_test_data(N_train = N_design_points, 
                                           N_test = N_test_points,
                                           prior_params = theta_prior_params, 
                                           joint = joint, 
                                           extrapolate = extrapolate, 
                                           ref_pars = ref_pars, 
                                           pars_cal_sel = pars_cal_sel, 
                                           data_obs = data_obs,
                                           PAR = PAR, 
                                           output_vars = output_vars, 
                                           scale_X = TRUE, 
                                           normalize_Y = TRUE, 
                                           log_SSR = any_log_SSR)
      
    # Loop over each model specification
    for(j in seq_len(nrow(emulator_settings))) {

      # Input train (design) and test points, either scaled or un-scaled.
      if(emulator_settings[j, "scale_X"]) {
        X_design <- train_test_data$X_train_preprocessed
        X_pred <- train_test_data$X_test_preprocessed
      } else {
        X_design <- train_test_data$X_train
        X_pred <- train_test_data$X_test
      }
      
      # Training output points, potentially normalized and/or log-transformed.  
      normalize_y <- emulator_settings[j, "normalize_y"]
      log_y <- emulator_settings[j, "target"] == "log_SSR"
      y_train_sel_string <- "Y_train"
      output_stats_sel_string <- "output_stats"
      if(normalize_y) {
        y_train_sel_string <- paste0(y_train_sel_string, "_preprocessed")
      }
      if(log_y) {
        y_train_sel_string <- paste0("log_", y_train_sel_string)
        output_stats_sel_string <- paste0("log_", output_stats_sel_string)
      }
      Y_design <- train_test_data[[y_train_sel_string]]
      
      # Fit GP and calculate error metrics
      gp_lib <- emulator_settings[j, "gp_lib"]
      gp_fit_list <- fit_independent_GPs(X_design, Y_design, gp_lib, emulator_settings[j, "kernel"])
      gp_pred_list <- predict_independent_GPs(X_pred, gp_fit_list$fits, gp_lib, cov_mat = FALSE,
                                              denormalize_predictions = normalize_y, 
                                              output_stats = train_test_data[[output_stats_sel_string]], 
                                              exponentiate_predictions = log_y)
      gp_err_list <- calc_independent_gp_pred_errs(gp_pred_list, train_test_data$Y_test)
      
      # Populate results data.frame
      test_results[idx, "iter"] <- iter
      test_results[idx, emulator_settings_cols] <- emulator_settings[j, emulator_settings_cols]
      test_results[idx, rmse_cols] <- sapply(gp_err_list, function(x) x$rmse)
      test_results[idx, rmse_scaled_cols] <- sapply(gp_err_list, function(x) x$rmse_scaled)
      test_results[idx, fit_time_cols] <- gp_fit_list$times
      idx <- idx + 1
    }

    
  }
  
  return(test_results)

}


calc_independent_gp_pred_errs <- function(gp_pred_list, Y_true) {
  
  lapply(seq_along(gp_pred_list), function(j) calc_gp_pred_err(gp_pred_list[[j]]$mean, gp_pred_list[[j]]$sd2, Y_true[,j]))
  
}


# TODO: look into variance vs. nugget variance returned in pred object 
calc_gp_pred_err <- function(gp_pred_mean, gp_pred_var, y_true) {
  
  sq_diff <- (y_true - gp_pred_mean)^2
  N_obs <- length(y_true)
  
  return(list(rmse = sum(sq_diff) / N_obs, 
              rmse_scaled = sum(sq_diff / gp_pred_var) / N_obs))
  
}


# plot_lnp_fit_1d <- function(X_test, y_test, X_train, y_train, lnp_log_mean_pred, lnp_log_var_pred, 
#                            vertical_line = NULL, xlab = "", ylab = "", main_title = "", CI_prob = 0.9) {
# 
#   order_pred <- order(X_test)
#   order_train <- order(X_train)
#   lnp_log_sd_pred <- sqrt(lnp_log_var_pred)
#   
#   # Confidence Intervals
#   CI_tail_prob <- 0.5 * (1 - CI_prob)
#   CI_plot_label <- paste0(CI_prob * 100, "% CI")
#   CI_lower <- qlnorm(CI_tail_prob, lnp_log_mean_pred, lnp_log_sd_pred, lower.tail = TRUE)
#   CI_upper <- qlnorm(CI_tail_prob, lnp_log_mean_pred, lnp_log_sd_pred, lower.tail = FALSE)
#   
#   # ggplot 
#   df <- data.frame(x_test = X_test[order_pred,1], 
#                    y_test = y_test[order_pred], 
#                    y_test_pred = lnp_log_mean_pred[order_pred],
#                    CI_lower = CI_lower[order_pred], 
#                    CI_upper = CI_upper[order_pred], 
#                    x_train = X_train[order_train,1], 
#                    y_train = y_train[order_train])
#   
#   lnp_plot <- ggplot(data=df, aes(x = x_test, y = y_test_pred)) + 
#                 geom_line(color = "blue") + 
#                 geom_line(aes(y = y_test), color = "red") + 
#                 geom_line(aes(y = CI_lower), color = "gray") + 
#                 geom_line(aes(y = CI_upper), color = "gray") + 
#                 geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), fill = "gray", alpha = 0.5) + 
#                 geom_point(aes(x = x_train, y = y_train), color = "red") + 
#                 xlab(xlab) + 
#                 ylab(paste0(ylab, " (", CI_plot_label, ")")) + 
#                 ggtitle(main_title)
#   
#   # Add vertical line
#   if(!is.null(vertical_line)) {
#     gp_plot <- gp_plot + geom_vline(xintercept = vertical_line, linetype = 2, color = "pink1")
#   }
#   
# }  


generate_vsem_test_data <- function(random_seed, N_time_steps, Sig_eps, pars_cal_names, pars_cal_vals, 
                                    ref_pars, output_vars, output_frequencies, obs_start_day) {
  # A helper function to generate data associated with a specific VSEM test example. Returns a list of all of the information
  # necessary to perform a VSEM emulation test. 
  #
  # Args:
  #    random_seed: integer, random seed. `generate_vsem_test_data()` should be the first function called that utilizes the random 
  #                 seed after the random seed is set in the main script. 
  #    N_time_steps: integer, number of time steps to integrate the VSEM ODE system. 
  #    Sig_eps: matrix of shape (p, p), where p is the number of output variables. Must have row and column names set to the 
  #             names of the associated output variables. 
  #    pars_cal_names: character(), vector of names of the calibration parameters. 
  #    pars_cal_vals: numeric(), vector of equal length as `pars_cal_names` with the true values of the calibration parameters
  #                   to use when generating the synthetic ground truth data. 
  #    ref_pars: data.frame, containing column "best" and with rownames set to the names of the VSEM parameters. The parameters 
  #              that are not specified in `pars_cal_names` will be set to their values in the "best" column when 
  #              generating the synthetic ground truth data.
  #    output_vars: character(), vector of names of the output variables to consider. Options to include are "NEE", "Cv", "Cs", "CR".
  #    output_frequencies: integer(), vector of length equal to length of `output_vars`. The ith element of this vector is the 
  #                        frequency with which observations will be generated for the respective output 
  #                        (1 = daily, 7 = weekly, 30 = monthly, etc.). Daily is the smallest allowed frequency. 
  #    obs_start_day: integer(), vector of length equal to length of `output_vars`. The ith element of this vector is the time step 
  #                   at which the observations will begin. For example, 1 indicates the observations will begin immediately, while 
  #                   10 will start observations at the 10th time step. The observation frequency will be unaffected by the start 
  #                   day. For example, a start day of 213 and observation frequency of 365 will imply annual observations 
  #                   on August 1st (assuming non leap-year). 
  #
  # Returns:
  #    List containing simulated ground truth data, simulated observed data, data summarizing the parameter values used, etc. 
  
  N_outputs <- length(output_vars)
  N_pars <- length(pars_cal_names)
  
  # Generate time series of photosynthetically active radiation (PAR), which drives the model. 
  PAR <- VSEMcreatePAR(seq_len(N_time_steps))
  
  # Specifying which parameters to calibrate. The remaining are fixed at their best values specified in 
  # 'ref_pars'. The values of the calibration parameters used to generate the ground truth data 
  # (the "true parameters") are given by `pars_cal_vals`, which may or may not be the same as the 
  # value of the calibration parameters in the "best" column of `ref_pars`. 
  pars_cal_sel <- sapply(pars_cal_names, function(par_name) which(rownames(ref_pars) == par_name))
  ref_pars[["calibrate"]] <- rownames(ref_pars) %in% pars_cal_names
  ref_pars[pars_cal_sel, "true_value"] <- pars_cal_vals
  ref_pars[ref_pars$calibrate == FALSE, "true_value"] <- ref_pars[ref_pars$calibrate == FALSE, "best"]
  
  # Run the model to generate the reference data, the ground truth. 
  data_ref <- run_VSEM(pars_cal_vals, ref_pars, pars_cal_sel, PAR, output_vars)

  # Add observational noise. Assumes Gaussian noise, potentially correlated across output variables. 
  Lt <- chol(Sig_eps[output_vars, output_vars]) # Upper triangular Cholesky factor of output covariance
  Z <- matrix(rnorm(N_time_steps*N_outputs), N_time_steps, N_outputs) 
  data_obs_complete <- data_ref + Z %*% Lt
  
  # Account for observation frequency
  data_obs <- data_obs_complete
  observation_selector <- matrix(nrow = N_time_steps, ncol = N_outputs)
  for(j in seq(1, N_outputs)) {
    obs_idx <- seq(obs_start_day[j], N_time_steps, by = output_frequencies[j])
    obs_sel <- rep(0, N_time_steps)
    obs_sel[obs_idx] <- 1
    observation_selector[,j] <- obs_sel
    data_obs[!obs_sel, j] <- NA
  }
  colnames(observation_selector) <- output_vars
  names(output_frequencies) <- output_vars
  names(obs_start_day) <- output_vars
  
  return(list(ref_pars = ref_pars, 
              PAR_data = PAR, 
              data_ref = data_ref, 
              pars_cal_sel = pars_cal_sel, 
              data_obs_complete = data_obs_complete, # Includes all daily data
              data_obs = data_obs, # Might have lower frequency data and/or missing values
              obs_selector = observation_selector,
              n_obs = colSums(observation_selector),
              Sig_eps = Sig_eps, 
              random_seed = random_seed, 
              N_time_steps = N_time_steps, 
              pars_cal_names = pars_cal_names, 
              output_vars = output_vars, 
              output_frequencies = output_frequencies, 
              obs_start_day = obs_start_day))
  
}


generate_vsem_test_case <- function(test_case_number) {
  # A convenience function that loads the VSEM test case given the number of the test case. 
  # First sets the random seed to the test case number and then loads the test case data.
  #
  # Args:
  #    integer(1), the test case number. 
  #
  # Returns:
  #    list, containing VSEM data for the test case. See `generate_vsem_test()` for details on 
  #    this list. 
  
  random_seed <- test_case_number
  f <- get(paste0("generate_vsem_test_", test_case_number))
  
  return(f(random_seed))
}


generate_vsem_test_1 <- function(random_seed) {
  # A convenience function to generate the VSEM data for "test case 1". This is the simplest 
  # test case, with a single calibration parameter, all outputs observed daily with no missing 
  # values, and no output correlation. 
  #
  # Args:
  #    random_seed: integer(1), random seed used to generate observation noise. 
  #
  # Returns:
  #    list, the list returned by `generate_vsem_test_data()`. 
  
  # Number of days.
  N_time_steps <- 1000
  
  # Diagonal covariance across outputs.
  Sig_eps <- diag(c(4.0, 1.0, 4.0, 1.0))
  rownames(Sig_eps) <- c("NEE", "Cv", "Cs", "CR")
  colnames(Sig_eps) <- c("NEE", "Cv", "Cs", "CR")
  
  # Single calibration parameter; extinction coefficient in Beer-Lambert law.
  pars_cal_names <- c("KEXT")
  pars_cal_vals <- 0.5
  ref_pars <- VSEMgetDefaults()
  
  # All outputs are observed daily, with no missing values. 
  output_vars <- c("NEE", "Cv", "Cs", "CR")
  output_frequencies <- rep(1, 4)
  obs_start_day <- rep(1, 4)
  
  test_list <- generate_vsem_test_data(random_seed, N_time_steps, Sig_eps, pars_cal_names, pars_cal_vals,
                                       ref_pars, output_vars, output_frequencies, obs_start_day)
  return(test_list)
  
}


generate_vsem_test_2 <- function(random_seed) {
  # A convenience function to generate the VSEM data for "test case 2". This is another 
  # single calibration parameter test case, but adds complexity by varying the frequency 
  # of output observations. It also varies the observations variances and increases the 
  # number of time steps (days). 
  #
  # Args:
  #    random_seed: integer(1), random seed used to generate observation noise. 
  #
  # Returns:
  #    list, the list returned by `generate_vsem_test_data()`. 
  
  # Number of days.
  N_time_steps <- 2048
  
  # Diagonal covariance across outputs.
  Sig_eps <- diag(c(4.0, 1.0, 2.0, 2.0))
  rownames(Sig_eps) <- c("NEE", "Cv", "Cs", "CR")
  colnames(Sig_eps) <- c("NEE", "Cv", "Cs", "CR")
  
  # Single calibration parameters; 1.) extinction coefficient in Beer-Lambert law 
  # and 2.) residence time of above-ground vegetation. 
  pars_cal_names <- c("LUE")
  pars_cal_vals <- c(0.002)
  ref_pars <- VSEMgetDefaults()
  
  # NEE is observed daily with no missing values. Soil and roots are observed annually on 
  # August 1st (assuming non leap-year). Above-ground vegetation is measured 6 times a year 
  # starting on March 1st (assuming non leap-year). 
  output_vars <- c("NEE", "Cv", "Cs", "CR")
  output_frequencies <- c(1, 60, 365, 365)
  obs_start_day <- c(1, 60, 213, 213)
  
  test_list <- generate_vsem_test_data(random_seed, N_time_steps, Sig_eps, pars_cal_names, pars_cal_vals,
                                       ref_pars, output_vars, output_frequencies, obs_start_day)
  return(test_list)
  
}


generate_vsem_test_3 <- function(random_seed) {
  # A convenience function to generate the VSEM data for "test case 3". This is another 
  # single calibration parameter test case, but adds complexity by varying the frequency 
  # of output observations. It also varies the observations variances and increases the 
  # number of time steps (days). 
  #
  # Args:
  #    random_seed: integer(1), random seed used to generate observation noise. 
  #
  # Returns:
  #    list, the list returned by `generate_vsem_test_data()`. 
  
  # Number of days.
  N_time_steps <- 3650
  
  # Diagonal covariance across outputs.
  Sig_eps <- diag(c(4.0, 0.36))
  rownames(Sig_eps) <- c("NEE", "LAI")
  colnames(Sig_eps) <- c("NEE", "LAI")
  
  # Single calibration parameter: Leaf Area Ratio (LAR)
  pars_cal_names <- c("LAR")
  pars_cal_vals <- c(1.500)
  ref_pars <- VSEMgetDefaults()
  
  # NEE is observed daily with no missing values. LAI is observed every three days.   
  output_vars <- c("NEE", "LAI")
  output_frequencies <- c(1, 3)
  obs_start_day <- c(1, 1)
  
  test_list <- generate_vsem_test_data(random_seed, N_time_steps, Sig_eps, pars_cal_names, pars_cal_vals,
                                       ref_pars, output_vars, output_frequencies, obs_start_day)
  return(test_list)
  
}




generate_vsem_test_4 <- function(random_seed) {
  # A convenience function to generate the VSEM data for "test case 4". This builds adds complexity 
  # to test case 1 by adding an additional calibration parameter, and varying the observation 
  # frequencies. 
  #
  # Args:
  #    random_seed: integer(1), random seed used to generate observation noise. 
  #
  # Returns:
  #    list, the list returned by `generate_vsem_test_data()`. 
  
  # Number of days.
  N_time_steps <- 2048
  
  # Diagonal covariance across outputs.
  Sig_eps <- diag(c(4.0, 2.0, 2.0, 2.0))
  rownames(Sig_eps) <- c("NEE", "Cv", "Cs", "CR")
  colnames(Sig_eps) <- c("NEE", "Cv", "Cs", "CR")
  
  # Single calibration parameters; 1.) extinction coefficient in Beer-Lambert law 
  # and 2.) residence time of above-ground vegetation. 
  pars_cal_names <- c("KEXT", "tauV")
  pars_cal_vals <- c(0.5, 1440)
  ref_pars <- VSEMgetDefaults()
  
  # All outputs are observed daily, with no missing values. 
  output_vars <- c("NEE", "Cv", "Cs", "CR")
  output_frequencies <- c(1, 300, 365, 30)
  obs_start_day <- c(1, 1, 1, 1)
  
  test_list <- generate_vsem_test_data(random_seed, N_time_steps, Sig_eps, pars_cal_names, pars_cal_vals,
                                       ref_pars, output_vars, output_frequencies, obs_start_day)
  return(test_list)
  
}


generate_vsem_test_5 <- function(random_seed) {
  # A convenience function to generate the VSEM data for "test case 5". This is another two calibration 
  # parameter test, which varies observation frequency, and now also includes time-independent correlation 
  # between observation noise in the output variables. 
  #
  # Args:
  #    random_seed: integer(1), random seed used to generate observation noise. 
  #
  # Returns:
  #    list, the list returned by `generate_vsem_test_data()`. 
  
  # Number of days.
  N_time_steps <- 3065
  
  # Diagonal covariance across outputs.
  Sig_eps <- diag(c(5.0, 2.0, 2.0, 2.0))
  rownames(Sig_eps) <- c("NEE", "Cv", "Cs", "CR")
  colnames(Sig_eps) <- c("NEE", "Cv", "Cs", "CR")
  Sig_eps["Cs", "CR"] <- -0.2
  Sig_eps["CR", "Cs"] <- -0.2
  
  # Single calibration parameters; 1.) extinction coefficient in Beer-Lambert law 
  # and 2.) residence time of above-ground vegetation. 
  pars_cal_names <- c("Cv", "tauR")
  pars_cal_vals <- c(3.0, 1440)
  ref_pars <- VSEMgetDefaults()
  
  # All outputs are observed daily, with no missing values. 
  output_vars <- c("NEE", "Cv", "Cs", "CR")
  output_frequencies <- c(1, 300, 365, 30)
  obs_start_day <- c(1, 1, 1, 1)
  
  test_list <- generate_vsem_test_data(random_seed, N_time_steps, Sig_eps, pars_cal_names, pars_cal_vals,
                                       ref_pars, output_vars, output_frequencies, obs_start_day <- c(1, 1, 1, 1))
  return(test_list)
  
}


GP_pointwise_loss <- function(val, gp_mean, gp_var) {
  # This is the negative log-predictive density of a GP evaluated at a single input point, 
  # up to an additive constant. Can evaluate the loss at a vector of values, in which 
  # case the three arguments `val`, `gp_mean`, and `gp_var` should all be of the same 
  # length. 
  #
  # Args:
  #    val: numeric, value at which to evaluate the loss/predictive density. 
  #    gp_mean: numeric, the GP predictive mean. 
  #    gp_var: numeric, the GP predictive variance. 
  #
  # Returns the negative log-predictive density, plus 0.5 * log(2*pi) to remove the constant. 
  
  log(gp_var) + 0.5 * (val - gp_mean)^2 / gp_var
  
}


get_gp_approx_posterior_LNP_params<- function(sig2_outputs, lprior_theta, gp_mean_pred, gp_var_pred, n_obs) {
  # This function assumes a product Gaussian likelihood, where independent GPs have been used to emulate
  # the sum of squared errors for each output. Under these assumptions, the approximation to the true 
  # posterior resulting from the GP approximations of the squared L2 error has a log-normal distribution. 
  # This function calculates the log mean and log variance of this log-normal process evaluated at some 
  # test inputs theta_1, ..., theta_M. 
  #
  # Args:
  #    sig2_outputs: numeric(p), vector of the observation variances for each output, where p is the number of 
  #                  outputs. 
  #    lprior_theta: numeric(M), vector of log prior evaluations at test points theta_1, ..., theta_M. 
  #    gp_mean_pred: matrix, of dimension M x p. The ith row contains the GP posterior mean evaluations 
  #                  mu^*_1(theta_i), ..., mu^*_p(theta_i) at test input theta_i. 
  #    gp_var_pred: matrix, of dimension M x p. The ith row contains the GP posterior variance evaluations 
  #                 k^*_1(theta_i), ..., k^*_p(theta_i) at test input theta_i. 
  #    n_obs: integer(p), the number of observations for each output. 
  #
  # Returns:
  #    list, with elements "mean_log" and "var_log". Each are numeric vectors of length M containing 
  #    the mean and variance, respectively, of the log of the GP posterior approximation at the 
  #    test inputs theta_1, ..., theta_M. 
  
  p <- length(sig2_outputs)
  log_C <- -0.5 * log(2*pi) * sum(n_obs) - 0.5 * sum(n_obs * log(sig2_outputs)) 
  
  gp_scaled_means <- gp_mean_pred %*% diag(1/sig2_outputs)
  gp_scaled_vars <- gp_var_pred %*% diag(1/sig2_outputs^2)
  
  mean_log_lnorm <- log_C + lprior_theta - 0.5 * rowSums(gp_scaled_means)
  var_log_lnorm <- 0.25 * rowSums(gp_scaled_vars)
  
  return(list(mean_log = mean_log_lnorm, var_log = var_log_lnorm))
  
}


gp_approx_posterior_pred_log_density <- function(log_vals, sig2_outputs = NULL, lprior_theta = NULL, gp_mean_pred = NULL, gp_var_pred = NULL, n_obs = NULL, 
                                                 lnorm_log_mean = NULL, lnorm_log_var = NULL) {
  # This function assumes a product Gaussian likelihood, where independent GPs have been used to emulate
  # the sum of squared errors for each output. The density evaluated is the predictive density of the posterior 
  # approximation ("pi star"). Under the stated assumptions, this is a log-normal density. Either the two parameters 
  # of this density may be passed (log mean and log var), or this function can calculate these paramneters first, then 
  # evaluate the density. 
  # The first argument is the log of the values at which to evaluate the density. 
  # Exponentiating these values will typically result in numerical underflow, so for numerical stability the data is 
  # shifted so that the density is always evaluated at the point 1.0, with the log-normal distribution appropriately 
  # scaled to account for this. Important note: this density is defined for the normalized 
  #
  # Args:
  #    log_vals: numeric(M), vector of (the log of) the points at which to evaluate the predictive density.
  #              This vector typically looks like [log pi(theta1), ..., log pi(thetaM)], where theta1, ..., thetaM
  #              are the test/validation inputs. This function takes the log of the evaluation points instead of the 
  #              values itself so that it can appropriately shift the data before exponentiating to avoid 
  #              numerical overflow/underflow. 
  #    sig2_outputs: numeric(p), vector of the observation variances for each output, where p is the number of 
  #                  outputs. 
  #    lprior_theta: numeric(M), vector of log prior evaluations at test points theta_1, ..., thetaM. 
  #    gp_mean_pred: matrix, of dimension M x p. The ith row contains the GP posterior mean evaluations 
  #                  mu^*_1(theta_i), ..., mu^*_p(theta_i) at test input theta_i. 
  #    gp_var_pred: matrix, of dimension M x p. The ith row contains the GP posterior variance evaluations 
  #                  k^*_1(theta_i), ..., k^*_p(theta_i) at test input theta_i. 
  #    n_obs: integer(p), the number of observations for each output.
  #    lnorm_log_mean: numeric(M), the first parameter (the mean of the log) for the lognormal distribution evaluated at the M 
  #                    points. If NULL, the log mean and variance will be computed. Default is NULL. 
  #    lnorm_log_var: numeric(M), the first parameter (the variance of the log) for the lognormal distribution evaluated at the M 
  #                   points. If NULL, the log mean and variance will be computed. Default is NULL.
  #
  # Returns:
  #    numeric(M), the log predictive density evaluations log p(pi(theta_i)|pi_hat(theta_i)) for the M input points
  #    theta_i. Note that the randomness is only coming from the p independent GPs used to approximate the posterior. 
  #    The GP posterior Gaussian distributions on the L2 error of each output induces a log-normal distribution on 
  #    the approximate posterior. 
  
  # Avoiding overflow/underflow 
  scaling_factors <- -log_vals
  
  # Calculate log mean and log variance, if not passed as arguments. 
  if(is.null(lnorm_log_mean) || is.null(lnorm_log_var)) {
    lnorm_params <- get_gp_approx_posterior_LNP_params(sig2_outputs, lprior_theta, gp_mean_pred, gp_var_pred, n_obs)
    lnorm_log_mean <- lnorm_params$mean_log
    lnorm_log_var <- lnorm_params$var_log
  }
  
  return(dlnorm(1.0, lnorm_log_mean + scaling_factors, sqrt(lnorm_log_var), log = TRUE))
  
}


# TODO: I think the range of the y-axis is just too large 
integrate_loss_1d_log_weights <- function(loss_vals, log_weights, thetas, max_range = log(.Machine$double.xmax), equally_spaced) {
  
  min_idx <- 1
  max_idx <- length(thetas)

  if(diff(range(log_weights)) < max_range) {
    return(loss_vals, exp(log_weights), thetas)
  }
  
  integral_approx <- 0
  max_idx <- floor(max_idx / 2)
  while(max_idx < length(thetas)) {
    m <- min(log_weights[min_idx:max_idx])
    M <- max(log_weights[min_idx:max_idx])
    
    if(M - m > max_range) {
      max_idx <- max(floor(max_idx / 2), min_idx + 1)
    } else {
      shifted_log_weights <- log_weights[min_idx:max_idx] - m
      val <- integrate_loss_1d(loss_vals[min_idx:max_idx], exp(shifted_log_weights), ) 
    }
    
  }  

  
}


integrate_loss_1d <- function(loss_vals, weights, thetas, equally_spaced) {
  # Approximates the integral ell(theta)* w(theta) dtheta over a one-dimensional region, 
  # where ell() is some loss function and w() is some density function (typically the posterior). 
  # Uses a trapezoidal numerical approximation. This function does not create the grid of theta
  # values for the approximation, but rather assumes the arguments are vectors of function 
  # evaluations of ell() and w() at some evenly spaced grid of theta values. 
  # 
  # Args:
  #    loss_vals: numeric(N), vector of ell() evaluations at the grid of evenly spaced points, where 
  #          N is the number of grid points. 
  #    weights: numeric(N), vector of weights for integral; this is often the post() evaluations 
  #             at the grid of evenly spaced points.
  #    thetas: numeric(N), the grid points. 
  #    equally_spaced: logical(1), is the grid `thetas` equally spaced or not? 
  #
  # Returns:
  #    numeric, approximation of the integral given by the trapezoidal rule. 
  
  if(equally_spaced) {
    d_theta <- thetas[2] - thetas[1]
  } else {
    stop("Non equally spaced grid not currently supported.")
  }
  
  N_grid <- length(loss_vals)
  0.5 * d_theta * (loss_vals[1]*weights[1] + loss_vals[N_grid]*weights[N_grid]) + d_theta * sum(loss_vals[2:(N_grid-1)] * weights[2:(N_grid-1)])
  
}


calculate_ground_truth_quantities <- function(theta_train_test_data, computer_model_data, theta_prior_params) {
  # A convenience function used in emulator comparison tests that computes the ground truth log prior, likelihood, and posterior
  # to be used for comparison with the emulated analogs. 
  #
  # Args:
  #    theta_train_test_data: list, as returned by `get_train_test_data()`. 
  #    computer_model_data: list, as returned by `generate_vsem_test_case()`. 
  #    theta_prior_params: data.frame, containing prior information for calibration parameters. 
  # 
  # Returns:
  #    list, containing the calcualted ground truth log density evaluations (log prior, log likelihood, log posterior up to normalizing constant). 

  # True log likelihood evaluations at test inputs
  vars_obs <- diag(computer_model_data$Sig_eps)
  llik_test_per_output <- llik_product_Gaussian(SSR = theta_train_test_data$Y_test, vars_obs = vars_obs, n_obs = computer_model_data$n_obs, 
                                                normalize = TRUE, na.rm = TRUE, sum_output_lliks = FALSE) 
  llik_test <- rowSums(llik_test_per_output)
                                                  
  # True log likelihood evaluations at design points
  llik_train_per_output <- llik_product_Gaussian(SSR = theta_train_test_data$Y_train, vars_obs = vars_obs, n_obs = computer_model_data$n_obs, 
                                                 normalize = TRUE, na.rm = TRUE, sum_output_lliks = FALSE) 
  llik_train <- rowSums(llik_train_per_output)
  
  # Log prior evaluations
  lprior_theta_test <- sapply(seq(1, nrow(theta_train_test_data$X_test)), function(i) calc_lprior_theta(theta_train_test_data$X_test[i,], theta_prior_params))
  lprior_theta_train <- sapply(seq(1, nrow(theta_train_test_data$X_train)), function(i) calc_lprior_theta(theta_train_test_data$X_train[i,], theta_prior_params))
  
  # Log posterior evaluations
  lpost_theta_test <- llik_test + lprior_theta_test
  lpost_theta_train <- llik_train + lprior_theta_train
  
  return(list(llik_test = llik_test, 
              llik_train = llik_train, 
              lprior_test = lprior_theta_test, 
              lprior_train = lprior_theta_train, 
              lpost_test = lpost_theta_test, 
              lpost_train = lpost_theta_train, 
              llik_test_per_output = llik_test_per_output, 
              llik_train_per_output = llik_train_per_output))
  
}


run_emulator_comparison <- function(computer_model_data, emulator_settings, theta_prior_params, train_test_data_settings, true_theta_samples = NULL) {
  # Currently just allow emulator_settings and train_test_data_settings to have multiple options (rows). 
  
  output_list <- list()
  output_list[["emulator_settings"]] <- emulator_settings
  output_list[["theta_prior_params"]] <- theta_prior_params
  
  # Synthetic data generation for VSEM model
  output_list[["computer_model_data"]] <- computer_model_data
  
  # Design data and test data
  output_list[["train_test_data_settings"]] <- train_test_data_settings
  output_list[["theta_train_test_data"]] <- list()
  any_log_SSR = any(emulator_settings["target"] == "log_SSR")
  for(j in seq(1, nrow(train_test_data_settings))) {
    data_setting_name <- paste0("data_setting_", j)
    theta_train_test_data <- get_train_test_data(N_train = train_test_data_settings[j, "N_train"], 
                                                 N_test = train_test_data_settings[j, "N_test"], 
                                                 prior_params = theta_prior_params,
                                                 extrapolate = train_test_data_settings[j, "extrapolate"], 
                                                 ref_pars = computer_model_data$ref_pars,
                                                 pars_cal_sel = computer_model_data$pars_cal_sel,
                                                 data_obs = computer_model_data$data_obs, 
                                                 PAR = computer_model_data$PAR_data,
                                                 output_vars = computer_model_data$output_vars,
                                                 scale_X = train_test_data_settings[j, "scale_X"], 
                                                 normalize_Y = train_test_data_settings[j, "normalize_Y"],
                                                 log_SSR = any_log_SSR, 
                                                 joint_LHS = train_test_data_settings[j, "joint_LHS"], 
                                                 method = train_test_data_settings[j, "method"], 
                                                 order_1d = TRUE, 
                                                 true_theta_samples = true_theta_samples)
    output_list[["theta_train_test_data"]][[j]] <- theta_train_test_data
  }
  
  # Run tests; each test is defined by an emulator setting/train_test_data setting combination.
  output_list[["test_results"]] <- list()
  for(i in seq(1, nrow(emulator_settings))) {
    emulator_setting <- emulator_settings[i, ]
    log_SSR = (emulator_setting["target"] == "log_SSR")
    output_list[["test_results"]][[i]] <- list()
    
    for(j in seq(1, nrow(train_test_data_settings))) {
      data_setting_name <- paste0("data_setting_", j)
      output_list[["test_results"]][[i]][[j]] <- run_emulator_test(computer_model_data, emulator_setting, theta_prior_params,
                                                                   output_list[["theta_train_test_data"]][[j]])
    }
  }
  
  return(output_list)
  
}


# TODO: generalize to LNP
run_emulator_test <- function(computer_model_data, emulator_settings, theta_prior_params, theta_train_test_data) {
  
  Sig_eps <- computer_model_data$Sig_eps
  if(!all(Sig_eps[!diag(nrow(Sig_eps))] == 0)) stop("Non-diagonal output covariances not supported yet.")
  sig2_outputs <- diag(Sig_eps)
  
  # Calculate the ground truth quantities (likelihood, posterior, etc.)
  ground_truth_quantities <- calculate_ground_truth_quantities(theta_train_test_data, computer_model_data, theta_prior_params)
  
  # Fit GPs
  fit_LNP <- (emulator_settings$target == "log_SSR")
  if(fit_LNP) stop("LNP not supported yet.")
  gp_fits_list <- fit_independent_GPs(theta_train_test_data$X_train_preprocessed, theta_train_test_data$Y_train_preprocessed, 
                                      emulator_settings$gp_lib, emulator_settings$kernel)
  gp_fit_times <- gp_fits_list$times
  gp_fits <- gp_fits_list$fits
  
  # Predict with GPs: returns predictions from GP without performing transformations (exponentiating, truncating, rectifying). These transformations 
  # are intended to be performed later, e.g. when plotting with plot_gp_1d(). 
  gp_pred_test <- predict_independent_GPs(theta_train_test_data$X_test_preprocessed, gp_fits, emulator_settings$gp_lib, cov_mat = FALSE, denormalize_predictions = TRUE,
                                          output_stats = theta_train_test_data$output_stats, exponentiate_predictions = FALSE, 
                                          transformation_method = "default")
  gp_pred_train <- predict_independent_GPs(theta_train_test_data$X_train_preprocessed, gp_fits, emulator_settings$gp_lib, cov_mat = FALSE, denormalize_predictions = TRUE,
                                           output_stats = theta_train_test_data$output_stats, exponentiate_predictions = FALSE, 
                                           transformation_method = "default")
  
  # TODO: left off here:
  #    - Need to update the below lines of code since need to perform the transformations before calculating SSR. 
  #    - Should probably have a calc_SSR() function that allows the transformations to be performed. Or I should 
  #      just perform the transformations and store them as separate elements of the list. Should probably have a 
  #      "transform_gp_predictions()" function that accounts for all of the transformations.
  
  SSR_pred_mean_test <- sapply(gp_pred_test, function(pred) pred$mean)
  SSR_pred_var_test <- sapply(gp_pred_test, function(pred) pred$sd2)
  SSR_pred_mean_train <- sapply(gp_pred_train, function(pred) pred$mean)
  SSR_pred_var_train <- sapply(gp_pred_train, function(pred) pred$sd2)
  
  # Compute GP predictive density of posterior approximation at test points
  # post_approx_LNP_params_test <- get_gp_approx_posterior_LNP_params(diag(Sig_eps), ground_truth_quantities$lprior_test, 
  #                                                                   SSR_pred_mean_test, SSR_pred_var_test, computer_model_data$n_obs)
  # lpred_pi_test <- gp_approx_posterior_pred_log_density(ground_truth_quantities$lpost_test, 
  #                                                       lnorm_log_mean = post_approx_LNP_params_test$mean_log, lnorm_log_var = post_approx_LNP_params_test$var_log)
  # lpred_pi_train <- gp_approx_posterior_pred_log_density(ground_truth_quantities$lpost_train, sig2_outputs, 
  #                                                        ground_truth_quantities$lprior_train, SSR_pred_mean_train, SSR_pred_var_train, computer_model_data$n_obs)
  
  return(list(ground_truth_quantities = ground_truth_quantities, 
              gp_fits = gp_fits, 
              gp_times = gp_fit_times, 
              gp_pred_list_test = gp_pred_test, 
              gp_pred_list_train = gp_pred_train, 
              SSR_pred_mean_test = SSR_pred_mean_test, 
              SSR_pred_var_test = SSR_pred_var_test, 
              SSR_pred_mean_train = SSR_pred_mean_train, 
              SSR_pred_var_train = SSR_pred_var_train))
              
              # post_pred_mean_log_test = post_approx_LNP_params_test$mean_log, 
              # post_pred_var_log_test = post_approx_LNP_params_test$var_log, 
              # lpost_pred_density_test = lpred_pi_test, 
              # lpost_pred_density_train = lpred_pi_train))  
  
}






