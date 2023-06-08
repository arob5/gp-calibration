#
# mcmc_calibration_functions.r
# Functions defining likelihoods and priors, and implementions of MCMC for parameter
# calibration of computer model. Used for tests with Very Simple Ecosystem (VSEM)
# model. 
#
# Andrew Roberts
#

library(dplyr)
library(parallel)
library(mvtnorm)
library(tmvtnorm)
library(truncnorm)
library(LaplacesDemon)

# TODO: 
#    - Updae llik_product_Gaussian() and functions (like the posterior density functions) that rely on this 
#      to fit within the new generic computer model framework. 
#    - Update comments for `get_input_output_design`.


llik_product_Gaussian <- function(computer_model_data, vars_obs, theta_vals = NULL, SSR = NULL, 
                                  normalize = TRUE, na.rm = FALSE, sum_output_lliks = TRUE) {
  # Evaluate Gaussian log-likelihood assuming independence across time and across
  # output variables. In order to calculate the log-likelihood, the forward model must be run and then 
  # the L2 difference between the forward model outputs and the observed data evaluated. This function can 
  # perform all of those steps, or the L2 error can be passed directly using the `SSR` argument to avoid 
  # running the forward model when it is not necessary. For the forward model to be run, the argument
  # `theta_vals` must be passed. Otherwise, this is not needed if `SSR` is provided. This function is 
  # vectorized in that it can evaluate  the log-likelihood at multiple input parameter values, `theta_vals`.
  #
  # Args:
  #    computer_model_data: the standard computer model data list. 
  #    vars_obs: numeric(p), vector of observation/noise variances for each output. 
  #    theta_vals: matrix, of dimension M x d, where M is the number of calibration parameter values
  #                and d is the number of calibration parameters. These are the parameter inputs at which 
  #                to run the forward to evaluate log-likelihood. 
  #    SSR: matrix, of dimension M x p. The (m, j) entry is the sum of squared errors for the 
  #         jth output variable and mth input parameter. 
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
    SSR <- get_computer_model_SSR(computer_model_data = computer_model_data,
                                  theta_vals = theta_vals, 
                                  na.rm = na.rm)
  }
  
  p <- length(vars_obs)
  n_obs <- computer_model_data$n_obs
  llik_outputs <- matrix(nrow = nrow(SSR), ncol = ncol(SSR))
  
  for(j in seq(1, p)) {
    llik_outputs[,j] <- -0.5 * n_obs[j] * log(vars_obs[j]) - 0.5 * SSR[,j] / vars_obs[j]
    if(normalize) llik_outputs[,j] <- llik_outputs[,j] - 0.5 * n_obs[j] * log(2*pi)
  }
  
  # Either return likelihood for each output separately, or return single likelihood across all outputs. 
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


llik_Gaussian_err <- function(model_errs, Sig_eps = NULL, L = NULL, output_vars = NA, normalize = TRUE) {
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
  #     L: matrix of dimension p x p, the lower Cholesky factor of Sig_eps. If not provided, the Cholesky 
  #        factor will be computed. 
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
  
  if(is.null(L)) {
    L <- t(chol(Sig_eps))
  }
  
  log_quadratic_form <- sum(forwardsolve(L, t(model_errs))^2)

  if(normalize) {
    return(-0.5 * log_quadratic_form - 0.5 * prod(dim(model_errs)) * log(2*pi) - sum(log(diag(L))))
  }
  
  return(-0.5 * log_quadratic_form)
  
}


run_VSEM <- function(theta_vals, computer_model_data = NULL, ref_pars = NULL, pars_cal_sel = NULL, PAR_data = NULL, output_vars = NULL) {
  # Runs the VSEM model using the specified parameter settings, thus returning the outputs of the forward model. 
  # If multiple input parameters are passed, then the model is run at each input and the results are returned in a 
  # list, one for each input. If a single input parameter is passed, this function reduces to `run_VSEM_single_input()`.
  # Note that `run_computer_model()` is a generic interface that should typically be used to execute a forward model, instead 
  # of directly calling model-specific functions. 
  #
  # Args:
  #    theta_vals: matrix of dimension M x d, the values of calibration parameters used to run the model. Each row is 
  #                a setting in the d-dimensional space of calibration parameters. 
  #    computer_model_data: list, must have named elements "par_ref", "par_cal_sel", "PAR_data", "output_vars".
  #    ref_pars: data.frame, rownames should correspond to parameters of computer model. 
  #             Must contain column named "true_value". Parameters that are fixed at nominal 
  #             values (not calibrated) are set to their values given in the "true_value" column. 
  #    pars_cal_sel: integer vector, selects the rows of 'par_ref' that correspond to 
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
  
  if(!is.matrix(theta_vals) || (nrow(theta_vals) == 1)) {
    return(run_VSEM_single_input(par_val = theta_vals, computer_model_data = computer_model_data, 
                                 ref_pars = ref_pars, pars_cal_sel = pars_cal_sel, PAR_data = PAR_data, 
                                 output_vars = output_vars))
  } else {
    return(apply(theta_vals, 1, function(theta) run_VSEM_single_input(par_val = theta, computer_model_data = computer_model_data, 
                                                                      ref_pars = ref_pars, pars_cal_sel = pars_cal_sel, PAR_data = PAR_data, 
                                                                      output_vars = output_vars), simplify = FALSE))
  }
  
}


run_VSEM_single_input <- function(par_val, computer_model_data = NULL, ref_pars = NULL, pars_cal_sel = NULL, 
                                  PAR_data = NULL, output_vars = NULL) {
  # Runs the VSEM model using specified parameter setting, returning the outputs of the model. This runs
  # the computer model as a single input, intended to work as a generic computer model`f()` mapping 
  # as used in the `run_computer_model()` function. 
  #
  # Args:
  #    computer_model_data: list, must have named elements "par_ref", "par_cal_sel", "PAR_data", "output_vars". 
  #    par_val: numeric vector, values of calibration parameters used to run the model. 
  #    The remaining parameters provide an alternative to passing `computer_model_data`; these are primarily 
  #    intended to be used when initially creating data, before the `computer_model_list` is generated. 
  #
  # Returns:
  #   matrix of dimensions n x p where n is the length of the time series and p is the number
  #   of output variables. 
  
  if(!is.null(computer_model_data)) {
    ref_pars <- computer_model_data$ref_pars
    pars_cal_sel <- computer_model_data$pars_cal_sel
    PAR_data <- computer_model_data$PAR_data
    output_vars <- computer_model_data$output_vars
  }
  
  # Parameters not calibrated are fixed at default values
  theta <- ref_pars$true_value
  theta[pars_cal_sel] <- par_val
  
  # Run forward model, re-scale NEE.
  VSEM_output <- as.matrix(VSEM(theta, PAR_data))
  if("NEE" %in% output_vars) {
    VSEM_output[, "NEE"] <- VSEM_output[, "NEE"] * 1000
  }

  # Compute LAI, if included in output variables. LAI is simply LAR times the above-ground vegetation pool
  # at time t. 
  if("LAI" %in% output_vars) {
    VSEM_output <- cbind(VSEM_output, theta[rownames(ref_pars) == "LAR"] * VSEM_output[, "Cv"])
    colnames(VSEM_output)[ncol(VSEM_output)] <- "LAI"
  }
  
  # Select only the output variables 
  VSEM_output <- VSEM_output[, output_vars]
  
  return(VSEM_output)
  
}


# Deprecated: replaced by new generic run_computer_model() framework. 
run_VSEM_single_input_old <- function(par_val, ref_pars = NULL, pars_cal_sel = NULL, PAR_data = NULL, output_vars = c("NEE", "Cv", "Cs", "CR"), 
                                      computer_model_data = NULL) {
  # Runs the VSEM model using specified parameter setting, returning the outputs of the model. Can either pass in the necessary data 
  # to run the forward model as individual arguments (par_ref, par_cal_sel, PAR_data, output_vars), or can pass in a list `computer_model_data`
  # will all of these arguments are named elements. 
  #
  # Args:
  #    par_val: numeric vector, values of calibration parameters used to run the model. 
  #    ref_pars: data.frame, rownames should correspond to parameters of computer model. 
  #             Must contain column named "true_value". Parameters that are fixed at nominal 
  #             values (not calibrated) are set to their values given in the "true_value" column. 
  #    pars_cal_sel: integer vector, selects the rows of 'par_ref' that correspond to 
  #                 parameters that will be calibrated. 
  #    PAR_data: numeric vector, time series of photosynthetically active radiation used as forcing
  #         term in VSEM. 
  #    output_vars: character vector, selects columns of output from VSEM model. Default returns
  #                 all four columns. 
  #    computer_model_data: list, an alternative to having to pass in the computer model data one argument at a time. Must have named 
  #                         elements "par_ref", "par_cal_sel", "PAR_data", "output_vars". 
  #
  # Returns:
  #   matrix of dimensions n x p where n is the length of the time series and p is the number
  #   of output variables. 
  
  if(!is.null(computer_model_data)) {
    ref_pars <- computer_model_data$ref_pars
    pars_cal_sel <- computer_model_data$pars_cal_sel
    PAR_data <- computer_model_data$PAR_data
    output_vars <- computer_model_data$output_vars
  }
  
  # Parameters not calibrated are fixed at default values
  theta <- ref_pars$true_value
  theta[pars_cal_sel] <- par_val
  
  # Run forward model, re-scale NEE.
  VSEM_output <- as.matrix(VSEM(theta, PAR_data))
  if("NEE" %in% output_vars) {
    VSEM_output[, "NEE"] <- VSEM_output[, "NEE"] * 1000
  }
  
  # Compute LAI, if included in output variables. LAI is simply LAR times the above-ground vegetation pool
  # at time t. 
  if("LAI" %in% output_vars) {
    VSEM_output <- cbind(VSEM_output, theta[rownames(ref_pars) == "LAR"] * VSEM_output[, "Cv"])
    colnames(VSEM_output)[ncol(VSEM_output)] <- "LAI"
  }
  
  # Select only the output variables 
  VSEM_output <- VSEM_output[, output_vars]
  
  return(VSEM_output)
  
}


calc_lprior_theta <- function(theta, theta_prior_params, check_bounds = FALSE) {
  # Evaluates the log prior density on calibration functions at specific values of the settings.  
  #
  # Args:
  #    theta: numeric vector, the value of the calibration parameters at which to evaluate the prior density. 
  #    theta_prior_params: data.frame, with columns "dist", "param1", and "param2", and potentially columns 
  #                        "bound_lower" and "bound_upper". The ith row of the data.frame
  #                        should correspond to the ith entry of 'theta'. Currently, accepted values of "dist" are 
  #                        "Gaussian" (param1 = mean, param2 = std dev) and "Uniform" (param1 = lower, param2 = upper)
  #                         and "Truncatd_Gaussian" (param1 = mean, param2 = std dev, bound_lower = lower truncation value, 
  #                         bound_upper = upper truncation value). 
  #    check_bounds: logical(1), if TRUE checks if the parameter `theta` lies within the bounds defined in the 
  #                  columns of `theta_prior_params` called "bound_lower" and "bound_upper". If these columns do not 
  #                  exist there is not effect. Similarly, if the columns exist but some rows have NA values these rows 
  #                  will also be ignored. Otherwise, if `theta` does not lie within the bounds, then the function returns 
  #                  -Inf. This prevents having to evaluate all of the prior densities, which will result in the same return 
  #                  value. If `check_bounds` is FALSE, no check is performed. Default is FALSE 
  #
  # Returns:
  #    The prior density evaluation log p(theta). Assumes prior independence, so the log-prior is the sum of the log-prior
  #    evaluations for each entry of 'theta'. Note that in certain cases the log prior can be negative infinity; e.g. for 
  #    a uniform prior where `theta` is not contained within the upper and lower bound. 
  
  if(check_bounds) {
    if(any(theta < theta_prior_params[["bound_lower"]], na.rm = TRUE) ||
       any(theta > theta_prior_params[["bound_upper"]], na.rm = TRUE)) {
      return(-Inf)
    }
  }
  
  lprior <- 0
  
  theta_prior_params[["val"]] <- theta
  Gaussian_priors <- theta_prior_params[theta_prior_params$dist == "Gaussian",]
  Uniform_priors <- theta_prior_params[theta_prior_params$dist == "Uniform",]
  Truncated_Gaussian_priors <- theta_prior_params[theta_prior_params$dist == "Truncated_Gaussian",]
  
  if(nrow(Gaussian_priors) > 0) {
    lprior <- lprior + sum(dnorm(Gaussian_priors$val, Gaussian_priors$param1, Gaussian_priors$param2, log = TRUE))
  }
  
  if(nrow(Uniform_priors) > 0) {
    lprior <- lprior + sum(dunif(Uniform_priors$val, Uniform_priors$param1, Uniform_priors$param2, log = TRUE))
  }
  
  if(nrow(Truncated_Gaussian_priors) > 0) {
    lprior <- lprior + sum(sapply(seq(1, nrow(Truncated_Gaussian_priors)), function(i) log(dtruncnorm(Truncated_Gaussian_priors$val[i], 
                                                                                           a = Truncated_Gaussian_priors$bound_lower[i], 
                                                                                           b = Truncated_Gaussian_priors$bound_upper[i], 
                                                                                           mean = Truncated_Gaussian_priors$param1[i],
                                                                                           sd = Truncated_Gaussian_priors$param2[i]))))
  }
  
  return(lprior)  

}


calc_lpost_theta_product_lik <- function(computer_model_data, lprior_vals = NULL, llik_vals = NULL, theta_vals = NULL, 
                                         SSR = NULL, vars_obs = NULL, normalize_lik = TRUE, na.rm = FALSE,
                                         theta_prior_params = NULL, return_list = TRUE) {
  # Evaluates the exact log-posterior density, up to the normalizing constant (model evidence). 
  # This functions provides the option to calculate the log prior and log likelihood from scratch, at the given parameter values `theta_vals` 
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
    lprior_vals <- apply(theta_vals, 1, function(theta) calc_lprior_theta(theta, theta_prior_params))
  }
  
  # Likelihood
  if(is.null(llik_vals)) {
    llik_vals <- llik_product_Gaussian(computer_model_data = computer_model_data, 
                                       vars_obs = vars_obs, 
                                       theta_vals = theta_vals, 
                                       SSR = SSR, 
                                       normalize = normalize_lik, 
                                       na.rm = na.rm, 
                                       sum_output_lliks = TRUE)
  }
  
  if(!return_list) return(lprior_vals + llik_vals)
  
  return(list(lpost = lprior_vals + llik_vals, 
              lprior = lprior_vals, 
              llik = llik_vals))
  
}


run_computer_model <- function(theta_vals, computer_model_data) {
  # This is the generic interface for evaluating a computer model at a single or set of multiple input 
  # parameters and returning the resulting outputs. `theta_vals` can either be a single input 
  # parameter vector (a numeric vector or a matrix with a single row), or it can be a 
  # matrix of dimension N_param x d, where d is the dimension of the parameter input space. In the former 
  # case the model output is returned directly by evaluating at the single parameter vector. In the latter 
  # case, a list is returned in which element i corresponds to the output resulting from the computer 
  # model evaluation at the parameter vector in the ith row of `theta_vals`. The computer model output 
  # is a matrix of dimension N x p, where N is typically the number of time steps and p is the 
  # number of outputs/data constraints. 
  #
  # Args:
  #    theta_vals: For a single parameter, a numeric(d) vector or matrix of dimension 1xd. For multiple 
  #                parameters a matrix of dimension N_param x d. 
  #    computer_model_data: the standard computer model data list. 
  #
  # Returns:
  #    Either Nxp matrix of list of Nxp matrices of length equal to the number of rows in `theta_vals`. 
  #    See above description for details. 
  
  if(!is.matrix(theta_vals) || (nrow(theta_vals) == 1)) {
    return(computer_model_data$f(theta_vals, computer_model_data = computer_model_data))
  } else {
    return(apply(theta_vals, 1, function(theta) computer_model_data$f(theta, computer_model_data = computer_model_data), simplify = FALSE))
  }
  
}


get_computer_model_errs <- function(theta_vals, computer_model_data) {
  # A convenience function that first runs the computer model to obtain outputs f(theta)
  # and then returns the errors Y - f(theta) with respect to observed data Y. 
  # Both Y and f(thets) are assumed to be matrices of shape N x p, where
  # N is the number of observations and p is the number of model outputs. This function 
  # works for a single input parameter vector or a matrix of multiple parameters. 
  # In the latter case, the errors are returned in a list, with each element corresponding 
  # to a different parameter. See `run_computer_model()` for details on the single-parameter
  # vs. multiple-parameter argument requirements. 
  #
  # Args:
  #    theta_vals: For a single parameter, a numeric(d) vector or matrix of dimension 1xd. For multiple 
  #                parameters a matrix of dimension N_param x d. 
  #    computer_model_data: the standard computer model data list. 
  #
  # Returns:
  #    Either Nxp matrix or list of Nxp matrices of length equal to the number of rows in `theta_vals`. 
  #    See above description for details. 
  
  output_vars <- computer_model_data$output_vars
  computer_model_output <- run_computer_model(theta_vals, computer_model_data)
  data_obs <- computer_model_data$data_obs[, output_vars, drop=FALSE]
  
  if(is.list(computer_model_output)) {
    model_errs <- lapply(computer_model_output, function(output) data_obs - output)
  } else {
    model_errs <- data_obs - computer_model_output
  }
  
  return(model_errs)

}


get_computer_model_SSR <- function(computer_model_data, model_outputs_list = NULL, theta_vals = NULL, na.rm = TRUE) {
  # Computes the sum of squared residuals (SSR) between model runs and observed 
  # data on a per-output basis. Can handle multiple model runs (e.g. one per 
  # design point for emulation) or outputs from single model run (e.g. as required
  # during MCMC).
  #
  # Args:
  #    computer_model_data: list, the standard computer model data list. 
  #    model_outputs_list: either a matrix of dimension n x p corresponding to the 
  #                        model outputs f(theta) from a single model run. Or a list 
  #                        {f(theta_1), ..., f(theta_N)} of such matrices, collecting
  #                        the outputs from multiple model runs. 
  #    theta_vals: numeric vector or matrix of input values, must be provided if `model_outputs_list` 
  #                is NULL, in which case the forward model  will be run at these inputs to obtain 
  #                `model_outputs_list`. 
  #    na.rm: logical(1), whether or not to remove NA values from the sum of squares 
  #           calculation, passed to the `colSums()` functions. Default is TRUE. 
  # 
  # Returns:
  #    matrix of dimension N x p, where N is the number of model runs and p is the 
  #    number of output variables. The (i, j) entry of the matrix is the SSR
  #    for the jth output of the ith model run, i.e. ||Y_j - f(j, theta_i)||^2.
  
  # Observed data. 
  data_obs <- computer_model_data$data_obs[, computer_model_data$output_vars, drop = FALSE]
  
  # Run forward model if it has not been run yet. 
  if(is.null(model_outputs_list)) {
    model_outputs_list <- run_computer_model(theta_vals = theta_vals, computer_model_data = computer_model_data)
  }
  
  # If model only run at single set of calibration parameter values. 
  if(is.matrix(model_outputs_list)) {
    model_outputs_list <- list(model_outputs_list)
  }
  
  N_param_runs <- length(model_outputs_list)
  N_outputs <- ncol(model_outputs_list[[1]])
  SSR <- matrix(nrow = N_param_runs, ncol = N_outputs)
  
  for(i in seq_along(model_outputs_list)) {
    SSR[i,] <- colSums((data_obs - model_outputs_list[[i]])^2, na.rm = na.rm)
  }
  colnames(SSR) <- colnames(data_obs)
  
  return(SSR)
  
}


# TODO: need to update this so that it can handle missing data. Main issue is updating the Gaussian likelihood function 
# that operates on model_errs. 
# TODO: Update description and argument comments; have changes this function to be exclusively for the case where the 
#       covariance matrix Sig_eps is non-diagonal. 
mcmc_calibrate <- function(computer_model_data, theta_prior_params,
                           theta_init = NULL, Sig_eps_init = NULL, learn_Sig_eps = FALSE, Sig_eps_prior_params = NULL, 
                           N_mcmc = 50000, adapt_frequency = 1000, adapt_min_scale = 0.1, accept_rate_target = 0.24, 
                           proposal_scale_decay = 0.7, Cov_prop_init_diag = 0.1, adapt_cov_method = "AM", 
                           adapt_scale_method = "MH_ratio", adapt_init_threshold = 3) {
  # MCMC implementation for VSEM carbon model. Accommodates Gaussian likelihood, possibly with correlations 
  # between the different output variables, but assumes independence across time. Samples from posterior over both  
  # calibration parameters (theta) and observation covariance (Sig_eps), or just over theta if `learn_Sig_eps` is FALSE.
  # Allows for arbitrary prior over theta, but assumes Inverse Wishart prior on Sig_eps (if treated as random).
  # MCMC scheme is adaptive random walk Metropolis. This MCMC algorithm does not involve any model emulation/likelihood 
  # approximation. 
  #
  # Args:
  #    computer_model_data:
  #    theta_prior_params: data.frame, with columns "dist", "param1", and "param2". The ith row of the data.frame
  #                        should correspond to the ith entry of 'theta'. Currently, accepted values of "dist" are 
  #                        "Gaussian" (param1 = mean, param2 = std dev) and "Uniform" (param1 = lower, param2 = upper).
  #    diag_cov: logical, if TRUE constrains Sig_eps to be diagonal. If the prior distribution is specified to be product Inverse Gamma then this 
  #              is automatically set to TRUE. Default is FALSE. 
  #    theta_init: numeric vector of length p, the initial value of the calibration parameters to use in MCMC. If NA, samples
  #                the initial value from the prior. 
  #    learn_Sig_eps: logical, if TRUE treats observation covariance matrix as random and MCMC samples from joint 
  #                   posterior over Sig_eps and theta. Otherwise, fixes Sig_eps at value `Sig_eps_init`.
  #    Sig_eps_init: matrix, p x p covariance matrix capturing dependence between output variables. If `learn_Sig_eps` is 
  #                  TRUE then `Sig_eps_init` can either be set to the initial value used for MCMC, or set to NA in which 
  #                  case the initial value will be sampled from the prior. If `learn_Sig_eps` is FALSE, then a non-NA value 
  #                  is required, and treated as the fixed nominal value of Sig_eps_init. 
  #    Sig_eps_prior_params: list, defining prior on Sig_eps. See `sample_prior_Sig_eps()` for details. Only required if `learn_Sig_eps` is TRUE.
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

  # Number observations in time series, number output variables, and dimension of parameter space.
  p <- length(computer_model_data$output_vars)
  d <- length(computer_model_data$pars_cal_names)
  
  # Ensure parameters in `theta_prior_params` are sorted correctly based on ordering in `computer_model_data`. 
  theta_prior_params <- theta_prior_params[computer_model_data$pars_cal_names,]
  
  # Objects to store samples.
  theta_samp <- matrix(nrow = N_mcmc, ncol = d)
  colnames(theta_samp) <- paste("theta", computer_model_data$pars_cal_names, sep = "_")
  Sig_eps_samp <- matrix(nrow = N_mcmc, ncol = 0.5*p*(p+1)) # Each row stores lower triangle of Sig_eps
  output_var_sel <- matrix(1:p^2, nrow = p)[lower.tri(matrix(1:p^2, nrow = p), diag = TRUE)]
  colnames(Sig_eps_samp) <- paste("sig_eps", output_var_sel, sep = "_") # Should replace this with actual output var names (e.g. "sig_eps_y1_y2").

  # Set initial conditions.
  if(is.null(theta_init)) {
    theta_init <- sample_prior_theta(theta_prior_params)
  }
  if(learn_Sig_eps) {
    if(is.null(Sig_eps_init)) {
      Sig_eps_init <- sample_prior_Sig_eps(Sig_eps_prior_params, return_matrix = TRUE)
    }
  } else {
    if(is.null(Sig_eps_init)) stop("Value for `Sig_eps_init` must be provided when `learn_Sig_eps` is FALSE.")
  }
  
  theta_samp[1,] <- theta_init
  Sig_eps_samp[1,] <- Sig_eps_init[lower.tri(Sig_eps_init, diag = TRUE)]

  Sig_eps_curr <- Sig_eps_init
  L_Sig_eps_curr <- t(chol(Sig_eps_init))
  model_errs_curr <- get_computer_model_errs(theta_init, computer_model_data)
  lprior_theta_curr <- calc_lprior_theta(theta_init, theta_prior_params)
  
  # Proposal covariance
  Cov_prop <- diag(Cov_prop_init_diag, nrow = d)
  log_scale_prop <- 0
  L_prop <- t(chol(Cov_prop))
  accept_count <- 0
  samp_mean <- theta_init
  
  for(itr in seq(2, N_mcmc)) {
    #
    # Metropolis step for theta
    #

    # theta proposals.
    theta_prop <- theta_samp[itr-1,] + (exp(log_scale_prop) * L_prop %*% matrix(rnorm(d), ncol = 1))[,1]
    
    # Calculate log-likelihoods for current and proposed theta.
    model_errs_prop <- get_computer_model_errs(theta_prop, computer_model_data)
    llik_curr <- llik_Gaussian_err(model_errs_curr, L = L_Sig_eps_curr) 
    llik_prop <- llik_Gaussian_err(model_errs_prop, L = L_Sig_eps_curr)
    
    # Metropolis-Hastings Accept-Reject Step.
    lprior_theta_prop <- calc_lprior_theta(theta_prop, theta_prior_params)
    lpost_theta_curr <- llik_curr + lprior_theta_curr
    lpost_theta_prop <- llik_prop + lprior_theta_prop
    alpha <- min(1.0, exp(lpost_theta_prop - lpost_theta_curr))
    
    if(runif(1) <= alpha) {
      theta_samp[itr,] <- theta_prop
      lprior_theta_curr <- lprior_theta_prop
      model_errs_curr <- model_errs_prop
      accept_count <- accept_count + 1 
    } else {
      theta_samp[itr,] <- theta_samp[itr-1,]
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
    # Gibbs step for Sig_eps
    #
    if(learn_Sig_eps) {
      Sig_eps_curr <- sample_cond_post_Sig_eps(model_errs = model_errs_curr, 
                                               Sig_eps_prior_params = Sig_eps_prior_params, 
                                               n_obs = computer_model_data$n_obs, return_matrix = TRUE)
      L_Sig_eps_curr <- t(chol(Sig_eps_curr))
      Sig_eps_samp[itr,] <- Sig_eps_curr[lower.tri(Sig_eps_curr, diag = TRUE)]
    }
    
  }
  
  return(list(samp = cbind(theta_samp, Sig_eps_samp), Cov_prop = Cov_prop, scale_prop = exp(log_scale_prop)))
  
}


# TODO: update comments. 
mcmc_calibrate_product_lik <- function(computer_model_data, theta_prior_params,
                                       theta_init = NULL, sig_eps_init = NULL, learn_sig_eps = FALSE, sig_eps_prior_params = NULL, 
                                       N_mcmc = 50000, adapt_frequency = 1000, adapt_min_scale = 0.1, accept_rate_target = 0.24, 
                                       proposal_scale_decay = 0.7, Cov_prop_init_diag = 0.1, adapt_cov_method = "AM", 
                                       adapt_scale_method = "MH_ratio", adapt_init_threshold = 3) {
  # MCMC implementation for VSEM carbon model. Accommodates Gaussian likelihood, possibly with correlations 
  # between the different output variables, but assumes independence across time. Samples from posterior over both  
  # calibration parameters (theta) and observation covariance (Sig_eps), or just over theta if `learn_Sig_eps` is FALSE.
  # Allows for arbitrary prior over theta, but assumes Inverse Wishart prior on Sig_eps (if treated as random).
  # MCMC scheme is adaptive random walk Metropolis. This MCMC algorithm does not involve any model emulation/likelihood 
  # approximation. 
  #
  # Args:
  #    computer_model_data: list, standard computer model data list. 
  #    theta_prior_params: data.frame, with columns "dist", "param1", and "param2". The ith row of the data.frame
  #                        should correspond to the ith entry of 'theta'. Currently, accepted values of "dist" are 
  #                        "Gaussian" (param1 = mean, param2 = std dev) and "Uniform" (param1 = lower, param2 = upper).
  #    diag_cov: logical, if TRUE constrains Sig_eps to be diagonal. If the prior distribution is specified to be product Inverse Gamma then this 
  #              is automatically set to TRUE. Default is FALSE. 
  #    theta_init: numeric vector of length p, the initial value of the calibration parameters to use in MCMC. If NA, samples
  #                the initial value from the prior. 
  #    learn_Sig_eps: logical, if TRUE treats observation covariance matrix as random and MCMC samples from joint 
  #                   posterior over Sig_eps and theta. Otherwise, fixes Sig_eps at value `Sig_eps_init`.
  #    Sig_eps_init: matrix, p x p covariance matrix capturing dependence between output variables. If `learn_Sig_eps` is 
  #                  TRUE then `Sig_eps_init` can either be set to the initial value used for MCMC, or set to NA in which 
  #                  case the initial value will be sampled from the prior. If `learn_Sig_eps` is FALSE, then a non-NA value 
  #                  is required, and treated as the fixed nominal value of Sig_eps_init. 
  #    Sig_eps_prior_params: list, defining prior on Sig_eps. See `sample_prior_Sig_eps()` for details. Only required if `learn_Sig_eps` is TRUE.
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
  
  if(learn_sig_eps && sig_eps_prior_params$dist != "IG") {
    stop("`mcmc_calibrate_product_lik()` requires inverse gamma priors on observation variances.")
  }
  
  # Number observations in time series, number output variables, and dimension of parameter space.
  p <- length(computer_model_data$output_vars)
  d <- length(computer_model_data$pars_cal_names)
  
  # Ensure parameters in `theta_prior_params` are sorted correctly based on ordering in `computer_model_data`. 
  theta_prior_params <- theta_prior_params[computer_model_data$pars_cal_names,]
  
  # Objects to store samples.
  theta_samp <- matrix(nrow = N_mcmc, ncol = d)
  colnames(theta_samp) <- paste("theta", computer_model_data$pars_cal_names, sep = "_")
  sig_eps_samp <- matrix(nrow = N_mcmc, ncol = p)
  colnames(sig_eps_samp) <- paste("sig_eps", computer_model_data$output_vars, sep = "_")

  # Set initial conditions.
  if(is.null(theta_init)) {
    theta_init <- sample_prior_theta(theta_prior_params)
  }
  if(learn_sig_eps) {
    if(is.null(sig_eps_init)) {
      sig_eps_init <- sample_prior_Sig_eps(sig_eps_prior_params)
    }
  } else {
    if(is.null(sig_eps_init)) stop("Value for `sig_eps_init` must be provided when `learn_sig_eps` is FALSE.")
  }
  
  theta_samp[1,] <- theta_init
  sig_eps_samp[1,] <- sig_eps_init

  sig_eps_curr <- sig_eps_init
  SSR_curr <- get_computer_model_SSR(computer_model_data, theta_vals = theta_init, na.rm = TRUE)
  lprior_theta_curr <- calc_lprior_theta(theta_init, theta_prior_params)
  
  # Proposal covariance
  Cov_prop <- diag(Cov_prop_init_diag, nrow = d)
  log_scale_prop <- 0
  L_prop <- t(chol(Cov_prop))
  accept_count <- 0
  samp_mean <- theta_init
  
  for(itr in seq(2, N_mcmc)) {
    #
    # Metropolis step for theta
    #
    
    # theta proposals.
    theta_prop <- theta_samp[itr-1,] + (exp(log_scale_prop) * L_prop %*% matrix(rnorm(d), ncol = 1))[,1]
    
    # Calculate log-likelihoods for current and proposed theta.
    SSR_prop <- get_computer_model_SSR(computer_model_data, theta_vals = theta_prop, na.rm = TRUE)
    llik_prop <- llik_product_Gaussian(computer_model_data, sig_eps_curr, SSR = SSR_prop, normalize = FALSE, na.rm = TRUE)
    llik_curr <- llik_product_Gaussian(computer_model_data, sig_eps_curr, SSR = SSR_curr, normalize = FALSE, na.rm = TRUE) 

    # Metropolis-Hastings Accept-Reject Step.
    lprior_theta_prop <- calc_lprior_theta(theta_prop, theta_prior_params)
    lpost_theta_curr <- llik_curr + lprior_theta_curr
    lpost_theta_prop <- llik_prop + lprior_theta_prop
    alpha <- min(1.0, exp(lpost_theta_prop - lpost_theta_curr))
    
    if(runif(1) <= alpha) {
      theta_samp[itr,] <- theta_prop
      lprior_theta_curr <- lprior_theta_prop
      SSR_curr <- SSR_prop
      accept_count <- accept_count + 1 
    } else {
      theta_samp[itr,] <- theta_samp[itr-1,]
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
    # Gibbs step for sig_eps
    #
    if(learn_sig_eps) {
      sig_eps_curr <- sample_cond_post_Sig_eps(SSR = SSR_curr, Sig_eps_prior_params = sig_eps_prior_params, n_obs = computer_model_data$n_obs)
      sig_eps_samp[itr,] <- sig_eps_curr
    }
    
  }
  
  return(list(samp = cbind(theta_samp, sig_eps_samp), Cov_prop = Cov_prop, scale_prop = exp(log_scale_prop)))
  
}


# emulator_info: list with "gp_fits", "output_stats", "settings", and "input_bounds"
# TODO: allow joint sampling, incorporating covariance between current and proposal. "output_stats" must be 
# on the correct scale (e.g. it should be on log scale for LNP).
mcmc_calibrate_ind_GP <- function(computer_model_data, theta_prior_params, emulator_info,
                                  theta_init = NULL, sig_eps_init = NULL, learn_sig_eps = FALSE, sig_eps_prior_params = NULL, 
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
  colnames(theta_samp) <- paste("theta", computer_model_data$pars_cal_names, sep = "_")
  sig_eps_samp <- matrix(nrow = N_mcmc, ncol = p)
  colnames(sig_eps_samp) <- paste("sig_eps", computer_model_data$output_vars, sep = "_")
  
  # Set initial conditions. 
  if(is.null(theta_init)) {
    theta_init <- sample_prior_theta(theta_prior_params)
  }

  if(learn_sig_eps) {
    if(is.null(sig_eps_init)) {
      sig_eps_init <- sample_prior_Sig_eps(sig_eps_prior_params)
    }
  } else {
    if(is.null(sig_eps_init)) stop("Value for `sig_eps_init` must be provided when `learn_sig_eps` is FALSE.")
  }
  
  theta_samp[1,] <- theta_init
  sig_eps_samp[1,] <- sig_eps_init
  
  sig_eps_curr <- sig_eps_init
  lprior_theta_curr <- calc_lprior_theta(theta_init, theta_prior_params)
  
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
    
      # Approximate SSR by sampling from GP. 
      # TODO: shouldn't have to re-scale theta_curr; should have `theta_curr_scaled` variable or something. 
      # Could also re-use mean prediction, but if including GP covariance then this will change. 
      thetas_scaled <- scale_input_data(rbind(theta_samp[itr-1,], theta_prop), input_bounds = emulator_info$input_bounds)
      gp_pred_list <- predict_independent_GPs(X_pred = thetas_scaled, gp_obj_list = emulator_info$gp_fits, 
                                              gp_lib = emulator_info$settings$gp_lib, include_cov_mat = FALSE, denormalize_predictions = TRUE,
                                              output_stats = emulator_info$output_stats)
      SSR_samples <- sample_independent_GPs_pointwise(gp_pred_list, transformation_methods = emulator_info$settings$transformation_method, include_nugget = TRUE)
      
      # Accept-Reject step. 
      lpost_theta_curr <- calc_lpost_theta_product_lik(computer_model_data, lprior_vals = lprior_theta_curr, SSR = SSR_samples[1,,drop=FALSE], 
                                                       vars_obs = sig_eps_curr, normalize_lik = FALSE, na.rm = TRUE, return_list = FALSE)
      lpost_theta_prop_list <- calc_lpost_theta_product_lik(computer_model_data, theta_vals = matrix(theta_prop, nrow=1), SSR = SSR_samples[2,,drop=FALSE], 
                                                            vars_obs = sig_eps_curr, normalize_lik = FALSE, na.rm = TRUE, theta_prior_params = theta_prior_params, 
                                                            return_list = TRUE)
      alpha <- min(1.0, exp(lpost_theta_prop_list$lpost - lpost_theta_curr))
    
      if(runif(1) <= alpha) {
        theta_samp[itr,] <- theta_prop
        lprior_theta_curr <- lpost_theta_prop_list$lprior
        accept_count <- accept_count + 1 
        pred_idx_curr <- 2
      } else {
        theta_samp[itr,] <- theta_samp[itr-1,]
        pred_idx_curr <- 1
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
    # Gibbs step for Sig_eps.
    #
    
    if(learn_sig_eps) {
      SSR_sample <- sample_independent_GPs_pointwise(gp_pred_list, transformation_methods = emulator_info$settings$transformation_method,
                                                     idx_selector = pred_idx_curr, include_nugget = TRUE)
      sig_eps_curr <- sample_cond_post_Sig_eps(SSR = SSR_sample, Sig_eps_prior_params = sig_eps_prior_params, n_obs = computer_model_data$n_obs)
      sig_eps_samp[itr,] <- sig_eps_curr
    } else {
      sig_eps_samp[itr,] <- sig_eps_init
    }
    
  }
  
  return(list(samp = cbind(theta_samp, sig_eps_samp), Cov_prop = Cov_prop, scale_prop = exp(log_scale_prop)))
  
}


# emulator_info: list with "gp_fits", "output_stats", "settings", and "input_bounds"
# TODO: allow joint sampling, incorporating covariance between current and proposal. "output_stats" must be 
# on the correct scale (e.g. it should be on log scale for LNP).
mcmc_calibrate_ind_GP_trunc_proposal <- function(computer_model_data, theta_prior_params, emulator_info,
                                                 theta_init = NULL, sig_eps_init = NULL, learn_sig_eps = FALSE, sig_eps_prior_params = NULL, 
                                                 N_mcmc = 50000, adapt_frequency = 1000, adapt_min_scale = 0.1, accept_rate_target = 0.24, 
                                                 proposal_scale_decay = 0.7, Cov_prop_init_diag = 0.1, adapt_cov_method = "AM", 
                                                 adapt_scale_method = "MH_ratio", adapt_init_threshold = 3) {

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
  colnames(theta_samp) <- paste("theta", computer_model_data$pars_cal_names, sep = "_")
  sig_eps_samp <- matrix(nrow = N_mcmc, ncol = p)
  colnames(sig_eps_samp) <- paste("sig_eps", computer_model_data$output_vars, sep = "_")

  # Set initial conditions, ensuring the initial value of the calibration parameters satisfies the 
  # input bounds constraint determined by the extent of the design points. 
  if(is.null(theta_init)) {
    init_samp_max <- 100
    t <- 1
    while(t < init_samp_max) {
      theta_init <- sample_prior_theta(theta_prior_params)
      if(all(theta_init >= emulator_info$input_bounds[1,]) && all(theta_init <= emulator_info$input_bounds[2,])) break
      if(t == init_samp_max) stop("Failed to generate initial sample from theta prior that satisfies input bounds.")
      t <- t + 1
    }
  }
  
  if(learn_sig_eps) {
    if(is.null(sig_eps_init)) {
      sig_eps_init <- sample_prior_Sig_eps(sig_eps_prior_params)
    }
  } else {
    if(is.null(sig_eps_init)) stop("Value for `sig_eps_init` must be provided when `learn_sig_eps` is FALSE.")
  }
  
  theta_samp[1,] <- theta_init
  sig_eps_samp[1,] <- sig_eps_init
  
  sig_eps_curr <- sig_eps_init
  lprior_theta_curr <- calc_lprior_theta(theta_init, theta_prior_params)
  
  # Proposal covariance.
  Cov_prop_curr <- diag(Cov_prop_init_diag, nrow = d)
  Cov_prop <- diag(Cov_prop_init_diag, nrow = d)
  log_scale_prop <- 0
  accept_count <- 0
  samp_mean <- theta_init
  
  for(itr in seq(2, N_mcmc)) {
    
    #
    # Metropolis step for theta.
    #
    
    # theta proposals.
    theta_prop <- tmvtnorm::rtmvnorm(1, mean = theta_samp[itr-1,], sigma = exp(2*log_scale_prop)*Cov_prop, 
                                     lower = emulator_info$input_bounds[1,], upper = emulator_info$input_bounds[2,])[1,]
    # theta_prop <- theta_samp[itr-1,] + (exp(log_scale_prop) * L_prop %*% matrix(rnorm(d), ncol = 1))[,1]
    
    # Approximate SSR by sampling from GP. 
    # TODO: shouldn't have to re-scale theta_curr; should have `theta_curr_scaled` variable or something. 
    # Could also re-use mean prediction, but if including GP covariance then this will change. 
    thetas_scaled <- scale_input_data(rbind(theta_samp[itr-1,], theta_prop), input_bounds = emulator_info$input_bounds)
    gp_pred_list <- predict_independent_GPs(X_pred = thetas_scaled, gp_obj_list = emulator_info$gp_fits, 
                                            gp_lib = emulator_info$settings$gp_lib, include_cov_mat = FALSE, denormalize_predictions = TRUE,
                                            output_stats = emulator_info$output_stats)
    SSR_samples <- sample_independent_GPs_pointwise(gp_pred_list, transformation_methods = emulator_info$settings$transformation_method, include_nugget = TRUE)
    
    # Accept-Reject step. 
    lpost_theta_curr <- calc_lpost_theta_product_lik(computer_model_data, lprior_vals = lprior_theta_curr, SSR = SSR_samples[1,,drop=FALSE], 
                                                     vars_obs = sig_eps_curr, normalize_lik = FALSE, na.rm = TRUE, return_list = FALSE)
    lpost_theta_prop_list <- calc_lpost_theta_product_lik(computer_model_data, theta_vals = matrix(theta_prop, nrow=1), SSR = SSR_samples[2,,drop=FALSE], 
                                                          vars_obs = sig_eps_curr, normalize_lik = FALSE, na.rm = TRUE, theta_prior_params = theta_prior_params, 
                                                          return_list = TRUE)
    q_curr_prop <- tmvtnorm::dtmvnorm(theta_prop, mean = theta_samp[itr-1,], sigma = exp(2*log_scale_prop)*Cov_prop,
                                     lower = emulator_info$input_bounds[1,], upper = emulator_info$input_bounds[2,], log = TRUE)
    q_prop_curr <- tmvtnorm::dtmvnorm(theta_samp[itr-1,], mean = theta_prop, sigma = exp(2*log_scale_prop)*Cov_prop,
                                     lower = emulator_info$input_bounds[1,], upper = emulator_info$input_bounds[2,], log = TRUE)
    alpha <- min(1.0, exp(lpost_theta_prop_list$lpost - lpost_theta_curr + q_prop_curr - q_curr_prop))

    if(runif(1) <= alpha) {
      theta_samp[itr,] <- theta_prop
      lprior_theta_curr <- lpost_theta_prop_list$lprior
      accept_count <- accept_count + 1 
      pred_idx_curr <- 2
    } else {
      theta_samp[itr,] <- theta_samp[itr-1,]
      pred_idx_curr <- 1
    }
    
    # Adapt proposal covariance matrix and scaling term.
    adapt_list <- adapt_cov_proposal(Cov_prop_curr, log_scale_prop, theta_samp, itr, accept_count, alpha, samp_mean, adapt_frequency, 
                                     adapt_cov_method, adapt_scale_method, accept_rate_target, adapt_min_scale,  
                                     proposal_scale_decay, adapt_init_threshold, C_L = Cov_prop)
    Cov_prop_curr <- adapt_list$C
    Cov_prop <- adapt_list$C_L
    log_scale_prop <- adapt_list$log_scale
    samp_mean <- adapt_list$samp_mean
    accept_count <- adapt_list$accept_count
    
    #
    # Gibbs step for Sig_eps.
    #
    
    if(learn_sig_eps) {
      SSR_sample <- sample_independent_GPs_pointwise(gp_pred_list, transformation_methods = emulator_info$settings$transformation_method,
                                                     idx_selector = pred_idx_curr, include_nugget = TRUE)
      sig_eps_curr <- sample_cond_post_Sig_eps(SSR = SSR_sample, Sig_eps_prior_params = sig_eps_prior_params, n_obs = computer_model_data$n_obs)
      sig_eps_samp[itr,] <- sig_eps_curr
    } else {
      sig_eps_samp[itr,] <- sig_eps_init
    }
    
  }
  
  return(list(samp = cbind(theta_samp, sig_eps_samp), Cov_prop = Cov_prop, scale_prop = exp(log_scale_prop)))
  
}


adapt_cov_proposal <- function(C, log_scale, sample_history, itr, accept_count, alpha, samp_mean, 
                               adapt_frequency = 1000, cov_method = "AM", scale_method = "MH_ratio", 
                               accept_rate_target = 0.24, min_scale = 0.1, tau = 0.7, init_threshold = 3, 
                               L = NULL, C_L = NULL) {
  # Returns an adapted proposal covariance matrix. The covariance matrix is assumed to be of the form `scale * C`, where 
  # `scale` is a scaling factor. This function supports different algorithms to adapt C and to adapt `scale`, and the methods 
  # can be mixed and matched. Only a portion of the function arguments are needed for certain methods, so take care 
  # to ensure the correct arguments are being passed to the function. The function updates C and also computes the lower 
  # triangular Cholesky factor L of C. For certain methods (e.g. the AM algorithm) C will be updated every iteration, but 
  # L may only be updated intermittently depending on `adapt_cov_frequency`. The matrix `C_L` is defined as the covariance 
  # matrix that corresponds to the Cholesky factor `L`; i.e. C_L = LL^T. The covariance `C` will always be the most up-to-date, 
  # but one may instead choose to utilize `L` or `C_L` for proposals if the proposal covariance is not to be updated 
  # every MCMC iteration; this saves on computation (not having to compute a Cholesky decomposition every iteration) and 
  # can help prevent singular proposal covariances. 
  #
  # Args:
  #    C: matrix, the current d x d positive definite covariance matrix (not multiplied by the scaling factor). 
  #    log_scale: numeric, the log of the scaling factor. 
  #    L: matrix, the lower triangular Cholesky factor of C. This is what is actually used to generate proposals. In general, this 
  #       can be "out of alignment" with C; for example, C can be updated every iteration, but the Cholesky factor can be updated 
  #       less frequently to save computation time. 
  #    sample_history: matrix, itr_curr x p, where itr_curr is the current MCMC iteration. Note that this must be the entire 
  #                    sample history, even for methods like the AP algorithm, which only require a recent subset of the history. 
  #    itr_curr: integer, the current MCMC iteration. 
  #    adapt_cov_frequency: integer, number of iterations that specifies how often C or L will be updated. The AP algorithm 
  #                         updates both C and L together; the AM algorithm updates C every iteration, but only updates L according 
  #                         to `update_cov_frequency`. 
  #    adapt_scale_frequency: integer, number of iterations that specifies how often log_scale will be updated. 
  #    cov_method: character(1), currently either "AP" (Adaptive Proposal, Haario 1999) or "AM" (Adaptive Metropolis, Haario 2001). 
  #    scale_method: character(1), currently either "pecan" (method currently used in PEcAn) or "MH_ratio" (based on Metropolis-Hastings 
  #                  acceptance probability). Both methods rely on `accept_rate_target` as a target. 
  #    accept_rate_target: numeric(1), the target acceptance rate.
  #    min_scale: numeric(1), used as a floor for the scaling factor for the "pecan" method. 
  #    accept_count: integer(1), the number of accepted MH proposals over the subset of the sample history that will be used in the updates (not the 
  #                  acceptance count over the whole history!). This is only used for the "AP" and "pecan" methods.
  #    samp_mean: numeric(d), used only by "AM". This is the cumulative sample mean over the MCMC samples. It is updated every iteration. 
  #    alpha: numeric(1), used only by "MH_ratio". This is the most recent Metropolis-Hastings acceptance probability. The "MH_ratio" method 
  #           updates log_sd2 every iteration. 
  #    tau: numeric(1), used only by "MH_ratio". This is the extinction coefficient used to determine the rate at which the log_scale updating occurs. 
  #
  # Returns:
  #    list, with elements "C", "C_L", L", "log_sd2", and "samp_mean". The first three are as described above. The fourth is only used for the "AM" method 
  #    and is the mean of the MCMC samples up through the current iteration. 
  
  accept_rate <- accept_count / adapt_frequency
  
  # Adapt proposal covariance. 
  if(cov_method == "AM" ) { 
    
    samp_mean <- samp_mean + (1 / itr) * (sample_history[itr,] - samp_mean)
    
    if(itr == 3) { # Sample covariance from first 3 samples. 
      samp_centered <- t(sample_history[1:3,]) - samp_mean
      C <- 0.5 * tcrossprod(samp_centered)
    } else if(itr > 3) { # Begin covariance updates every iteration. 
      samp_centered <- sample_history[itr,] - samp_mean
      w <- 1/(itr-1)
      C <- C + w * (itr*w*outer(samp_centered, samp_centered) - C)
    }
    
    if((itr >= init_threshold) && (itr %% adapt_frequency == 0)) {
      L <- t(chol(C))
      C_L <- C
    }
  
  } else if(cov_method == "AP") {
    
    if((itr >= init_threshold) && (itr %% adapt_frequency == 0)) {
      sample_history <- sample_history[(itr - adapt_frequency + 1):itr,,drop=FALSE]
      if(accept_rate == 0) {
        C <- min_scale * C
        L <- sqrt(min_scale) * L
      } else {
        C <- stats::cov(sample_history)
        L <- t(chol(C + diag(sqrt(.Machine$double.eps), nrow=nrow(C))))
      }
      C_L <- C
    }
    
  }
  
  # Adapt scaling factor for proposal covariance matrix. 
  if(scale_method == "pecan") {
    if((itr >= init_threshold) && (itr %% adapt_frequency == 0)) {
      log_scale <- 2 * log(max(accept_rate / accept_rate_target, min_scale))
    }
  } else if(scale_method == "MH_ratio") {
    if(itr >= init_threshold) log_scale <- log_scale + (1 / itr^tau) * (alpha - accept_rate_target) 
  }
  
  if((itr >= init_threshold) && (itr %% adapt_frequency == 0)) accept_count <- 0
  
  return(list(C = C, C_L = C_L, L = L, log_scale = log_scale, samp_mean = samp_mean, accept_count = accept_count))
  
}


sample_prior_theta <- function(theta_prior_params) {
  # Return sample from prior distribution on the calibration parameters (theta). 
  #
  # Args:
  #    theta_prior_params: data.frame, with columns "dist", "param1", and "param2"; and potentially also columns 
  #                        "bound_lower" and "bound_upper". The ith row of the data.frame
  #                        should correspond to the ith entry of 'theta'. Currently, accepted values of "dist" are 
  #                        "Gaussian" (param1 = mean, param2 = std dev) and "Uniform" (param1 = lower, param2 = upper), 
  #                        and "Truncatd_Gaussian" (param1 = mean, param2 = std dev, bound_lower = lower truncation value, 
  #                        bound_upper = upper truncation value).  
  #
  # Returns:
  #    numeric vector of length equal to number of rows of `theta_prior_params`, the prior sample. 

  theta_samp <- vector(mode = "numeric", length = nrow(theta_prior_params))
  
  for(i in seq_along(theta_samp)) {
    if(theta_prior_params[i, "dist"] == "Gaussian") {
      theta_samp[i] <- rnorm(1, theta_prior_params[i, "param1"], theta_prior_params[i, "param2"])
    } else if(theta_prior_params[i, "dist"] == "Uniform") {
      theta_samp[i] <- runif(1, theta_prior_params[i, "param1"], theta_prior_params[i, "param2"])
    } else if(theta_prior_params[i, "dist"] == "Truncated_Gaussian") {
      theta_samp[i] <- rtruncnorm(1, a = theta_prior_params[i, "bound_lower"], b = theta_prior_params[i, "bound_upper"], 
                                  mean = theta_prior_params[i, "param1"], sd = theta_prior_params[i, "param2"])
    } else {
      stop("Prior distribution ", theta_prior_params[i, "dist"], " not supported.")
    }
  }

  return(theta_samp)  
  
}


truncate_prior_theta <- function(theta_prior_params, input_bounds) {
  # Converts the prior parameters on the calibration parameters (theta) so that they are truncated
  # in that they assign 0 probability mass beyond the bounds specified in `input_bounds`. 
  # `input_bounds` is typically determined by the extent of the design points, so this
  # function modifies the prior so that posterior evaluations are 0 outside of the 
  # extent of the design points (where the GP would have to interpolate). This is an 
  # alternative to allowing an unbounded prior and instead truncating the MCMC 
  # proposals. The truncation applied to uniform priors simply sets the bounds on the 
  # uniform prior to the bounds provided in `input_bounds`. Applied to Gaussian priors, 
  # the Gaussian distributions are converted to truncated Gaussian distributions, with 
  # the truncation bounds again determined by `input_bounds`. For compactly supported priors 
  # (e.g. uniform or truncated Gaussian), the bounds are only updated if they fall outside of 
  # the bounds in `input_bounds` (e.g. an existing lower bound will not be made any lower). 
  #
  # Args:
  #    theta_prior_params: data.frame, with columns "dist", "param1", and "param2". The ith row of the data.frame
  #                        should correspond to the ith entry of 'theta'. Currently, accepted values of "dist" are 
  #                        "Gaussian" (param1 = mean, param2 = std dev) and "Uniform" (param1 = lower, param2 = upper).
  #    input_bounds: matrix of dimension 2 x d, where d is the number of parameters. The first row contains lower bounds 
  #                  on the parameters and the second row contains upper bounds (these are often determined by the 
  #                  extent of the design points). If provided, the prior sample will be required to satisfy the lower 
  #                  and upper bounds, so sampled parameters that do not satisfy the constraints will be rejected until 
  #                  the constraint is satisfied. Note that this implies that the resulting samples are from a truncated 
  #                  version of the prior distribution, rather than the prior itself. The columns of `input_bounds` must be 
  #                  sorted in the same order as the rows of `theta_prior_params.` Default is NULL, which imposes 
  #                  no constraints. 
  #
  # Returns:
  #    data.frame, the updated version of `theta_prior_params`. 
  
  for(i in seq(1, nrow(theta_prior_params))) {
    l <- input_bounds[1, i]
    u <- input_bounds[2, i]
    
    if(theta_prior_params[i, "dist"] == "Gaussian") {
      theta_prior_params[i, "dist"] <- "Truncated_Gaussian"
      theta_prior_params[i, "bound_lower"] <- l
      theta_prior_params[i, "bound_upper"] <- u
    } else if(theta_prior_params[i, "dist"] == "Uniform") {
      if(theta_prior_params[i, "param1"] < l) theta_prior_params[i, "param1"] <- l
      if(theta_prior_params[i, "param2"] > u) theta_prior_params[i, "param2"] <- u
    } else if(theta_prior_params[i, "dist"] == "Truncated_Gaussian") { 
      if(theta_prior_params[i, "bound_lower"] < l) theta_prior_params[i, "bound_lower"] <- l
      if(theta_prior_params[i, "bound_upper"] > u) theta_prior_params[i, "bound_upper"] <- u
    } else {
      stop("Prior distribution ", theta_prior_params[i, "dist"], " not supported.")
    }
  }
  
  return(theta_prior_params)
  
}


sample_prior_Sig_eps <- function(Sig_eps_prior_params, return_matrix = FALSE) {
  # Returns sample from prior distribution on the p x p observation covariance matrix Sig_eps.
  #
  # Args:
  #    Sig_eps_prior_params: list, must contain element named "dist" specifying the prior distribution. Current accepted values are 
  #                          "IW" (Inverse Wishart) or "IG" (independent Inverse Gamma priors). Depending on value of "dist", 
  #                          must also contain either either 1.) names "scale_matrix" and "dof" 
  #                          that are the arguments of the Inverse Wishart distribution, or 2.) names "IG_shape" and "IG_scale" which 
  #                          each correspond to p-length vectors storing the parameters for the independent Inverse Gamma priors on each
  #                          variance parameter. Note that dist "IG" constrains Sig_eps to be diagonal, while "IW" does not. 
  #    return_matrix: logical(1), only relevant for inverse Gamma prior. When prior is inverse Wishart, the return value is always 
  #                   a matrix. If `return_matrix` is TRUE, then the return value under the inverse Gamma prior will also be a matrix. 
  #                   Otherwise it is a numeric vector. 
  #
  # Returns:
  #    matrix or numeric vector. In prior is inverse Wishart, returns p x p positive definite matrix sampled from the prior. If prior 
  #    instead returns numeric vector of length p containing the sampled variances. Setting `return_matrix` to TRUE will force the 
  #    return value to be a matrix in either case (this will return a diagonal matrix in the inverse gamma case). 
  
  if(Sig_eps_prior_params$dist == "IW") {
    return(LaplacesDemon::rinvwishart(nu = Sig_eps_prior_params$dof, S = Sig_eps_prior_params$scale_matrix))
  } else if(Sig_eps_prior_params$dist == "IG") {
    p <- length(Sig_eps_prior_params$IG_shape)
    sig2_eps_vars <- vector(mode = "numeric", length = p)
    for(j in seq_len(p)) {
      sig2_eps_vars[j] <- 1/rgamma(1, shape = Sig_eps_prior_params$IG_shape[j], rate = Sig_eps_prior_params$IG_scale[j])
    }
    
    if(return_matrix) return(diag(sig2_eps_vars, nrow = p))
    return(sig2_eps_vars)
    
  }
  
}


sample_cond_post_Sig_eps <- function(model_errs = NULL, SSR = NULL, Sig_eps_prior_params, n_obs, return_matrix = FALSE) {
  # Return sample of the observation covariance matrix Sig_eps, drawn from the conditional  
  # posterior p(Sig_eps|theta, Y). Under the model assumptions, this conditional posterior 
  # has an inverse Wishart or Inverse Gamma product distribution. This function does not explicitly take theta as an 
  # argument; rather, it assumes the forward model f(theta) has already been run and 
  # the error matrix Y - f(theta) is provided by the argument `model_errs`, or alternatively the sum of squared errors 
  # are provided in the Inverse Gamma case. 
  #
  # Args:
  #    model_errs: matrix of dimensions n x p, where n = number observations in time series and p = number output variables.
  #                This is the model error matrix Y - f(theta). Can be null for inverse gamma prior, in which case `SSR` must be provided. 
  #    SSR: numeric(), vector of length equal to the number of outputs p, containing the squared L2 errors for each output. This is only 
  #         used in the case of inverse gamma prior. If NULL, `model_errs` must be provided instead. Alternatively, this can be a 
  #         matrix of dimension 1 x p. 
  #    Sig_eps_prior_params: list, must contain element named "dist" specifying the prior distribution. Current accepted values are 
  #                          "IW" (Inverse Wishart) or "IG" (independent Inverse Gamma priors). Depending on value of "dist", 
  #                          must also contain either either 1.) names "scale_matrix" and "dof" 
  #                          that are the arguments of the Inverse Wishart distribution, or 2.) names "IG_shape" and "IG_scale" which 
  #                          each correspond to p-length vectors storing the parameters for the independent Inverse Gamma priors on each
  #                          variance parameter. Note that dist "IG" constrains Sig_eps to be diagonal, while "IW" does not. 
  #    n_obs: numeric(), vector of length equal to the number of outputs p. The number of observations for each output. 
  #    return_matrix: logical(1), only relevant for inverse Gamma prior. When prior is inverse Wishart, the return value is always 
  #                   a matrix. If `return_matrix` is TRUE, then the return value under the inverse Gamma prior will also be a matrix. 
  #                   Otherwise it is a numeric vector. 
  #
  # Returns:
  #    If prior is inverse wishart, returns a p x p positive definite matrix sampled from the conditional posterior distribution p(Sig_eps|theta, Y).
  #    If prior is inverse Gamma, returns a numeric vector of length p of the sampled variances from the conditional posterior. Setting `return_matrix`
  #    to TRUE forces the return value to be a matrix in this case as well. 
  
  if(Sig_eps_prior_params$dist == "IW") {
    stop("Need to update IW prior for case where there are missing observations.")
    inv_wishart_scale <- crossprod(model_errs_curr) + Sig_eps_prior_params$scale_matrix
    inv_wishart_dof <- n + Sig_eps_prior_params$dof
    return(LaplacesDemon::rinvwishart(nu = inv_wisharat_dof, S = inv_wishart_scale))
  } else if(Sig_eps_prior_params$dist == "IG") {

    if(is.null(SSR)) {
      SSR <- colSums(model_errs^2)
    }
    
    p <- length(n_obs)
    sig2_eps_vars <- vector(mode = "numeric", length = p)
    for(j in seq_len(p)) {
      sig2_eps_vars[j] <- 1/rgamma(1, shape = 0.5*n_obs[j] + Sig_eps_prior_params$IG_shape[j], rate = 0.5*SSR[j] + Sig_eps_prior_params$IG_scale[j])
    }
    
    if(return_matrix) return(diag(sig2_eps_vars, nrow = p))
    return(sig2_eps_vars)
    
  }
  
}


generate_linear_Gaussian_test_data <- function(random_seed, N_obs, D, Sig_theta, G, sig2_eps = NULL, sig_eps_frac = 0.1, pars_cal_sel = NULL) {
  # Sets up a test example (including computer model data and prior distributions) in which the forward model is linear 
  # (represented by matrix `G`), the likelihood is Gaussian (with variance `sig2_eps`), and the prior on the calibration 
  # parameters is zero-mean Gaussian (with diagonal covariance matrix `Sig_theta`). Note that currently this function 
  # only creates an example with a single output variable (p = 1). Also, for the time being `Sig_theta` must be diagonal, 
  # until the prior code is updated to allow for correlated priors. This linear Gaussian setup admits a closed form 
  # posterior so is useful for validating MCMC schemes, etc. This function also returns the mean and covariance matrix 
  # of the true posterior. 
  #
  # Args:
  #    random_seed: integer(1), the seed for the random number generator. 
  #    N_obs: integer(1), the number of observed data points that will be generated. 
  #    D: integer(1), the dimension of the input parameter space; if `pars_cal_sel` is NULL or corresponds to all parameters, 
  #       then D is the dimension of the calibration parameter space. However, if `pars_cal_sel` only selects a subset of parameters, 
  #       then D will be larger than the dimension of the parameter calibration space. 
  #    Sig_theta: matrix of dimension D x D. Note that if only a subset of parameters are calibrated, then the prior covariance on the 
  #               calibration parameters with be a sub-matrix of `Sig_theta`.
  #    sig2_eps: numeric(1), the observation variance. If not provided, then the obsevation variance will be set using `coef_var`, or 
  #              if `coef_var`.
  #    sig_eps_frac: numeric(1), If `sig2_eps` is provided directly then `coef_var` will not be used.
  #                  Otherwise, the noise variance  will be set so that the standard deviation sqrt(sig2_eps) equals 
  #                  the range of the data * `sig_eps_frac`. 
  #    G: matrix of dimension N_obs x D, the linear forward model.
  #    pars_cal_sel: integer(), vector containing the indices used to select the parameters which will be calibrated. The remaining parameters
  #                  will be fixed. If NULL, calibrates all parameters. 
  #    
  # Returns:
  #    list with 3 elements: computer_model_data, theta_prior_params, and true_posterior. 
  
  if(!all.equal(dim(G), c(N_obs, D))) {
    stop("Forward model G must be matrix of dimension N_obs x D.")
  }
  
  if(!isTRUE(all.equal(Sig_theta, diag(diag(Sig_theta))))) {
    stop("Code does not currently support correlated prior parameters. `Sig_theta` should be diagonal.")
  }
  
  # Sample from model to generate observed data. 
  L_theta <- t(chol(Sig_theta))
  theta <- L_theta %*% matrix(rnorm(D), ncol=1)
  data_ref <- G %*% theta

  set.seed(random_seed)
  if(is.null(sig2_eps)) {
    sig2_eps <- (diff(range(data_ref)) * sig_eps_frac)^2
  }
  Sig_t <- diag(sig2_eps, nrow = N_obs)
  L_t <- t(chol(Sig_t))
  data_obs <- data_ref + L_t %*% matrix(rnorm(N_obs), ncol=1)
  output_vars <- "y"
  colnames(data_obs) <- output_vars
  
  # Select parameters to calibrate. 
  if(is.null(pars_cal_sel)) pars_cal_sel <- seq(1,D)
  theta_names <- paste0("theta", seq(1,D))
  pars_cal_names <- theta_names[pars_cal_sel]
  if(length(pars_cal_names) > N_obs) {
    stop("For now number of calibration parameters must be <= number of observations.")
  }
  theta_true <- theta[pars_cal_sel]
  names(theta_true) <- pars_cal_names
  
  # Forward map. 
  f <- function(par_val, computer_model_data) {
    theta <- computer_model_data$ref_pars$true_value
    theta[computer_model_data$pars_cal_sel] <- par_val
    
    return(computer_model_data$G %*% matrix(theta, ncol=1))
  }
  
  # Computer model data.
  computer_model_data <- list(f = f, 
                              data_ref = data_ref,
                              data_obs = data_obs, 
                              theta_true = theta_true,
                              n_obs = N_obs, 
                              output_vars = output_vars, 
                              pars_cal_names = pars_cal_names,
                              pars_cal_sel = pars_cal_sel,
                              forward_model = "custom_likelihood", 
                              G = G, 
                              Sig_eps = matrix(sig2_eps), 
                              ref_pars = data.frame(true_value = theta, 
                                                    row.names = pars_cal_names))
  
  # Prior Parameters. 
  theta_prior_params <- data.frame(dist = rep("Gaussian", length(pars_cal_names)), 
                                   param1 = rep(0, length(pars_cal_names)), 
                                   param2 = sqrt(diag(Sig_theta)[pars_cal_sel]))
  rownames(theta_prior_params) <- pars_cal_names
  
  # True posterior (note that we need to adjust for the case where only a subset of the parameters 
  # are calibrated). The posterior moments are computed using the SVD and Woodbury identity. This 
  # assumes the number of calibration parameters is <= N_obs. 
  pars_fixed_sel <- setdiff(seq(1,D), pars_cal_sel)
  G_cal <- G[, pars_cal_sel]
  theta_cal <- theta[pars_cal_sel]
  Sig_theta_cal <- Sig_theta[pars_cal_sel, pars_cal_sel]
  
  if(length(pars_fixed_sel) == 0) {
    G_fixed <- matrix(0, nrow = N_obs, ncol = 1)
    theta_fixed <- 0
  } else {
    G_fixed <- G[, pars_fixed_sel]
    theta_fixed <- theta[pars_fixed_sel]
  }
  y_adjusted <- matrix(data_obs, ncol=1) - G_fixed %*% theta_fixed
  
  # TODO: check SVD calculation. For now just directly taking inverse. 
  # svd_list <- svd(G_cal)
  # Cov_post <- Sig_theta_cal - Sig_theta_cal %*% diag(1 / (sig2_eps * (svd_list$d^(-2)) + diag(Sig_theta_cal))) %*% Sig_theta_cal
  
  Cov_post <- sig2_eps * solve(crossprod(G_cal) + sig2_eps * diag(1/diag(Sig_theta_cal)))
  mean_post <- (1/sig2_eps) * tcrossprod(Cov_post, G) %*% y_adjusted
  true_posterior <- list(mean = mean_post, Cov = Cov_post)
  
  return(list(computer_model_data = computer_model_data, theta_prior_params = theta_prior_params, 
              true_posterior = true_posterior))
  
}


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
  theta_true <- ref_pars[pars_cal_sel, "true_value"]
  names(theta_true) <- pars_cal_names
  
  # Run the model to generate the reference data, the ground truth. 
  data_ref <- run_VSEM(theta_vals = pars_cal_vals, ref_pars = ref_pars, pars_cal_sel = pars_cal_sel, 
                       PAR_data = PAR, output_vars = output_vars)

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
              obs_start_day = obs_start_day, 
              forward_model = "VSEM", 
              f = run_VSEM_single_input, 
              theta_true = theta_true))
  
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
  # A convenience function to generate the VSEM data for "test case 4". This test is 
  # identical to test case 3, with the addition of a second parameter `KEXT` (the 
  # extinction coefficient in the GPP calculation). KEXT and LAR are not identifiable 
  # in the ODE; they completely trade off. Thus, this test is intended to explore 
  # calibration behavior in an extreme case of non-identifiability. 
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
  pars_cal_names <- c("LAR", "KEXT")
  pars_cal_vals <- c(1.5, 0.5)
  ref_pars <- VSEMgetDefaults()
  
  # NEE is observed daily with no missing values. LAI is observed every three days.   
  output_vars <- c("NEE", "LAI")
  output_frequencies <- c(1, 3)
  obs_start_day <- c(1, 1)
  
  test_list <- generate_vsem_test_data(random_seed = random_seed,
                                       N_time_steps = N_time_steps , 
                                       Sig_eps = Sig_eps, 
                                       pars_cal_names = pars_cal_names, 
                                       pars_cal_vals = pars_cal_vals,
                                       ref_pars = ref_pars, 
                                       output_vars = output_vars, 
                                       output_frequencies = output_frequencies, 
                                       obs_start_day = obs_start_day)
  return(test_list)
  
}


generate_vsem_test_5 <- function(random_seed) {
  # A convenience function to generate the VSEM data for "test case 5". This test is 
  # identical to test case 3, with the addition of a second parameter `Cs` (the 
  # initial condition for the soil pool). This is also a counterpart ot test case 4, 
  # which calibrates LAR and KEXT, which are completely unidentifiable. In this case
  # the parameters LAR and Cs are less related, so comparing test cases 4 and 5 allow 
  # for different degrees of non-identifiability to be tested. 
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
  pars_cal_names <- c("LAR", "Cs")
  pars_cal_vals <- c(1.5, 15.0)
  ref_pars <- VSEMgetDefaults()
  
  # NEE is observed daily with no missing values. LAI is observed every three days.   
  output_vars <- c("NEE", "LAI")
  output_frequencies <- c(1, 3)
  obs_start_day <- c(1, 1)
  
  test_list <- generate_vsem_test_data(random_seed, N_time_steps, Sig_eps, pars_cal_names, pars_cal_vals,
                                       ref_pars, output_vars, output_frequencies, obs_start_day)
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



# ------------------------------------------------------------------------------
# MCMC Plotting Functions. 
# ------------------------------------------------------------------------------


get_hist_plot <- function(samples_list, col_sel = 1, bins = 30, vertical_line = NULL, xlab = "samples", ylab = "density", 
                          main_title = "Histogram", data_names = NULL) {
  # Generates a single plot with one or more histograms. The input data `samples_list` must be a list of matrices. 
  # For example, the first element of the list could be a matrix of samples from a reference distribution, and the 
  # second element a matrix of samples from an approximate distribution. This function will select one column from 
  # each matrix based on the index `col_sel`. Each selected column will be turned into a histogram. 
  # Note that the matrices may have different lengths (numbers of samples).
  #
  # Args:
  #    samples_list: list of matrices, each of dimension (num samples, num parameters). 
  #    col_sel: integer(1), the column index to select from each matrix. i.e. this selects a particular parameter, whose
  #             samples will be transformed into a histogram. 
  #    bins: integer(1), number of bins to use in histogram. 
  #    vertical_line: If not NULL, the x-intercept to be used for a plotted vertical line. 
  #    xlab, ylab, main_title: x and y axis labels and plot title. 
  #    data_names: character(), vector equal to the length of the list `samples_list`. These are the names used in 
  #                the legend to identify the different histograms. 
  #
  # Returns:
  #    ggplot2 object. 
  
  # Pad with NAs in the case that the number of samples is different across different parameters. 
  N_samp <- sapply(samples_list, nrow)
  if(length(unique(N_samp)) > 1) {
    N_max <- max(N_samp)
    for(j in seq_along(samples_list)) {
      if(nrow(samples_list[[j]]) < N_max) {
        samples_list[[j]] <- rbind(samples_list[[j]], matrix(NA, nrow = N_max - nrow(samples_list[[j]]), ncol = ncol(samples_list[[j]])))
      }
    }
  }
  
  dt <- as.data.table(lapply(samples_list, function(mat) mat[,col_sel]))
  if(!is.null(data_names)) setnames(dt, colnames(dt), data_names)
  dt <- melt(dt, measure.vars = colnames(dt), na.rm = TRUE)
  
  plt <- ggplot(data = dt, aes(x = value, color = variable)) + 
          geom_histogram(aes(y = ..density..), bins = bins, fill = "white", alpha = 0.2, position = "identity") + 
          xlab(xlab) + 
          ylab(ylab) + 
          ggtitle(main_title)
  
  if(!is.null(vertical_line)) {
    plt <- plt + geom_vline(xintercept = vertical_line, color = "red")
  }
  
  return(plt)
  
}


get_2d_density_contour_plot <- function(samples_list, col_sel = c(1,2), xlab = "theta1", ylab = "theta2", main_titles = NULL) {
  # Plots the contours of a 2D kernel density estimate. If `samples_list` contains multiple elements, then one plot will be 
  # returned per element. Each element of `samples_list` is matrix of dimension N_samples x N_params. The vector `col_sel`
  # determines which 2 columns will be used for the 2D KDE plot in each matrix. 
  #
  # Args:
  #    samples_list: list of matrices, each of dimension (num samples, 2). 
  #    col_sel: integer(1), the two column indees to select from each matrix. i.e. this selects two particular parameters, whose
  #             samples will be used to compute the KDE. 
  #    xlab, ylab, main_title: x and y axis labels and plot title. 
  #
  # Returns:
  #    list, of equal length to `samples_list`. Each element is a ggplot object. 
                                        
  if(is.null(main_titles)) main_titles <- paste0("2D KDE Countours: ", seq_along(samples_list))
  plts <- vector(mode = "list", length = length(samples_list))
  
  for(j in seq_along(samples_list)) {
    df <- data.frame(samples_list[[j]])
    x <- colnames(df)[col_sel[1]]
    y <- colnames(df)[col_sel[2]]
    plts[[j]] <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) + 
                  geom_density_2d_filled() + 
                  xlab(xlab) + 
                  ylab(ylab) + 
                  ggtitle(main_titles[[j]])
  }
  
  return(plts)
  
}


get_2d_heatmap_plot <- function(X, y, param_names, samples_kde = NULL, samples_points = NULL,  
                                raster = FALSE, point_coords = NULL, main_title = "Heatmap") {
  # Plots a 2d heatmap of a scalar quantity `y`. Optionally overlays contours of a 2d 
  # kernel density estimate from `samples`. The input locations are given by the 
  # M x 2 matrix `X`. If these input locations correspond to an evenly-spaced grid, 
  # then `raster = TRUE` may be set to create a classic heatmap produced over a dense 
  # grid. Altenatively, if `X` consists of more sparsely sampled or non-evenly-spaced 
  # locations, `raster = FALSE` should be set and the resulting plot will plot the 
  # individual points, which will still be colored in heatmap fashion. 
  #
  # Args:
  #    X: matrix, of dimension M x D with colnames specified. `param_names` will be used to 
  #       select the two input dimensions of `X`. These define the input locations used in the 
  #       heatmap. 
  #    y: numeric(M), the scalar output value used to determine the colors in the heatmap. 
  #    samples_kde: matrix, with colnames specified. `param_names` will be used to 
  #                 select the two input dimensions of `samples_kde`. These input points 
  #                 will not be part of the heatmap. Instead, they will be used to construct 
  #                 a 2D KDE estimate and the contours of this KDE will be overlaid on the 
  #                 heatmap. 
  #    samples_points: matrix, with colnames specified. `param_names` will be used to 
  #                    select the two input dimensions of `samples_points`. These input points 
  #                    will be directly plotted as points on the plot and colored red. This 
  #                    argument is typically used to plot design points. 
  #    raster: logical(1), see above description. Set to TRUE when `X` is a dense grid of evenly-spaced  
  #            points FALSE when the points in `X` are not evenly-spaced or are sparse. 
  #    point_coords: numeric(2), coordinates to plot a single point as a red triangle. This typically 
  #                  corresponds to the location of the true parameter value. 
  #    main_title: character(1), the title of the plot. 
  
  if(length(param_names) != 2) stop("<param_names> must have length 2.")
  
  df <- as.data.frame(cbind(X[, param_names], y))
  colnames(df) <- c("theta1", "theta2", "y")
  
  # Heatmap. 
  if(raster) {
    plt <- ggplot() + 
            geom_tile(data = df, aes(x = theta1, y = theta2, fill = y)) + 
            scale_fill_viridis(discrete=FALSE)
  } else {
    plt <- ggplot() + 
            geom_point(data = df, aes(x = theta1, y = theta2, color = y)) + 
            scale_color_viridis(discrete=FALSE)
  }
  
  # Title and axis labels. 
  plt <- plt + ggtitle(main_title) + xlab(param_names[1]) + ylab(param_names[2])
  
  # Density contours from samples. 
  if(!is.null(samples_kde)) {
    if(!all(param_names %in% colnames(samples_kde))) stop("<param_names> must be column names of <samples_kde>.")
  
    samples_kde <- as.data.frame(samples_kde[, param_names])
    colnames(samples_kde) <- c("theta1", "theta2")
    
    plt <- plt + geom_density_2d(data = samples_kde, mapping = aes(x = theta1, y = theta2))
  }
  
  # Plot points. 
  if(!is.null(samples_points)) {
    if(!all(param_names %in% colnames(samples_points))) stop("<param_names> must be column names of <samples_points>.")
    samples_points <- as.data.frame(samples_points[, param_names])
    colnames(samples_points) <- c("theta1", "theta2")
    
    plt <- plt + geom_point(data = samples_points, mapping = aes(x = theta1, y = theta2), color = "red")
  }
  
  # Mark specific point in plot. 
  if(!is.null(point_coords)) {
    plt <- plt + geom_point(data = data.frame(theta1 = point_coords[1], theta2 = point_coords[2]), 
                            aes(x = theta1, y = theta2), color = "red", shape = 17)
  }
  
  return(plt)
  
}


get_2d_response_surface_plot <- function(computer_model_data, theta_vals, param_names, response_surface, 
                                         SSR_vals = NULL, theta_prior_params = NULL, output_variables = NULL, 
                                         combine_outputs = TRUE, raster = FALSE, point_coords = NULL, 
                                         samples_kde = NULL, samples_points = NULL, scale_inputs = FALSE, 
                                         input_bounds = NULL, ...) {
  # A wrapper for `get_2d_heatmap_plot()` that specializes in plotting quantities of interest 
  # related to the calibration problem. This function can be used when the dimension of the 
  # calibration parameter space is 2, or to plot 2-dimensional projections when the input 
  # dimension is larger than 2. This function can plot a heatmap 
  # of 1.) the SSR (sum of squared residual) surface, 2.) the likelihood surface,  
  # 3.) the posterior surface, or 4.) the prior surface. It also allows the user to control 
  # whether individual plots are generated for each output variable, or if all the output variables 
  # should  be combined to plot the overall likelihood, posterior, etc. It also allows points 
  # to be overlaid on the plot, which typically represent design points or samples 
  # from a ground truth distribution. Additionally a single point can be added to 
  # mark the true value of the parameters. 
  #
  # Args:
  #    response_surface: character(1), either "SSR", "likelihood", "prior", or "posterior". 
  #    ...: Other arguments passed to ggplot()
  #
  #
  
  theta_vals_unscaled <- as.matrix(theta_vals)
  if(scale_inputs) {
    theta_vals <- scale_input_data(theta_vals, input_bounds)
    if(!is.null(point_coords)) point_coords <- scale_input_data(matrix(point_coords, nrow=1), input_bounds)
    if(!is.null(samples_kde)) samples_kde <- scale_input_data(samples_kde, input_bounds)
    if(!is.null(samples_points)) samples_points <- scale_input_data(samples_points, input_bounds)
  } 
  
  if(response_surface == "SSR") {
    
    if(is.null(SSR_vals)) SSR_vals <- get_computer_model_SSR(computer_model_data, theta_vals = theta_vals_unscaled)
    
    plts <- get_2d_response_surface_plot_SSR(theta_vals = theta_vals, 
                                             SSR_vals = SSR_vals, 
                                             param_names = param_names, 
                                             output_variables = output_variables, 
                                             raster = raster, 
                                             point_coords = point_coords, 
                                             samples_kde = samples_kde, 
                                             samples_points = samples_points)
  }
  
  return(plts)
  
}


get_2d_response_surface_plot_SSR <- function(theta_vals, SSR_vals, param_names, output_variables = NULL, 
                                             raster = FALSE, point_coords = NULL, samples_kde = NULL,
                                             samples_points = NULL) {
  
  if(!is.null(output_variables)) {
    SSR_vals <- SSR_vals[, output_variables]
  } else {
    output_variables <- colnames(SSR_vals)
  }
  
  plts <- vector(mode = "list", length = length(output_variables))
  for(j in seq_along(plts)) {
    plts[[j]] <- get_2d_heatmap_plot(X = theta_vals, 
                                     y = SSR_vals[,j], 
                                     param_names = param_names,
                                     samples_kde = samples_kde, 
                                     samples_points = samples_points, 
                                     raster = raster, 
                                     main_title = paste0("SSR: ", output_variables[j]), 
                                     point_coords = point_coords)
  }
  
  return(plts)
  
}
                                  

get_trace_plots <- function(samp_df, burn_in_start = 1, ...) {
  # See `select_mcmc_samp_cols()` for details on how to specify columns which to plot. 

  # Select columns to plot.  
  n <- nrow(samp_df)
  cols_sel <- select_mcmc_samp_cols(...)
  samp_df_plot <- samp_df %>% select(matches(cols_sel))
  col_names <- colnames(samp_df_plot)
  
  # Get starting iteration number for each parameter. 
  start_itrs <- get_mcmc_burn_in_start_itrs(burn_in_start, colnames(samp_df_plot))
  
  # Generate one trace plot per column. 
  plts <- vector(mode = "list", length = length(col_names))
  
  for(j in seq_along(plts)) {
    y <- col_names[j]
    itr_start <- start_itrs[j]
    df <- samp_df_plot[itr_start:n, y, drop = FALSE]
    colnames(df) <- "param"
    df$itr <- itr_start:n
    
    plts[[j]] <- ggplot(data = df, aes(x = itr, y = param)) + 
                 geom_line() + 
                 ggtitle(paste0("Trace Plot: ", y)) + 
                 xlab("Iteration")
  }
  
  return(plts)
  
}


get_mcmc_marginal_hist_plot <- function(samp_df, param_names, burn_in_start = 1, bins = 30, vertical_lines = NULL, ...) {
  # Generates one plot per parameter name. 
  
  # Select columns to plot.  
  n <- nrow(samp_df)
  cols_sel <- select_mcmc_samp_cols(param_names = param_names, ...)
  samp_df_plot <- samp_df %>% select(matches(cols_sel))
  col_names <- colnames(samp_df_plot)
  
  # Get starting iteration number for each parameter. 
  start_itrs <- get_mcmc_burn_in_start_itrs(burn_in_start, colnames(samp_df_plot))
  
  # Convert data.frame to long format. 
  samp_df_plot <- melt(as.data.table(samp_df_plot), measure.vars = colnames(samp_df_plot), na.rm = TRUE)
  
  # Produce one plot per parameter name. 
  plts <- vector(mode = "list", length = length(param_names))
  for(j in seq_along(plts)) {
    df_param_cols <- grep(param_names[j], col_names, value = TRUE)
    plts[[j]] <- ggplot(data = samp_df_plot[variable %in% df_param_cols], aes(x = value, color = variable)) + 
                        geom_histogram(aes(y = ..density..), bins = bins, fill = "white", alpha = 0.2, position = "identity") + 
                        xlab(param_names[j]) + 
                        ylab("Frequency") + 
                        ggtitle(paste0("Marginal Distribution: ", param_names[j]))
    
    if(!is.null(vertical_lines)) {
      plts[[j]] <- plts[[j]] + geom_vline(xintercept = vertical_lines[param_names[j]], color = "red")
    }
    
  }
  
  return(plts)
  
}


# get_2d_gp_pred_heatmap_plots <- function(emulator_info_list, X_test, raster = TRUE, emulator_pred_list = NULL) {
#   
#   emulator_pred_list <- vector(mode = "list", length = length(emulator_info_list))
#   for(j in seq_along(emulator_info_list)) {
#     predict_independent_GPs(X_pred = X_test, 
#                             gp_obj_list = emulator_info_list[[j]]$gp_fits, 
#                             gp_lib = emulator_info_list[[j]]$settings$gp_lib, 
#                             include_cov_mat = FALSE, 
#                             denormalize_predictions = TRUE, 
#                             output_stats = emulator_info_list[[j]]$output_stats, 
#                             transformation_method = emulator_info_list[[j]]$settings$transformation_method)
#   }
#                             
#   
# }
# 
# 
# get_2d_gp_heatmap_plot <- function() {
#   
# }


select_mcmc_samp_cols <- function(test_labels = NULL, param_types = NULL, param_names = NULL, col_names = NULL, ...) {
  # This function generates a vector of column names intended to select columns from a data.frame, in which each 
  # column of the data.frame stores a set of samples from a particular parameter. The assumed column name structure 
  # is "<test_label>_<parameter_label>_<parameter_types>", where <test_label> is used to identify different 
  # MCMC schemes (e.g. for comparing different algorithms), <parameter_types> is currently either "theta" or 
  # "sig_eps", and <param_names> provides the specific name for identifying a parameter. Specifying these values 
  # in the arguments `test_labels`, `param_types`, `params_names` will concatenate every combination and return 
  # the resulting column names. Specifying `col_names` overrides all other specifications, and simply returns 
  # `col_names` itself; this is used when the user wants to manually select certain columns. Specifying only one 
  # or two of `test_labels`, `param_types`, `params_names` will generate regular expressions to select all 
  # of the variable(s) that are missing. For example, one could specify only `test_labels` = c("test1", "test2)
  # and `param_types` = c("theta"), which would then generate a character vector used to select all columns 
  # corresponding to tests 1 and 2 with parameter theta (and will include all theta parameters). Calling 
  # this function will all arguments equal to NULL will return "*_*_*", which is a regular expression that 
  # will select all columns under the assumed column name format. 
  #
  # Args:
  #    See above description. 
  #
  # Returns:
  #    character vector, containing either column names or regular expressions to select column names. 
  
  if(!is.null(col_names)) {
    return(col_names)
  } else {
    
    # Deal with case that one or two of the inputs are missing. 
    if(is.null(test_labels)) test_labels <- "*"
    if(is.null(param_types)) param_types <- "*"
    if(is.null(param_names)) param_names <- "*"
    
    # Generate column selector. 
    col_combs <- paste(rep(test_labels, each = length(param_types)), param_types, sep = "_")
    col_sel <- paste(rep(col_combs, each = length(param_names)), param_names, sep = "_")
    return(col_sel)
    
  } 

}


get_mcmc_burn_in_start_itrs <- function(burn_in_start, samp_df_col_names) {
  # Returns an integer vector of equal length as `samp_df_col_names` which contains 
  # starting MCMC iteration for each respective variable in `samp_df_col_names` (all 
  # iterations prior to the start iteration are considered burn-in, i.e. warm-up). 
  # The logic for `burn_in_start` is as follows. If it is an unnamed vector of length 1  
  # then it is assumed this is the starting iteration that is desired for all plots. 
  # If it is an unnamed vector of multiple elements, then it is assumed that its length 
  # and order corresponds to the plots that will be produced, in which case the respective 
  # iterations will be used for each plot. An error will be thrown if the length does 
  # not agree with the number of plots. Otherwise, `burn_in_start` must be a named 
  # integer vector with the names equal to the values of the `test_labels`. `burn_in_start`
  # is assumed to only be a function of the test labels, as typically the same burn-in is 
  # used for all variables within the same MCMC run. In this case, the starting iterations 
  # are mapped onto the appropriate plots using the names of the vector elements. If plots 
  # of a certain test label are selected but not burn-in values passed, then the 
  # burn-in start iteration defaults to 1. There is no issue if test labels are included 
  # in `burn_in_start` but not in `samp_df_col_names`. This allows the user to set 
  # `burn_in_start` for all tests at once and then simply pass in the whole thing 
  # regardless of which subset of plots are being created without having to worry. 
  
  N_cols <- length(samp_df_col_names)
  
  # Set burn-in starting iteration values. 
  if(is.null(names(burn_in_start)) && (length(burn_in_start) == 1)) {
    burn_in_start <- rep(burn_in_start, N_cols)
  } else if(is.null(names(burn_in_start))) {
    if(length(burn_in_start) != N_cols) {
      stop("<burn_in_start> is unnamed and its length does not agree with number of plots.")
    }
  } else {
    burn_in_itrs <- vector(mode = "integer", length = N_cols)
    for(i in 1:N_cols) {
      test_label <- names(burn_in_start)[i]
      idx_sel <- grep(test_label, samp_df_col_names)
      burn_in_itrs[idx_sel] <- burn_in_start[i]
    }
    
    # Assign remaining columns burn-in start iteration of 1. 
    burn_in_itrs[burn_in_itrs == 0L] <- 1
  }
  
  return(burn_in_itrs)
  
}


# ------------------------------------------------------------------------------
# Validation Functions
# ------------------------------------------------------------------------------

MCMC_Gibbs_linear_Gaussian_model <- function(computer_model_data, theta_prior_params, Sig_eps_prior_params, N_mcmc) {
  # Implementation of a Gibbs sampler to sample posterior over `theta` and noise variance `sig2_eps`
  # for the linear Gaussian model (see `generate_linear_Gaussian_test_data()`). Current assumptions: 
  #    - theta has a zero-mean prior with diagonal prior covariance.
  #    - Single output (p=1), with iid Gaussian noise model, meaning there is only a single variance parameter `sig2_eps`. 
  #    - `sig2_eps` has an inverse Gamma prior. 
  #    - All thetas are calibrated (none are fixed). See `generate_linear_Gaussian_test_data()` for comments on this. 
  #      It will not be hard to remove this assumption. 
  #
  # Args:
  #    `computer_model_data`, `theta_prior_params`, `Sig_eps_prior_params` are all the standard objects fed into 
  #     the calibration procedure. The first two are returned by `generate_linear_Gaussian_test_data()`
  #    N_mcmc: integer(1), number of MCMC iterations. 
  #
  # Returns:
  #    list, with two elements `theta` and `sig2_eps`, each of which are matrices storing the MCMC samples. 

  d <- length(computer_model_data$pars_cal_sel)
  Sig_theta <- diag(theta_prior_params$param2^2)
  
  theta_samp_Gibbs <- matrix(nrow = N_mcmc, ncol = d)
  sig2_eps_samp_Gibbs <- matrix(nrow = N_mcmc, ncol = 1)
  
  theta_samp_Gibbs[1,] <- computer_model_data$theta_true
  sig2_eps_samp_Gibbs[1,] <- diag(computer_model_data$Sig_eps)
  theta_curr <- theta_samp_Gibbs[1,]
  sig2_eps_curr <- sig2_eps_samp_Gibbs[1,] 
  
  
  for(itr in seq(2, N_mcmc)) {
    
    # Sample theta conditional posterior. 
    Cov_post <- sig2_eps_curr * solve(crossprod(G) + sig2_eps_curr * diag(1/diag(Sig_theta)))
    mean_post <- (1/sig2_eps_curr) * tcrossprod(Cov_post, G) %*% computer_model_data$data_obs
    L_post <- t(chol(Cov_post))
    theta_curr <- mean_post + L_post %*% matrix(rnorm(nrow(L_post)), ncol=1)
    theta_samp_Gibbs[itr,] <- theta_curr
    
    # Sample sig2_eps conditional posterior. 
    a_post <- Sig_eps_prior_params$IG_shape + 0.5 * computer_model_data$n_obs + 1
    b_post <- Sig_eps_prior_params$IG_scale + 0.5 * sum((computer_model_data$data_obs - G %*% matrix(theta_curr, ncol=1))^2)
    sig2_eps_curr <- 1/rgamma(1, shape = a_post, rate = b_post)
    sig2_eps_samp_Gibbs[itr,] <- sig2_eps_curr
    
  }
  
  return(list(theta = theta_samp_Gibbs, sig2_eps = sig2_eps_samp_Gibbs))
  
}


# ------------------------------------------------------------------------------
# Calibration Test Functions: automating testing of different calibration schemes
# ------------------------------------------------------------------------------

run_calibration_emulator_comparison <- function(computer_model_data, theta_prior_params, emulator_settings, train_data, 
                                                learn_sig_eps = FALSE, sig_eps_prior_params = NULL, 
                                                test_data = NULL, N_mcmc_exact = 50000, N_mcmc_approx = 50000, 
                                                run_exact_mcmc = TRUE,
                                                test_labels = as.character(seq(1, nrow(emulator_settings))), ...) {
  # Automates the comparison of emulator-based MCMC calibration across different emulator settings. 
  # In particular, runs `run_calibration_single_emulator_setting()` for each emulator setting, where 
  # distinct emulator settings are provided in the rows of `emulator_settings`. Optionally also runs
  # exact (no emulator) calibration as a basis for comparison. Note that aside from the emulator settings, 
  # everything else (e.g. train data, priors, computer model data) are held constant across runs. 
  #
  # Args:
  #    computer_model_data: list, the standard computer model data list. 
  #    theta_prior_params: data.frame, with columns "dist", "param1", and "param2". The ith row of the data.frame
  #                        should correspond to the ith entry of 'theta'. Currently, accepted values of "dist" are 
  #                        "Gaussian" (param1 = mean, param2 = std dev) and "Uniform" (param1 = lower, param2 = upper).
  #    emulator_settings: data.frame, with columns "gp_lib", "kernel", "transformation_method", "scale_X", "normalize_y" used to specify 
  #                       GP emulator settings. Each row is interpreted as one emulator setting, and a calibration will be executed 
  #                       for each setting. 
  #    train_data: matrix of dimension N x d, where N is the number of design points and d the dimension of the input space. The training 
  #                inputs used to fit the GP emulator. 
  #    learn_sig_eps: logical, if TRUE treats observation covariance matrix as random and MCMC samples from joint 
  #                   posterior over Sig_eps and theta. Otherwise, fixes Sig_eps at value `Sig_eps_init`.
  #    sig_eps_prior_params: list, defining prior on Sig_eps. See `sample_prior_Sig_eps()` for details. Only required if `learn_Sig_eps` is TRUE.
  #    test_data: matrix, of dimension M x d, where M is the number of test points and d the number of calibration parameters. If provided, GP 
  #               predictions are computed at this set of points. Default is NULL. 
  #    N_mcmc_exact: integer, the number of MCMC iterations used for the exact MCMC scheme (used only if `run_exact_mcmc` is TRUE). 
  #    N_mcmc_approx: integer, the number of MCMC iterations used for each approximate emulator-assisted MCMC scheme. 
  #    run_exact_mcmc: logical(1), if TRUE runs exact non-emulator assisted MCMC and includes results in return list. Default is TRUE. 
  #    test_labels: character vector, of equal length to the number of rows in `emulator_settings`. These are labels used to identify each 
  #                 test run for each emulator setting. Default is just to label them by the row numbers in `emulator_settings`. 
  #    ...: additional MCMC arguments, see `mcmc_calibrate_ind_GP()` for options.
  #
  # Returns:
  #    list, including element "mcmc_approx_list" which is a list of length equal to the number of rows in `emualator_settings`. This
  #    sub-list contains the return value of `run_calibration_single_emulator_setting()` applied to each emulator setting. Also contains
  #    element "test_labels" containing the test labels, and element "samp_mcmc_exact" storing the MCMC results for the exact 
  #    calibration (NULL if `run_exact_mcmc` is FALSE). 
  
  # Run exact MCMC sampler (no emulator). 
  if(run_exact_mcmc) {
    time_start <- proc.time()
    samp_mcmc_exact <- mcmc_calibrate_product_lik(computer_model_data = computer_model_data, 
                                                  theta_prior_params = theta_prior_params, 
                                                  learn_sig_eps = learn_sig_eps,
                                                  sig_eps_prior_params = sig_eps_prior_params,
                                                  N_mcmc = N_mcmc_exact, ...)
    mcmc_exact_runtime <- (proc.time() - time_start)[["elapsed"]]
  } else {
    samp_mcmc_exact <- NULL
    mcmc_exact_runtime <- NULL
  }
  
  # Run separate emulator-assisted MCMC for each emulator setting.
  # TODO: get parallelization to work.
  FUN <- function(j) run_calibration_single_emulator_setting(computer_model_data = computer_model_data,
                                                             theta_prior_params = theta_prior_params,
                                                             emulator_setting = emulator_settings[j,],
                                                             train_data = train_data,
                                                             test_data = test_data,
                                                             learn_sig_eps = learn_sig_eps,
                                                             sig_eps_prior_params = sig_eps_prior_params,
                                                             N_mcmc = N_mcmc_approx, ...)
  # N_cores <- min(detectCores(), nrow(emulator_settings))
  # approx_list <- mclapply(seq(1, nrow(emulator_settings)), FUN, mc.cores = N_cores)
  approx_list <- lapply(seq(1, nrow(emulator_settings)), FUN)
  names(approx_list) <- test_labels

  return(list(samp_mcmc_exact = samp_mcmc_exact, 
              mcmc_approx_list = approx_list, 
              test_labels = test_labels, 
              mcmc_exact_runtime = mcmc_exact_runtime))
  
}


run_calibration_single_emulator_setting <- function(computer_model_data, theta_prior_params, emulator_setting, train_data, 
                                                    learn_sig_eps = FALSE, sig_eps_prior_params = NULL, 
                                                    test_data = NULL, N_mcmc = 50000, ...) {
  # A convenience function that fits a GP emulator and runs emulator-assisted MCMC. Also optionally computes emulator predictions 
  # at a set of test/validation points. Returns a list of all the relevant quantities. 
  #
  # Args:
  #    computer_model_data: list, the standard computer model data list. 
  #    theta_prior_params: data.frame, with columns "dist", "param1", and "param2". The ith row of the data.frame
  #                        should correspond to the ith entry of 'theta'. Currently, accepted values of "dist" are 
  #                        "Gaussian" (param1 = mean, param2 = std dev) and "Uniform" (param1 = lower, param2 = upper).
  #    emulator_setting: list, with elements "gp_lib", "kernel", "transformation_method", "scale_X", "normalize_y" used to specify 
  #                      GP emulator settings. 
  #    train_data: matrix of dimension N x d, where N is the number of design points and d the dimension of the input space. The training 
  #                inputs used to fit the GP emulator. 
  #    learn_sig_eps: logical, if TRUE treats observation covariance matrix as random and MCMC samples from joint 
  #                   posterior over Sig_eps and theta. Otherwise, fixes Sig_eps at value `Sig_eps_init`.
  #    sig_eps_prior_params: list, defining prior on Sig_eps. See `sample_prior_Sig_eps()` for details. Only required if `learn_Sig_eps` is TRUE.
  #    test_data: matrix, of dimension M x d, where M is the number of test points and d the number of calibration parameters. If provided, GP 
  #               predictions are computed at this set of points. Default is NULL. 
  #    N_mcmc: integer, the number of MCMC iterations. 
  #    ...: additional MCMC arguments, see `mcmc_calibrate_ind_GP()` for options. 
  #
  # Returns:
  #    list, with elements "emulator_info", "samp_mcmc", and "gp_pred_list". The third is NULL if `test_data` is NULL. 
  
  # Select correct "output_stats" used to normalize response variable, depending on whether emulator is log-normal process or not. 
  if(emulator_setting$transformation_method == "LNP") {
    output_stats <- train_data$log_output_stats
  } else {
    output_stats <- train_data$output_stats
  }
  
  # Fit GPs. 
  gp_fits <- fit_independent_GPs(train_data$inputs_scaled, train_data$outputs_normalized, emulator_setting$gp_lib, emulator_setting$kernel)$fits

  # Predict at test points.
  if(!is.null(test_data)) {
    gp_pred_list <- predict_independent_GPs(X_pred = test_data$inputs_scaled, gp_obj_list = gp_fits, gp_lib = emulator_setting$gp_lib, 
                                            denormalize_predictions = TRUE, output_stats = output_stats)
  } else {
    gp_pred_list <- NULL
  }
  
  # Truncate priors on calibration parameters. 
  if(isTRUE(emulator_setting$truncate_theta_prior)) {
    theta_prior_params <- truncate_prior_theta(theta_prior_params, train_data$input_bounds)
  }
  
  # Emulator info list. 
  emulator_info <- list(gp_fits = gp_fits, 
                        input_bounds = train_data$input_bounds,
                        output_stats = output_stats, 
                        settings = emulator_setting)
  
  # Emulator-assisted MCMC. 
  time_start <- proc.time()
  samp_mcmc <- mcmc_calibrate_ind_GP(computer_model_data = computer_model_data, 
                                     theta_prior_params = theta_prior_params,
                                     emulator_info = emulator_info, 
                                     learn_sig_eps = learn_sig_eps, 
                                     sig_eps_prior_params = sig_eps_prior_params, 
                                     N_mcmc = N_mcmc, ...)
  mcmc_runtime <- (proc.time() - time_start)[["elapsed"]]
                                     
  return(list(emulator_info = emulator_info, 
              samp_mcmc = samp_mcmc, 
              gp_pred_list = gp_pred_list, 
              mcmc_runtime = mcmc_runtime, 
              theta_prior_params_trunc = theta_prior_params)) 
                                               
} 


get_mcmc_samp_df <- function(calibration_results_list) {
  # A convenience function that extracts the MCMC samples from `run_calibration_emulator_comparison()`
  # and arranges them in a data.frame. The columns will contain all "theta" and "sig_eps" samples from 
  # each approximate MCMC scheme. If exact MCMC results are also present in `calibration_results_list` 
  # then the exact samples will also be included. This function is robust to the case when different 
  # MCMC schemes have different numbers of samples, in which case those with fewer samples will 
  # be padded with NAs. Column names for the approximate sample schemes are determined by 
  # `calibration_results_list$test_labels` as well as the name of the sampled parameter. The name formatting
  # is "<test_label>_<parameter_label>_<parameter_name>", where "<parameter_label>" is either "sig_eps" or 
  # "theta". 
  # 
  # Samples from exact MCMC are labeled as "exact". Every MCMC algorithm is assumed to have "theta" 
  # samples, but  may or may not have "sig_eps" samples, depending on whether the observation variance  
  # was fixed or not. 
  #
  # Args:
  #    calibration_results_list: list, the returned value from `run_calibration_emulator_comparison()`. 
  #
  # Returns:
  #    data.frame, containing all MCMC samples as described above. 
  
  # Approximate samples list. 
  mcmc_samp_list <- lapply(calibration_results_list$mcmc_approx_list, function(l) l$samp_mcmc$samp)
  names(mcmc_samp_list) <- names(calibration_results_list$mcmc_approx_list)
  
  # If exact samples are included, add to list. 
  samp_mcmc_exact <- calibration_results_list$samp_mcmc_exact$samp
  if(!is.null(samp_mcmc_exact)) {
    mcmc_samp_list[["exact"]] <- samp_mcmc_exact
  }

  N_samp <- sapply(mcmc_samp_list, nrow)
  N_max <- max(N_samp)
  for(j in seq_along(mcmc_samp_list)) {
    
    # Pad with NAs in the case that the number of samples is different across different parameters.
    if(nrow(mcmc_samp_list[[j]]) < N_max) {
      mcmc_samp_list[[j]] <- rbind(mcmc_samp_list[[j]], matrix(NA, nrow = N_max - nrow(mcmc_samp_list[[j]]),
                                                                   ncol = ncol(mcmc_samp_list[[j]])))         
    }
      
    # Set column names. 
    colnames(mcmc_samp_list[[j]]) <- paste(names(mcmc_samp_list)[j], colnames(mcmc_samp_list[[j]]), sep = "_")
                                                        
  }

  # Combine approximate samples into single data.frame. 
  df_samples <- as.data.frame(do.call("cbind", mcmc_samp_list))
  
  return(df_samples)
  
}










