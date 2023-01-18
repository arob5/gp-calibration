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
  #                 series and p = number output variables. This is the model error matrix
  #                 Y - f(theta). 
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


mcmc_calibrate <- function(par_ref, par_cal_sel, data_obs, output_vars, PAR,
                           theta_init = NA, theta_prior_params, 
                           learn_Sig_eps = FALSE, Sig_eps_init = NA, Sig_eps_prior_params = NA, diag_cov = FALSE, 
                           N_mcmc, adapt_frequency, adapt_min_scale, accept_rate_target) {
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
  #    data_obs: data.table, dimension n x p (n = length of time series, p = number outputs).
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
  #    min_scale: numeric scalar, used as a floor for the scaling factor in covariance adaptation, see `adapt_cov_proposal()`.
  #    accept_rate_target: numeric scalar, the desired acceptance rate, see `adapt_cov_proposal()`.
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
  if(isTRUE(Sig_eps_prior_params$dist == "IG")) {
    diag_cov <- TRUE
  }
  if(diag_cov) {
    Sig_eps_samp <- matrix(nrow = N_mcmc, ncol = p)
  } else {
    Sig_eps_samp <- matrix(nrow = N_mcmc, ncol = 0.5*p*(p+1)) # Each row stores lower triangle of Sig_eps
  }

  # Set initial values
  if(is.na(theta_init)) {
    theta_init <- sample_prior_theta(theta_prior_params)
  }
  if(is.na(Sig_eps_init)) {
    Sig_eps_init <- sample_prior_Sig_eps(Sig_eps_prior_params)
  }
  
  theta_samp[1,] <- theta_init
  if(diag_cov) {
    Sig_eps_samp[1,] <- diag(Sig_eps_init)
  } else {
    Sig_eps_samp[1,] <- lower.tri(Sig_eps_init, diag = TRUE)
  }
  
  Sig_eps_curr <- Sig_eps_init
  model_errs_curr <- data_obs - run_VSEM(theta_init, par_ref, par_cal_sel, PAR)
  log_theta_post_curr <- llik_Gaussian_err(model_errs_curr, Sig_eps_curr, output_vars) + calc_lprior_theta(theta_init, theta_prior_params)
  
  # Proposal covariance
  Cov_prop <- diag(1, nrow = d)
  L_prop <- t(chol(Cov_prop))
  accept_count <- 0

  for(itr in seq(2, N_mcmc)) {

    #
    # Metropolis step for theta
    #

    # Adapt proposal covariance matrix
    if((itr > 3) && ((itr - 1) %% adapt_frequency) == 0) {
      Cov_prop <- adapt.cov.proposal(Cov_prop, theta_samp[(itr - adapt_frequency):(itr - 1),,drop=FALSE], 
                                     adapt_min_scale, accept_count / adapt_frequency, accept_rate_target)
      L_prop <- t(chol(Cov_prop))
      accept_count <- 0
    }
    
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
      accept_count <- accept_count + 1 
    } else {
      theta_samp[itr,] <- theta_samp[itr-1,]
    }

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
    output_l2_errs <- colsums(model_errs^2)
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
    scale.factor <- max(accept_rate / accept_rate_target, min_scale)
    stdev <- apply(sample_history, 2, stats::sd)
    scale_mat <- scale_factor * diag(stdev, nrow = length(stdev))
    return(scale_mat %*% cor_estimate %*% scale_mat)
  }
}











