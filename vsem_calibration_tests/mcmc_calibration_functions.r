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


llik_Gaussian_err <- function(model_errs, Sig_eps, output_vars = NA, normalize = FALSE) {
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
  #               log density (where Sig_eps is treated as a constant). Default is FALSE. 
  #
  # Returns:
  #    numeric, the unnormalized log-likelihood. 
  
  if(!is.na(output_vars)) {
    Sig_eps <- Sig_eps[output_vars, output_vars]
    model_errs <- model_errs[, ..output_vars]
  }
  L <- t(chol(Sig_eps))
  log_quadratic_form <- sum(forwardsolve(L, t(model_errs))^2)
  
  if(normalize) {
    return(-0.5 * log_quadratic_form - 0.5 * prod(dim(model_errs)) * log(2*pi) - nrow(model_errs) * sum(log(diag(L))))
  }
  
  return(-0.5 * log_quadratic_form)
  
}


run_VSEM <- function(par, par_ref, par_cal_sel, PAR, output_vars = c("NEE", "Cv", "Cs", "CR")) {
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
  #    output_vars: character vector, selects columns of output from VSEM model. Default returns
  #                 all four columns. 
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
  model_errs_curr <- data_obs[, ..output_vars] - run_VSEM(theta_init, par_ref, par_cal_sel, PAR, output_vars)
  lprior_theta_curr <- calc_lprior_theta(theta_init, theta_prior_params)
  
  # Proposal covariance
  Cov_prop <- diag(proposal_scale_init, nrow = d)
  log_scale_prop <- 0
  L_prop <- t(chol(Cov_prop))
  accept_count <- 0
  
  # TEMP
  prop_vals <- matrix(nrow = N_mcmc, ncol = par_cal_sel)
  prop_sd <- vector(length = N_mcmc, mode = "numeric")

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
    
    # theta proposal
    # print(sqrt(exp(log_scale_prop)) * as.numeric(L_prop))
    theta_prop <- (theta_samp[itr-1,] + sqrt(exp(log_scale_prop)) * L_prop %*% matrix(rnorm(d), ncol = 1))[1,]
    
    # TEMP
    prop_vals[itr,] <- theta_prop
    prop_sd[itr] <- as.numeric(L_prop * sqrt(exp(log_scale_prop)))

    # Calculate log-likelihoods for current and proposed theta
    model_errs_prop <- data_obs[, ..output_vars] - run_VSEM(theta_prop, par_ref, par_cal_sel, PAR, output_vars)
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
  
  return(list(theta = theta_samp, Sig_eps = Sig_eps_samp, prop_vals = prop_vals, prop_sd = prop_sd))

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


calc_SSR <- function(data_obs, model_outputs_list) {
  # Computes the sum of squared residuals (SSR) between model runs and observed 
  # data on a per-output basis. Can handle multiple model runs (e.g. one per 
  # design point for emulation) or outputs from single model run (e.g. as required
  # during MCMC).
  #
  # Args:
  #    data_obs: data.table of dimension n x p, where n is the length of the time 
  #              series outputs from the model, and p is the number of output variables.
  #              This is Y in my typical notation. 
  #    model_outputs_list: either a data.table of dimension n x p corresponding to the 
  #                        model outputs f(theta) from a single model run. Or a list 
  #                        {f(theta_1), ..., f(theta_N)} of such data.tables, collecting
  #                        the outputs from multiple model runs. 
  # 
  # Returns:
  #    matrix of dimension N x p, where N is the number of model runs and p is the 
  #    number of output variables. The (i, j) entry of the matrix is the SSR
  #    for the jth output of the ith model run, i.e. ||Y_j - f(j, theta_i)||^2.
  
  # If model only run at single set of calibration parameter values
  if(is.data.table(model_outputs_list)) {
    model_outputs_list <- list(model_outputs_list)
  }
  
  N_param_runs <- length(model_outputs_list)
  N_outputs <- ncol(model_outputs_list[[1]])
  SSR <- matrix(nrow = N_param_runs, ncol = N_outputs)
  
  for(i in seq_along(model_outputs_list)) {
    SSR[i,] <- colSums(data_obs - model_outputs_list[[i]])^2
  }
  
  return(SSR)
  
}

# TODO: Parallelize independent GP fitting
fit_independent_GPs <- function(X_train, Y_train, gp_lib, gp_kernel) {
  # Builds on top of `fit_GP()` to generalize to multivariate GP regression. Simply 
  # fits independent GPs for each output (where the same input points are used for each GP).
  #
  # Args:
  #    X_train: matrix of shape N x d, where N is the number of design (training) points
  #             and d is the dimension of the input space. 
  #    y_train: matrix of shape N x p, with jth column containing the training outputs 
  #             for the jth output variable corresponding to the training inputs. 
  #    gp_lib: character, string specifying the GP package to use. Currently supports 
  #            "mlegp" and "hetGP". 
  #    gp_kernel: character, string specifying the kernel family to use. Potential options are 
  #               limited by the specified GP library. Currently only supports "Gaussian" (i.e. the 
  #               squared exponential/exponentiated quadratic/radial basis function kernel).
  #
  # Returns:
  #    list, with named elements "fits" and "times". The first is itself a list, 
  #    with length `ncol(Y_train)`; i.e. length equal to the number of outputs. The jth 
  #    element stores the fit GP object returned by `fit_GP()`, containing the fit GP for the jth output. 
  #    The "times" element is a numeric vector storing the times required to fit each GP. 

  p <- ncol(Y_train)
  GP_objects <- vector(mode = "list", length = p)
  GP_fit_times <- vector(mode = "numeric", length = p)
  
  for(j in seq_len(p)) {
    GP_fit_list <- fit_GP(X_train, Y_train[,j, drop = FALSE], gp_lib, gp_kernel)
    GP_objects[[j]] <- GP_fit_list[["fit"]]
    GP_fit_times[j] <- GP_fit_list[["time"]]
  }
  
  return(list(fits = GP_objects, times = GP_fit_times))
  
}


fit_GP <- function(X_train, y_train, gp_lib, gp_kernel) {
  # Estimate kernel and mean function hyperparameters for a univariate GP regression 
  # using specified Gaussian Process (GP) package. 
  #
  # Args:
  #    X_train: matrix of shape N x d, where N is the number of design (training) points
  #             and d is the dimension of the input space. 
  #    y_train: matrix of shape N x 1, the training outputs corresponding to the training inputs. 
  #    gp_lib: character, string specifying the GP package to use. Currently supports 
  #            "mlegp" and "hetGP". 
  #    gp_kernel: character, string specifying the kernel family to use. Potential options are 
  #               limited by the specified GP library. Currently only supports "Gaussian" (i.e. the 
  #               squared exponential/exponentiated quadratic/radial basis function kernel). 
  #
  # Returns:
  #    Returns, list with named elements "fit" and "time". The former stores the fit GP object 
  #    returned by the fitting function in the specified GP package. The latter contains the 
  #    time taken by the fitting procedure. 
  
  eps_nug <- sqrt(.Machine$double.eps)
  
  tic <- proc.time()[3]
  if(gp_lib == "mlegp") {
    gp_fit <- mlegp(X_train, y_train, nugget.known = 1, nugget = eps_nug, constantMean = 1)
  } else if(gp_lib == "hetGP") {
    gp_fit <- mleHomGP(X_train, y_train, covtype = gp_kernel, known = list(g = eps_nug))
  } else {
    stop("Invalid GP library: ", gp_lib)
  }

  toc <- proc.time()[3]
  gp_fit_time <- toc - tic
  
  return(list(fit = gp_fit, time = gp_fit_time))
  
}


LHS_train_test <- function(N_train, N_test, prior_params, joint = TRUE, extrapolate = TRUE) {
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
  
  return(list(X_train = X_train, X_test = X_test))
  
}


prep_GP_training_data <- function(X = NULL, Y = NULL, scale_X = FALSE, normalize_Y = FALSE) {
  # Preps training data inputs X and outputs Y for fitting Gaussian Process (GP) model by 
  # scaling each input variable to the unit interval, and normalizing the output variable 
  # by subtracting its mean and dividing by its standard deviation. The bounds used to 
  # standardize X are computed using the range of each input variable in the training set X.
  # If Y has multiple columns, each column is treated as a different output and each column
  # is normalized independently. 
  #
  # Source for matrix operations code to normalize X: find_reps() function of hetGP package. 
  #
  # Args:
  #    X: matrix of dimension N x d, where N is the number of training points and d is the dimension 
  #       of the input space. 
  #    Y: matrix of dimension N x p where p is the number of output variables.
  #    scale_X: logical, if TRUE, maps X to d-dimensional unit hypercube. 
  #    normalize_Y: logical, if TRUE, transforms y via Z-score. 
  #
  # Returns:
  #    list, with named elements "X", "Y", "input_bounds", and "output_stats". X and Y are the 
  #    training inputs and outputs, and may or may not be standardized/normalized depending on 
  #    the arguments `scale_X` and `normalize_Y`. `input_bounds` is a 2 x d matrix summarizing the 
  #    range of each input variable used to standardize X. The first row contains the minimum 
  #    value of each of the respective input variables in the training set, and similarly the 
  #    second row stores the maxima. `output_stats` is a named vector with names 
  #    "mean_Y" and "var_Y" storing the mean and variance of Y used to compute the Z-scores. 
  #    If `scale_X` is FALSE then "input_bounds" will be NULL and likewise with 
  #    `normalize_Y` and "output_stats". 
  
  if(!is.null(X) && scale_X) {
    input_bounds <- apply(X, 2, range)
    X <- scale_input_data(X, input_bounds)
  } else {
    input_bounds <- NULL
  }
  
  if(normalize_Y) {
    output_stats <- rbind(apply(Y, 2, mean), apply(Y, 2, var))
    rownames(output_stats) <- c("mean_Y", "var_Y")
    Y <- normalize_output_data(Y, output_stats)
  } else {
    output_stats <- NULL
  }
  
  return(list(X = X, Y = Y, input_bounds = input_bounds, output_stats = output_stats))
  
}


scale_input_data <- function(X, input_bounds, inverse = FALSE) {
  # Transforms input data X (typically points at which to predict) by performing a linear scaling encoded
  # in `input_bounds`. Can also perform the inverse transformation. 
  #
  # Args:
  #    X: matrix of dimension N_pred x d, where N_pred is the number of input points and d is the dimension of the 
  #       input space. 
  #    input_bounds: 2 x d matrix summarizing the range of each input variable used to standardize the training inputs. 
  #                  The first row contains the minimum value of each of the respective input variables in the training set, 
  #                  and similarly the second row stores the maxima. This object is returned by prep_GP_training_data(). 
  #    inverse: logical, if TRUE, treats X as already scaled and reverses the scaling. Otherwise performs the forward 
  #             transformation. Default is FALSE. 
  #
  # Returns:
  #    matrix of dimension N_pred x d; the matrix X whose columns have been scaled according to `input_bounds`. 
  #
  
  if(inverse) {
    X <- X %*% diag(input_bounds[2,] - input_bounds[1,], ncol(X)) + matrix(input_bounds[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  } else {
    X <- (X - matrix(input_bounds[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% diag(1/(input_bounds[2,] - input_bounds[1,]), ncol(X))
  }
  
  return(X)
}


normalize_output_data <- function(Y, output_stats, inverse = FALSE) {
  # Transforms output data by transforming Y to Z-scores, or performing the reverse operation. 
  # Can handle multiple output variables, where each column of `Y` is treated as an output 
  # variable and is normalized independently. 
  #
  # Args:
  #    Y: matrix of dimension N x p where p is the number of output variables. 
  #    output_stats: a matrix of dimensions 2xp, where p is the number of output variables. The matrix 
  #                  must have rownames "mean_Y" and "var_Y" storing the mean and 
  #                  variance of each output variable used to compute the Z-scores. This object is 
  #                  returned by prep_GP_training_data(). 
  #
  # Returns: 
  #    matrix of dimension N x p, the transformed version of Y. 
  
  if(inverse) {
    Y <- Y * matrix(sqrt(output_stats["var_Y",]), nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE) + 
             matrix(output_stats["mean_Y",], nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
  } else {
    Y <- (Y - matrix(output_stats["mean_Y",], nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)) / 
          matrix(sqrt(output_stats["var_Y",]), nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
  }
  
  return(Y)
  
}


predict_independent_GPs <- function(X_pred, gp_obj_list, gp_lib, cov_mat = FALSE, 
                                    transform_predictions = FALSE, output_stats = NULL, inverse = TRUE) {
  # A wrapper function for predict_GP() that generalizes the latter to generating predictions for 
  # multi-output GP regression using independent GPs. 
  #
  # Args:
  #    X_pred: matrix of dimension N_pred x d, where d is the dimension of the input space. Each row is an input 
  #            point at which to predict. 
  #    gp_obj_list: A list of GP objects, each of which represents a GP fit to one of the outputs. The objects will 
  #                differ based on the specific GP library used for fitting. All of the objects in the list must have 
  #                been fit using the same library. 
  #    gp_lib: character(1), the library used to fit the GP. Currently supports "mlegp" or "hetGP". 
  #    cov_mat: logical(1), if TRUE, calculates and returns the N_pred x N_pred predictive covariance matrix 
  #             over the set of input points. Otherwise, only calculates the pointwise predictive variances. 
  #    transform_predictions: logical(1), if TRUE, applies linear transformation to predictions, 
  #                           inverting the Z-score transformation.
  #    output_stats: If not NULL, then a matrix of dimensions 2xp, where p is the number of output variables. The matrix 
  #                  must have rownames "mean_Y" and "var_Y" storing the mean and 
  #                  variance of each output variable used to compute the Z-scores. This object is 
  #                  returned by prep_GP_training_data(). Only required if `transform_predictions` is TRUE. 
  # 
  # Returns:
  #    list, with length equal to the length of `gp_obj_list`. Each element of this list is itself a list, 
  #    with named elements "mean", "sd2", "sd2_nug", "cov" (the output of the function `predict_GP()` applied 
  #    to each GP in `gp_obj_list`). 
  
  lapply(seq_along(gp_obj_list), function(j) predict_GP(X_pred, gp_obj_list[[j]], gp_lib, cov_mat, transform_predictions, 
                                              output_stats[,j,drop=FALSE]))
  
}


# TODO: look into rebuild(robust=TRUE) for hetGP prediction (using ginv rather than Cholesky for matrix inverse).
predict_GP <- function(X_pred, gp_obj, gp_lib, cov_mat = FALSE, 
                       transform_predictions = FALSE, output_stats = NULL) {
  # Calculate GP predictive mean, variance, and optionally covariance matrix at specified set of 
  # input points. 
  #
  # Args:
  #    X_pred: matrix of dimension N_pred x d, where d is the dimension of the input space. Each row is an input 
  #            point at which to predict. 
  #    gp_obj: An object representing a GP fit, which will differ based on the library used to fit the GP. 
  #    gp_lib: character(1), the library used to fit the GP. Currently supports "mlegp" or "hetGP". 
  #    cov_mat: logical(1), if TRUE, calculates and returns the N_pred x N_pred predictive covariance matrix 
  #             over the set of input points. Otherwise, only calculates the pointwise predictive variances. 
  #    transform_predictions: logical(1), if TRUE, applies linear transformation to predictions, 
  #                           inverting the Z-score transformation.
  #    output_stats: If not NULL, then a matrix of dimensions 2x1. The matrix 
  #                  must have rownames "mean_Y" and "var_Y" storing the mean and 
  #                  variance of the output variable used to compute the Z-scores. This object is 
  #                  returned by prep_GP_training_data(). Only required if `transform_predictions` is TRUE. 
  #
  # Returns:
  #    list, with named elements "mean", "sd2", "sd2_nug", "cov". 
  
  pred_list_names <- c("mean", "sd2", "sd2_nug", "cov")
  pred_list <- vector(mode = "list", length = length(pred_list_names))
  names(pred_list) <- pred_list_names
  
  if(gp_lib == "mlegp") {
    mlegp_pred <- predict(gp_obj, newData = X_pred, se.fit = TRUE)
    pred_list[["mean"]] <- mlegp_pred[["fit"]]
    pred_list[["sd2"]] <- mlegp_pred[["se.fit"]]^2
  } else if(gp_lib == "hetGP") {
    # Second matrix for computing predictive covariance
    if(cov_mat) {
      X_prime <- X_pred
    } else {
      X_prime <- NULL
    }
    hetGP_pred <- predict(gp_obj, X_pred, xprime = X_prime)
    pred_list[pred_list_names] <- hetGP_pred[c("mean", "sd2", "nugs", "cov")]
  }
  
  # Map back to original scale
  if(transform_predictions) {
    pred_list[["mean"]] <- output_stats["mean_Y",1] + sqrt(output_stats["var_Y",1]) * pred_list[["mean"]]
    pred_list[["sd2"]] <- output_stats["var_Y",1] * pred_list[["sd2"]]
    pred_list[["nugs"]] <- output_stats["var_Y",1] * pred_list[["nugs"]]
    pred_list[["cov"]] <- output_stats["var_Y",1] * pred_list[["cov"]]
  }
  
  return(pred_list)
  
}


get_train_test_data <- function(N_train, N_test, prior_params, joint, extrapolate, 
                                ref_pars, pars_cal_sel, data_obs, PAR, output_vars, scale_X, normalize_Y) {
  # Generates training (design) dataset and test dataset via Latin Hypercube sampling. 
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
  #    joint: logical, if TRUE jointly samples train/test. Default is TRUE. 
  #    extrapolate: logical, if TRUE no truncation is performed for the test samples. Otherwise 
  #                 truncation is performed to guarantee the samples are within the extent of the 
  #                 training set. Default is TRUE. 
  #    ref_pars: data.frame, rownames should correspond to parameters of computer model. 
  #             Must contain column named "best". Parameters that are fixed at nominal 
  #             values (not calibrated) are set to their values given in the "best" column. 
  #    par_cal_sel: integer vector, selects the rows of 'ref_pars' that correspond to 
  #                 parameters that will be calibrated. 
  #    data_obs: data.table, dimension n x p (n = length of time series, p = number outputs).
  #              Colnames set to output variable names. 
  #    output_vars: character vector, used to the select the outputs to be considered in 
  #                 the likelihood; e.g. selects the correct sub-matrix of 'Sig_eps' and the 
  #                 correct columns of 'data_obs'. 
  #    PAR: numeric vector, time series of photosynthetically active radiation used as forcing
  #         term in VSEM.
  #
  # Returns:
  #    list, with names "X_train", "X_test", "X_train_preprocessed", "X_test_preprocessed", 
  #    "Y_train", "Y_train_preprocessed", "Y_test", "input_bounds", and "output_stats".
  #    The latter two encode the transformations applied to the inputs and outputs, respectively 
  #    (see prep_GP_training_data()). 
  
  # Generate train and test input datasets via Latin Hypercube Sampling
  X_train_test <- LHS_train_test(N_train, N_test, prior_params, joint, extrapolate)
  X_train <- X_train_test$X_train
  X_test <- X_train_test$X_test
  
  # Run VSEM on train and test sets to obtain outputs
  model_outputs_list_train <- lapply(X_train, function(theta) run_VSEM(theta, ref_pars, pars_cal_sel, PAR, output_vars))
  model_outputs_list_test <- lapply(X_test, function(theta) run_VSEM(theta, ref_pars, pars_cal_sel, PAR, output_vars))
  Y_train <- calc_SSR(data_obs, model_outputs_list_train)
  Y_test <- calc_SSR(data_obs, model_outputs_list_test)
  
  # Pre-process training data: scale inputs and normalize outputs
  # TODO: normalizing Y in the log(SSR) case must be done separately
  GP_data_preprocessed <- prep_GP_training_data(X_train, Y_train, scale_X, normalize_Y)
  X_train_preprocessed <- GP_data_preprocessed$X
  Y_train_preprocessed <- GP_data_preprocessed$Y
  
  # Pre-process testing input data
  X_test_preprocessed <- scale_input_data(X_test, GP_data_preprocessed$input_bounds)
  
  return(list(X_train = X_train, X_test = X_test, X_train_preprocessed = X_train_preprocessed, 
              X_test_preprocessed = X_test_preprocessed, Y_train = Y_train, 
              Y_train_preprocessed = Y_train_preprocessed, Y_test = Y_test, 
              input_bounds = GP_data_preprocessed$input_bounds, 
              output_stats = GP_data_preprocessed$output_stats))
  
}


evaluate_GP_emulators <- function(emulator_settings, N_iter, N_design_points, N_test_points, 
                                  theta_prior_params, ref_pars, pars_cal_sel, data_obs, PAR, output_vars, joint = TRUE, extrapolate = FALSE) {
  
  # TODO: allow storing separate error metrics for each output
  # Create list to store results
  nrow_test_results <- N_iter * nrow(emulator_settings)
  colnames_test_results <- c("gp_libs", "target", "kernel", "scale_X", "normalize_y",
                             paste0("rmse_", output_vars), paste0("rmse_scaled_", output_vars),
                             paste0("fit_time_", output_vars))
  test_results <- data.frame(matrix(nrow = nrow_test_results, ncol = length(colnames_test_results)))
  names(test_results) <- colnames_test_results
  
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
                                           output_vars = output_vars)
      
    # Loop over each model specification
    # TODO: modify for log(SSR).
    for(j in seq_along(emulator_settings)) {
      
      if(emulator_settings[j, "scale_X"]) X <- train_test_data$X_train_preprocessed
      else X <- train_test_data$X_train
      
      if(emulator_settings[j, "normalize_y"]) Y <- train_test_data$Y_train_preprocessed
      else Y <- train_test_data$Y_train
      
      GP_fit_list <- fit_independent_GPs(X, Y, emulator_settings[j, "gp_lib"], emulator_settings[j, "kernel"])
      GP_pred_list <- predict_independent_GPs(X_pred, gp_obj_list, gp_lib)
      GP_pred_errs <- calc_GP_pred_errs(GP_pred_list, train_test_data$Y_test)
      
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











