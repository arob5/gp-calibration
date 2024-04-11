#
# statistical_helper_functions.r
# General helper functions for statistical inverse problems.  
#
# Andrew Roberts
# 

convert_Gaussian_to_LN <- function(mean_Gaussian, var_Gaussian=NULL, cov_Gaussian=NULL, 
                                   return_mean=TRUE, return_var=TRUE, return_cov=FALSE, log_scale=FALSE) {
  # Given the mean and either variance or covariance matrix of a Gaussian random vector X, computes
  # the mean, variance, and covariance of Y := exp(X), which is log-normally (LN) distributed. 
  # Optionally, returns the log of these quantities, which is often recommended for numerical stability. 
  # The log scale option is not allowed when `return_cov` is TRUE, since the covariance matrix can contain
  # negative values. In the below descriptions, let N denote the length of X. 
  #
  # Args:
  #    mean_Gaussian: numeric(N), or Nx1 matrix; the mean vector of X. 
  #    var_Gaussian: numeric(N), or Nx1 matrix; the variance of each component of X. Required
  #                  if `cov_Gaussian` is NULL. 
  #    cov_Gaussian: NxN matrix, the covariance matrix of X. Required if `var_Gaussian` is NULL.
  #                  Required if `return_cov` is NULL. 
  #    return_mean: logical(1), whether or not to return the (log) mean of Y. 
  #    return_var: logical(1), whether or not to return the (log) variances of Y.
  #    return_cov: logical(1), whether or not to return the covariance matrix of Y. 
  #    log_scale: logical(1), whether to return the log of the mean/variance of Y. Must be 
  #               FALSE if `return_cov` is TRUE. 
  #
  # Returns: 
  #    list with the LN computations. Potential list arguments include "mean", "var", "cov", 
  #    "log_mean", and "log_var". 
                                   
  assert_that(!is.null(var_Gaussian) || !is.null(cov_Gaussian), 
              msg="Gaussian variance or covariance matrix required to compute log-normal moments.")
  assert_that(!(log_scale && return_cov), msg="Cannot return LN moments on log scale if `return_cov` is TRUE.")
  if(is.null(var_Gaussian)) var_Gaussian <- diag(cov_Gaussian)

  # Compute the log of the LN mean and variance. Note that this can't be done for the 
  # covariance given the that covariance matrix can have negative values. 
  log_mean <- drop(mean_Gaussian) + 0.5*drop(var_Gaussian)
  if(return_var) log_var <- log_exp_minus_1(var_Gaussian) + 2*log_mean
  if(log_scale) {
    return_list <- list(log_mean=log_mean)
    if(return_var) return_list$log_var <- log_var
    return(return_list)
  }
  
  # Continue if LN moments are requested on the original scale. The mean is always returned.
  return_list <- list()
  if(return_mean) return_list$mean <- exp(log_mean)
  if(return_var) return_list$var <- exp(log_var)
  if(return_cov) {
    assert_that(!is.null(return_cov), msg="Computing cov requires non-NULL `cov_Gaussian`.")
    return_list$cov <- exp(outer(log_mean, log_mean, FUN="+")) * (exp(cov_Gaussian)-1)
  }
  
  return(return_list)
                                
}


gen_lin_Gaus_NIW_test_data <- function(G_list, Sig_eps=NULL, mu0=NULL, Sig0=NULL, 
                                       IW_scale=NULL, IW_dof=NULL, N_missing_obs=NULL) {
  # TODO: write up blog post on Normal Inverse Wishart model and then write this function.  
  
  .NotYetImplemented()
  
}


generate_linear_Gaussian_test_data <- function(random_seed, N_obs, D, Sig_theta, G, 
                                               sig2_eps=NULL, sig_eps_frac=0.1, pars_cal_sel=NULL) {
  # Sets up a test example (including computer model data and prior distributions) in which the forward model is linear 
  # (represented by matrix `G`), the likelihood is Gaussian (with variance `sig2_eps`), and the prior on the calibration 
  # parameters is zero-mean Gaussian (with diagonal covariance matrix `Sig_theta`). Note that currently this function 
  # only creates an example with a single output variable (p = 1). Also, for the time being `Sig_theta` must be diagonal, 
  # until the prior code is updated to allow for correlated priors. Currently, a zero-mean prior is generated, which 
  # should also be generalized in the future. This linear Gaussian setup admits a closed form 
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
  #    sig2_eps: numeric(1), the observation variance. If not provided, then the observation variance will be set using `coef_var`, or 
  #              if `coef_var`.
  #    sig_eps_frac: numeric(1), If `sig2_eps` is provided directly then `sig_eps_frac` will not be used.
  #                  Otherwise, the noise variance  will be set so that the standard deviation sqrt(sig2_eps) equals 
  #                  the range of the data * `sig_eps_frac`. 
  #    G: matrix of dimension N_obs x D, the linear forward model.
  #    pars_cal_sel: integer(), vector containing the indices used to select the parameters which will be calibrated. The remaining parameters
  #                  will be fixed. If NULL, calibrates all parameters. 
  #    
  # Returns:
  #    list with 3 elements: computer_model_data, theta_prior_params, and true_posterior. 
  
  set.seed(random_seed)
  
  if(!all.equal(dim(G), c(N_obs, D))) {
    stop("Forward model G must be matrix of dimension N_obs x D.")
  }
  
  diag_prior_cov <- TRUE
  if(!isTRUE(all.equal(Sig_theta, diag(diag(Sig_theta))))) {
    diag_prior_cov <- FALSE
    if(!is.null(pars_cal_sel)) {
      stop("Currently linear Gaussian data function does not support factor fixing with non-diagonal prior covariance.")
    }
  }
  
  # Sample from model to generate observed data. 
  L_theta <- t(chol(Sig_theta))
  theta <- L_theta %*% matrix(rnorm(D), ncol=1)
  data_ref <- G %*% theta
  
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
  if(diag_prior_cov) {
    theta_prior_params <- data.frame(dist = rep("Gaussian", length(pars_cal_names)), 
                                     param1 = rep(0, length(pars_cal_names)), 
                                     param2 = sqrt(diag(Sig_theta)[pars_cal_sel]))
    rownames(theta_prior_params) <- pars_cal_names
  } else {
    theta_prior_params <- NULL
  }
  
  # True posterior (note that we need to adjust for the case where only a subset of the parameters 
  # are calibrated). The posterior moments are computed using the SVD and Woodbury identity. This 
  # assumes the number of calibration parameters is <= N_obs. 
  pars_fixed_sel <- setdiff(seq(1,D), pars_cal_sel)
  G_cal <- G[, pars_cal_sel]
  theta_cal <- theta[pars_cal_sel]
  Sig_theta_cal <- Sig_theta[pars_cal_sel, pars_cal_sel]
  
  if(length(pars_fixed_sel) == 0) {
    G_fixed <- matrix(0, nrow=N_obs, ncol=1)
    theta_fixed <- 0
  } else {
    G_fixed <- G[,pars_fixed_sel]
    theta_fixed <- theta[pars_fixed_sel]
  }
  y_adjusted <- matrix(data_obs, ncol=1) - G_fixed %*% theta_fixed
  
  # TODO: check SVD calculation. For now just directly taking inverse. 
  # svd_list <- svd(G_cal)
  # Cov_post <- Sig_theta_cal - Sig_theta_cal %*% diag(1 / (sig2_eps * (svd_list$d^(-2)) + diag(Sig_theta_cal))) %*% Sig_theta_cal
  
  # Cov_post <- sig2_eps * solve(crossprod(G_cal) + sig2_eps * diag(1/diag(Sig_theta_cal))) # old
  
  Sig_theta_cal_inv <- chol2inv(chol(Sig_theta_cal))
  Cov_post <- chol2inv(chol(crossprod(G_cal)/sig2_eps + Sig_theta_cal_inv))
  mean_post <- (1/sig2_eps) * tcrossprod(Cov_post, G) %*% y_adjusted
  
  rownames(mean_post) <- pars_cal_names
  rownames(Cov_post) <- pars_cal_names
  colnames(Cov_post) <- pars_cal_names
  true_posterior <- list(mean=mean_post, Cov=Cov_post)
  
  return(list(computer_model_data=computer_model_data, theta_prior_params=theta_prior_params, 
              true_posterior=true_posterior, random_seed=random_seed))
  
}

