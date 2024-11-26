#
# statistical_helper_functions.r
# General helper functions for statistical inverse problems.  
#
# Andrew Roberts
# 

get_lprior_dens <- function(par_prior, check_bounds=FALSE) {
  # Returns a function representing the log-prior density. The returned function 
  # is vectorized so that it accepts matrix inputs where each row is a different 
  # parameter vector at which to evaluate the prior. In this case, the log-prior
  # density function returns a vector of length equal to the number of rows in 
  # the input matrix. See `calc_lprior_denssingle_input()` for requirements on 
  # `par_prior` and `check_bounds`.  
  
  function(par) calc_lprior_dens(par, par_prior, check_bounds)
}


calc_lprior_dens <- function(par, par_prior, check_bounds=FALSE) {
  # A wrapper around `calc_lprior_dens_single_input()` that allows computation 
  # of the log prior density at multiple input values. Here, `par` is a matrix
  # with each row being an input parameter vector, and each column representing
  # a dimension of the parameter space. The returned value is a numeric vector 
  # of length equal to the number of rows in `par`, with each entry being the 
  # corresponding log-prior density evaluation. `par` can also be a numeric 
  # vector corresponding to a single input vector. See 
  # `calc_lprior_dens_single_input()` for requirements on `par_prior` and 
  # `check_bounds`.  

  # If single input is passed as a numeric vector. 
  if(is.null(nrow(par))) par <- matrix(par, nrow=1)
  
  return(apply(par, 1, function(u) calc_lprior_dens_single_input(u, par_prior, check_bounds)))
  
}


calc_lprior_dens_single_input <- function(par, par_prior, check_bounds=FALSE) {
  # Evaluates the log prior density at a single input parameter vectir `par`. 
  # This function currently only works for independent priors on each parameter.
  # The definition of the prior for each input parameter is encoded in the 
  # data.frame `par_prior`. Currently accepted distributions:
  #    "Gaussian", param1=mean, param2=standard deviation
  #    "Uniform", param1=lower, param2=upper
  #    "Truncated_Gaussian", param1=mean, param2=standard deviation (the moments
  #     of the Gaussian distribution inducing the truncated Gaussian), 
  #     bound_lower=lower truncation value, bound_upper=upper truncation value.
  #
  # Args:
  #    par: numeric, the vector at which to evaluate the log prior density.
  #    par_prior: data.frame, with columns "dist", "param1", and "param2". The 
  #               ith row of the data.frame should correspond ot the ith entry 
  #               of `par`. `par_prior` can also optionally include the columns 
  #               "bound_lower" and "bound_upper". Some distributions (e.g., 
  #               truncated Gaussian) require these parameters, but they can 
  #               be provided for other distributions as well. See `check_bounds`
  #               for more details.
  #    check_bounds: logical(1), if TRUE checks if the parameter `par` lies 
  #                  within the bounds defined in the columns of `par_prior`
  #                  called "bound_lower" and "bound_upper". If these columns do 
  #                  not exist or contain NA values then there is not effect.
  #                  Otherwise, if `par` does not lie within the bounds, then 
  #                  the function returns -Inf. If `check_bounds` is FALSE, no 
  #                  check is performed.
  #
  # Returns:
  #    The prior density evaluation log p(par). Note that in certain cases the 
  #    log prior can be negative infinity; e.g. for a uniform prior where `par` 
  #    is not contained within the upper and lower bound, or if `check_bounds`
  #    is enforced.
  
  if(check_bounds) {
    if(any(par < par_prior[["bound_lower"]], na.rm=TRUE) ||
       any(par > par_prior[["bound_upper"]], na.rm=TRUE)) {
      return(-Inf)
    }
  }
  
  lprior <- 0
  par_prior[, "val"] <- par
  Gaussian_priors <- par_prior[par_prior$dist == "Gaussian",]
  Uniform_priors <- par_prior[par_prior$dist == "Uniform",]
  Truncated_Gaussian_priors <- par_prior[par_prior$dist == "Truncated_Gaussian",]
  
  if(nrow(Gaussian_priors) > 0) {
    lprior <- lprior + sum(dnorm(Gaussian_priors$val, Gaussian_priors$param1, 
                                 Gaussian_priors$param2, log=TRUE))
  }
  
  if(nrow(Uniform_priors) > 0) {
    lprior <- lprior + sum(dunif(Uniform_priors$val, Uniform_priors$param1, 
                                 Uniform_priors$param2, log=TRUE))
  }
  
  if(nrow(Truncated_Gaussian_priors) > 0) {
    lprior <- lprior + sum(sapply(seq(1, nrow(Truncated_Gaussian_priors)), 
                                  function(i) log(dtruncnorm(Truncated_Gaussian_priors$val[i],                                                                                                    sd = Truncated_Gaussian_priors$param2[i]))))
  }
  
  return(lprior)  
}


get_prior_sampler <- function(par_prior) {
  # Returns a function with no arguments that returns a single sample from the 
  # prior distribution encoded by `par_prior`.
  function(n=1L, ...) sample_prior(par_prior, n=n)
}


sample_prior <- function(par_prior, n=1L) {
  # Return independent samples from prior distribution defined by `par_prior`. 
  # Currently only supports priors that assume prior independence across 
  # parameters. 
  #
  # Args:
  #    par_prior: See `calc_lprior_dens_single_input()` for the requirements on 
  #               `par_prior`.
  #    n: integer, the number of samples to draw from the prior.
  #
  # Returns:
  # matrix, with number of rows equal to `n`, the number of samples, and 
  # number of columns equal to the number of parameters. The rows contain 
  # independent samples from the prior. Column names are set to parameter 
  # names.
  
  n_pars <- nrow(par_prior)
  par_samp <- matrix(nrow=n, ncol=n_pars)

  for(j in 1:n_pars) {
    if(par_prior[j, "dist"]=="Gaussian") {
      par_samp[,j] <- rnorm(n, par_prior[j, "param1"], par_prior[j, "param2"])
    } else if(par_prior[j, "dist"]=="Uniform") {
      par_samp[,j] <- runif(n, par_prior[j, "param1"], par_prior[j, "param2"])
    } else if(par_prior[j, "dist"]=="Truncated_Gaussian") {
      par_samp[,j] <- truncnorm::rtruncnorm(n, a=par_prior[j, "bound_lower"], 
                                            b=par_prior[j, "bound_upper"], 
                                            mean=par_prior[j, "param1"], 
                                            sd=par_prior[j, "param2"])
    } else {
      stop("Prior distribution ", par_prior[j, "dist"], " not supported.")
    }
  }
  
  colnames(par_samp) <- par_prior$par_name
  return(par_samp)  
}


plot_prior_samp <- function(par_prior, n=10000L) {
  # Samples from the prior and then returns a list of histograms, summarizing 
  # the marginal distribution of each parameter. Currently only supports 
  # independent priors for each parameter.
  #
  # Args:
  #    par_prior: See `calc_lprior_dens_single_input()` for the requirements on 
  #               `par_prior`.
  #    n: integer, the number of samples to draw from the prior.
  #
  # Returns:
  #    list of ggplot objects, of length equal to the number of parameters.
  
  samp <- sample_prior(par_prior, n=n)
  plot_list <- setNames(vector(mode="list", length=ncol(samp)),
                        colnames(samp))
  
  for(par_name in names(plot_list)) {
    x <- sym(par_name)
    plot_list[[par_name]] <- ggplot(data.frame(x=samp[,par_name])) + 
                             geom_histogram(aes(x=x)) + 
                             labs(title=paste0("Prior Samples: ", par_name),
                                  y=par_name)
  }
  
  return(plot_list)
}


truncate_prior <- function(par_prior, input_bounds) {
  # Converts the prior parameters on the calibration parameters (par) so that they are truncated
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
  #    par_prior: data.frame, with columns "dist", "param1", and "param2". The ith row of the data.frame
  #                        should correspond to the ith entry of 'par'. Currently, accepted values of "dist" are 
  #                        "Gaussian" (param1 = mean, param2 = std dev) and "Uniform" (param1 = lower, param2 = upper).
  #    input_bounds: matrix of dimension 2 x d, where d is the number of parameters. The first row contains lower bounds 
  #                  on the parameters and the second row contains upper bounds (these are often determined by the 
  #                  extent of the design points). If provided, the prior sample will be required to satisfy the lower 
  #                  and upper bounds, so sampled parameters that do not satisfy the constraints will be rejected until 
  #                  the constraint is satisfied. Note that this implies that the resulting samples are from a truncated 
  #                  version of the prior distribution, rather than the prior itself. The columns of `input_bounds` must be 
  #                  sorted in the same order as the rows of `par_prior.` Default is NULL, which imposes 
  #                  no constraints. 
  #
  # Returns:
  #    data.frame, the updated version of `par_prior`. 
  
  for(i in seq(1, nrow(par_prior))) {
    l <- input_bounds[1, i]
    u <- input_bounds[2, i]
    
    if(par_prior[i, "dist"] == "Gaussian") {
      par_prior[i, "dist"] <- "Truncated_Gaussian"
      par_prior[i, "bound_lower"] <- l
      par_prior[i, "bound_upper"] <- u
    } else if(par_prior[i, "dist"] == "Uniform") {
      if(par_prior[i, "param1"] < l) par_prior[i, "param1"] <- l
      if(par_prior[i, "param2"] > u) par_prior[i, "param2"] <- u
    } else if(par_prior[i, "dist"] == "Truncated_Gaussian") { 
      if(par_prior[i, "bound_lower"] < l) par_prior[i, "bound_lower"] <- l
      if(par_prior[i, "bound_upper"] > u) par_prior[i, "bound_upper"] <- u
    } else {
      stop("Prior distribution ", par_prior[i, "dist"], " not supported.")
    }
  }
  
  return(par_prior)
}


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

