#
# statistical_helper_functions.r
# General helper functions for statistical inverse problems.  
#
# Andrew Roberts
# 
# Depends: general_helper_functions.r

library(truncnorm)

# ------------------------------------------------------------------------------
# Prior Distribution object:
# At present, prior distributions are encoded by data.frames in which each 
# row is interpreted as a distribution for an independent scalar parameter.
# This data.frame is required to have column names "par_name", "dist", 
# "param1", "param2", "bound_lower", and "bound_upper". "param1" and "param2" 
# store the parameters defining a distributional family (e.g., mean and 
# standard deviation for Gaussian). In the case that these parameters also 
# encode the bounds for the parameter (e.g., uniform distribution) then they 
# must align with the values of "bound_lower" and "bound_upper". For example, 
# for a uniform prior, it must be true that "param1 == bound_lower" and 
# "param2 == bound_upper". In some cases, the lower/upper bounds are actually 
# required to define the distribution. For example, for a truncated normal prior 
# "param1" is the mean, "param2" is the standard deviation, and "param_lower"/
# "param_upper" provide the bounds defining the truncation. 
#
# Note: 
# This is clearly not a flexible approach and ought to be changed. When time 
# allows, this should be updated so that the prior distribution is a list (or 
# class) with each element corresponding to a parameter or group of parameters 
# (in the case of non-independent priors). Each element should store parameter
# bounds, optionally a function representing the log density, optionally a 
# function that draws samples, and a list of hyperparameter names/values.
# ------------------------------------------------------------------------------

get_lprior_dens <- function(par_prior) {
  # Returns a function representing the log-prior density. The returned function 
  # is vectorized so that it accepts matrix inputs where each row is a different 
  # parameter vector at which to evaluate the prior. In this case, the log-prior
  # density function returns a vector of length equal to the number of rows in 
  # the input matrix. See `calc_lprior_denssingle_input()` for requirements on 
  # `par_prior`.
  
  function(par) calc_lprior_dens(par, par_prior)
}


calc_lprior_dens <- function(par, par_prior) {
  # A wrapper around `calc_lprior_dens_single_input()` that allows computation 
  # of the log prior density at multiple input values. Here, `par` is a matrix
  # with each row being an input parameter vector, and each column representing
  # a dimension of the parameter space. The returned value is a numeric vector 
  # of length equal to the number of rows in `par`, with each entry being the 
  # corresponding log-prior density evaluation. `par` can also be a numeric 
  # vector corresponding to a single input vector. See 
  # `calc_lprior_dens_single_input()` for requirements on `par_prior`.  

  # Special case: flat prior.
  if(isTRUE(par_prior == "flat")) return(0)
  
  # If single input is passed as a numeric vector. 
  if(is.null(nrow(par))) par <- matrix(par, nrow=1)
  
  return(apply(par, 1, function(u) calc_lprior_dens_single_input(u, par_prior)))
}


calc_lprior_dens_single_input <- function(par, par_prior) {
  # Evaluates the log prior density at a single input parameter vector `par`. 
  # This function currently only works for independent priors on each parameter.
  # The definition of the prior for each input parameter is encoded in the 
  # data.frame `par_prior`. Currently accepted distributions:
  #    "Gaussian", param1=mean, param2=standard deviation
  #    "Uniform", param1=lower, param2=upper
  #    "Gamma", param1=shape, param2=rate
  #    "Beta", param1=shape1, param2=shape2
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
  #               be provided for other distributions as well.
  #
  # Returns:
  #    The prior density evaluation log p(par). Note that in certain cases the 
  #    log prior can be negative infinity; e.g. for a uniform prior where `par` 
  #    is not contained within the upper and lower bound, or if `par` falls 
  #    outside of the parameter bounds.
  
  # Special case: flat prior.
  if(isTRUE(par_prior == "flat")) return(0)
  
  # Return -infinity if parameter falls outside parameter bounds.
  if(any(par < par_prior[["bound_lower"]], na.rm=TRUE) ||
     any(par > par_prior[["bound_upper"]], na.rm=TRUE)) {
      return(-Inf)
  }

  lprior <- 0
  par_prior[, "val"] <- par
  Gaussian_priors <- par_prior[par_prior$dist == "Gaussian",]
  Uniform_priors <- par_prior[par_prior$dist == "Uniform",]
  Gamma_priors <- par_prior[par_prior$dist == "Gamma",]
  Beta_priors <- par_prior[par_prior$dist == "Beta",]
  Truncated_Gaussian_priors <- par_prior[par_prior$dist == "Truncated_Gaussian",]
  
  if(nrow(Gaussian_priors) > 0) {
    lprior <- lprior + sum(dnorm(Gaussian_priors$val, Gaussian_priors$param1, 
                                 Gaussian_priors$param2, log=TRUE))
  }
  
  if(nrow(Uniform_priors) > 0) {
    lprior <- lprior + sum(dunif(Uniform_priors$val, Uniform_priors$param1, 
                                 Uniform_priors$param2, log=TRUE))
  }
  
  if(nrow(Gamma_priors) > 0) {
    lprior <- lprior + sum(dgamma(Gamma_priors$val, shape=Gamma_priors$param1, 
                                  rate=Gamma_priors$param2, log=TRUE))
  }
  
  if(nrow(Beta_priors) > 0) {
    lprior <- lprior + sum(dbeta(Beta_priors$val, shape1=Beta_priors$param1, 
                                 shape2=Beta_priors$param2, log=TRUE))
  }
  
  if(nrow(Truncated_Gaussian_priors) > 0) {
    lprior <- lprior + sum(sapply(seq(1, nrow(Truncated_Gaussian_priors)), 
                                  function(i) log(dtruncnorm(Truncated_Gaussian_priors$val[i], sd = Truncated_Gaussian_priors$param2[i]))))
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
    } else if(par_prior[j, "dist"] =="Gamma") {
      par_samp[,j] <- rgamma(n, shape=par_prior[j, "param1"], rate=par_prior[j, "param2"])
    } else if(par_prior[j, "dist"] == "Beta") {
      par_samp[,j] <- rbeta(n, shape1=par_prior[j, "param1"], shape2=par_prior[j, "param2"])
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

get_prior_bounds <- function(par_prior, tail_prob_excluded=0.01, 
                             set_hard_bounds=FALSE) {
  # Assembles a matrix of bounds based on the prior on a parameter-by-parameter
  # basis. In general, the bounds are set to lower and upper quantiles of the 
  # distribution. If `set_hard_bounds = TRUE`, then the bounds will be
  # set to "hard bounds" of the distribution when relevant; e.g., the lower
  # bound will be set to 0 for distributions supported on the non-negative reals.
  # For the quantile case, the quantiles are determined by `tail_prob_excluded`, 
  # which is defined so that `1 - tail_prob_excluded` is the prior mass 
  # (for the univariate parameter) that is contained within the returned bounds. 
  # The returned matrix is 2 x d, where d is the number of parameters. The 
  # first and second rows correspond to the lower and upper bounds, respectively.
  
  bounds <- matrix(NA, nrow=2, ncol=nrow(par_prior))
  p <- 0.5 * tail_prob_excluded
  
  for(j in 1:ncol(bounds)) {
    dist_name <- par_prior[j, "dist"]
    if(dist_name == "Uniform") {
      if(set_hard_bounds) {
        bounds[1,j] <- par_prior[j, "param1"]
        bounds[2,j] <- par_prior[j, "param2"]
      } else {
        bounds[1,j] <- qunif(p, par_prior[j,"param1"], par_prior[j,"param2"])
        bounds[2,j] <- qunif(p, par_prior[j,"param1"], par_prior[j,"param2"], lower.tail=FALSE)
      }
    } else if(dist_name == "Gaussian") {
      bounds[1,j] <- qnorm(p, par_prior[j,"param1"], par_prior[j,"param2"])
      bounds[2,j] <- qnorm(p, par_prior[j,"param1"], par_prior[j,"param2"], lower.tail=FALSE)
    } else if(dist_name == "Gamma") {
      if(set_hard_bounds) {
        bounds[1,j] <- 0
      } else {
        bounds[1,j] <- qgamma(p, shape=par_prior[j,"param1"], rate=par_prior[j,"param2"])
      }
      bounds[2,j] <- qgamma(p, shape=par_prior[j,"param1"], rate=par_prior[j,"param2"], lower.tail=FALSE)
    } else if(dist_name == "Beta") {
      if(set_hard_bounds) {
        bounds[1,j] <- 0
        bounds[2,j] <- 1
      } else {
        bounds[1,j] <- qbeta(p, par_prior[j,"param1"], par_prior[j,"param2"])
        bounds[2,j] <- qbeta(p, par_prior[j,"param1"], par_prior[j,"param2"], lower.tail=FALSE)
      }
    } else if(dist_name == "Truncated_Gaussian") {
      if(set_hard_bounds) {
        bounds[1,j] <- par_prior[j,"bound_lower"]
        bounds[2,j] <- par_prior[j,"bound_upper"]
      } else {
        bounds[1,j] <- truncnorm::qtruncnorm(p, mean=par_prior[j,"param1"], 
                                             sd=par_prior[j,"param2"],
                                             a=par_prior[j,"bound_lower"],
                                             b=par_prior[j,"bound_upper"])
        bounds[2,j] <- truncnorm::qtruncnorm(1-p, mean=par_prior[j,"param1"], 
                                             sd=par_prior[j,"param2"],
                                             a=par_prior[j,"bound_lower"],
                                             b=par_prior[j,"bound_upper"])
      }
    } else {
      stop("Unsupported prior distribution: ", dist_name)
    }
  }
  
  colnames(bounds) <- par_prior$par_name
  return(bounds)
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
  # Converts the prior parameters on the calibration parameters (par) so that 
  # they are truncated in that they assign 0 probability mass beyond the bounds 
  # specified in `input_bounds`. Note that this function does more than simply 
  # setting lower/upper bounds on an existing prior object; it actually modifies
  # the parameter values of the prior distributions so that they become 
  # well-defined truncated versions of the distributions.
  # The truncation applied to uniform priors simply sets the bounds on the 
  # uniform prior to the bounds provided in `input_bounds`. Applied to Gaussian
  # priors, the Gaussian distributions are converted to truncated Gaussian 
  # distributions, with the truncation bounds again determined by 
  # `input_bounds`. For compactly supported priors (e.g. uniform or truncated 
  # Gaussian), the bounds are only updated if they fall outside of the bounds in 
  # `input_bounds` (e.g. an existing lower bound will not be made any lower). 
  #
  # Args:
  #    par_prior: prior distribution data.frame object.
  #    input_bounds: matrix of dimension 2 x d, where d is the number of 
  #                  parameters. The first row contains lower bounds on the 
  #                  parameters and the second row contains upper bounds.
  #                  Column names must be set to the parameter names.
  #
  # Returns:
  #    data.frame, the updated version of `par_prior`. 
  
  # Ensure the parameters in `input_bounds` are ordered as in `par_prior`.
  par_names <- par_prior[,"par_name"]
  input_bounds <- input_bounds[,par_names, drop=FALSE]
  
  for(i in seq_along(par_names)) {
    l <- input_bounds[1, i]
    u <- input_bounds[2, i]
    update_l <- isTRUE(par_prior[i, "bound_lower"] < l)
    update_u <- isTRUE(par_prior[i, "bound_upper"] > u)
    
    # Update distributions and distribution parameters.
    if(par_prior[i, "dist"] == "Gaussian") {
      par_prior[i, "dist"] <- "Truncated_Gaussian"
    } else if(par_prior[i, "dist"] == "Uniform") {
      if(update_l) par_prior[i, "param1"] <- l
      if(update_u) par_prior[i, "param2"] <- u
    } else if(par_prior[i, "dist"] == "Truncated_Gaussian") { 
      # No parameters to update, only the bounds.
    } else {
      stop("Prior distribution ", par_prior[i, "dist"], " not supported.")
    }
    
    # Update distribution bounds.
    if(update_l) par_prior[i, "bound_lower"] <- l
    if(update_u) par_prior[i, "bound_upper"] <- u
  }
  
  return(par_prior)
}


convert_Gaussian_to_LN <- function(mean_Gaussian, var_Gaussian=NULL, 
                                   cov_Gaussian=NULL, return_mean=TRUE, 
                                   return_var=TRUE, return_cov=FALSE, 
                                   log_scale=FALSE, adjustment=NULL,
                                   lower=-Inf, upper=Inf) {
  # Given the mean and either variance or covariance matrix of a Gaussian random 
  # vector X, computes the mean, variance, and covariance of Y := exp(X), which 
  # is log-normally (LN) distributed. Optionally, returns the log of these 
  # quantities, which is often recommended for numerical stability. The log 
  # scale option is not allowed when `return_cov` is TRUE, since the covariance 
  # matrix can contain negative values. In the below descriptions, let N denote 
  # the length of X. 
  #
  # Args:
  #    mean_Gaussian: numeric(N), or Nx1 matrix; the mean vector of X. 
  #    var_Gaussian: numeric(N), or Nx1 matrix; the variance of each component of 
  #                  X. Required if `cov_Gaussian` is NULL. 
  #    cov_Gaussian: NxN matrix, the covariance matrix of X. Required if 
  #                  `var_Gaussian` is NULL. Required if `return_cov` is NULL. 
  #    return_mean: logical(1), whether or not to return the (log) mean of Y. 
  #    return_var: logical(1), whether or not to return the (log) variances of Y.
  #    return_cov: logical(1), whether or not to return the covariance matrix of Y. 
  #    log_scale: logical(1), whether to return the log of the mean/variance of  
  #               Y. Must be FALSE if `return_cov` is TRUE.
  #    lower: vector of lower bounds (scalar assumed to mean constant bound).
  #    upper: vector of upper bounds (scalar assumed to mean constant bound).
  #
  # Returns: 
  #    list with the LN computations. Potential list arguments include "mean", 
  #    "var", "cov",  "log_mean", and "log_var". 
                         

  if(!is.null(adjustment) && !(adjustment %in% c("truncated", "rectified"))) {
    stop("`convert_Gaussian_to_LN` only supported truncated or rectified adjustments.")
  }
  
  if(all(is.infinite(lower)) && all(is.infinite(upper))) adjustment <- NULL

  # Truncated or rectified log-normal. Does not support multivariate - only
  # pointwise means/variances.
  if(!is.null(adjustment)) {
    if(adjustment == "truncated") {
      return(convert_Gaussian_to_trunc_LN(mean_Gaussian, var_Gaussian,
                                          return_var=return_var, lower=lower,
                                          upper=upper, log_scale=log_scale))
    } else if(adjustment == "rectified") {
      return(convert_Gaussian_to_rect_LN(mean_Gaussian, var_Gaussian,
                                         return_var=return_var, lower=lower,
                                         upper=upper, log_scale=log_scale))
    }
  }
  
  # Otherwise compute log-normal moments.
  assert_that(!is.null(var_Gaussian) || !is.null(cov_Gaussian), 
              msg=paste0("Gaussian variance or covariance matrix required to",
                         "compute log-normal moments."))
  assert_that(!(log_scale && return_cov), 
              msg="Cannot return LN moments on log scale if `return_cov` is TRUE.")
  if(is.null(var_Gaussian)) var_Gaussian <- diag(cov_Gaussian)

  # Compute the log of the LN mean and variance. Note that this can't be done 
  # for the  covariance given the that covariance matrix can have negative values. 
  log_mean <- drop(mean_Gaussian) + 0.5*drop(var_Gaussian)
  if(return_var) log_var <- log_exp_minus_1(var_Gaussian) + 2*log_mean
  if(log_scale) {
    return_list <- list(log_mean=log_mean)
    if(return_var) return_list$log_var <- log_var
    return(return_list)
  }
  
  # Continue if LN moments are requested on the original scale. 
  # The mean is always returned.
  return_list <- list()
  if(return_mean) return_list$mean <- exp(log_mean)
  if(return_var) return_list$var <- exp(log_var)
  if(return_cov) {
    assert_that(!is.null(return_cov), 
                msg="Computing cov requires non-NULL `cov_Gaussian`.")
    return_list$cov <- exp(outer(log_mean, log_mean, FUN="+")) * (exp(cov_Gaussian)-1)
  }
  
  return(return_list)
}


convert_Gaussian_to_trunc_LN <- function(mean_Gaussian, var_Gaussian,
                                         return_var=TRUE, 
                                         lower=-Inf, upper=Inf, 
                                         log_scale=FALSE,
                                         return_intermediate_calcs=FALSE) {
  # Converts (univariate) Gaussian means and variances to rectified log-normal 
  # means and  variances. If x ~ N(m, v) then y ~ rect-LN(m, v; b1, b2) with
  # respect to bounds (b1, b2) if y = exp(b1) when x < b1, y = exp(b2) when 
  # y = exp(b2) when x > b2, and y = exp(x) when b1 <= x <= b2. Computing the 
  # rectified  moments first requires computing the truncated log-normal moments. 
  # Note that these moments are not the same as first computing the 
  # rectified Gaussian expectations and then exponentiating. Also supports
  # one-sided bounds, by setting one of the elements in `bounds` to -Inf
  # of Inf. Vectorized so that `mean_Gaussian` and `var_Gaussian` can be
  # vectors of equal length; the same bounds will be applied to all entries.
  # `return_intermediate_calcs` returns intermediate computations that are
  # used when computing the rectified LN moments; see 
  # `convert_Gaussian_to_rect_LN`. The bounds can be vectors if the bounds
  # differ by entry, or they can be scalars for uniform bounds.

  mean_Gaussian <- drop(mean_Gaussian)
  var_Gaussian <- drop(var_Gaussian)
  assert_that(length(mean_Gaussian) == length(var_Gaussian))
  n <- length(mean_Gaussian)
  
  return_list <- list()
  
  # First compute the (untruncated) log-normal mean.
  ln_moments <- convert_Gaussian_to_LN(mean_Gaussian, var_Gaussian, 
                                       adjustment=NULL, log_scale=TRUE,
                                       return_mean=TRUE, return_var=FALSE)
  
  # Upper and lower bounds.
  b1 <- drop(lower)
  b2 <- drop(upper)
  if(length(b1) == 1L) b1 <- rep(b1, n)
  if(length(b2) == 1L) b2 <- rep(b2, n)
  if(length(b1) != n) stop("Lower bound must be length of `mean_Gaussian` or length 1.")
  if(length(b2) != n) stop("Upper bound must be length of `mean_Gaussian` or length 1.")
  
  constrain_lower <- is.finite(b1)
  constrain_upper <- is.finite(b2)
  
  if(return_intermediate_calcs) {
    return_list$constrain_lower <- constrain_lower
    return_list$constrain_upper <- constrain_upper
  }
  
  # Probabilities involved in moment expressions.
  lprob_leq_upper <- lprob_num_upper <- rep(0, n)
  lprob_leq_upper[constrain_upper] <- pnorm(b2[constrain_upper], 
                                            mean_Gaussian[constrain_upper], 
                                            sqrt(var_Gaussian[constrain_upper]), log.p=TRUE)
  lprob_num_upper[constrain_upper] <- pnorm(b2[constrain_upper], 
                                            mean_Gaussian[constrain_upper] + var_Gaussian[constrain_upper], 
                                            sqrt(var_Gaussian[constrain_upper]), log.p=TRUE)
  
  lprob_leq_lower <- lprob_num_lower <- rep(-Inf, n)
  lprob_leq_lower[constrain_lower] <- pnorm(b1[constrain_lower], 
                                            mean_Gaussian[constrain_lower], 
                                            sqrt(var_Gaussian[constrain_lower]), log.p=TRUE)
  lprob_num_lower[constrain_lower] <- pnorm(b1[constrain_lower], 
                                            mean_Gaussian[constrain_lower] + var_Gaussian[constrain_lower], 
                                            sqrt(var_Gaussian[constrain_lower]), log.p=TRUE)

  if(return_intermediate_calcs) {
    return_list$lprob_leq_upper <- lprob_leq_upper
    return_list$lprob_leq_lower <- lprob_leq_lower
  }
  
  # Compute truncated LN mean.
  log_diff_num <- rep(NA_real_, n)
  lprob_num_lower_inf <- is.infinite(lprob_num_lower)
  log_diff_num[lprob_num_lower_inf] <- lprob_num_upper[lprob_num_lower_inf]
  log_diff_num[!lprob_num_lower_inf] <- log_diff_exp(lprob_num_upper[!lprob_num_lower_inf], 
                                                     lprob_num_lower[!lprob_num_lower_inf])

  log_diff_denom <- rep(NA_real_, n)
  lprob_leq_lower_inf <- is.infinite(lprob_leq_lower)
  log_diff_denom[lprob_leq_lower_inf] <- lprob_leq_upper[lprob_leq_lower_inf]
  log_diff_denom[!lprob_leq_lower_inf] <- log_diff_exp(lprob_leq_upper[!lprob_leq_lower_inf], 
                                                       lprob_leq_lower[!lprob_leq_lower_inf])
  
  ln_trunc_mean <- ln_moments$log_mean + log_diff_num - log_diff_denom
                   
  if(log_scale) return_list$log_mean <- ln_trunc_mean
  else return_list$mean <- exp(ln_trunc_mean)
  
  if(!return_var) return(return_list)
  
  # Computing (log of) truncated LN second moment: E[y^2].
  lprob_num_upper_y2 <- rep(0, n)
  lprob_num_upper_y2[constrain_upper] <- pnorm(b2[constrain_upper], 
                                               mean_Gaussian[constrain_upper] + 2*var_Gaussian[constrain_upper], 
                                               sqrt(var_Gaussian[constrain_upper]), log.p=TRUE)

  lprob_num_lower_y2 <- rep(-Inf, n)
  lprob_num_lower_y2[constrain_lower] <- pnorm(b1[constrain_lower], 
                                               mean_Gaussian[constrain_lower] + 2*var_Gaussian[constrain_lower], 
                                               sqrt(var_Gaussian[constrain_lower]), log.p=TRUE)

  log_diff_num_y2 <- rep(NA_real_, n)
  lprob_num_lower_y2_inf <- is.infinite(lprob_num_lower_y2)
  log_diff_num_y2[lprob_num_lower_y2_inf] <- lprob_num_upper_y2[lprob_num_lower_y2_inf]
  log_diff_num_y2[!lprob_num_lower_y2_inf] <- log_diff_exp(lprob_num_upper_y2[!lprob_num_lower_y2_inf], 
                                                           lprob_num_lower_y2[!lprob_num_lower_y2_inf])

  log_Ey2 <- log_diff_num_y2 + 2 * (mean_Gaussian + var_Gaussian) - log_diff_denom
  
  if(return_intermediate_calcs) return_list$log_Ey2 <- log_Ey2
  
  # Use mean and second moment to compute variance.
  ln_trunc_var <- log_diff_exp(log_Ey2, 2*ln_trunc_mean)

  if(log_scale) return_list$log_var <- ln_trunc_var
  else return_list$var <- exp(ln_trunc_var)
  
  return(return_list)
}


convert_Gaussian_to_rect_LN <- function(mean_Gaussian, var_Gaussian,
                                        return_var=TRUE, lower=-Inf, upper=Inf,
                                        log_scale=FALSE) {
  # Converts Gaussian means and variances to truncated log-normal means and 
  # variances. If x ~ N(m, v) then y ~ trunc-LN(m, v; b1, b2) with respect to
  # bounds (b1, b2) is the random variable equal in distribution to
  # exp(x) | b1 <= x <= b2. Note that the truncated LN moments are not the same 
  # as exponentiating the truncated Gaussian moments. Also supports
  # one-sided bounds, by setting one of the elements in `bounds` to -Inf
  # of Inf. Vectorized so that `mean_Gaussian` and `var_Gaussian` can be
  # vectors of equal length; the same bounds will be applied to all entries.
  #
  # TODO: this function needs more testing.

  return_list <- list()
  
  # First compute truncated moments.
  ln_trunc <- convert_Gaussian_to_trunc_LN(mean_Gaussian, var_Gaussian,
                                           return_var=return_var,
                                           lower=lower, upper=upper,
                                           log_scale=TRUE,
                                           return_intermediate_calcs=TRUE)
  n <- length(ln_trunc$log_mean)
  
  # Lower and upper bounds (potentially infinite).
  b1 <- drop(lower)
  b2 <- drop(upper)
  if(length(b1) == 1L) b1 <- rep(b1, n)
  if(length(b2) == 1L) b2 <- rep(b2, n)
  if(length(b1) != n) stop("Lower bound must be length of `mean_Gaussian` or length 1.")
  if(length(b2) != n) stop("Upper bound must be length of `mean_Gaussian` or length 1.")

  # Probabilities for the three pieces of the rectified random variable. Some
  # values may be -Inf here, which will imply zero probability.
  # prob_upper = 1 - prob_leq_upper
  # prob_middle = 1 - prob_leq_upper - prob_leq_lower
  eps <- 1e-16
  upper_prob_zero <- (abs(ln_trunc$lprob_leq_upper) < eps) # i.e., prob <= upper bound is 1.
  upper_prob_one <- is.infinite(ln_trunc$lprob_leq_upper)
  upper_prob_not_zero_one <- !(upper_prob_zero | upper_prob_one)

  # log P(y > b2)
  lprob_upper <- rep(NA_real_, n)
  lprob_upper[upper_prob_zero] <- -Inf
  lprob_upper[upper_prob_one] <- 0
  lprob_upper[upper_prob_not_zero_one] <- log_1_minus_exp(ln_trunc$lprob_leq_upper[upper_prob_not_zero_one])
    
  # log P(y < b1)
  lprob_lower <- ln_trunc$lprob_leq_lower
  lower_prob_zero <- is.infinite(lprob_lower)
  lower_prob_one <- (abs(lprob_lower) < eps)
  lower_prob_not_zero_one <- !(lower_prob_zero | lower_prob_one)
  
  # log P(b1 <= y <= b2)
  lprob_middle <- rep(NA_real_, n)
  lprob_middle[lower_prob_one | upper_prob_one] <- -Inf
  lprob_middle[lower_prob_zero & upper_prob_zero] <- 0
  lprob_middle[lower_prob_zero & upper_prob_not_zero_one] <- log_1_minus_exp(lprob_upper[lower_prob_zero & upper_prob_not_zero_one])
  lprob_middle[upper_prob_zero & lower_prob_not_zero_one] <- log_1_minus_exp(lprob_lower[upper_prob_zero & lower_prob_not_zero_one])
  idx_temp <- (lower_prob_not_zero_one & upper_prob_not_zero_one)
  lprob_middle[idx_temp] <- log_1_minus_exp(matrixStats::logSumExp(lprob_lower[idx_temp], lprob_upper[idx_temp]))

  # Compute rectified LN mean. Note lprob_middle can contain -Inf.
  log_summands <- cbind(ln_trunc$log_mean + lprob_middle)
  constrain_lower <- ln_trunc$constrain_lower
  constrain_upper <- ln_trunc$constrain_upper
  add_lower <- add_upper <- rep(0, n)
  add_lower[constrain_lower] <- b1[constrain_lower] + lprob_lower[constrain_lower]
  add_upper[constrain_upper] <- b2[constrain_upper] + lprob_upper[constrain_upper]
  log_summands <- cbind(log_summands, add_lower)
  log_summands <- cbind(log_summands, add_upper)
  
  # Note: `matrixStats::rowLogSumExps` can handle the presence of -Inf values.
  log_rect_ln_mean <- matrixStats::rowLogSumExps(log_summands)
  if(log_scale) return_list$log_mean <- log_rect_ln_mean
  else return_list$mean <- exp(log_rect_ln_mean)
  
  if(!return_var) return(return_list)
  
  # Computing rectified LN second moment.
  log_summands_y2 <- cbind(ln_trunc$log_Ey2 + lprob_middle)
  add_lower_y2 <- add_upper_y2 <- rep(0, n)
  add_lower_y2[constrain_lower] <- 2*b1[constrain_lower] + lprob_lower[constrain_lower]
  add_upper_y2[constrain_upper] <- 2*b2[constrain_upper] + lprob_upper[constrain_upper]
  log_summands_y2 <- cbind(log_summands_y2, add_lower_y2)
  log_summands_y2 <- cbind(log_summands_y2, add_upper_y2)
  
  # Very small values can cause issues in `log_diff_exp`. Currently just 
  # thresholding as a hack.
  log_rect_Ey2 <- matrixStats::rowLogSumExps(log_summands_y2)
  log_rect_Ey2[log_rect_Ey2 < 1e-8] <- 1e-8
  neg_var_idx <- (log_rect_Ey2 < 2*log_rect_ln_mean)
  log_rect_ln_var <- rep(NA, n)
  log_rect_ln_var[!neg_var_idx] <- log_diff_exp(log_rect_Ey2[!neg_var_idx], 
                                                2*log_rect_ln_mean[!neg_var_idx])
  log_rect_ln_var[neg_var_idx] <- -70

  if(log_scale) return_list$log_var <- log_rect_ln_var
  else return_list$var <- exp(log_rect_ln_var)

  return(return_list)
}


rect_norm_quantile <- function(q, mean=0, sd=1, lower=-Inf, upper=Inf, lower.tail=TRUE) {
  # NOTE: this function was written with assistance from ChatGPT. It has been
  # numerically tested for correctness.
  #
  # Computes a quantile of a rectified normal distribution; i.e., a Gaussian
  # that has been "clipped" at potentially lower and upper bounds, thus
  # creating point masses at the bound(s) (not to be confused with the truncated
  # Gaussian). `mean`, `sd` are the moments of the underlying Gaussian.
  # `mean`, `sd`, `lower`, `upper` can be vectors of the same length, or 
  # a mixture of vectors and scalars. `q` is a probability in [0,1].
  
  # Ensure q is a vector
  q <- as.numeric(q)
  if(!lower.tail) lower.tail <- 1 - q
  
  # Determine target length
  arg_lengths <- c(length(mean), length(sd), length(lower), length(upper))
  target_len <- max(arg_lengths)
  
  # Recycle scalars
  recycle_if_scalar <- function(x) {
    if (length(x) == 1) rep(x, target_len)
    else if (length(x) == target_len) x
    else stop("mean, sd, lower, and upper meanst be scalars or vectors of the same length.")
  }
  
  mean    <- recycle_if_scalar(mean)
  sd <- recycle_if_scalar(sd)
  lower <- recycle_if_scalar(lower)
  upper <- recycle_if_scalar(upper)
  
  if (length(q) != target_len) {
    if (length(q) == 1) {
      q <- rep(q, target_len)
    } else if (length(q) != target_len) {
      stop("q meanst be a scalar or have the same length as other vector arguments.")
    }
  }
  
  # Compute cumulative probabilities
  p_lower <- pnorm(lower, mean=mean, sd=sd)
  p_upper <- pnorm(upper, mean=mean, sd=sd)
  
  mass_lower <- p_lower
  mass_upper <- 1 - p_upper
  mass_middle <- p_upper - p_lower
  
  result <- numeric(target_len)
  
  # Case 1: q <= mass_lower → return lower
  result[q <= mass_lower] <- lower[q <= mass_lower]
  
  # Case 2: q >= 1 - mass_upper → return upper
  result[q >= 1 - mass_upper] <- upper[q >= 1 - mass_upper]
  
  # Case 3: middle region → inverse CDF
  in_middle <- (q > mass_lower) & (q < 1 - mass_upper)
  q_adj <- (q[in_middle] - mass_lower[in_middle]) / mass_middle[in_middle]
  target_p <- p_lower[in_middle] + q_adj * mass_middle[in_middle]
  result[in_middle] <- qnorm(target_p, mean=mean[in_middle], sd=sd[in_middle])
  
  return(result)
}


rect_lnorm_quantile <- function(q, meanlog, sdlog, lower, upper, log_scale=TRUE) {
  # NOTE: this function was written with assistance from ChatGPT. It has been
  # numerically tested for correctness.
  
  # Ensure q, meanlog, sdlog, lower, upper are same length (or scalars)
  args <- list(q, meanlog, sdlog, lower, upper)
  n <- max(lengths(args))
  args <- lapply(args, function(x) if (length(x) == 1) rep(x, n) else x)
  if (!all(lengths(args) == n)) stop("Arguments must be scalars or vectors of the same length.")
  
  q <- args[[1]]; meanlog <- args[[2]]; sdlog <- args[[3]]
  lower <- args[[4]]; upper <- args[[5]]
  
  # Compute bounds in quantile space
  p_lower <- pnorm(lower, meanlog, sdlog)
  p_upper <- pnorm(upper, meanlog, sdlog)
  
  # Output vector
  result <- numeric(n)
  
  for (i in seq_len(n)) {
    if (q[i] <= p_lower[i]) {
      result[i] <- exp(lower[i])
    } else if (q[i] >= p_upper[i]) {
      result[i] <- exp(upper[i])
    } else {
      # Rescale quantile to truncated normal quantile range
      q_rescaled <- (q[i] - p_lower[i]) / (p_upper[i] - p_lower[i])
      result[i] <- EnvStats::qlnormTrunc(q, meanlog=meanlog[i],
                                         sdlog=sdlog[i],
                                         min=exp(lower[i]),
                                         max=exp(upper[i]))
    }
  }
  
  if(log_scale) return(log(result))
  return(result)
}


weighted_quantile_logw <- function(x, log_w, probs = c(0.25, 0.5, 0.75)) {
  # Computes weighted empirical quantiles, where the weights are potentially
  # unnormalized, and on the log scale. It is assumed that the weights may
  # have a large range, so the LogSumExp trick is used to exponentiate and 
  # normalize them in a numerically stable fashion.

  stopifnot(length(x) == length(log_w))
  
  # Use log-sum-exp trick for stability
  log_w_max <- max(log_w)
  w <- exp(log_w - log_w_max)
  w <- w / sum(w)  # now normalized
  
  # Sort x and weights
  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- w[ord]
  cum_w <- cumsum(w_sorted)
  
  # Find quantiles
  sapply(probs, function(p) {
    x_sorted[which(cum_w >= p)[1]]
  })
}


gen_lin_Gaus_NIW_test_data <- function(G_list, Sig_eps=NULL, mu0=NULL, Sig0=NULL, 
                                       IW_scale=NULL, IW_dof=NULL, N_missing_obs=NULL) {
  # TODO: write up blog post on Normal Inverse Wishart model and then write this function.  
  .NotYetImplemented()
}


generate_linear_Gaussian_test_data <- function(random_seed, N_obs, D, Sig_theta, G, 
                                               sig2_eps=NULL, sig_eps_frac=0.1, 
                                               pars_cal_sel=NULL) {
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

