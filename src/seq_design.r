#
# seq_design.r
#
# General functions related to sequential design and optimization. While 
# `seq_design_gp.r` concerns sequential design for Gaussian processes (GPs)
# within the acquisition function minimization framework, and 
# `seq_design_for_post_approx.r` considers the same framework but specialized 
# to the used of log-likelihood emulators for posterior approximation, this 
# file contains functions related to experimental design which may not 
# fall into the acquisition function framework and may not involve GPs at all. 
# In particular, this includes "model free" space-filling methods such as 
# Latin Hypercube sampling (LHS) which are often used to generate an initial 
# design that will subsequently be refined by sequential methods. Also 
# included are functions tailored to solving Bayesian inverse problems 
# such as Ensemble Kalman Inversion (EKI). This also includes a class 
# of methods which construct a design by sub-sampling a set of candidate 
# points (e.g. support points and Stein points). 
# 
# Andrew Roberts
# 


get_batch_design <- function(method, N_batch, prior_params=NULL, design_candidates=NULL, 
                             design_candidate_weights=NULL, ...) {
  
  if(method == "LHS") return(get_LHS_sample(N_batch, prior_dist_info=prior_params, ...))
  else if(method == "EKI_finite_time") return(run_EKI_finite_time(N_steps))
  else if(method == "tensor_product_grid") return(get_tensor_product_grid(N_batch, prior_dist_info=prior_params, ...))
  else stop("Design method ", method, " not supported.")
}


get_LHS_sample <- function(N_batch, prior_dist_info=NULL, bounds=NULL, order_1d=FALSE, ...) {
  # Produces a Latin Hypercube Sample (LHS) with `N_batch` points. The LHS is first sampled
  # uniformly in the unit hypercube using `lhs` package, then an inverse cumulative 
  # distribution function (CDF) tranformation is applied corresponding to the probability 
  # distributions specified in `prior_dist_info`. The argument `bounds` also allows 
  # bounds to be enforced on the sample space. `bounds` is interpreted in the context of 
  # the specific probability distributions specified; see details below. 
  #
  # Args:
  #    N_batch: integer(1), the number of points to sample. 
  #    prior_dist_info: data.frame with columns "dist", "param1", "param2" containing one
  #                     row per input dimension. Each dimension is considered independent. 
  #                     Currently supports "Uniform" with the parameters interpreted as the
  #                     lower/upper bounds, respectively; and "Gaussian", with the parameters
  #                     interpreted as the mean/sd. `prior_dist_info` must be non-NULL if 
  #                     `bounds` is NULL.  
  #    bounds: matrix of dimension (2,input dimension). The first row contains lower
  #            bounds for each parameter, and the second row contains the upper bounds. 
  #            If `prior_dist_info` is NULL, then `bounds` must be provided and is 
  #            interpreted as defining the bounds for a hyperrectangle which will be sampled
  #            uniformly. If `prior_dist_info` is provided, then `bounds` is applied 
  #            differently depending on the distribution. For "Uniform", `bounds` may 
  #            truncate the lower/upper parameters, but is not allowed to enlarge them. 
  #            For "Gaussian", `bounds` is used to convert the Gaussian distribution into 
  #            a truncated Gaussian with the specified bounds (-Inf/Inf is allowed for 
  #            one-sided truncation). 
  #    order_1d: If TRUE, sorts the LHS in increasing order. This is only used when the 
  #              input dimension is 1. 
  #    
  # Returns: 
  #    matrix of dimension N_batch x input dimension. Each row contains a sampled point. 

  # Either `prior_dist_info` or `bounds` must be provided. 
  assert_that(!is.null(prior_dist_info) || !is.null(bounds))

  # The dimension of the input space.
  dim_input <- nrow(prior_dist_info)
  
  # If `prior_dist_info` is NULL, samples points uniformly in the region 
  # defined by `bounds`.
  if(is.null(prior_dist_info)) {
    prior_dist_info <- data.frame(dist="Uniform", param1=bounds[1,], param2=bounds[2,])
  }
  
  # If `bounds` is NULL, set to NA. No bounds will be enforced beyond those
  # implied by the prior distributions. For distributions with unbounded support. 
  if(is.null(bounds)) bounds <- matrix(NA, nrow=2, ncol=dim_input)
  
  # Generate LHS design on unit hypercube. 
  X_lhs <- lhs::randomLHS(N_batch, dim_input)
  
  # Apply inverse CDF transform using prior distributions.
  for(j in seq_len(dim_input)) {
    dist_name <- prior_dist_info[j,"dist"]
    if(dist_name == "Uniform") {
      lower_bound <- ifelse(is.na(bounds[1,j]), prior_dist_info[j,"param1"], 
                            max(bounds[1,j], prior_dist_info[j,"param1"]))
      upper_bound <- ifelse(is.na(bounds[2,j]), prior_dist_info[j,"param2"], 
                            min(bounds[2,j], prior_dist_info[j,"param2"]))
      X_lhs[,j] <- qunif(X_lhs[,j], lower_bound, upper_bound)
    } else if(dist_name == "Gaussian") {
      lower_bound <- ifelse(is.na(bounds[1,j]), -Inf, bounds[1,j])
      upper_bound <- ifelse(is.na(bounds[2,j]), Inf, bounds[2,j])
      X_lhs[,j] <- truncnorm::qtruncnorm(X_lhs[,j], a=lower_bound, b=upper_bound, 
                                         mean=prior_dist_info[j,"param1"], sd=prior_dist_info[j, "param2"])
    } else {
      stop("Unsupported prior distribution: ", dist_name)
    }
  }
  
  # For 1 dimensional data, optionally order samples in increasing order. 
  if(order_1d && (dim_input == 1)) X_lhs <- X_lhs[order(X_lhs),,drop=FALSE]
    
  colnames(X_lhs) <- rownames(prior_dist_info)
  return(X_lhs)
}


get_tensor_product_grid <- function(N_batch, prior_dist_info=NULL, bounds=NULL, 
                                    tail_prob_excluded=0.01, ...) {
  # Generates a grid of points in d-dimensions, constructed as the Cartesian 
  # product of d one-dimensional evenly-spaced grids. The boundaries for the 
  # grid are specified either directly by `bounds`, or are inferred from 
  # probability distributions specified in `prior_dist_info`. The latter 
  # currently supports "Uniform" or "Gaussian". The former sets the bounds
  # based on the upper/lower bounds of the Uniform distribution. The latter 
  # sets the bounds at quantiles of the Gaussian distribution specified by 
  # `tail_prob_excluded`. Exactly one of `prior_dist_info` or `bounds` should 
  # be non-NULL. 
  #
  # Args
  #    N_batch: integer(1), the number of points to sample. 
  #    prior_dist_info: data.frame with columns "dist", "param1", "param2" containing one
  #                     row per input dimension (d rows). Each dimension is considered  
  #                     independent. See above description for how the grid bounds are 
  #                     derived from this argument. 
  #    bounds: matrix of dimension (2,d). The first row contains lower
  #            bounds for each parameter, and the second row contains the upper bounds. 
  #
  # Returns:
  #    matrix, of dimension N_batch x d containing the grid points as rows. 
  
  # Exactly one of `prior_dist_info` or `bounds` must be non-NULL.  
  assert_that(xor(is.null(prior_dist_info), is.null(bounds)))
  
  # If `prior_dist_info` is provided, use it to compute bounds.
  if(is.null(bounds)) {
    bounds <- matrix(NA, nrow=2, ncol=nrow(prior_dist_info))
    for(j in 1:ncol(bounds)) {
      dist_name <- prior_dist_info[j, "dist"]
      if(dist_name == "Uniform") {
        bounds[1,j] <- prior_dist_info[j, "param1"]
        bounds[2,j] <- prior_dist_info[j, "param2"]
      } else if(dist_name == "Gaussian") {
        bounds[1,j] <- qnorm(tail_prob_excluded/2, prior_dist_info[j,"param1"], 
                             prior_dist_info[j,"param2"])
        bounds[2,j] <- qnorm(tail_prob_excluded/2, prior_dist_info[j,"param1"], 
                             prior_dist_info[j,"param2"], lower.tail=FALSE)
      } else {
        stop("Unsupported prior distribution: ", dist_name)
      }
    }
  }
  
  # The dimension of the input space.
  dim_input <- ncol(bounds)

  # Number of points marginally for each dimension.
  N_per_dim <- N_batch^(1/dim_input)
  if(round(N_per_dim) != N_per_dim) stop("N_batch must have an integer d^th root; d = ", dim_input)
  
  # Create the 1d grids.  
  X_marginals <- matrix(nrow=N_per_dim, ncol=dim_input)
  for(j in seq_len(dim_input)) X_marginals[,j] <- seq(bounds[1,j], bounds[2,j], length.out=N_per_dim)
  
  # Create grid by taking tensor product of 1d grids. 
  X_grid <- as.matrix(expand.grid(lapply(seq(1, dim_input), function(j) X_marginals[,j])))
  colnames(X_grid) <- rownames(prior_dist_info)
  
  return(X_grid)
}


run_EKI_finite_time <- function(N_steps=1, init_ensemble=NULL) {
  .NotYetImplemented()
}


run_EKI_unit_step <- function(total_steps=1) {
  .NotYetImplemented()
}












