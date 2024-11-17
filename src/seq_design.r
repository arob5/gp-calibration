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
# Dependencies:
#    gpWrapper.r, llikEmulator.r, statistical_helper_functions.r

# -----------------------------------------------------------------------------
# Acquisition Function Framework 
# The acquisition function framework is generally split into two pieces: 
# acquisition functions that target Gaussian processes (GPs) (seq_design_gp.r)
# and those that opetate on log-likelihood emulators 
# (seq_design_for_post_approx.r). The functions here are high-level functions 
# that interface between the two. This primarily means that the functions 
# here will identify whether the passed model object is a `gpWrapper` or 
# `llikEmulator` and then call the appropriate function. 
# -----------------------------------------------------------------------------

identify_model_type <- function(model) {
  # Identify whether a model object `model` is a Gaussian process (GP)
  # emulator, a log-likelihood emulator, or neither. In particular, 
  # identifies whether `model` inherits from the `gpWrapper` or 
  # `llikEmulator` class. `model` is allowed to be NULL (due to the 
  # fact that some acquisition functions may not require a model), 
  # but any other class besides these three throws an exception.
  
  if(inherits(model, "gpWrapper")) return("gp")
  else if(inherits(model, "llikEmulator")) return("llik_em")
  else if(is.null(model)) return(NULL)
  
  stop("Unrecognized model ", model)
}

is_gp <- function(model) {
  inherits(model, "gpWrapper")
}


evaluate_acq_func_vectorized <- function(acq_func, input_mat, model=NULL, ...) {
  # Evaluates the objective function `acq_func` at a set of inputs `input_mat` using 
  # the model `model`. Calls either `evaluate_acq_gp_func_vectorized()` or 
  # `evaluate_acq_llik_func_vectorized()` depending on whether `model` is a 
  # gpWrapper or llikEmulator. 
  #
  # Args:
  #    acq_func: An acquisition function for llikEmulator objects. 
  #    input_mat: matrix, of dimension (N,D) containing N D-dimensional inputs 
  #               stacked in the rows of the matrix. The acquisition functions 
  #               will be evaluated at these points. 
  #    model: object that inherits from the `gpWrapper` or `llikEmulator` class. 
  #    ...: other arguments passed to `acq_func`. 
  #
  # Returns: 
  #    numeric vector of length `nrow(input_mat)` containing the acqusition function 
  #    evaluations. 
  
  if(is_gp(model) || is.null(model)) evaluate_gp_acq_func_vectorized(acq_func, input_mat, model, ...)
  else if(is_llik_em(model)) evaluate_llik_acq_func_vectorized(acq_func, input_mat, model, ...)
  else stop("Unrecognized model ", model)
}

compare_acq_funcs <- function(input, acq_func_names, model=NULL, ...) {
  # A convenience function that takes a character string of acquisition 
  # function names, and evaluates the respective acquisition functions 
  # at a set of points.
  #
  # Args:
  #    input: matrix, of dimension (N,D) containing N D-dimensional inputs 
  #           stacked in the rows of the matrix. The acquisition functions 
  #           will be evaluated at these points. 
  #    acq_func_names: character, vector of acquisition function names, 
  #                    excluding the "acq_" prefix. 
  #    model: object that inherits from the `gpWrapper` or `llikEmulator` class.
  #    ...: Other arguments passed to the acquisition functions. 
  #
  # Returns: 
  #    Matrix of acqusition function evaluations. Has dimension 
  #    `(nrow(input), length(acq_func_names))`. The colnames attribute
  #    is set to `acq_func_names`. 
  
  acq_func_str <- paste0("acq_", acq_func_names)
  acq_vals <- sapply(acq_func_str, function(x) evaluate_acq_func_vectorized(get(x), input, model, ...))
  acq_vals <- matrix(acq_vals, ncol=length(acq_func_names))
  colnames(acq_vals) <- acq_func_names
  return(acq_vals)
}


compare_acq_funcs_by_model <- function(input, acq_func_names, model_list, ...) {
  # A convenience function that acts as a wrapper around `compare_acq_funcs`,
  # which calls this function for each Gaussian process in the list `gp_list`. 
  #
  # Args:
  #    input: matrix, of dimension (N,D) containing N D-dimensional inputs 
  #           stacked in the rows of the matrix. The acquisition functions 
  #           will be evaluated at these points. 
  #    acq_func_names: character, vector of acquisition function names, 
  #                    excluding the "acq_" prefix. 
  #    model_list: list of objects, each of which inherit from `gpWrapper`
  #                or `llikEmulator`. 
  #    ...: other arguments passed to `compare_acq_funcs()`. 
  #
  # Returns:
  #    list of length `length(gp_list)`. Each element contains a matrix, 
  #    the output of `compare_acq_funcs()` evaluated using the respective model. 
  #    The names attribute of the list is set to `names(model_list)`. 
  
  assert_that(is.list(model_list))
  lapply(model_list, function(model) compare_acq_funcs(input, acq_func_names, model, ...))
}


# -----------------------------------------------------------------------------
# General functions for batch design:
# These are typically intended for "one-shot" design, which often implies
# that they are used to generate an initial design, which will then be refined 
# by sequential design methods. 
# -----------------------------------------------------------------------------

get_batch_design <- function(method, N_batch, prior_params=NULL, design_candidates=NULL, 
                             design_candidate_weights=NULL, ...) {
  # LHS = Latin Hypercube Sample
  # tensor_product_grid = tensor product of equally-spaced 1d grids in each dimension. 
  # simple = simple iid random sample.
  
  if(method == "LHS") return(get_LHS_sample(N_batch, prior_dist_info=prior_params, ...))
  else if(method == "EKI_finite_time") return(run_EKI_finite_time(N_steps))
  else if(method == "tensor_product_grid") return(get_tensor_product_grid(N_batch, prior_dist_info=prior_params, ...))
  else if(method == "simple") return(sample_prior(prior_params, n=N_batch))
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

  # If `prior_dist_info` is NULL, samples points uniformly in the region 
  # defined by `bounds`.
  if(is.null(prior_dist_info)) {
    prior_dist_info <- data.frame(dist="Uniform", param1=bounds[1,], param2=bounds[2,])
  }
  
  # The dimension of the input space.
  dim_input <- nrow(prior_dist_info)
  
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
      } else if(dist_name == "Truncated_Gaussian") {
        bounds[1,j] <- prior_dist_info[j,"bound_lower"]
        bounds[2,j] <- prior_dist_info[j,"bound_upper"]
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


gen_extrapolation_test_inputs <- function(d, N_points_per_dim, max_scaler=2.0, scale="linear", 
                                          target_bounds=NULL, return_long=FALSE) {
  # By default, this function implicitly treats data as lying in the hypercube [-1,1]^d.
  # One can think of [-1,1]^d as the hypercube associated with the bounds of
  # the design points, so test points generated outside of this cube can be thought of
  # as testing the extrapolation of the model.
  # The output test points can be appropriately scaled with respect to a different
  # hypercube by specifiying the `target_bounds` argument.
  # Returns an array `X_test` of shape (d, N_points, d), where `X_test[j,:,:]` gives
  # the `N_points` test inputs spread along the jth standard basis vector. These
  # points will be centered relative to the middle of the target hyperrectangle, 
  # and extend equally in both directions along each standard coordinate direction 
  # outwards from this center point. If `return_long = TRUE` then instead 
  # of returning the 3-dimensional array, a list with elements "inputs" and 
  # "idx" is returned. The former contains a matrix of shape (d*N_points_per_dim, d),
  # which is the inputs from `X_test` stacked row-wise into a single matrix. 
  # The `idx` element is a list of length `d` storing the beginning and end 
  # indices of `inputs` that correspond to each respective dimension. This allows
  # easy extraction of inputs along the jth direction. 

  # These scalers are wrt the hypercube [-1,1]^d.
  reference_bounds <- rbind(rep(-1,d), rep(1,d))
  if(scale=="linear") {
    scalers <- seq(-max_scaler, max_scaler, length.out=N_points_per_dim)
  } else if(scale=="log") {
    .NotYetImplemented()
  }

  # First dimension of `X_test` corresponds to the basis vector/direction being 
  # considered.
  X_test <- array(dim=c(d,N_points_per_dim,d))
  for(j in seq_len(d)) {
    basis_vec <- rep(0,d)
    basis_vec[j] <- 1.0
    X_test[j,,] <- matrix(basis_vec, nrow=N_points_per_dim, ncol=d, byrow=TRUE)
    X_test[j,,] <- mult_vec_with_mat_cols(scalers, X_test[j,,])

    if(!is.null(target_bounds)) {
      X_test[j,,] <- scale_inputs(X_test[j,,], source_bounds=reference_bounds, 
                                  target_bounds=target_bounds)
    }
  }
  
  # Return 3-dimensional array. 
  if(!return_long) return(X_test)
  
  # Alternatively, return list. 
  X_test_long <- do.call(rbind, lapply(1:d, function(j) X_test[j,,]))
  break_points <- vector(mode="list", length=d)
  for(j in seq_len(d)) {
    start_idx <- 1L + (j-1)*N_points_per_dim
    end_idx <- start_idx + N_points_per_dim - 1L
    break_points[[j]] <- c(start_idx, end_idx)
  }
  
  return(list(inputs=X_test_long, idx=break_points))
}


run_EKI_finite_time <- function(N_steps=1, init_ensemble=NULL) {
  .NotYetImplemented()
}


run_EKI_unit_step <- function(total_steps=1) {
  .NotYetImplemented()
}












