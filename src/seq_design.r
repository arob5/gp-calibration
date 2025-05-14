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
# and those that operate on log-likelihood emulators 
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
  
  if(is_gp(model)) return("gp")
  else if(is_llik_em(model)) return("llik_em")
  else if(is.null(model)) return(NULL)
  
  stop("Unrecognized model ", model)
}

is_gp <- function(model) {
  inherits(model, "gpWrapper")
}

update_model <- function(model, input_new, response_new, 
                         reoptimize_hyperpar=FALSE, ...) {
  
  if(!is.matrix(input_new)) input_new <- matrix(input_new, nrow=1L)
  n_new_inputs <- nrow(input_new)
  
  if(is_gp(model)) {
    if(!is.matrix(response_new)) response_new <- matrix(response_new, nrow=n_new_inputs)
    model$update(input_new, response_new, 
                 update_hyperpar=reoptimize_hyperpar, ...)
  } else if(is_llik_em(model)) {
    model$update_emulator(input_new, response_new, 
                          update_hyperpar=reoptimize_hyperpar, ...)
  } else if(is.null(model)) {
    return(model)
  } else {
    stop("Unrecognized `model`; currently only supports gpWrapper, llikEmulator, and NULL.")
  }
  
  return(model)
}


get_model_input_dim <- function(model) {
  if(is_gp(model)) {
    return(model$X_dim)
  } else if(is_llik_em(model)) {
    return(model$dim_input)
  } else if(is.null(model)) {
    return(NULL)
  } else {
    stop("Unrecognized `model`; currently only supports gpWrapper, llikEmulator, and NULL.")
  }
}

get_model_input_names <- function(model) {
  if(is_gp(model)) {
    return(model$X_names)
  } else if(is_llik_em(model)) {
    return(model$input_names)
  } else if(is.null(model)) {
    return(NULL)
  } else {
    stop("Unrecognized `model`; currently only supports gpWrapper, llikEmulator, and NULL.")
  }
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


run_seq_design <- function(model, acq_func_name, n_batch, opt_method,
                           response_heuristic=NULL, true_func=NULL, 
                           reoptimize_hyperpar=FALSE, tracking_settings=NULL, 
                           candidate_grid=NULL, ...) {
  # At each iteration, the minimum value of acquisition function is stored 
  # in the "tracking list". Optionally, users can provide additional settings
  # in `tracking_settings` that can compute additionally quantities that will be
  # tracked. These quantities will be computed every `tracking_settings$interval`
  # iterations. `candidate_grid` is only used when `opt_method = "grid"`; it is 
  # a matrix containing the finite set of inputs over which the discrete 
  # optimization is performed.
 
  # Model must inherit from gpWrapper, llikEmulator, or be NULL.
  identify_model_type(model)
  if(is.null(model)) stop("At present, NULL model is not supported.")
  
  if((response_heuristic != "none") && (reoptimize_hyperpar)) {
    message("`reoptimize_hyperpar` is TRUE but `model_response_heuristic` is",
            " not none. It is not recommended to reoptimize hyperparameters ",
            " using pseudo-observations.")
  }
  
  # If using a constant liar heuristic, then the "lie" is set here and will
  # be constant throughout the whole batch selection process.
  response_lie <- NULL
  if(isTRUE(response_heuristic %in% c("cl_optimist", "cl_pessimist"))) {
    response_lie <- get_pseudo_response(input=NA, 
                                        response_heuristic=response_heuristic, 
                                        model=model, ...)
  }

  # List that stores optimal acquisition values and other computed quantities
  # as the design loop progresses.
  tracking_list <- list(acq_val=rep(NA_real_, n_batch),
                        computed_quantities=list())
  
  # Make a copy to avoid modifying the model provided in argument. 
  model_copy <- model$copy(shallow=FALSE)
  
  # Objects to store the acquired inputs and the associated function 
  # (perhaps pseudo) responses. 
  input_dim <- get_model_input_dim(model_copy)
  input_names <- get_model_input_names(model_copy)
  inputs <- matrix(nrow=n_batch, ncol=input_dim, dimnames=list(NULL, input_names))
  responses <- rep(NA_real_, n_batch)

  for(i in 1:n_batch) {
    print(paste0("Iteration: ", i))

    # Acquire new input point.
    acq_info <- optimize_acq(acq_func_name, model_copy, opt_method, 
                             candidate_grid=candidate_grid, ...)
    input_new <- acq_info$input
    candidate_grid <- acq_info$candidate_grid # NULL if not using grid-based opt.
    inputs[i,] <- input_new
    
    # Acquire function response or pseudo response at acquired input.
    response_new <- get_pseudo_response(input_new, response_heuristic,
                                        model_copy, true_func, 
                                        response_lie=response_lie, ...)
    responses[i] <- response_new

    # Update model. 
    model_copy <- update_model(model_copy, input_new, response_new, 
                               reoptimize_hyperpar, ...)
    
    # Update tracking list.
    tracking_list <- update_tracking_info(tracking_list, model_copy, 
                                          acq_info$acq_val, i, tracking_settings,
                                          ...)
  }
  
  return(list(inputs=inputs, responses=responses, tracking_list=tracking_list,
              updated_model=model_copy, response_lie=response_lie))
}


run_batch_seq_design <- function() {
  return(.NotYetImplemented())
}


optimize_acq <- function(acq_func_name, model, opt_method, 
                         candidate_grid=NULL, ...) {
  # Optimizes an acquisition function. Currently only supports discrete 
  # optimization over a grid of candidate values.
  
  # Define objective function for the optimization. 
  acq_func <- get(paste0("acq_", acq_func_name))
  
  # Dispatch to the correct optimization algorithm. 
  if(opt_method == "grid") {
    acq_info <- minimize_objective_grid(acq_func, model=model, 
                                        candidate_grid=candidate_grid, ...)
  } else {
    stop("`opt_method` ", opt_method, " not supported.")
  }
  
  return(acq_info)
}


minimize_objective_grid <- function(acq_func, model, candidate_grid, 
                                    remove_acquired_point=TRUE, ...) {
  # Evaluates the acquisition function at each input and then returns 
  # the input with the minimum acquisition function value. If 
  # `remove_acquired_point` is TRUE, then the optimal input is removed from 
  # the candidate set and the updated candidate set is returned. This is 
  # not strictly necessary for most GP acquisition criteria, as the criteria
  # will naturally not favor points that the GP is already conditioning on.
  # However, it is often useful for a couple reasons:
  #  - It reduces computation in future rounds, as the acquisition function 
  #    will be evaluated at a smaller candidate set.
  #  - In the case that an acquisition function does not penalize already 
  #    acquired points, then this prevents acquiring the same input again.

  if(!is.matrix(candidate_grid)) {
    stop("`candidate_grid` must be a matrix, with each row containing an input.")
  }
  
  if(nrow(candidate_grid) < 1L) {
    stop("No points in `candidate_grid`.")
  }
  
  acq_func_evals <- evaluate_acq_func_vectorized(acq_func, 
                                                 input_mat=candidate_grid, 
                                                 model=model, ...)
  argmin_idx <- which.min(acq_func_evals)
  
  return_list <- list(input = candidate_grid[argmin_idx,],
                      acq_val = acq_func_evals[argmin_idx])
  
  if(remove_acquired_point) {
    candidate_grid <- candidate_grid[-argmin_idx,, drop=FALSE]
  }
  
  return_list$candidate_grid <- candidate_grid
  
  return(return_list)
}


get_pseudo_response <- function(input, response_heuristic, model=NULL, 
                                true_func=NULL, response_lie=NULL, ...) {
  # This function is used to either evaluate the true function 
  # `true_func(input)`, or to return a "lie"; that is, a pseudo-response
  # instead of the true response. The latter option is used in greedy algorithms
  # for batch design, where the true function is not evaluated until the whole
  # batch has been selected. Note that in greedy batch design, the
  # `response_heuristic` values of NULL and "kb" (kriging believer) imply
  # that the response value should be updated every iteration. On the other
  # hand, the constant liar ("cl") methods maintain a constant lie throughout
  # the batch selection process, in which case this function should only be
  # called once at the beginning of the process. This is handled by the 
  # `response_lie` argument; if non-NULL, then this value will be used instead
  # of computing a new lie.
  
  if(is.null(response_heuristic)) {
    return(true_func(input))
  } else if(response_heuristic == "kb") {
    .NotYetImplemented()
  } else if(response_heuristic %in% c("cl_optimist", "cl_pessimist")) {
    return(get_constant_liar_response(input, response_heuristic, model, 
                                      response_lie=response_lie, ...))
  }
}


get_constant_liar_response <- function(input, response_heuristic, model, 
                                       response_lie=NULL, ...) {
  # Given a gpWrapper or llikEmulator model, either returns the maximum or 
  # minimum value of the response found in the current design used by the model.
  # This is controlled by the value of `response_heuristic`:
  #   "cl_pessimist": returns the minimum value.
  #   "cl_optimist": returns the maximum value.
  # The use of "optimist" and "pessimist" align with the convention of 
  # maximizing functions, and the fact that larger log-likelihood values are 
  # considered better. Note that in greedy batch sequential design, the "lie" 
  # used by the constant liar methods is intended to remain constant throughout
  # the whole batch selection process; i.e., the lie should be based on current
  # true responses, not pseudoresponses that are added throughout the batch
  # selection. This is handled by the `response_lie` argument; if non-NULL, 
  # then this value will be used instead of computing a new lie.
  
  # Non-NULL value is interpreted as a previously fixed lie which should not
  # be updated.
  if(!is.null(response_lie)) return(response_lie)
  
  # Responses in the current design.
  if(is_gp(model)) {
    current_responses <- drop(model$Y)
  } else if(is_llik_em(model)) {
    current_responses <- model$get_design_llik(...)
  } else {
    stop("`get_constant_liar_response` requires `model` to be gpWrapper or llikEmulator object.")
  }
  
  # Select the "lie". 
  if(response_heuristic == "cl_optimist") {
    return(max(current_responses))
  } else if(response_heuristic == "cl_pessimist") {
    return(min(current_responses))
  } else {
    stop("Invalid constant liar heuristic: ", response_heuristic)
  }
  
}


update_tracking_info <- function(tracking_list, model, acq_val, itr, 
                                 tracking_settings=NULL, ...) {
  # Note that when specifying `tracking_settings$interval`, the first iteration
  # will be tracked. For example, if the interval is 10 then tracking will
  # occur at iterations 1, 11, 21, etc.
  
  # Store acquisition function value for current iteration.
  tracking_list$acq_val[itr] <- acq_val
  
  # If `tracking_settings` is provided, compute additional quantities.
  if(is.null(tracking_settings)) return(tracking_list)
  if(is.null(tracking_settings$func_list)) return(tracking_list)
  
  interval <- tracking_settings$interval
  if(is.null(interval)) interval <- 1L
  
  if((itr+interval-1) %% interval == 0) {
    itr_lbl <- paste0("itr_", itr)
    tracking_list$computed_quantities[[itr_lbl]] <- lapply(tracking_settings$func_list, function(f) f(model))
  }
  
  return(tracking_list)
}


# -----------------------------------------------------------------------------
# General functions for batch design:
# These are typically intended for "one-shot" design, which often implies
# that they are used to generate an initial design, which will then be refined 
# by sequential design methods. 
# -----------------------------------------------------------------------------

get_batch_design <- function(method, N_batch, prior_params=NULL, 
                             design_candidates=NULL,  
                             design_candidate_weights=NULL, ...) {
  # LHS = Latin Hypercube Sample
  # tensor_product_grid = tensor product of equally-spaced 1d grids in each dimension. 
  # simple = simple iid random sample.
  # subsample = subsample from `design_candidates` weighted by `design_candidate_weights`.
  
  if(method == "LHS") return(get_LHS_sample(N_batch, prior_dist_info=prior_params, ...))
  else if(method == "EKI_finite_time") return(run_EKI_finite_time(N_steps))
  else if(method == "tensor_product_grid") return(get_tensor_product_grid(N_batch, prior_dist_info=prior_params, ...))
  else if(method == "simple") {
    samp <- sample_prior(prior_params, n=N_batch)
    colnames(samp) <- rownames(prior_params)
    return(samp)
  } else if(method == "subsample") return(subsample_design(design_candidates, n=N_batch, 
                                                           design_candidate_weights, ...))
  else stop("Design method ", method, " not supported.")
}


get_init_design_list <- function(inv_prob, design_method, N_design, inputs=NULL, ...) {
  # A convenience function that constructs an initial design for an inverse 
  # problem and returns a list containing the design inputs, as well as 
  # associated log-likelihood, log-prior, and forward model evaluations.
  # 
  # Args:
  #    inv_prob: list, defining the inverse problem, following the conventions
  #              described in `inv_prob_test_functions.r`. Must have elements
  #              "llik_obj" and "par_prior". 
  #    design_method: character, passed to the "method" argument of 
  #                   `get_batch_design()`.
  #    N_design: integer, number of design points.
  #    inputs: If design inputs have already been sampled, then they can be provided 
  #            via this argument. `design_method` and `N_design` are still required 
  #            so that they can be included in the design info list.
  #    ...: additional arguments passed to `get_batch_design()`.
  #
  # Returns:
  # list, with elements:
  #    "input": matrix with design inputs in the rows.
  #    "llik": numeric vector of log-likelihood evaluations at design inputs.
  #    "lprior": numeric vector of log-prior evaluations at design inputs.
  #    "bounds": matrix storing bounding box for design inputs.
  #    "fwd": if applicable, matrix of forward model outputs at design points.
  
  design_info <- list(design_method=design_method, N_design=N_design)
  
  # Design inputs, ensuring correct parameter order.
  if(is.null(inputs)) {
    design_info$input <- get_batch_design(design_method, N_design, 
                                          prior_params=inv_prob$par_prior, ...)
  } else {
    assert_that(setequal(colnames(inputs), inv_prob$par_names))
    design_info$input <- inputs
  }
  
  design_info$input <- design_info$input[,inv_prob$par_names,drop=FALSE]
  
  # Associated log-likelihood and log prior evaluations.
  design_info$llik <- inv_prob$llik_obj$assemble_llik(design_info$input)
  design_info$lprior <- calc_lprior_dens(design_info$input, inv_prob$par_prior)
  design_info$lpost <- design_info$llik + design_info$lprior
  design_info$bounds <- get_bounds(design_info$input)
  
  # Not all log-likelihood objects have the `run_fwd_model()` method.
  design_info$fwd <- try(inv_prob$llik_obj$run_fwd_model(design_info$input))
  
  return(design_info)
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
  if(!is.null(prior_dist_info)) dim_input <- nrow(prior_dist_info)
  else if(!is.null(bounds)) dim_input <- ncol(bounds)

  # Generate LHS design on unit hypercube. 
  X_lhs <- lhs::randomLHS(N_batch, dim_input)
  
  # Apply inverse CDF transform using prior distributions.
  X_lhs <- map_from_uniform(X_lhs, prior_dist_info, bounds)
  
  # For 1 dimensional data, optionally order samples in increasing order. 
  if(order_1d && (dim_input == 1)) X_lhs <- X_lhs[order(X_lhs),,drop=FALSE]
    
  colnames(X_lhs) <- rownames(prior_dist_info)
  return(X_lhs)
}


update_LHS_sample <- function(X, n_batch, prior_dist_info=NULL, bounds=NULL) {
  # The function lhs::augmentLHS updates an LHS sample in the unit hypercube.
  # This wrapper function assumes that `X` is a transformed LHS sample, where
  # the marginals have been transformed to align with the priors in 
  # `prior_dist_info`. So this function simply undoes this transform, then 
  # uses lhs::augmentLHS to augment the design, then transforms back. This 
  # will return a matrix with number of rows equal to `nrow(X) + n_batch`. The
  # first `nrow(X)` rows of the returned matrix will correspond to `X`, while
  # the remainder of the rows will constitute the newly sampled points.
  
  # First map so that marginals are U(0,1).
  X_unif <- map_to_uniform(X, prior_dist_info, bounds)
  
  # Augment uniform Latin hypercube design.
  X_unif_new <- lhs::augmentLHS(X_unif, m=n_batch)
  
  # Map back to desired marginals.
  map_from_uniform(X_unif_new, prior_dist_info, bounds)
}


map_from_uniform <- function(X, prior_dist_info=NULL, bounds=NULL) {
  # Note that this function relies on the ordering of the columns of X, while
  # `map_to_uniform` requires the parameter names to be set to the column 
  # names of `X`.

  # Either `prior_dist_info` or `bounds` must be provided. 
  assert_that(!is.null(prior_dist_info) || !is.null(bounds))
  
  # If `prior_dist_info` is NULL, assumes uniform priors with support determined
  # by `bounds`.
  if(is.null(prior_dist_info)) {
    prior_dist_info <- data.frame(dist="Uniform", param1=bounds[1,], param2=bounds[2,])
  }
  
  # The dimension of the input space.
  dim_input <- nrow(prior_dist_info)
  
  # If `bounds` is NULL, set to NA. No bounds will be enforced beyond those
  # implied by the prior distributions. For distributions with unbounded support. 
  if(is.null(bounds)) bounds <- matrix(NA, nrow=2, ncol=dim_input)
  
  # Apply inverse CDF transform using prior distributions.
  for(j in seq_len(dim_input)) {
    dist_name <- prior_dist_info[j,"dist"]
    if(dist_name == "Uniform") {
      lower_bound <- ifelse(is.na(bounds[1,j]), prior_dist_info[j,"param1"], 
                            max(bounds[1,j], prior_dist_info[j,"param1"]))
      upper_bound <- ifelse(is.na(bounds[2,j]), prior_dist_info[j,"param2"], 
                            min(bounds[2,j], prior_dist_info[j,"param2"]))
      X[,j] <- qunif(X[,j], lower_bound, upper_bound)
    } else if(dist_name == "Gaussian") {
      lower_bound <- ifelse(is.na(bounds[1,j]), -Inf, bounds[1,j])
      upper_bound <- ifelse(is.na(bounds[2,j]), Inf, bounds[2,j])
      X[,j] <- truncnorm::qtruncnorm(X[,j], a=lower_bound, b=upper_bound, 
                                     mean=prior_dist_info[j,"param1"], 
                                     sd=prior_dist_info[j, "param2"])
    } else if(dist_name == "Gamma") {
      if(!is.na(bounds[1,j]) || !is.na(bounds[1,j])) {
        stop("`map_from_uniform()`: Bounds not supported for Gamma dist.")
      }
      X[,j] <- qgamma(X[,j], shape=prior_dist_info[j,"param1"], rate=prior_dist_info[j,"param2"])
    } else if(dist_name == "Beta") {
      if(!is.na(bounds[1,j]) || !is.na(bounds[1,j])) {
        stop("`map_from_uniform()`: Bounds not supported for Beta dist.")
      }
      X[,j] <- qbeta(X[,j], shape1=prior_dist_info[j,"param1"], shape2=prior_dist_info[j,"param2"])
    } else {
      stop("Unsupported prior distribution: ", dist_name)
    }
  }
  
  return(X)
}


map_to_uniform <- function(X, prior_dist_info=NULL, bounds=NULL) {
  # Given a matrix with rows corresponding to samples, and columns to 
  # parameters, maps the samples so that each transformed sample has 
  # uniform(0,1) marginals. The original samples are assumed to be independent
  # across dimensions. See `get_LHS_sample` for description of how `bounds` is
  # used. These functions are essentially all defunct, and will be replaced
  # by a more robust class for encoding probability distributions.

  # Either `prior_dist_info` or `bounds` must be provided. 
  assert_that(!is.null(prior_dist_info) || !is.null(bounds))
  
  # If `prior_dist_info` is NULL, assumes uniform priors with support determined
  # by `bounds`.
  if(is.null(prior_dist_info)) {
    prior_dist_info <- data.frame(dist="Uniform", param1=bounds[1,], param2=bounds[2,])
  }
  
  # The dimension of the input space.
  dim_input <- ncol(X)
  
  # If `bounds` is NULL, set to NA. No bounds will be enforced beyond those
  # implied by the prior distributions. For distributions with unbounded support. 
  if(is.null(bounds)) bounds <- matrix(NA, nrow=2, ncol=dim_input)
  
  # Apply CDF transform to map to uniform.
  for(j in seq_len(dim_input)) {
    par_name <- rownames(prior_dist_info)[j]
    dist_name <- prior_dist_info[j,"dist"]
    if(dist_name == "Uniform") {
      lower_bound <- ifelse(is.na(bounds[1,j]), prior_dist_info[j,"param1"], 
                            max(bounds[1,j], prior_dist_info[j,"param1"]))
      upper_bound <- ifelse(is.na(bounds[2,j]), prior_dist_info[j,"param2"], 
                            min(bounds[2,j], prior_dist_info[j,"param2"]))
      X[,par_name] <- punif(X[,par_name], lower_bound, upper_bound)
    } else if(dist_name == "Gaussian") {
      lower_bound <- ifelse(is.na(bounds[1,j]), -Inf, bounds[1,j])
      upper_bound <- ifelse(is.na(bounds[2,j]), Inf, bounds[2,j])
      X[,par_name] <- truncnorm::ptruncnorm(X[,par_name], a=lower_bound, b=upper_bound, 
                                            mean=prior_dist_info[j,"param1"], 
                                            sd=prior_dist_info[j, "param2"])
    } else if(dist_name == "Gamma") {
      if(!is.na(bounds[1,j]) || !is.na(bounds[1,j])) {
        stop("`map_from_uniform()`: Bounds not supported for Gamma dist.")
      }
      X[,par_name] <- pgamma(X[,par_name], shape=prior_dist_info[j,"param1"], rate=prior_dist_info[j,"param2"])
    } else if(dist_name == "Beta") {
      if(!is.na(bounds[1,j]) || !is.na(bounds[1,j])) {
        stop("`map_from_uniform()`: Bounds not supported for Beta dist.")
      }
      X[,par_name] <- pbeta(X[,par_name], shape1=prior_dist_info[j,"param1"], shape2=prior_dist_info[j,"param2"])
    } else {
      stop("Unsupported prior distribution: ", dist_name)
    }
  }
  
  return(X)
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
    bounds <- get_prior_bounds(prior_dist_info, 
                               tail_prob_excluded=tail_prob_excluded)
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


subsample_design <- function(design_candidates, n, 
                             design_candidate_weights=NULL, ...) {
  # Randomly samples `n` rows (without replacement) from the rows of 
  # `design_candidates`. If a vector of weights is provided by 
  # `design_candidate_weights` then a weighted sample is taken.
  #
  # Args:
  #    design_candidates: matrix, with rows interpreted as points.
  #    n: integer, number of subsamples to extract.
  #    design_candidate_weights: numeric, nonnegative, of length equal to the 
  #                              number of rows of `design_candidates`. The 
  #                              weights used in sampling.
  #
  # Returns:
  # matrix of min(`n`, nrow(design_candidates)) rows, the row-subsetted 
  # version of `design_candidates`.
  
  assert_that(is.matrix(design_candidates))
  n_rows <- nrow(design_candidates)
  
  if(n >= n_rows) {
    message("Warning: `n` is >= number of rows in `design_candidates`. Returning 
            full matrix without subsampling.")
    return(design_candidates)
  }
  
  # Randomly sample row indices.
  idcs <- sample(1:n_rows, size=n, replace=FALSE, prob=design_candidate_weights)
  design_candidates[idcs,,drop=FALSE]
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













