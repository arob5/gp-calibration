#
# seq_design_for_post_approx.r
#
# Implements a framework for sequential design and optimization using log-likelihood
# emulators. The goal is to refine a current log-likelihood emulator by selecting 
# new design points (input parameters values) in such a way that improves 
# the approximation to the posterior distribution induced by the log-likelihood 
# emulator. Note that the target function (the posterior density) is univariate, 
# but the underlying emulator may be multivariate. For example, a multi-output
# Gaussian process (GP) may emulate different quantities which then all 
# feed into the log-likelihood, thus combining to produce a scalar output value. 
#
# Many of the methods implemented here fall under the framework of optimizing 
# some objective function to select the design points. We refer to this objective
# function as an "acquisition function" ("acq" for short), but it is also sometimes
# called a design criterion or refinement function. Other design methods 
# (e.g. Latin Hypercube Design) fall outside of this framework.
#
# The file "seq_design_gp.r" contains similar functions that are designed to act
# on objects that inherit from `gpWrapper`. The functions here are instead 
# defined to act on objects that inherit from `llikEmulator`. 
#
# Depends: These functions are designed to work with objects from the 
#          llikEmulator class, which also involve the llikSumEmulator and 
#          gpWrapper classes. 
#
# Author: Andrew Roberts
#

library(matrixStats)


# -----------------------------------------------------------------------------
# Acquisition Function Framework 
#
# Acquisition functions must be called as `acq_llik_<acq_func_name>(input, ...)`, 
# where `input` is a numeric D-length vector or 1xD matrix representing a 
# parameter value. This section contains functions that operate on acquisition 
# functions. The next section actually defines specific acqusition functions. 
# -----------------------------------------------------------------------------

evaluate_llik_acq_func_vectorized <- function(acq_func, input_mat, llik_em=NULL, ...) {
  # Evaluates the objective function `acq_func` at a set of inputs `input_mat` using 
  # the log-likelihood emulator `llik_em`. 
  #
  # Args:
  #    acq_func: An acquisition function for llikEmulator objects. 
  #    input_mat: matrix, of dimension (N,D) containing N D-dimensional inputs 
  #               stacked in the rows of the matrix. The acquisition functions 
  #               will be evaluated at these points. 
  #    llik_em: object that inherits from the `llikEmulator` class. 
  #    ...: other arguments passed to `acq_func`. 
  #
  # Returns: 
  #    numeric vector of length `nrow(input_mat)` containing the acqusition function 
  #    evaluations. 
  
  input_names <- colnames(input_mat)
  apply(input_mat, 1, function(input) acq_func(matrix(input, nrow=1, 
                                                      dimnames=list(NULL,input_names)), 
                                               llik_em=llik_em, ...))
}


acquire_llik_batch_input_sequentially <- function(llik_emulator, acq_func_name, N_batch, 
                                                  model_response_heuristic, opt_method, 
                                                  llik_exact=NULL, reoptimize_hyperpar=FALSE, 
                                                  lik_par_val=NULL, ...) {
  # This function is the analog of `acquire_batch_input_sequentially()` but 
  # operates on objects of class llikEmulator, rather than gpWrapper. It is thus 
  # the entrypoint for an acquisition optimization-based sequential design loop 
  # which targets improvements to likelihood or posterior approximation. This 
  # function specifically handles the case of one-at-a-time sequential 
  # acquisition; setting `model_response_heuristic` is "none" implies a true 
  # sequential design, in which the exact likelihood is evaluated following 
  # each design point acquisition. Setting `model_response_heuristic` to "KB", 
  # "CL_pessimist", or "CL_optimist" alternatively causes this function to 
  # operate in "heuristic mode", in which the exact likelihood is never  
  # evaluated within this function. `opt_method` currently only allows "grid". 
  #
  # Args:
  #    llik_emulator: object that inherits from class llikEmulator. 
  #    acq_func_name: character(1), the acquisition function name. 
  #    N_batch: integer(1), the number of design points to acquire. 
  #    model_response_heuristic: character(1), currently accepts "none" 
  #                              (requires exact llik evaluations), "KB" 
  #                              (kriging believer), "CL_pessimist" 
  #                              (constant liar pessimist), "CL_optimist" 
  #                              (constant liar optimist). 
  #    opt_method: character(1), the optimization method to use to minimize 
  #                the acquisition function; currently only allows "grid", 
  #                which selects points from a finite set of candidates. 
  #    llik_exact: Only required if `model_response_heuristic == "none"`. 
  #                `llik_exact` can either be a function or an object that 
  #                inherits from `llikEmulator` with `exact_llik == TRUE`. 
  #                If a function must have two arguments; the first is 
  #                interpreted as the input parameter value and the second 
  #                optional argument is the likelihood parameter. The function 
  #                must return exact log-likelihood evaluations evaluated at 
  #                the parameter/likelihood parameter values (the user should 
  #                take care that this function uses the same normalization of 
  #                the likelihood as used by `llik_emulator`).
  #    reoptimize_hyperpar: If TRUE, re-optimizes emulator hyperparameters after 
  #                         acquiring each new design point. Otherwise the 
  #                         hyperparameters are not changed at all in this 
  #                         function. This should only be set to TRUE if 
  #                         `model_response_heuristic == "none"`. 
  #    lik_par_val: No optimization is currently done over the likelihood 
  #                 parameter value.
  
  
  # `model_func` is a function with single argument `input` and the output should be the response value
  # associated with `emulator_obj`. Note that this is often different from the output of the "forward
  # model". For example, if `emulator_obj` is a log-likelihood emulator then `model_func` should output
  # log-likelihood values. If `emulator_obj` implements independent GP emulators for P outputs, 
  # then `model_func` should output these P outputs. In short, the returned value of `model_func` must 
  # be compatible with `emulator_obj` as in `emulator_obj$update(input, model_func(input))`.
  #
  # Returns:
  #     Note that `emulator_obj_updated` will contain the pseudo model responses if a model response 
  #     heuristic is used. 
  
  validate_args_acquire_llik_batch_input_sequentially(llik_emulator, acq_func_name, N_batch, 
                                                      model_response_heuristic, opt_method, 
                                                      llik_exact=NULL, reoptimize_hyperpar=FALSE, 
                                                      lik_par_val=NULL, ...)
  
  if((model_response_heuristic != "none") && (reoptimize_hyperpar)) {
    message("`reoptimize_hyperpar` is TRUE but `model_response_heuristic` is not none.")
  }
  
  # Make a copy to avoid modifying the emulator provided in argument. 
  llik_emulator_copy <- llik_emulator$copy(shallow=FALSE)
  
  # Objects to store the acquired inputs and the associated model 
  # (perhaps pseudo) responses. 
  input_batch <- matrix(nrow=N_batch, ncol=emulator_obj$dim_input)
  response_batch <- matrix(nrow=N_batch, ncol=1L)
  
  for(i in 1:N_batch) {
    # Acquire new input point. 
    input_new <- optimize_acq_single_input(acq_func_name, llik_emulator_copy, 
                                           opt_method, ...)
    input_batch[i,] <- input_new
    
    # Acquire model response or pseudo model response at acquired input. 
    response_new <- get_acq_model_response(input_new, model_response_heuristic, 
                                           llik_emulator, model_func, ...)
    response_batch[i,] <- response_new
    
    # Update emulator. 
    llik_emulator_copy$update(matrix(input_new, nrow=1), response_new, 
                              update_hyperpar=reoptimize_hyperpar, ...)
  }
  
  return(list(input_batch=input_batch, response_batch=response_batch, 
              llik_emulator_updated=llik_emulator_copy))
  
}


get_constant_liar_value <- function(lpost_emulator, constant_liar_method) {
  
  get_design_llik = function(lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, 
                             return_list=FALSE, labels=llik_label, ...)
    
    if(constant_liar_method == "constant_liar_pessimist") {
      return(min(lpost_emulator$outputs_lpost))
    } else if(constant_liar_method == "constant_liar_optimist") {
      return(max(lpost_emulator$outputs_lpost))
    } else if(constant_liar_method == "constant_liar_mean") {
      return(mean(lpost_emulator$outputs_lpost))
    } else {
      stop("Invalid constant liar method: ", constant_liar_method)
    }
  
}


optimize_acq_single_input <- function(acq_func_name, emulator_obj, 
                                      opt_method, candidate_grid=NULL, ...) {
  
  # Define objective function for the optimization. 
  objective_func <- function(input) get(paste0("acq_", acq_func_name))(input, emulator_obj=emulator_obj, ...)
  
  # Dispatch to the correct optimization algorithm. 
  if(opt_method == "grid") input_new <- minimize_objective_grid(objective_func, candidate_grid, ...)
  else stop("`opt_method` ", opt_method, " not supported.")
  
  return(input_new)
}


get_acq_model_response <- function(input, model_response_heuristic, 
                                   emulator_obj=NULL, model_func=NULL, ...) {
  
  if(model_response_heuristic == "none") {
    return(model_func(input))
  } else if(model_response_heuristic == "kriging_believer") {
    return(emulator_obj$predict(input, return_mean=TRUE, return_var=FALSE, ...))
  }
  
}


# -----------------------------------------------------------------------------
# Acquisition functions. 
# These acquisitions are defined with respect to a log-likelihood emulator, 
# which means an object that inherits from the `llikEmulator` class. 
# These functions must all be named using the convention 
# "acq_llik_<name>", where `<name>` is the acquisition function name. 
# The first argument is required to be "input" and the second "llik_em", 
# specifying the input points and log-likelihood emulator, respectively. 
# For now, acquisition functions are defined to accept a single or multiple 
# inputs, with multiple inputs implying a "batch" of candidate points; note 
# that this is different vectorized evaluation at multiple points. For the 
# latter, use `evaluate_llik_acq_func_vectorized()`. 
#
# Note that the "acq_llik_" prefix indicates that these acquisition functions 
# are defined with respect to a log-likelihood emulator object. This does 
# NOT mean that the acquisition function targets the log-likelihood emulator 
# directly. This may be the case, but many of the acqusition functions 
# target the exponential of the log-likelihood emulator, or the induced 
# unnormalized posterior density emulator.
#
# For acquisition functions defined on the exponentiated (likelihood) 
# scale, numerical under/overflow often becomes an issue. Therefore, 
# these acquisition functions have a `log_scale` argument that 
# defaults to TRUE.
#
# TODO: how should we handle the prior here? Have it as an optional 
# argument so that the functions can operate on either the unduced 
# likelihood emulator or the induced unnormalized posterior density 
# emulator? 
# TODO: need to think of a consistent way to handle `lik_par_val`. 
# -----------------------------------------------------------------------------

acq_llik_neg_var_gp <- function(input, llik_em, ...) {
  # Defined for `llik_em` objects that depend on an underlying Gaussian 
  # process (GP); i.e., `is_gp(llik_em$emulator_model)` must be `TRUE`. 
  # Simply returns the negative predictive variance of this underlying 
  # GP at the input point `input`. The negative is due to the fact that 
  # the framework assumes that acquisition functions are always minimized,  
  # so to implement the "maximum variance" acquisition it must be 
  # negated here.
  
  assert_that(is_gp(llik_em$emulator_model))

  args <- c(list(...), 
            list(input=input, gp=llik_em$emulator_model, 
                 bounds=llik_em$get_llik_bounds()))
  
  do.call(acq_neg_var, args)
}


acq_llik_neg_entropy_gp <- function(input, llik_em, ...) {
  # Defined for `llik_em` objects that depend on an underlying Gaussian 
  # process (GP); i.e., `is_gp(llik_em$emulator_model)` must be `TRUE`. 
  # Returns the negative predictive entropy of this underlying 
  # GP at the input point `input`. The negative is due to the fact that 
  # the framework assumes that acquisition functions are always minimized,  
  # so to implement the "maximum entropy" acquisition it must be 
  # negated here.
  
  assert_that(is_gp(llik_em$emulator_model))
  
  args <- c(list(...), 
            list(input=input, gp=llik_em$emulator_model, 
                 bounds=llik_em$get_llik_bounds()))
  do.call(acq_neg_entropy, args)
}


acq_llik_neg_var_lik <- function(input, llik_em, log_scale=TRUE, ...) {
  # Defined for `llik_em` objects that have the method `predict_lik()` 
  # implemented with  the option to return the variance of the likelihood
  # emulator. This function returns the negative variance of the induced 
  # likelihood emulator (exponential of the log-likelihood emulator) at  
  # inputs `input`. If `log_scale` is TRUE (the default), then the 
  # negative of the log of the variance is returned. 
  -llik_em$predict_lik(input, return_mean=FALSE, return_var=TRUE, 
                       log_scale=log_scale, ...)$log_var  
}


acq_llik_IVAR_grid_gp <- function(input, llik_em, grid_points, weights=NULL, ...) {
  # Defined for `llik_em` objects that depend on an underlying Gaussian process 
  # (GP); i.e., `is_gp(llik_em$emulator_model)` must be `TRUE`. Approximates 
  # the standard integrated variance (i.e., integrated mean squared 
  # prediction error) GP criterion by approximating the integral with a discrete 
  # sum at inputs `grid_points`.
  
  assert_that(is_gp(llik_em$emulator_model))
  
  args <- c(list(...), 
            list(input=input, gp=llik_em$emulator_model, grid_points=grid_points, 
                 weights=weights, bounds=llik_em$get_llik_bounds()))
  do.call(acq_IVAR_grid, args)
}


acq_llik_IENT_grid_gp <- function(input, llik_em, grid_points, weights=NULL, ...) {
  # Defined for `llik_em` objects that depend on an underlying Gaussian process (GP); 
  # i.e., `is_gp(llik_em$emulator_model)` must be `TRUE`. Approximates the integrated 
  # entropy criterion by approximating the integral with a descrete sum at inouts `grid_points`. 
  
  assert_that(is_gp(llik_em$emulator_model))
  
  args <- c(list(...), 
            list(input=input, gp=llik_em$emulator_model, grid_points=grid_points, 
                 weights=weights, bounds=llik_em$get_llik_bounds()))
  do.call(acq_IENT_grid, args)
}


acq_llik_IEVAR_grid <- function(input, llik_em, grid_points, weights=NULL, 
                                log_scale=TRUE, plugin=FALSE, ...) {
  # Approximates the integrated expected variance criterion for a random likelihood induced 
  # by a log-likelihood emulator. The integral approximation is obtained via a discretization
  # at inputs `grid_points`. The argument `plugin` is passed to the `calc_expected_lik_cond_var`
  # method of `llik_em` (see the comments in this method for details). In brief, this argument 
  # controls whether the conditional expectation is considered directly for the random 
  # likelihood, or if the underlying emulator model is first conditioned and then "plugged in" 
  # to the likelihood before the expectation is taken. 
  
  # Defined for `llik_em` objects that depend on an underlying Gaussian process (GP); i.e., 
  # `is_gp(llik_em$emulator_model)` must be `TRUE`. Approximates the standard integrated variance 
  # (i.e., integrated mean squared prediction error) GP criterion by approximating the integral 
  # with a discrete sum at inputs `grid_points`.

  log_evar <- llik_em$calc_expected_lik_cond_var(grid_points, input, log_scale=TRUE, 
                                                 plugin=plugin, ...)
  
  if(is.null(weights)) {
    log_weights <- 0
  } else {
    log_weights <- log(weights)
    if(length(log_weights) == 1L) log_weights <- rep(log_weights, length(log_evar))
  }
  
  log_summands <- log_evar + log_weights
  log_IEVAR <- matrixStats::logSumExp(log_summands)
  
  if(log_scale) return(log_IEVAR)
  return(exp(log_IEVAR))                                   
}


acq_llik_IEVAR_grid_old <- function(input, llik_em, grid_points, weights=NULL, log_scale=TRUE, ...) {
  # TODO: how should lik_par be handled here? 
  # TODO: validate_args_acq_IEVAR_grid()
  
  N_grid <- nrow(grid_points)
  if(is.null(weights)) weights <- rep(1/N_grid, N_grid)
  llik_em_copy <- llik_em$copy(shallow=FALSE)
  
  # Emulator predictions at acquisition evaluation locations and at grid locations. 
  pred <- emulator_obj_copy$predict(input, return_mean=TRUE, return_var=TRUE, ...)
  pred_grid <- llik_em_copy$predict(grid_points, return_mean=FALSE, return_var=TRUE, ...)
  
  # Update the GP model, treating the predictive mean as the observed 
  # response at the acquisition evaluation locations. 
  llik_em_copy$update(input, pred$mean, update_hyperpar=FALSE, ...)
  
  # Predict with the conditional ("cond") GP (i.e., the updated GP) at the grid locations. 
  # Convert to the exponentiated scale to obtain log-normal predictive quantities. 
  pred_cond <- llik_em_copy$predict(grid_points, return_mean=TRUE, return_var=TRUE, ...)
  log_pred_cond_LN <- convert_Gaussian_to_LN(mean_Gaussian=pred_cond$mean, var_Gaussian=pred_cond$var,
                                             return_mean=FALSE, return_var=TRUE, log_scale=TRUE)
  
  # Compute the variance inflation term on the log scale. 
  log_var_inflation <- 2 * (pred_grid$var - pred_cond$var)
  
  # Compute log IEVAR. 
  log_summands <- log_pred_cond_LN$log_var + log_var_inflation + log(weights)
  log_IEVAR <- matrixStats::logSumExp(log_summands)
  
  if(log_scale) return(log_IEVAR)
  return(exp(log_IEVAR))
  
}








# -----------------------------------------------------------------------------
# Argument validation functions 
# -----------------------------------------------------------------------------

validate_args_acquire_llik_batch_input_sequentially <- function(llik_emulator, acq_func_name, N_batch, 
                                                                model_response_heuristic, opt_method, 
                                                                llik_exact=NULL, reoptimize_hyperpar=FALSE, 
                                                                lik_par_val=NULL, ...) {
  
  # It is not recommended to optimize hyperparameters based on pseudo model responses. 
  if((model_response_heuristic != "none") && (reoptimize_hyperpar)) {
    message("`reoptimize_hyperpar` is TRUE but `model_response_heuristic` is not none.")
  }
  
  # Exact log-likelihood function is required when no model response heuristic is to be used. 
  if(model_response_heuristic == "none") {
    assert_that(!is.null(llik_exact), 
                msg="`llik_exact` is required when `model_response_heuristic` is 'none'")
    using_exact_llik_obj <- inherits(llik_exact, "llikEmulator")
    assert_that(using_exact_llik_obj || is.function(llik_exact))
    if(using_exact_llik_obj) {
      assert_that(llik_exact$exact_llik)
      if(is.null(lik_par_val)) assert_that(all.equal(llik_exact$lik_par, llik_emulator$lik_par))
    }
  }
  
}








