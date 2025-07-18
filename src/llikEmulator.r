# 
# llikEmulator.r
#
# Class definitions for the `llikEmulator`, `llikSumEmulator`, and classes which 
# inherit from these classes. All classes defined using reference classes. The 
# operational classes in this file are:
#   llikEmulator: The base log-likelihood (llik) emulator class. 
#   llikEmulatorGP: llik emulator derived from direct Gaussian process (GP)
#                   emulation of the llik.
#   llikEmulatorExactGaussDiag: Exact Gaussian likelihood with diagonal 
#                               covariance matrix (no emulation). Mean of the 
#                               Gaussian given by the output of a forward model.
#   llikEmulatorGPFwdGaussDiag: The same likelihood structure as above, but the
#                               forward model is emulated by a (potentially 
#                               multi-output) GP.
#   llikEmulatorGPFwdGauss: The same likelihood structure as above, but the
#                           forward model is emulated by a (potentially 
#                           multi-output) GP.
# 
#
# The following classes are in development/currently not operational.
#   llikSumEmulator: sum of independent llik emulators.
#   llikEmulatorMultGausGP: Multiplicative Gaussian likelihood with a single 
#                           variance parameter. GP emulates the quadratic 
#                           model-data misfit term (i.e., the sum of squared 
#                           errors).
#   llikEmulatorExactGauss: Generalization of `llikEmulatorExactGaussDiag` that
#                           allows arbitrary (non-diagonal) covariance matrix.
#
# Andrew Roberts
#
# Depends on: general_helper_functions.r, statistical_helper_functions.r, 
#             gpWrapper.r
# 

library(assertthat)
library(ggmatplot)
library(abind)

# -----------------------------------------------------------------------------
# llikEmulator: Class encapsulating a (typically stochastic) approximation to 
#               the log-likelihood. 
# 
# This is the base class for a log likelihood approximation. It is intentionally 
# quite general and intended to be the parent class for more specialized 
# log likelihood approximation classes (e.g. multiplicative Gaussian llik 
# with sum of squares emulated by GPs). This class is primarily intended 
# to encapsulate a stochastic approximation to the log-likelihood, but 
# it can be set up to function with a deterministic approximation (or 
# just the exact likelihood, which is useful for testing code). 
# The llik emulator is assumed to be defined by some underlying 
# surrogate/emulator model, which is stored in the field `emulator_model` of the 
# class. Currently, this model is almost always a Gaussian process (GP), in 
# which case `emulator_model` is a `gpWrapper` object. `lik_par` stores the 
# likelihood parameters. These can be stored in advance if the likelihood 
# parameters will be fixed throughout an analysis; otherwise, this field can be 
# left NULL. Note that the likelihood approximation may be a function of the 
# likelihood parameters (e.g., the likelihood parameters are included as inputs 
# to a GP emulator). Note that many methods in this base class are left as 
# `.NotYetImplemented()`, but are implemented by classes which inherit from 
# the base class.
#
# This class organizes the pieces comprising the likelihood into 
# (1) "inputs" or "par": the primary parameters of interest; e.g. the parameters 
# characterizing the forward model in an inverse problem; and (2) "lik_par", 
# other likelihood parameters (e.g. variance parameters) that define 
# the observation model, may be unknown, but are not of primary interest. 
#
# llikEmulator fields:
# lik_description: Character description of the likelihood.
# llik_label: A short character label for the class. 
# lik_par: Stores any likelihood parameters (e.g. observation variances) 
#          beyond the core input parameters.
# input_names: character vector of input parameter names. The order of this 
#              vector sets the official parameter order.
# dim_input: integer, the length of `input_names`; i.e., the dimension of the 
#            input parameter space.
# default_conditional: Should the likelihood be viewed as a function of, or 
#                      conditional on, the likelihood parameters `lik_par` by 
#                      default? If the latter, then the likelihood may exclude 
#                      terms that only depend on `lik_par`. Only used to set 
#                      the default behavior, can always be adjusted when calling
#                      class methods.
# default_normalize: Should the likelihood drop constant terms by default? 
#                    If `default_conditional = TRUE`, then by default the 
#                    `lik_par` parameters are viewed as constant terms in 
#                    addition to any true constants. Also just sets the default
#                    behavior.
#                      
# emulator_model: The emulator underlying the likelihood approximation. At 
#                 present, this is typically a `gpWrapper` object.
# emulator_description: Character description of the emulator model.
# llik_pred_dist: character, for stochastic llik emulators this is the 
#                 distribution of the llik emulator at a fixed input parameter.
#                 e.g., is a GP directly emulates the log-likelihood, then this
#                 field will be "Gaussian", since the GP predictive distribution 
#                 is Gaussian at any input. This field will be empty if the 
#                 log-likelihood is deterministic or the predictive distribution
#                 is not known.
# exact_llik: logical, is this actually an exact llik without any emulation or 
#             approximation? Useful for testing.
# shift_func: function, a deterministic trend that is added to llikEmulator.
#             It thus affects the predictive mean but not the variance. This
#             function must be vectorized to accept multiple inputs (one per
#             row of a matrix).
# llik_bounds: either a length-2 vector or a matrix. If a vector, of the
#              form c(lower, upper) providing lower and upper bounds that 
#              are independent of the input. If a function, a vectorized 
#              function that accepts a matrix of inputs and returns a list
#              with elements "lower" and "upper, where each is 
#              a vector of length equal to the number of inputs, thus providing
#              input-dependent bounds. It is important to note that these 
#              bounds are applied to the llikEmulator predictive distribution
#              after applying the `shift_func`; i.e., the bounds apply to the 
#              final shifted predictive distribution.
# is_lpost_em: logical(1), if TRUE indicates that the object approximates
#              the log-posterior (i.e., the unnormalized log-posterior density).
#              If FALSE, indicates a log-likelihood approximation.
#
# Core Methods:
# predict: computes log-likelihood predictions at specific input parameter
#          values. Typically returns a list summarizing the llik emulator 
#          predictive distribution at the inputs.
# sample: samples one or more values from the llik emulator distribution at  
#         specific inputs. For deterministic llik emulators, this reduces to 
#         outputing the deterministic prediction.
# predict_emulator: return predictions from `emulator_model`, the underlying 
#                   emulator, at specific inputs.
# sample_emulator: return samples from `emulator_model` at specific inputs.
# assemble_llik: defines a map from the `emulator_model` output to the llik 
#                values. Specifically, `assemble_llik(em_val, lik_par, ...)`
#                computes the llik using the quantity `em_val` and likelihood 
#                parameters `lik_par`. `em_val` is commonly some sort of 
#                sufficient statistic. If `exact_llik = TRUE`, then `em_val`
#                is just the input parameter and the `assemble_llik` represents
#                a typical llik function. If `emulator_model` directly emulates
#                the llik, then `em_val` is the llik itself, and thus 
#                `assemble_llik` is the identity function. As another example, 
#                supposing the predictive distribution of `emulator_model` has 
#                a "mean" attribute, then the composition
#                `assemble_llik(predict_emulator(input,...)$mean, lik_par, ...)`
#                represents a llik prediction using a "plug-in" mean approach.
# predict_lik: returns likelihood (not log-likelihood predictions), which are 
#              typically log-transformed for numerical stability. If the 
#              llik emulator is deterministic, then the log of the likelihood 
#              predictions are equal to the llik predictions from the `predict()`
#              method. Otherwise, they may differ; for example, if 
#              `llik_pred_dist = "Gaussian"` then `predict_lik` summarizes the 
#              log-normal distribution that results from exponentiating the llik.
# calc_lik_approx: computes a deterministic approximation to the likelihood, 
#                  derived from the underlying emulator model. If the llik 
#                  emulator is deterministic, then there is only one such 
#                  approximation, and this function aligns with `predict_lik()`.
#                  Otherwise, there are many options. The simplest option is 
#                  the "plug-in" mean approach mentioned above. 
# There are other methods for plotting, computing quantiles/intervals, and 
# model validation. Need to add documentation for these.
# -----------------------------------------------------------------------------

llikEmulator <- setRefClass(
  Class = "llikEmulator", 
  fields = list(emulator_model="ANY", lik_description="character", 
                llik_label="character", emulator_description="character", 
                default_conditional="logical", default_normalize="logical", 
                lik_par="ANY", input_names="character", dim_input="integer", 
                llik_pred_dist="character", exact_llik="logical", 
                llik_bounds="ANY", shift_func="ANY", is_lpost_em="logical")
)

llikEmulator$lock("llik_label")
llikEmulator$lock("llik_pred_dist")
llikEmulator$lock("exact_llik")

llikEmulator$methods(
  
  initialize = function(llik_label, input_names, lik_description, 
                        emulator_description, dim_input,  emulator_model=NULL, 
                        default_conditional=FALSE, default_normalize=FALSE,
                        lik_par=NULL, llik_bounds=c(-Inf, Inf), 
                        llik_pred_dist="unspecified", exact_llik=FALSE, 
                        shift_func=NULL, is_lpost_em=FALSE, ...) {

    assert_that(is.function(llik_bounds) || isTRUE(length(llik_bounds)==2L))
    if(!is.function(llik_bounds)) assert_that(llik_bounds[2] >= llik_bounds[1])

    initFields(llik_label=llik_label, input_names=input_names, 
               lik_description=lik_description, dim_input=dim_input, 
               emulator_description=emulator_description, 
               emulator_model=emulator_model, 
               default_conditional=default_conditional,
               default_normalize=default_normalize, 
               llik_bounds=llik_bounds, is_lpost_em=is_lpost_em,
               lik_par=lik_par, llik_pred_dist=llik_pred_dist, 
               exact_llik=exact_llik, shift_func=shift_func)  
  }, 
  
  get_lik_par = function(lik_par_val=NULL, ...) {
    if(is.null(lik_par_val)) lik_par_val <- .self$lik_par
    # assert_that(!is.null(lik_par_val), 
    #             msg="Missing likelihood parameter.")
    return(lik_par_val)
  },
  
  get_llik_attr = function(attr_name, ...) {
    # Extracts a field and assigns a name attribute given by `llik_label`.
    setNames(.self$field(attr_name), llik_label)
  },
  
  get_llik_bounds = function(input, lik_par_val=NULL, shift_bounds=TRUE, 
                             apply_shift=TRUE, shift_func_new=NULL, add_shifts=TRUE, ...) {
    # Returns list of (potentially input-dependent) lower and upper bounds on
    # the predictive distribution. By default, the bounds are interpreted
    # as applying to the already shifted distribution; e.g., for upper bounds, 
    # L(u) + shift(u) <= bound(u). If `shift_bounds` is TRUE, the shifted
    # bound bound(u) - shift(u) is returned (which is a bound on L(u)). If 
    # FALSE, the bound is not shifted. Note that `apply_shift` must also 
    # be TRUE for the shift to occur.

    input <- .self$get_input(input)
    lb <- .self$llik_bounds
    
    # Get shift function.
    sf <- NULL
    if(apply_shift && shift_bounds) {
      sf <- .self$get_shift_func(shift_func_new, add_shifts)
    }
    
    if(is.null(sf)) shift_vals <- 0
    else shift_vals <- sf(input)
    
    # Evaluate bound function at inputs.
    if(is.function(lb)) {
      bounds <- lb(input, lik_par_val=lik_par_val, ...)
      bounds$lower <- bounds$lower - shift_vals
      bounds$upper <- bounds$upper - shift_vals
      return(bounds)
    }
    
    # Otherwise apply bounds which only potentially depend on the input
    # through the shift.
    list(lower=lb[1] - shift_vals, upper=lb[2] - shift_vals)
  },
  
  get_shift_func = function(shift_func_new=NULL, add_shifts=TRUE) {
    # Returns the deterministic shift function. The default behavior is to
    # use the sum of the shift function provided in the argument 
    # `shift_func_new` and the default `.self$shift_func`. If
    # exactly one of these is non-NULL, then the non-NULL shift function is
    # used. If both are non-NULL but `add_shifts = FALSE`, then only
    # `shift_func_new` is used (i.e., it overwrites the default).
    # The function returns NULL if both are NULL.
    
    f <- NULL
    
    if(is.null(shift_func_new)) {
      f <- .self$shift_func
    } else {
      if(add_shifts && !is.null(.self$shift_func)) {
        f <- function(input) shift_func_new(input) + .self$shift_func(input)
      } else {
        f <- shift_func_new
      }
    }

    return(f)
  },
  
  get_input = function(input, ...) {
    # Helper function to preprocess input values passed by the user. Allows 
    # for `input` to be a numeric vector or a matrix in which each row is 
    # a different input parameter value.
    
    # If not a matrix, assume that a single input point has been passed as 
    # a numeric vector.
    if(is.null(nrow(input))) input <- matrix(input, nrow=1L)
    
    assert_that(is.matrix(input) && (ncol(input)==dim_input), 
                msg="`input` must be matrix with ncol equal to `dim_input`.")
    
    # TODO: checking input names becomes problematic when this needs to run 
    # using external code; e.g., `mcmc_bt_wrapper()`. Need to think of the best 
    # way to handle this. For now, not checking names and assuming parameter 
    # ordering is correct.
    # assert_that(!is.null(colnames(input)) && all(is.element(colnames(input), input_names)),
    #             msg="`input` must have colnames set to subset of `input_names`.")
    # 
    # return(input[,input_names, drop=FALSE])
    
    if(is.null(colnames(input))) colnames(input) <- .self$input_names
    
    return(input)
  },
  
  sample_emulator = function(input, N_samp=1, ...) {
    .NotYetImplemented()
  },
  
  update_emulator = function() {
    .NotYetImplemented()
  },
  
  predict_emulator = function(input, lik_par_val=NULL, em_pred_list=NULL, 
                              return_mean=TRUE, return_var=TRUE, 
                              return_cov=FALSE, return_cross_cov=FALSE, 
                              input_cross=NULL, ...) {
    # Compute predictions for `emulator_model` (note: not necessarily log
    # likelihood predictions). This default method assumes the emulator inherits 
    # from `gpWrapper`, but different llikEmulator classes may want to override 
    # this default method. 
    
    assert_that(inherits(.self$emulator_model, "gpWrapper"))
    if(!is.null(em_pred_list)) return(em_pred_list)
    
    .self$emulator_model$predict(.self$get_input(input), return_mean=return_mean, 
                                 return_var=return_var, return_cov=return_cov, 
                                 return_cross_cov=return_cross_cov, 
                                 input_cross=input_cross, ...)
  },
  
  assemble_llik = function(em_val, lik_par_val=NULL, 
                           conditional=default_conditional, 
                           normalize=default_normalize, ...) {
    # Since a log-likelihood is scalar-valued, this should always return a 
    # numeric vector of length equal to the number of rows in `input`. This
    # method provides the main interface for computing the log-likelihood as
    # a function of the output quantity which the emulator predicts. It 
    # calls the child class internal method to do the bulk of the work, but
    # this parent method adds processing of the likelihood parameters and 
    # applies the constant shift function to the outputs. Note that the 
    # deterministic shift is not applied here, as this function is generally
    # not a function of the inputs and hence the shift cannot be computed.
    
    # Fetch the likelihood parameters. 
    lik_par_val <- .self$get_lik_par(lik_par_val)
    
    # Run child method.
    llik_vals <- .self$.assemble_llik(em_val, lik_par_val=lik_par_val,
                                      conditional=conditional,
                                      normalize=normalize, ...)

    return(llik_vals)
  },
  
  sample = function(input, lik_par_val=NULL, N_samp=1L,
                    em_pred_list=NULL, llik_pred_list=NULL,
                    conditional=default_conditional,
                    normalize=default_normalize, apply_shift=NULL, 
                    shift_func_new=NULL, add_shifts=TRUE, ...) {
    # Returns samples from the log-likelihood at `M` different inputs, where
    # `M = ncol(input)`. Returns (M, N_samp) matrix of samples.
    # `apply_shift` default to TRUE if `llik_pred_list` is NULL and FALSE
    # otherwise. This is because it is assumed that `llik_pred_list` has 
    # already been shifted.

    # Setup.
    input <- .self$get_input(input)
    lik_par_val <- .self$get_lik_par(lik_par_val)
    if(is.null(apply_shift)) apply_shift <- is.null(llik_pred_list)
    
    # Run child method.
    llik_samp <- .self$.sample(input, lik_par_val=lik_par_val, N_samp=N_samp,
                               em_pred_list=em_pred_list, 
                               llik_pred_list=llik_pred_list,
                               conditional=default_conditional,
                               normalize=default_normalize, ...)
    
    # Apply shift.
    if(apply_shift) {
      s <- get_shift_func(shift_func_new=shift_func_new, add_shifts=add_shifts)
      if(!is.null(s)) llik_samp <- llik_samp + s(input)
    }
    
    if(sum(is.na(llik_samp)) > 0L) {
      message("Some llikEmulator samples are NA.")
    }
    
    return(llik_samp)
  },
  
  predict = function(input, lik_par_val=NULL, return_mean=TRUE, return_var=TRUE, 
                     return_cov=FALSE, return_cross_cov=FALSE, input_cross=NULL,
                     conditional=default_conditional, normalize=default_normalize, 
                     em_pred_list=NULL, apply_shift=TRUE, shift_func_new=NULL,
                     add_shifts=TRUE, ...) {
    # Setup.
    input <- .self$get_input(input)
    lik_par_val <- .self$get_lik_par(lik_par_val)

    # Run child method.
    llik_pred_list <- .self$.predict(input, lik_par_val=lik_par_val, 
                                     return_mean=return_mean, return_var=return_var,
                                     return_cov=return_cov, 
                                     return_cross_cov=return_cross_cov,
                                     input_cross=input_cross, conditional=conditional,
                                     normalize=normalize, em_pred_list=em_pred_list, ...)
    
    # Apply deterministic shift (only affects mean).
    if(apply_shift) {
      s <- get_shift_func(shift_func_new=shift_func_new, add_shifts=add_shifts)
      if(!is.null(s)) llik_pred_list$mean <- llik_pred_list$mean + s(input)
    }
    
    return(llik_pred_list)
  },
  
  
  predict_lik = function(input, lik_par_val=NULL, em_pred_list=NULL, 
                         llik_pred_list=NULL, return_mean=TRUE, return_var=TRUE, 
                         return_cov=FALSE, return_cross_cov=FALSE, 
                         input_cross=NULL, conditional=default_conditional,  
                         normalize=default_normalize, log_scale=TRUE,
                         apply_shift=NULL, shift_func_new=NULL, add_shifts=TRUE,
                         ...) {
    # Returns a list of moments/predictive quantities that summarize the 
    # likelihood predictive distribution at inputs `input`. By default,
    # the quantities are returned on the log-scale and can be accessed
    # like `l$log_mean` and `l$log_var`. If `log_scale = FALSE` then the
    # quantities are on the likelihood scale, and can be accessed as
    # `l$mean` and `l$var`. `apply_shift` default to TRUE if `llik_pred_list` 
    # is NULL and FALSE otherwise. This is because it is assumed that 
    # `llik_pred_list` has already been shifted.

    # Setup.
    input <- .self$get_input(input)
    lik_par_val <- .self$get_lik_par(lik_par_val)
    if(is.null(apply_shift)) apply_shift <- is.null(llik_pred_list)
    
    # Run child method.
    lik_pred_list <- .self$.predict_lik(input, lik_par_val=lik_par_val, 
                                        return_mean=return_mean, return_var=return_var,
                                        return_cov=return_cov, log_scale=log_scale,
                                        return_cross_cov=return_cross_cov,
                                        input_cross=input_cross, conditional=conditional,
                                        normalize=normalize, em_pred_list=em_pred_list,
                                        llik_pred_list=llik_pred_list, ...)
    
    # Apply deterministic shift (affects mean and variance/cov).
    if(apply_shift) {
      s <- get_shift_func(shift_func_new=shift_func_new, add_shifts=add_shifts)
      if(!is.null(s)) {
        shift_vals <- s(input)
        
        # Shift mean.
        if(log_scale) lik_pred_list$log_mean <- lik_pred_list$log_mean + shift_vals
        else lik_pred_list$mean <- lik_pred_list$mean * exp(shift_vals)
        
        # Shift variance.
        if(return_var) {
          if(log_scale) lik_pred_list$log_var <- lik_pred_list$log_var + 2*shift_vals
          else lik_pred_list$var <- lik_pred_list$var * shift_vals^2
        }
        
        # Shift covariance.
        if(!is.null(lik_pred_list$cov)) {
          lik_pred_list$cov <- mult_vec_with_mat_cols(shift_vals, lik_pred_list$cov)
          lik_pred_list$cov <- mult_vec_with_mat_rows(shift_vals, lik_pred_list$cov)
        }
      }
    }
    
    return(lik_pred_list)
  },
  
  calc_quantiles = function(input, p, lik_par_val=NULL, 
                            conditional=default_conditional, 
                            normalize=default_normalize, em_pred_list=NULL,
                            lower_tail=TRUE, apply_shift=TRUE,
                            shift_func_new=NULL, add_shifts=TRUE, ...) {
    # Computes quantiles of the log-likelihood predictive distribution pointwise
    # over a set of inputs. `p` is the probability defining the quantile.
    
    input <- .self$get_input(input)
    
    # Call child method.
    q <- .self$.calc_quantiles(input, p, lik_par_val=lik_par_val, 
                               em_pred_list=em_pred_list,
                               conditional=conditional, normalize=normalize, 
                               lower_tail=lower_tail, ...)
    
    # Apply deterministic shift. 
    if(apply_shift) {
      s <- get_shift_func(shift_func_new=shift_func_new, add_shifts=add_shifts)
      if(!is.null(s)) q <- q + s(input) 
    }
    
    return(q)
  },
  
  calc_lik_quantiles = function(input, p, lik_par_val=NULL, 
                                 em_pred_list=NULL, log_scale=TRUE, 
                                 conditional=default_conditional, 
                                 normalize=default_normalize, 
                                 lower_tail=TRUE, apply_shift=TRUE,
                                 shift_func_new=NULL, add_shifts=TRUE, ...) {
    # Computes quantiles of the likelihood predictive distribution pointwise
    # over a set of inputs. `p` is the probability defining the quantile.
    
    input <- .self$get_input(input)
    
    # Call child method.
    q <- .self$.calc_lik_quantiles(input, p, lik_par_val=lik_par_val, 
                                   em_pred_list=em_pred_list, log_scale=log_scale, 
                                   conditional=conditional, 
                                   normalize=normalize, 
                                   lower_tail=lower_tail, ...)
    
    # Apply deterministic shift. 
    if(apply_shift) {
      s <- get_shift_func(shift_func_new=shift_func_new, add_shifts=add_shifts)
      if(!is.null(s)) {
        if(log_scale) q <- q + s(input)
        else q <- q * exp(s(input))
      }
    }
    
    return(q)
  },
  
  
  .assemble_llik = function(em_val, lik_par_val=NULL, 
                           conditional=default_conditional, 
                           normalize=default_normalize, ...) {
    .NotYetImplemented()
  },
  
  .sample = function(input, lik_par_val=NULL, N_samp=1L,
                     em_pred_list=NULL, llik_pred_list=NULL,
                     conditional=default_conditional,
                     normalize=default_normalize, ...) {
    .NotYetImplemented()  
  },
  
  .predict = function(input, lik_par_val=NULL, return_mean=TRUE, return_var=TRUE, 
                     return_cov=FALSE, return_cross_cov=FALSE, input_cross=NULL,
                     conditional=default_conditional, normalize=default_normalize, 
                     em_pred_list=NULL, ...) {
    .NotYetImplemented()
  },
  
  .predict_lik = function(input, lik_par_val=NULL, em_pred_list=NULL, return_mean=TRUE,  
                          return_var=TRUE, return_cov=FALSE, return_cross_cov=FALSE, 
                          input_cross=NULL, conditional=default_conditional,  
                          normalize=default_normalize, log_scale=FALSE, ...) {
    .NotYetImplemented()
  },
  
  
  calc_func = function(func, input=NULL, em_pred_list=NULL, lik_par_val=NULL,
                       conditional=default_conditional, 
                       normalize=default_normalize, ...) {
    # `func` can either be a function with signature `func(em_pred_list, ...)`
    # that returns a vector of function evaluations at each input associated 
    # with `em_pred_list`. Alternatively, `func` can be a string specifying 
    # a class method. `calc_lik_approx()` and `calc_err_metric()` are two 
    # higher-level wrappers around `calc_func()`. Note that in the case 
    # that `func` is a string specifying the method name, the method must
    # be included in the call to `usingMethods()` below or it won't be 
    # available in the scope of this function.
    #
    # See description on usingMethods here: 
    # https://docs.tibco.com/pub/enterprise-runtime-for-R/6.0.0/doc/html/Language_Reference/methods/setRefClass.html
    
    # Note that for some reason this doesn't work if called like `.self$usingMethods()`. 
    usingMethods(calc_lik_marginal_approx, calc_lik_mean_approx, 
                 calc_lik_quantile_approx, calc_lik_sample_approx)

    # Select method to use. If `func` is character, ensure that it is a valid 
    # method.
    assert_that(is.function(func) || is.character(func))
    if(is.character(func)) {
      assert_that(func %in% names(.self), 
                  msg=paste0("Method not found: ", func))
      func <- .self[[func]]
    }
    
    # Compute emulator predictions, if not provided. No adjustments are made
    # here. Methods assume that passed in emulator predictions are unadjusted.
    # if(is.null(em_pred_list)) {
    #   em_pred_list <- .self$predict_emulator(input, lik_par_val=lik_par_val, 
    #                                          adjustment=NULL, ...)
    # }
    
    # Evaluate the function at the inputs.
    func(em_pred_list=em_pred_list, input=input, lik_par_val=lik_par_val, 
         conditional=conditional, normalize=normalize, ...)
  },
  
  calc_multi_func = function(func_list, input=NULL, em_pred_list=NULL,  
                             lik_par_val=NULL, conditional=default_conditional, 
                             normalize=default_normalize,
                             return_type="list", func_names=NULL, ...) {
    # A wrapper around `calc_func()` that calls the latter method for each 
    # function in `func_list`. `return_type` can either be "list" or "matrix". 
    # "list" implies that each element will be the associated 
    # return from the `calc_func()` call. "matrix" will stack the output for 
    # each call into the columns of a matrix. `func_names` is an optional 
    # character vector of length equal to `length(func_list)` providing names 
    # to use in the returned list/matrix. Otherwise, default names are 
    # constructed.
    
    assert_that(return_type %in% c("matrix", "list"))
    if(!is.null(func_names)) {
      assert_that(is.character(func_names))
      assert_that(length(func_names) == length(func_list))
    }
    
    # Compute emulator predictions, if not provided.
    # TODO: removing this for now. The treatment of existing pred lists is currently inconsistent.
    # if(is.null(em_pred_list)) {
    #   em_pred_list <- .self$predict_emulator(input, lik_par_val=lik_par_val, ...)
    # }
    
    # Return function evaluations for each function in list.
    l <- lapply(func_list, function(f) .self$calc_func(f, input=input, 
                                                       lik_par_val=lik_par_val,
                                                       em_pred_list=NULL, # TODO: setting this to NULL for now (see above) 
                                                       conditional=conditional, 
                                                       normalize=normalize, ...))
    
    # Obtain names identifying each function.
    if(is.null(func_names)) {
      if(is.character(func_list)) func_names <- func_list
      else if(!is.null(names(func_list))) func_names <- names(func_list)
      else func_names <- paste0("f", seq_along(func_list))
    }

    # Return as list, if requested.
    if(return_type == "list") {
      names(l) <- func_names
      return(l)
    }
    
    # Otherwise stack into columns of matrix.
    mat <- do.call(cbind, l)
    colnames(mat) <- func_names
    
    return(mat)
  },
  
  calc_lik_approx = function(approx_type, input=NULL, em_pred_list=NULL,
                             lik_par_val=NULL, conditional=default_conditional, 
                             normalize=default_normalize, log_scale=TRUE, 
                             return_type="list", simplify=FALSE, ...) {
    # A wrapper around `calc_multi_func()` that computes functions of the  
    # emulator predictive distribution that are interpreted as deterministic 
    # approximations of the likelihood function (or the log of such 
    # deterministic predictions). The string `approx_type` is used to specify 
    # a method of the form `calc_lik_<approx_type>_approx()`. The argument 
    # `approx_type` may be a vector of multiple approximation types. In this 
    # case `return_type` controls the type of the returned object; see 
    # `calc_multi_func()` for details. If `length(approx_type)==1` and 
    # `simplify = TRUE` then the return value from the single call to 
    # `calc_func()` will be returned, without constructing a matrix or list.

    assert_that(is.character(approx_type))
    simplify_output <- (length(approx_type)==1L) && isTRUE(simplify)
    if(simplify_output) return_type <- "list"
    
    # Select methods to use.
    method_names <- paste0("calc_lik_", approx_type, "_approx")

    # Evaluate the function at the inputs.
    func_evals <- .self$calc_multi_func(func_list=method_names, input=input, 
                                        em_pred_list=em_pred_list, 
                                        lik_par_val=lik_par_val,
                                        conditional=conditional, 
                                        normalize=normalize,
                                        log_scale=log_scale, 
                                        return_type=return_type,
                                        func_names=approx_type, ...)
    
    if(simplify_output) return(func_evals[[1]])
    return(func_evals)
  },
  
  get_design_llik = function(lik_par_val=NULL, conditional=default_conditional, 
                             normalize=default_normalize, apply_shift=TRUE,
                             shift_func_new=NULL, add_shifts=TRUE, ...) {
    # Calls child class method `.get_design_llik()`, and then adds the shift
    # function. Note that it is possible that `.get_design_llik()` returns
    # NULL, since some log-likelihood emulators won't have design points.
    # Applying the shift requires access to the design inputs, so that the
    # shift can be computed at those inputs.
    
    lik_par_val <- .self$get_lik_par(lik_par_val)
    llik_vals <- .self$.get_design_llik(lik_par_val=lik_par_val, 
                                        conditional=conditional, 
                                        normalize=normalize)
    if(is.null(llik_vals)) return(NULL)
    
    # Applying the shift requires fetching design inputs.
    if(apply_shift) {
      s <- get_shift_func(shift_func_new=shift_func_new, add_shifts=add_shifts)
      if(!is.null(s)) {
        design_inputs <- .self$get_design_inputs()
        if(is.null(design_inputs)) {
          stop("`design_inputs` is NULL, so `get_design_llik()` cannot computed shift function.")
        }
        
        llik_vals <- llik_vals + s(design_inputs)
      }
    }
    
    return(llik_vals)
  },
  
  calc_lik_approx_pw_err = function(llik_true, approx_type="mean", err_type="mse",  
                                    em_pred_list=NULL, input=NULL, lik_par_val=NULL, 
                                    conditional=default_conditional,
                                    normalize=default_normalize, llik_pred=NULL,
                                    return_type="list", ...) {
    # A high-level convenience function to compute pointwise ("pw") errors 
    # between likelihood approximations and true baseline log-likelihood 
    # values `llik_true`. The errors are all computed on the log scale. 
    # Both `approx_type` and `err_type` are allowed to be vectors with 
    # multiple options. Currently `return_type` is allowed to be either "list"
    # or "data.table". The former returns a list with one element per value 
    # in `err_type`, with each element being a matrix storing the errors for 
    # each `approx_type` value. If "data.table", returns a data.table with 
    # columns "approx_type", "err_type", and "value". The "weighted" error 
    # measures "wmse" and "wmae" divide the pointwise absolute errors by 
    # the log-likelihood predictive standard deviations (for "wmse" the 
    # result is then squared). `llik_pred` allows passing in previously
    # computed likelihood approximations, which should be on the log 
    # scale. If provided, should be a matrix with one column per approximation
    # type.
    
    # Ensure valid error and return types.
    assert_that(all(err_type %in% c("mse", "mae", "wmse", "wmae")))
    assert_that(return_type %in% c("list", "data.table"))
    
    # If "weighted" error types are requested, then need to compute predictive 
    # variances of log-likelihood emulator.
    llik_pred_inv_sd <- NULL
    if(any(c("wmse", "wmae") %in% err_type)) {
      llik_pred_list <- .self$predict(input=input, return_mean=FALSE, 
                                      return_var=TRUE, conditional=conditional,
                                      normalize=normalize, 
                                      em_pred_list=em_pred_list)
      llik_pred_inv_sd <- 1 / drop(sqrt(llik_pred_list$var))
    }
    
    # Compute plug-in mean log-likelihood predictions.
    if(is.null(llik_pred)) {
      llik_pred <- .self$calc_lik_approx(approx_type, em_pred_list=em_pred_list,
                                         input=input, lik_par_val=lik_par_val,
                                         conditional=conditional, 
                                         normalize=normalize, return_type="matrix", 
                                         simplify=FALSE, ...)
    }
    assert_that(nrow(llik_pred) == length(drop(llik_true)))
    
    # Compute pointwise errors.
    err_list <- list()
    diff <- add_vec_to_mat_cols(-drop(llik_true), llik_pred)
    if("mse" %in% err_type) err_list$mse <- diff^2
    if("mae" %in% err_type) err_list$mae <- abs(diff)
    if("wmse" %in% err_type) err_list$wmse <- mult_vec_with_mat_cols(llik_pred_inv_sd, diff)^2
    if("wmae" %in% err_type) err_list$mae <- abs(mult_vec_with_mat_cols(llik_pred_inv_sd, diff))
    
    # If requested, return as list.
    if(return_type == "list") return(err_list)
    
    # Otherwise convert to data.table.
    for(i in seq_along(err_list)) {
      err_list[[i]] <- as.data.table(err_list[[i]])
      err_list[[i]][["err_type"]] <- names(err_list)[i]
    }
    
    dt <- data.table::rbindlist(err_list, use.names=TRUE)
    dt <- data.table::melt.data.table(dt, id.vars="err_type", 
                                      variable.name="approx_type", 
                                      value.name="value")
    return(dt)
  },
  
  get_llik_func = function(approx_type=NULL, lik_par_val=NULL, 
                           conditional=default_conditional, 
                           normalize=default_normalize, ...) {
    # Returns a function that represents a (deterministic) log-likelihood 
    # function. This is useful for interfacing with other software packages; 
    # e.g., MCMC functions that require a log-likelihood argument. The 
    # deterministic log-likelhoods are constructed via `calc_lik_approx()`.

    # If llikEmulator object encodes an exact log-likelihood, then return the 
    # true log-likelihood function.
    if(.self$exact_llik) {
      llik <- function(input) {
        if(is.null(dim(input))) input <- matrix(input, nrow=1L)
        .self$assemble_llik(input, conditional=conditional, 
                            normalize=normalize, ...)
      }
    } else {
      # Otherwise, returns an approximate log-likelihood.
      assert_that(!is.null(approx_type), 
                  msg="`approx_type` required when `exact_llik` is FALSE.")
      
      llik <- function(input) {
        if(is.null(dim(input))) input <- matrix(input, nrow=1)
        .self$calc_lik_approx(approx_type=approx_type, input=input, 
                              lik_par_val=lik_par_val, conditional=conditional, 
                              normalize=normalize, log_scale=TRUE, 
                              simplify=TRUE, ...)
      }
    }
    
    return(llik)
  },
  
  
  calc_lik_marginal_approx = function(em_pred_list=NULL, input=NULL, 
                                      lik_par_val=NULL, llik_pred_list=NULL,
                                      conditional=default_conditional,
                                      normalize=default_normalize, 
                                      include_noise=TRUE, log_scale=TRUE, ...) {
    # `Marginal approx` is defined to be the expectation of the likelihood 
    # surrogate (on a pointwise input-by-input basis). This implements a 
    # default method that will work for most llikEmulator classes. 
    # This default method can be written for llikEmulator classes with special 
    # structure that differs from this. 
    
    # TODO: setting previous pred lists to NULL for now until the assumptions on
    # existing prediction lists are solidified.
    em_pred_list <- NULL
    llik_pred_list <- NULL
    
    sel <- ifelse(log_scale, "log_mean", "mean")
    .self$predict_lik(input, lik_par_val=lik_par_val, em_pred_list=em_pred_list, 
                      llik_pred_list=llik_pred_list, return_mean=TRUE, 
                      return_var=FALSE, conditional=conditional,  
                      normalize=normalize, log_scale=log_scale, ...)[[sel]]
  },
  
  calc_lik_sample_approx = function(em_pred_list=NULL, input=NULL, 
                                    lik_par_val=NULL, 
                                    conditional=default_conditional,
                                    normalize=default_normalize, 
                                    include_noise=TRUE, log_scale=TRUE, ...) {
    .NotYetImplemented()  
  },
  
  calc_lik_quantile_approx = function(em_pred_list=NULL, input=NULL, alpha=0.9,  
                                      lik_par_val=NULL,
                                      conditional=default_conditional, 
                                      normalize=default_normalize,
                                      log_scale=TRUE, ...) {
    # Deterministic likelihood approximation that is given by the alpha
    # quantile of the likelihood surrogate.
    
    # TODO: setting `em_pred_list=NULL` for now until the assumptions on
    # existing prediction lists are solidified.
    .self$calc_lik_quantiles(alpha, input=input, lik_par_val=lik_par_val, 
                             em_pred_list=NULL, log_scale=log_scale, 
                             conditional=conditional, 
                             normalize=normalize, 
                             lower_tail=TRUE, include_noise=include_noise, ...)
  },
  
  calc_lik_mean_approx = function(em_pred_list=NULL, input=NULL, 
                                  lik_par_val=NULL,
                                  conditional=default_conditional,
                                  normalize=default_normalize, log_scale=TRUE, 
                                  apply_shift=TRUE, shift_func_new=NULL,
                                  add_shifts=TRUE, ...) {
    # `Mean approx` means that the `emulator_model` predictive mean is computed
    # at inputs `input`, then the predictive mean is passed to `assemble_llik()` 
    # and the result is exponentiated (if `log_scale` is FALSE). This implements 
    # the "plug-in emulator mean" approximation. This default method can be 
    # written for llikEmulator classes with special structure that differs 
    # from this. `input` is required to apply the shift, even if `em_pred_list`
    # is passed.
    
    # TODO: setting `em_pred_list=NULL` for now until the assumptions on
    # existing prediction lists are solidified.
    em_pred_list <- NULL
    
    # Compute emulator predictions, if not provided.
    if(is.null(em_pred_list)) {
      em_pred_list <- .self$predict_emulator(input, lik_par_val=lik_par_val, 
                                             return_var=FALSE, ...)
    } else {
      assert_that(!is.null(em_pred_list$mean))
    }
    
    # Evaluate log-likelihood by plugging in the emulator mean predictions.
    llik_pred <- .self$assemble_llik(em_pred_list$mean, lik_par_val=lik_par_val, 
                                     conditional=conditional, 
                                     normalize=normalize, ...)
    
    # Apply shift.
    if(apply_shift) {
      s <- get_shift_func(shift_func_new=shift_func_new, add_shifts=add_shifts)
      if(!is.null(s)) llik_pred <- llik_pred + s(input)
    }
    
    if(log_scale) return(llik_pred)
    return(exp(llik_pred))
  },
  
  
  calc_expected_lik_cond_var = function(input_eval, input_cond, 
                                        include_noise=TRUE, log_scale=TRUE, 
                                        plugin=FALSE, ...) {
    # Under the default setting `plugin = FALSE`, computes the log of 
    # E_l Var[exp(L(input_eval))|L(input_cond)=l]
    # where L is the log-likelihood emulator and the expectation is with respect 
    # to the current log-likelihood emulator predictive distribution (i.e., this 
    # is a posterior predictive quantity). If `plugin = TRUE`, then the function 
    # will instead compute E_g Var[exp(L(input_eval; G|{input_cond,g}))|g]
    # where `G` denotes the underlying emulator `emulator_model`. In words, the 
    # underlying emulator will first be conditioned at the new inputs 
    # `input_cond`- G|{input_cond,g} - where `g` is assumed to be distributed 
    # according to the current emulator model predictive distribution. This 
    # conditioned emulator is plugged into the log-likelihood and then the 
    # expectation of the exponentiated log-likelihood is computed with respect 
    # to `g`. Certain classes may implement this method with only one of the 
    # two options for `plugin`, in which case an error should be thrown if 
    # the other option is passed. Other classes may not have this method at all, 
    # or return an approximation to this expected conditional variance.
    #
    # Returns vector of length `length(X_eval)` containing the expected 
    # conditional variance at the evaluation points `input_eval`. 
    .NotYetImplemented()
  },
  
  get_pred_interval = function(input, lik_par_val=NULL, em_pred_list=NULL, 
                               target_pred_list=NULL, target="llik",  
                               log_scale=TRUE, method="CI", N_std_dev=2L, 
                               CI_prob=0.9, conditional=default_conditional, 
                               normalize=default_normalize, 
                               include_noise=TRUE, ...) {
    # Options for `target`: "llik" and "lik". 
    # Options for `method`: "pm_std_dev" (pm = "plus-minus"), "CI". The former 
    # computes the bounds by adding and subtracting `N_std_dev` standard 
    # deviations from the predictive mean. The latter computes a 100*`CI_prob`% 
    # confidence interval. `target_pred_list` is either the pred list returned 
    # by `predict` or `predict_lik`; which one it is should align with the value 
    # of `target`. 

    if(log_scale && (target=="lik") && (method=="CI")) {
      stop("Method 'CI' not yet implemented when `log_scale=TRUE` and `target='lik'`.")
    }
    
    if(target=="llik") log_scale <- TRUE
    
    assert_that(target %in% c("llik", "lik"))
    assert_that(method %in% c("pm_std_dev", "CI"))
    interval_list <- list()
    
    # Log likelihood or likelihood predictions. 
    mean_sel <- "mean"
    var_sel <- "var"
    if((target=="lik") && log_scale) {
      mean_sel <- "log_mean"
      var_sel <- "log_var"
    }
    
    if(is.null(target_pred_list)) {
      if(target == "llik") {
        target_pred_list <- .self$predict(input, lik_par_val=lik_par_val, 
                                          em_pred_list=em_pred_list, 
                                          return_mean=TRUE, return_var=TRUE, 
                                          conditional=conditional,
                                          normalize=normalize, ...)
      } else {
        target_pred_list <- .self$predict_lik(input, lik_par_val=lik_par_val, 
                                              em_pred_list=em_pred_list, 
                                              return_mean=TRUE, return_var=TRUE, 
                                              conditional=conditional,
                                              normalize=normalize, log_scale=log_scale, ...)
      }
    } else {
      assert_that(!is.null(target_pred_list[[mean_sel]]) && 
                  !is.null(target_pred_list[[var_sel]]))
    }
    
    # Plus/minus standard deviation method 
    if(method == "pm_std_dev") {
      m <- drop(target_pred_list[[mean_sel]])
      v <- drop(target_pred_list[[var_sel]])
      
      if((target=="lik") && log_scale) {
        log_s <- log(N_std_dev) +  0.5*v
        interval_list$lower <- log_diff_exp(m, log_s)
        interval_list$upper <- matrixStats::rowLogSumExps(cbind(m, log_s))
      } else {
        interval_list$lower <- m - N_std_dev * sqrt(v)
        interval_list$upper <- m + N_std_dev * sqrt(v)
      }
      
    } else {
      CI_tail_prob <- 0.5 * (1-CI_prob)
      interval_list$lower <- .self$calc_quantiles(p=CI_tail_prob, input=input, 
                                                  lik_par_val=lik_par_val, 
                                                  em_pred_list=em_pred_list, 
                                                  target=target, 
                                                  conditional=conditional, 
                                                  normalize=normalize, 
                                                  lower_tail=TRUE, 
                                                  include_noise=include_noise, ...)
      interval_list$upper <- .self$calc_quantiles(p=CI_tail_prob, input=input, 
                                                  lik_par_val=lik_par_val, 
                                                  em_pred_list=em_pred_list, 
                                                  target=target, 
                                                  conditional=conditional,
                                                  normalize=normalize, 
                                                  lower_tail=FALSE, 
                                                  include_noise=include_noise, ...)
    }
    
    return(interval_list)
  },
  
  get_design_inputs = function(...) {
    # By default returns NULL, since not all log-likelihood emulators may have 
    # design points. When desinging classes that inherit from `llikEmulator`, 
    # this method should be overwritten if the log-likelihood emulator has 
    # design points. 
    
    return(NULL)
  },
  
  
  .get_design_llik = function(lik_par_val=NULL, conditional=default_conditional, 
                              normalize=default_normalize, ...) {
    # Child class method defaults to returning NULL since not all log-likelihood
    # emulators have design points. Use `get_design_llik()`, which calls
    # `.get_design_llik()`.
    return(NULL)
  },
  
  plot_samp_1d = function(input, lik_par_val=NULL, em_pred_list=NULL, N_samp=1, 
                          plot_type="llik", conditional=default_conditional, 
                          normalize=default_normalize, true_llik=NULL, 
                          include_design=TRUE, use_cov=TRUE, 
                          ground_truth_col="black", design_col="black", 
                          include_bounds=FALSE, ...) {
    # `plot_type` options: "llik", "lik". If `include_bounds = TRUE`, then the llik/lik
    # bounds will be plotted as horizontal lines (if there are bounds).
    
    assert_that(dim_input==1, 
                msg=paste0("plot_llik_samp_1d() requires 1d input space. input_dim = ", dim_input))
    assert_that(plot_type %in% c("llik", "lik"))
    
    input <- get_input(input)
    samp <- .self$sample(input, lik_par=lik_par_val, em_pred_list=em_pred_list,
                         N_samp=N_samp, use_cov=use_cov, ...)
    
    # Adjustments in plotting likelihood (not log-likelihood). 
    if(plot_type == "lik") {
      samp <- exp(samp)
      base_plot_title <- "Likelihood"
      true_vals <- exp(true_llik)
    } else {
      base_plot_title <- "Log Likelihood"
      true_vals <- true_llik
    }
    
    plt <- ggmatplot(input, samp, plot_type="line", color="gray") + 
                  theme(legend.position = "none") + 
                  ggtitle(paste0(base_plot_title, " Samples")) + 
                  xlab(input_names) + ylab(paste0(base_plot_title, ": ", llik_label))

    if(!is.null(true_llik)) {
      df <- data.frame(x=input[,1], y=drop(true_vals))
      plt <- plt + geom_line(aes(x=x, y=y), df, inherit.aes=FALSE, 
                             color=ground_truth_col)
    }
    
    if(include_design) {
      design_response_vals <- drop(get_design_llik(lik_par_val, conditional, 
                                                   normalize))
      if(plot_type == "lik") design_response_vals <- exp(design_response_vals)
      
      design_df <- data.frame(x=drop(get_design_inputs()), 
                              y=design_response_vals)
      plt <- plt + geom_point(aes(x=x, y=y), design_df, inherit.aes=FALSE, 
                              color=design_col)
    }
    
    if(include_bounds) {
      plt <- .self$plot_llik_bounds(plt, plot_type)
    }
    
    return(plt)
  },
  
  plot_pred_1d = function(input, lik_par_val=NULL, em_pred_list=NULL, 
                          pred_list=NULL, plot_type="llik", log_scale=TRUE, 
                          conditional=default_conditional, 
                          normalize=default_normalize, include_interval=TRUE, 
                          interval_method="pm_std_dev", N_std_dev=1, 
                          CI_prob=0.9, true_llik=NULL, include_design=TRUE, 
                          xlab=NULL, ylab=NULL, plot_title=NULL, 
                          include_bounds=FALSE, ...) {
    # `pred_list` can be passed if llik or lik predictions have already been 
    # computed. In this case, it should be a llik prediction list (as returned
    # by `predict()`) if `plot_type = llik` and a likelihood prediction list
    # (as returned by `predict_lik()`) if `plot_type = lik`.
    
    if(plot_type=="llik") log_scale <- TRUE
    
    assert_that(dim_input==1, 
                msg=paste0("plot_llik_pred_1d() requires 1d input space. dim_input = ", dim_input))
    assert_that(plot_type %in% c("llik", "lik"))
    
    input <- .self$get_input(input)
    
    # Default axis labels. 
    if(is.null(xlab)) xlab <- input_names
    if(is.null(ylab)) ylab <- plot_type
    
    # Compute required predictive quantities if not already provided. 
    true_vals <- true_llik
    mean_sel <- "mean"
    var_sel <- "var"
    
    if(is.null(pred_list)) {
      if(plot_type == "llik") {
        pred_list <- .self$predict(input, lik_par_val=lik_par_val, 
                                   em_pred_list=em_pred_list, return_mean=TRUE, 
                                   return_var=include_interval, 
                                   conditional=conditional, 
                                   normalize=normalize, ...)
      } else {
        pred_list <- .self$predict_lik(input, lik_par_val=lik_par_val, 
                                       em_pred_list=em_pred_list,  
                                       return_mean=TRUE, 
                                       return_var=include_interval, 
                                       conditional=conditional, 
                                       normalize=normalize, log_scale=log_scale, ...)
      }
    }
    
    # Adjustments for likelihood plot.
    if(plot_type == "lik") {
      if(!is.null(true_vals) && !log_scale) true_vals <- exp(true_vals)
      if(log_scale) {
        mean_sel <- "log_mean"
        var_sel <- "log_var"
      }
    }
    
    # Plot title and labels.
    if(is.null(plot_title)) {
      plot_title <- "Likelihood Emulator Predictions"
      if(log_scale) plot_title <- paste0("Log ", plot_title)
      if(include_interval && (interval_method=="CI")) {
        plot_title <- paste0(plot_title, ", ", 100*CI_prob, "% CI")
      }
      if(include_interval && (interval_method=="pm_std_dev")) {
        plot_title <- paste0(plot_title, ", +/- ", N_std_dev, " std dev")
      }
    }
    
    if(!normalize) {
      ylab <- paste0(ylab, ", ", ifelse(conditional, "unnormalized/conditional", 
                                        "unnormalized"))
    }

    # Compute prediction interval.  
    if(include_interval) {
      interval_list <- .self$get_pred_interval(input, lik_par_val=lik_par_val, 
                                               target_pred_list=pred_list, 
                                               target=plot_type, 
                                               log_scale=log_scale,
                                               method=interval_method, 
                                               N_std_dev=N_std_dev,
                                               CI_prob=CI_prob, 
                                               conditional=conditional, 
                                               normalize=normalize, ...)
    } else {
      interval_list <- NULL
    }

    # Design points. 
    design_inputs <- NULL
    design_response_vals <- NULL 
    if(include_design) {
      design_inputs <- drop(.self$get_design_inputs(...))
      design_response_vals <- drop(get_design_llik(lik_par_val, conditional, 
                                                   normalize, ...))
      if(!log_scale) design_response_vals <- exp(design_response_vals)
    } 
    
    # Produce plot. 
    plt <- plot_pred_1d_helper(X_new=drop(input), pred_mean=drop(pred_list[[mean_sel]]), 
                               include_CI=include_interval, 
                               CI_lower=interval_list$lower,  
                               CI_upper=interval_list$upper, 
                               y_new=drop(true_vals), X_design=design_inputs, 
                               y_design=design_response_vals, 
                               plot_title=plot_title, xlab=xlab, ylab=ylab, ...)
    
    # Optionally plot (log) likelihood bounds.
    if(include_bounds) {
      plt <- .self$plot_llik_bounds(plt, plot_type)
    }
    
    return(plt)                     
  }, 
  
  
  plot_pred_validation = function(input, true_llik, lik_par_val=NULL, 
                                  em_pred_list=NULL, plot_type="llik", 
                                  conditional=default_conditional, 
                                  normalize=default_normalize, 
                                  include_interval=TRUE, 
                                  interval_method="CI", N_std_dev=2L, 
                                  CI_prob=0.9, include_design=TRUE, xlab=NULL, 
                                  ylab=NULL, plot_title=NULL, ...) {
    
    assert_that(plot_type %in% c("llik", "lik"))
    
    # Default plot title. 
    if(is.null(plot_title)) {
      if(plot_type == "llik") plot_title <- "Log-Likelihood Predictions"
      if(plot_type == "lik") plot_title <- "Likelihood Predictions"
    }
    
    input <- .self$get_input(input)
    
    # Compute required predictive quantities if not already provided. 
    true_vals <- true_llik
    if(plot_type == "llik") {
      pred_list <- .self$predict(input, lik_par_val=lik_par_val, 
                                 em_pred_list=em_pred_list, return_mean=TRUE, 
                                 return_var=include_interval, 
                                 conditional=conditional, 
                                 normalize=normalize, ...)
    } else {
      pred_list <- .self$predict_lik(input, lik_par_val=lik_par_val, 
                                     em_pred_list=em_pred_list, return_mean=TRUE, 
                                     return_var=include_interval, 
                                     conditional=conditional, 
                                     normalize=normalize, log_scale=FALSE, ...)
      if(!is.null(true_vals)) true_vals <- exp(true_vals)
    }
    
    # Compute prediction interval.  
    if(include_interval) {
      interval_list <- .self$get_pred_interval(input, lik_par_val=lik_par_val, 
                                               target_pred_list=pred_list, 
                                               target=plot_type, 
                                               method=interval_method, 
                                               N_std_dev=N_std_dev,
                                               CI_prob=CI_prob, 
                                               conditional=conditional, 
                                               normalize=normalize, ...)
    } else {
      interval_list <- NULL
    }
    
    # Design points. 
    design_inputs <- NULL
    design_response_vals <- NULL 
    if(include_design) {
      design_response_vals <- drop(.self$get_design_llik(lik_par_val, conditional, 
                                                         normalize, ...))
      if(plot_type == "lik") design_response_vals <- exp(design_response_vals)
    } 
    
    # Produce plot. 
    plt <- plot_true_pred_scatter(y_pred=drop(pred_list$mean), 
                                  y_true=drop(true_vals), 
                                  include_CI=include_interval, 
                                  CI_lower=interval_list$lower,  
                                  CI_upper=interval_list$upper, 
                                  y_design=design_response_vals, 
                                  plot_title=plot_title, xlab=xlab, 
                                  ylab=ylab, ...)
    
    return(plt)
  }, 
  
  plot_llik_bounds = function(plt, plot_type="llik") {
    # Takes an existing ggplot object `plt` and adds horizontal lines
    # corresponding to the (log)-likelihood bounds returned by 
    # `get_llik_bounds()`. Bounds that are -Inf or Inf are not 
    # plotted. If plot type is "lik" then the bounds are 
    # exponentiated.
  
    bounds <- .self$get_llik_bounds()
    if(plot_type == "lik") bounds <- exp(bounds)
      
    if(is.finite(bounds[1])) {
      plt <- plt + geom_hline(yintercept=bounds[1], 
                              color="green", linetype="dashed")
    }
      
    if(is.finite(bounds[2])) {
      plt <- plt + geom_hline(yintercept=bounds[2], 
                              color="green", linetype="dashed")
    }
    
    return(plt)
  }, 
  
  plot_1d_proj_approx_lik = function(input_names_proj=NULL, n_points=100L, 
                                     input_list_proj=NULL, input_fixed=NULL, 
                                     input_grids=NULL, lik_par_val=NULL, 
                                     llik_func_true=NULL,
                                     approx_type="mean", log_scale=TRUE, 
                                     input_bounds=NULL,  
                                     conditional=default_conditional, 
                                     normalize=default_normalize, ...) {
    # This method is very similar to the gpWrapper method of the same name (see 
    # that method for detailed comments). The major differences here include the 
    # fact that llikEmulator's are not required to have design points; if this 
    # is the case then default input bounds for the inputs cannot be defined so 
    # either `input_bounds` or `input_list_proj` and `input_fixed` must be 
    # explicitly passed. A list of length length(input_names_proj)` is returned 
    # containing the plots.
    #
    # If `input_grids` is passed, then no new grid points are created. `input_grids`
    # should be a list as returned by `get_input_grid_1d_projection()`.
    #
    # TODO: should include the option to use the log of the likelihood predictions.
    # Currently, `log_scale` only applies to the `lik_approx` type but it should 
    # also apply to `lik`.
    
    if(length(approx_type) > 1) .NotYetImplemented()

    # If provided, `llik_func_true` must either be a function (representing the 
    # true log-likelihood) or a `llikEmulator` class representing an exact 
    # log-likelihood.
    if(!is.null(llik_func_true)) {
      is_llik_em_obj <- is_llik_em(llik_func_true)
      assert_that(is.function(llik_func_true) || is_llik_em_obj)
      if(is_llik_em_obj) {
        assert_that(llik_func_true$exact_llik)
        llik_func_true <- llik_func_true$get_llik_func()
      }
    }
    
    # Determine if the y-axis will be on the likelihood or log-likelihood scale. 
    ylab <- ifelse(log_scale, "log likelihood", "likelihood")
    
    if(is.null(input_grids)) {
      # If provided, ensure `input_bounds` has colnames set to the input names. 
      if(!is.null(input_bounds)) {
        assert_that(setequal(colnames(input_bounds), .self$input_names))
        input_bounds <- input_bounds[, .self$input_names]
      }
      
      # Determine the inputs that will be varied. 
      if(is.null(input_names_proj)) {
        input_names_proj <- .self$input_names
      } else {
        input_names_proj <- unique(input_names_proj)
        assert_that(all(is.element(input_names_proj, .self$input_names)))
      }
  
      # Constructing input grids for varied parameters.
      input_grids <- get_input_grid_1d_projection(.self$input_names, 
                                                  x_vary=input_names_proj, 
                                                  X_list=input_list_proj, 
                                                  X_fixed=input_fixed,
                                                  X_bounds=input_bounds, 
                                                  n_points_default=n_points)
    }
    
    # Helper function for computing response values at the input points given
    # in `input_grids[[i]]`, which is a list containing one set of input 
    # points per set of fixed values.
    f <- function(i,j) {
      .self$calc_lik_approx(approx_type, input=input_grids[[i]][[j]],
                            lik_par_val=lik_par_val, conditional=conditional, 
                            normalize=normalize, log_scale=log_scale, 
                            return_type="matrix", ...)
    }
    
    compute_response <- function(i) {
      l <- lapply(seq_along(input_grids[[i]]), function(j) f(i,j))
      do.call(cbind, l)
    }
    
    # Create plots, loop over each input that is varied.
    plts <- list()
    n_fixed <- .self$dim_input - 1L
    for(i in seq_along(input_grids)) {
      input_name <- names(input_grids)[i]
      
      # Compute response values to be plotted, one per set of fixed values.
      vals <- compute_response(i)
      colnames(vals) <- paste0("fixed", 1:ncol(vals))

      # Optionally compute true baseline values.
      true_response <- NULL
      if(!is.null(llik_func_true)) {
        true_response <- llik_func_true(inputs)
        if(!log_scale) true_response <- exp(true_response)
      }
      
      # Construct plot.
      input_1d <- drop(input_grids[[input_name]][[1]][,input_name])
      plts[[input_name]] <- plot_curves_1d_helper(input_1d, pred=vals, 
                                                  y_new=true_response,
                                                  xlab=input_name, ylab=ylab)
    }
    
    return(plts)
  }, 
  
  
  plot_1d_projection = function(input_names_proj=NULL, n_points=100L, 
                                input_list_proj=NULL, input_fixed=NULL, 
                                lik_par_val=NULL, include_design=FALSE, 
                                llik_func_true=NULL, include_interval=TRUE,   
                                interval_method="pm_std_dev", N_std_dev=2, 
                                CI_prob=0.9, include_noise=TRUE, 
                                plot_type="llik", approx_type="mean", 
                                log_scale=TRUE, input_bounds=NULL, 
                                conditional=default_conditional, 
                                normalize=default_normalize, ...) {
    # This method is very similar to the gpWrapper method of the same name (see 
    # that method for detailed comments). The major differences here include the 
    # fact that llikEmulator's are not required to have design points; if this 
    # is the case then default input bounds for the inputs cannot be defined so 
    # either `input_bounds` or `input_list_proj` and `input_fixed` must be 
    # explicitly passed. A list of length `length(input_names_proj)` is returned 
    # containing the plots.
    #
    # `plot_type` can be "llik", "lik", or "lik_approx". If "lik_approx", the 
    # specific likelihood approximation is determined by "approx_type", and 
    # `log_scale` controls whether or not the log of the likelihood approximation 
    # is plotted. `approx_type` can be a vector of multiple approximation types. 
    #
    # TODO: should include the option to use the log of the likelihood predictions.
    # Currently, `log_scale` only applies to the `lik_approx` type but it should 
    # also apply to `lik`.
    
    assert_that(is.element(plot_type, c("llik", "lik", "lik_approx")))
    if(plot_type=="lik") .NotYetImplemented()
    
    # If provided, `llik_func_true` must either be a function (representing the 
    # true log-likelihood) or a `llikEmulator` class representing an exact 
    # log-likelihood.
    if(!is.null(llik_func_true)) {
      is_llik_em_obj <- is_llik_em(llik_func_true)
      assert_that(is.function(llik_func_true) || is_llik_em_obj)
      if(is_llik_em_obj) {
        assert_that(llik_func_true$exact_llik)
        llik_func_true <- llik_func_true$get_llik_func()
      }
    }
    
    # Determine if the y-axis will be on the likelihood or log-likelihood scale. 
    if(plot_type=="llik") log_scale <- TRUE
    ylab <- ifelse(log_scale, "log likelihood", "likelihood")
    
    # If provided, ensure `input_bounds` has colnames set to the input names. 
    if(!is.null(input_bounds)) {
      assert_that(setequal(colnames(input_bounds), .self$input_names))
      input_bounds <- input_bounds[, .self$input_names]
    }
    
    # Code currently only supports fixing fixed parameters at a single value;  
    # i.e., only one line plotted per plot. 
    if(!is.null(input_fixed)) {
      assert_that(nrow(input_fixed)==1, 
                  msg="Currently only supports the case where `input_fixed` has a single row.")
    }
    
    # Determine the inputs that will be varied. 
    if(is.null(input_names_proj)) {
      input_names_proj <- .self$input_names
    } else {
      input_names_proj <- unique(input_names_proj)
      assert_that(all(is.element(input_names_proj, .self$input_names)))
    }
    
    # Design points. These are used if `include_design` is TRUE or as a means 
    # to determine default bounds in the case that `input_list_proj` does not 
    # define grids for all variables and `input_bounds` is not provided.
    design_response <- NULL
    design_inputs <- .self$get_design_inputs(...)
    if(is.null(input_bounds) && !is.null(design_inputs)) {
      input_bounds <- get_bounds(design_inputs)
    }
    
    if(include_design) {
      design_response <- drop(.self$get_design_llik(lik_par_val, 
                                                    conditional, normalize, ...))
      if(is.element(plot_type, c("lik","lik_approx")) && !log_scale) {
        design_response <- exp(design_response)
      }
    } else {
      design_inputs <- NULL
    }
    
    # Constructing input grids for varied parameters.
    # TODO: use helper function here.
    n_proj <- length(input_names_proj)
    if(is.null(input_list_proj)) input_list_proj <- list()
    for(input_name in input_names_proj) {
      if(is.null(input_list_proj[[input_name]])) {
        if(is.null(input_bounds)) {
          stop("Failed to construct default `input_bounds`; must explicitly pass this argument.")
        }
        input_bounds_curr <- input_bounds[,input_name]
        input_list_proj[[input_name]] <- seq(input_bounds_curr[1], 
                                             input_bounds_curr[2], 
                                             length.out=n_points)
      }
    }
    
    # If not provided, set values of fixed parameters to midpoints with respect  
    # to design bounds.
    if(is.null(input_fixed)) {
      if(is.null(input_bounds)) {
        stop("Failed to construct default `input_bounds`; must explicitly pass this argument.")
      }
      input_fixed <- matrix(apply(input_bounds, 2, mean), nrow=1)
      colnames(input_fixed) <- colnames(input_bounds)
    } 
    
    # Create plots.
    plts <- list()
    n_fixed <- .self$dim_input - 1L
    for(i in seq_len(n_proj)) {
      # Input matrix for prediction.
      input_name <- input_names_proj[i]
      input_names_fixed <- setdiff(.self$input_names, input_name)
      input_1d <- matrix(input_list_proj[[input_name]], ncol=1, 
                         dimnames=list(NULL,input_name))
      n_grid <- nrow(input_1d)
      inputs <- cbind(input_1d, matrix(input_fixed[,input_names_fixed,drop=FALSE], 
                                       ncol=n_fixed, nrow=n_grid, byrow=TRUE, 
                                       dimnames=list(NULL,input_names_fixed)))
      inputs <- inputs[, .self$input_names]

      # Compute log-likelihood or likelihood predictions.
      if(plot_type=="llik") {
        pred_list <- .self$predict(inputs, lik_par_val=lik_par_val, 
                                   return_var=include_interval, 
                                   include_noise=include_noise, 
                                   conditional=conditional, 
                                   normalize=normalize, ...)
        list_sel <- "mean"
      } else if(plot_type=="lik") {
        pred_list <- .self$predict_lik(inputs, lik_par_val=lik_par_val, 
                                       return_var=include_interval,  
                                       conditional=conditional, 
                                       normalize=normalize, log_scale=FALSE, ...)
        list_sel <- ifelse(log_scale, "mean", "log_mean")
      } else if(plot_type=="lik_approx") {
        pred_list <- list()
        list_sel <- "approx"
        lik_approx_mat <- .self$calc_lik_approx_comparison(inputs, approx_type, 
                                                           lik_par_val=lik_par_val,
                                                           conditional=conditional, 
                                                           normalize=normalize,
                                                           log_scale=log_scale, ...)
      }
      
      # Optionally compute true baseline values.
      true_response <- NULL
      if(!is.null(llik_func_true)) {
        true_response <- llik_func_true(inputs)
        if(is.element(plot_type, c("lik","lik_approx")) && !log_scale) {
          true_response <- exp(true_response)
        }
      }
      
      # Compute prediction interval. Only relevant for "llik" and "lik" plot types.
      if(include_interval && (plot_type != "lik_approx")) {
        interval_list <- .self$get_pred_interval(inputs, lik_par_val=lik_par_val, 
                                                 target_pred_list=pred_list, 
                                                 target=plot_type, 
                                                 method=interval_method, 
                                                 N_std_dev=N_std_dev,
                                                 CI_prob=CI_prob, 
                                                 conditional=conditional, 
                                                 normalize=normalize, ...)
      } else {
        interval_list <- NULL
        include_interval <- FALSE
      }
      
      # Construct plots.
      if(plot_type != "lik_approx") {
        plts[[input_name]] <- plot_pred_1d_helper(drop(input_1d), 
                                                  pred_mean=drop(pred_list[[list_sel]]), 
                                                  include_CI=include_interval, 
                                                  CI_lower=drop(interval_list$lower),  
                                                  CI_upper=drop(interval_list$upper), 
                                                  X_design=design_inputs[,input_name],
                                                  y_design=design_response, 
                                                  y_new=true_response, 
                                                  xlab=input_name, ylab=ylab, ...)
      } else {
        plts[[input_name]] <- plot_curves_1d_helper(drop(input_1d), lik_approx_mat,
                                                    X_design=design_inputs[,input_name], 
                                                    y_design=design_response,
                                                    y_new=true_response, 
                                                    xlab=input_name, ylab=ylab, 
                                                    legend=TRUE, ...)
                                                    
      }
    }
    
    return(plts)
  } 
)


is_llik_em <- function(model) {
  inherits(model, "llikEmulator")
}


# -----------------------------------------------------------------------------
# llikEmulatorGP class: 
# 
# Direct GP emulation of the log-likelihood. Currently only supports fixed 
# `lik_par`, but this should be generalized to allow GPs that predict the 
# llik as a function of both `par` and `lik_par`. This class is likelihood 
# agnostic - any log likelihood can be emulated by a GP and this class does
# not store any information about the underlying likelihood structure. For 
# this reason, there is no way to change the settings `normalize` and 
# `conditional` on the fly, as in most of the other llikEmulator classes. 
# These must be set when creating the class; e.g., if the GP emulator was
# fit to normalized log likelihood values, then the class should be 
# initialized with `default_conditional=TRUE`, `default_normalize=TRUE`. 
# Unlike in other llikEmulator classes, if the user calls a method and 
# tries to pass a value for `normalize` or `conditional` that differs 
# from the default, then an error is thrown. While `lik_par` is not 
# required for any computations in this class, it is still required 
# for reference and validation purposes. Similar to  `normalize` and 
# `conditional`, an error is thrown if the user tries to pass a likelihood 
# parameter that differs from the fixed value set when initializing the class.
#
# Since this class does not know anything about the specific form of the 
# likelihood, the user may optionally provide bounds on the log-likelihood 
# that will be enforced by truncating the GP predictive distribution. This is 
# done through the `llik_bounds` attribute, which is simply a vector of length
# two storing lower and upper global bounds on the log-likelihood.
# -----------------------------------------------------------------------------

llikEmulatorGP <- setRefClass(
  Class = "llikEmulatorGP", 
  contains = "llikEmulator"
)

llikEmulatorGP$methods(
  
  initialize = function(llik_lbl, gp_model, default_conditional, 
                        default_normalize, lik_par=NULL, ...) {
    # There are issues related to the use of `$copy()` method regarding missing 
    # arguments. The below line is somewhat of a hack that deals with this case.
    if(any(missing(llik_lbl), missing(gp_model), missing(default_conditional),
           missing(default_normalize))) return(NULL)
    
    # Argument checking.
    assert_that(inherits(gp_model, "gpWrapper"),
                msg="`gp_model` must inherit from `gpWrapper` class.")
    assert_that(gp_model$Y_dim==1, 
                msg="`llikEmulatorGP` only supports single-output GP emulator.")
    assert_that(!is.null(gp_model$X_names) && noNA(gp_model$X_names), 
                msg=paste0("`llikEmulatorGP` requires that `gp_model`",
                           "has `X_names` field set."))

    callSuper(emulator_model=gp_model, llik_label=llik_lbl, lik_par=lik_par, 
              input_names=gp_model$X_names, dim_input=gp_model$X_dim, 
              default_conditional=default_conditional, 
              default_normalize=default_normalize,
              lik_description="Generic log likelihood",
              emulator_description="GP directly emulating the log-likelihood.", 
              llik_pred_dist="Gaussian", exact_llik=FALSE, ...)
  }, 
  
  check_fixed_quantities = function(conditional=NULL, normalize=NULL, 
                                    lik_par_val=NULL) {
    
    # For now, removing this.
    
    # if(!is.null(conditional)) {
    #   assert_that(conditional==default_conditional, 
    #               msg=paste0("`llikEmulatorGP` class requires `conditional` to",
    #                          "agree with `default_conditional`."))
    # }
    # 
    # if(!is.null(normalize)) {
    #   assert_that(normalize==default_normalize, 
    #               msg=paste0("`llikEmulatorGP` class requires `normalize` to",
    #                          "agree with `default_normalize`."))
    # }
    # 
    # if(!is.null(lik_par_val)) {
    #   assert_that(lik_par_val==.self$get_lik_par(lik_par), 
    #               msg=paste0("`llikEmulatorGP` class requires `lik_par_val` to",
    #                          "agree with `lik_par`."))
    # }
  },
  
  .assemble_llik = function(llik, lik_par_val=NULL, 
                            conditional=default_conditional, 
                            normalize=default_normalize, ...) {
    # llik should be N_input x 1. Since the llik is emulated directly, then this 
    # function simply returns the argument `llik` (flattened to a vector) after 
    # performing argument validation. 
    .self$check_fixed_quantities(conditional, normalize, lik_par_val)
    return(drop(llik))
  },
  
  get_design_inputs = function(...) {
    .self$emulator_model$X
  }, 
  
  .get_design_llik = function(lik_par_val=NULL, conditional=default_conditional,
                              normalize=default_normalize, ...) {
    # Returns the response values in the design of the GP, since the response 
    # is the llik in this case. 
    .self$check_fixed_quantities(conditional, normalize, lik_par_val)
    return(drop(emulator_model$Y))
  },
  
  sample_emulator = function(input, em_pred_list=NULL, N_samp=1L, use_cov=FALSE, 
                             include_noise=TRUE, adjustment="rectified", ...) {
    # Default behavior is to sample from a truncated Gaussian to satisfy 
    # the likelihood bounds.

    # The default behavior here is to subtract the shift from the log-likelihood
    # bounds. This is because the default is to treat the bounds as applying
    # to the shifted distribution, so we undo the shift here in order to 
    # apply the bounds to the log-likelihood emulator.
    bounds <- .self$get_llik_bounds(input, ...)
    
    # Supports the same adjustments as the underlying GP object, since the
    # log-likelihood emulator is a GP.
    adjustment <- .self$emulator_model$get_dist_adjustment(adjustment, 
                                                           lower_bound=bounds$lower,
                                                           upper_bound=bounds$upper)
    
    .self$emulator_model$sample(.self$get_input(input), use_cov=use_cov, 
                                include_noise=include_noise, 
                                N_samp=N_samp, 
                                pred_list=em_pred_list, 
                                adjustment=adjustment, 
                                lower_bound=bounds$lower, 
                                upper_bound=bounds$upper, ...)[,,1,drop=FALSE]             
  },
  
  update_emulator = function(input_new, llik_new, update_hyperpar=FALSE, ...) {
    # Condition the GP emulator on the new input pair (input_new, llik_new).
    .self$emulator_model$update(.self$get_input(input_new), matrix(llik_new, ncol=1), 
                                update_hyperpar=update_hyperpar, ...)
  },
  
  .sample = function(input, lik_par_val=NULL, em_pred_list=NULL, N_samp=1L, 
                     use_cov=FALSE, include_noise=TRUE, adjustment="rectified",
                     conditional=default_conditional, 
                     normalize=default_normalize, ...) {
    # Note that the bounds, etc. are dealt with in `sample_emulator` so no
    # need to deal with it here.

    # Directly returns the emulator samples, since these are llik samples. 
    .self$check_fixed_quantities(conditional, normalize, lik_par_val)
    .self$sample_emulator(input, em_pred_list, N_samp, use_cov, 
                          include_noise, adjustment=adjustment, ...)[,,1]
  }, 
  
  .predict = function(input, lik_par_val=NULL, em_pred_list=NULL, return_mean=TRUE, 
                      return_var=TRUE, return_cov=FALSE, return_cross_cov=FALSE, 
                      input_cross=NULL, conditional=default_conditional, 
                      normalize=default_normalize, include_noise=TRUE, 
                      adjustment="rectified", ...) {
    # Log-likelihood emulator mean/var/cov predictions. Since the GP directly  
    # emulates the log-likelihood, then simply return the GP predictions 
    # directly. The only modification is to flatten the mean/var/cov predictions
    # since the gpWrapper predictions can be multi-output, but the predictions
    # here will always be univariate.
    #
    # Note: for now it is assumed that if `em_pred_list` is explicitly passed,
    #       then it has already been adjusted as desired; i.e. `adjustment` is
    #       not applied to the existing predictions. This behavior may change.
    
    .self$check_fixed_quantities(conditional, normalize, lik_par_val)
    
    # Note that by default the shift is subtracted from the bounds. Default
    # behavior can be overwritten via `...` arguments.
    bounds <- .self$get_llik_bounds(input, ...)
    
    # Supports the same adjustments as the underlying GP object, since the
    # log-likelihood emulator is a GP.
    adjustment <- .self$emulator_model$get_dist_adjustment(adjustment, 
                                                           lower_bound=bounds$lower,
                                                           upper_bound=bounds$upper)
    
    if(is.null(em_pred_list)) {
      em_pred_list <- .self$predict_emulator(input, 
                                             return_mean=return_mean, 
                                             return_var=return_var,
                                             return_cov=return_cov, 
                                             return_cross_cov=return_cross_cov,
                                             X_cross=input_cross, 
                                             include_noise=include_noise, 
                                             adjustment=adjustment,
                                             lower_bound=bounds$lower,
                                             upper_bound=bounds$upper, ...)
    }
    
    # Flatten predictions.
    em_pred_list$mean <- drop(em_pred_list$mean)
    em_pred_list$var <- drop(em_pred_list$var)
    em_pred_list$trend <- drop(em_pred_list$trend)
    em_pred_list$cov <- em_pred_list$cov[,,1]
    em_pred_list$cross_cov <- em_pred_list$cross_cov[,,1]
    
    return(em_pred_list)
  }, 
  
  
  .predict_lik = function(input, lik_par_val=NULL, em_pred_list=NULL, 
                          llik_pred_list=NULL, return_mean=TRUE, return_var=TRUE, 
                          return_cov=FALSE, return_cross_cov=FALSE, 
                          input_cross=NULL, conditional=default_conditional, 
                          normalize=default_normalize, include_noise=TRUE, 
                          log_scale=FALSE, adjustment="rectified", ...) {
    # Likelihood emulator mean/var/cov predictions. For this class (under direct 
    # GP emulation of the log-likelihood), the likelihood emulator is a 
    # log-normal process. Thus, the GP log-likelihood predictions can simply be 
    # transformed to obtain log-normal likelihood predictions. Currently 
    # `return_cross_cov` is not supported. Note that both the GP mean and 
    # variance are required to compute the log-normal mean and/or variance. 

    # Note that by default the shift is subtracted from the bounds. Default
    # behavior can be overwritten via `...` arguments.
    bounds <- .self$get_llik_bounds(input, ...)
    
    # Supports the same adjustments as the underlying GP object. In this case,
    # rectified/truncated refer to rectified/truncated log-normal, as opposed
    # to rectified/truncated Gaussian in the GP case.
    adjustment <- .self$emulator_model$get_dist_adjustment(adjustment, 
                                                           lower_bound=bounds$lower,
                                                           upper_bound=bounds$upper)

    if(!is.null(adjustment) && return_cov) {
      stop("Adjustments to `llikEmulatorGP$predict_lik()` lognormal distribution ",
           " are not currently supported for predictive covariances. Truncated and ",
           " rectified lognormal adjustments are supported on a pointwise basis.")
    }
    
    if(return_cross_cov) {
      stop("`return_cross_cov` is not yet supported for `llikEmulatorGP$predict_lik()`.")
    }
    
    if(is.null(llik_pred_list)) {
      llik_pred_list <- .self$predict(input=input, lik_par_val=lik_par_val, 
                                      em_pred_list=em_pred_list, 
                                      return_mean=TRUE, return_var=TRUE, 
                                      return_cov=return_cov,  
                                      return_cross_cov=return_cross_cov, 
                                      input_cross=input_cross, 
                                      conditional=conditional, 
                                      normalize=normalize, 
                                      include_noise=include_noise,
                                      adjustment=NULL, apply_shift=FALSE, ...)
    }
                               
    convert_Gaussian_to_LN(mean_Gaussian=llik_pred_list$mean, 
                           var_Gaussian=llik_pred_list$var, 
                           cov_Gaussian=llik_pred_list$cov, 
                           return_mean=return_mean, return_var=return_var, 
                           return_cov=return_cov, log_scale=log_scale,
                           adjustment=adjustment, lower=bounds$lower,
                           upper=bounds$upper)
  }, 
  
  
  calc_lik_mean_approx = function(em_pred_list=NULL, input=NULL, 
                                  lik_par_val=NULL,
                                  conditional=default_conditional,
                                  normalize=default_normalize, 
                                  include_noise=TRUE, log_scale=TRUE, 
                                  adjustment="rectified", ...) {
    # This method overrides the llikEmulator default for the sole purpose of 
    # allowing for the predictive mean to be truncated/rectified. The plug-in
    # mean truncated/rectified method means that the GP predictive dist 
    # is first truncated rectified, and the mean of this dist is then 
    # exponentiated to obtain the likelihood approximation.

    llik_pred_list <- .self$predict(input, em_pred_list=em_pred_list,
                                    return_mean=TRUE, conditional=conditional,
                                    normalize=normalize, adjustment=adjustment, ...)
    
    if(log_scale) return(llik_pred_list$mean)
    return(exp(llik_pred_list$mean))
  },
  

  calc_lik_marginal_approx = function(em_pred_list=NULL, input=NULL, 
                                      lik_par_val=NULL, llik_pred_list=NULL,
                                      conditional=default_conditional, 
                                      normalize=default_normalize,
                                      include_noise=TRUE, log_scale=TRUE, 
                                      adjustment="rectified", ...) {
    # For log-likelihood distribution l(u) ~ N(m(u), k(u)) the marginal 
    # approximation is E[exp{l(u)}]. The truncated adjustment is 
    # E[exp{l(u)} | b1 <= l(u) <= b2], where (b1, b2) are the bounds provided 
    # by get_llik_bounds(). The rectified adjustment is E[T(u)], where T(u) is 
    # defined to be exp(b1) when l(u) < b1, exp(b2) when l(u) > b2, and 
    # exp(l(u)) when b1 <= l(u) <= b2. Computing the rectified expectation 
    # first requires computing the truncated expectation. Note that these 
    # expectations are not the same as first computing the truncated/rectified
    # expectations and then exponentiating.

    # Note that the bounds for the adjustment are dealt with within `predict_lik`,
    # so no need to deal with them here.
    log_mean <- .self$predict_lik(input, lik_par_val=lik_par_val, 
                                  em_pred_list=em_pred_list,
                                  llik_pred_list=llik_pred_list,
                                  return_mean=TRUE, return_var=FALSE, 
                                  conditional=default_conditional, 
                                  normalize=default_normalize, 
                                  include_noise=include_noise, 
                                  log_scale=TRUE,
                                  adjustment=adjustment, ...)$log_mean
    
    if(log_scale) return(log_mean)
    return(exp(log_mean))
  },
  
  calc_lik_quantile_approx = function(em_pred_list=NULL, input=NULL, alpha=0.9,  
                                      lik_par_val=NULL,
                                      conditional=default_conditional, 
                                      normalize=default_normalize,
                                      include_noise=TRUE, log_scale=TRUE,
                                      adjustment=NULL, ...) {
    # Deterministic likelihood approximation that is given by the alpha
    # quantile of the likelihood surrogate.
    
    .self$calc_lik_quantiles(input=input, p=alpha, lik_par_val=lik_par_val, 
                             em_pred_list=em_pred_list, log_scale=log_scale, 
                             conditional=conditional, normalize=normalize, 
                             lower_tail=TRUE, include_noise=include_noise, 
                             adjustment=adjustment, ...)
  },
  
  .calc_quantiles = function(input, p, lik_par_val=NULL, em_pred_list=NULL,
                             conditional=default_conditional, 
                             normalize=default_normalize, 
                             lower_tail=TRUE, include_noise=TRUE, 
                             adjustment=NULL, ...) {
    # The log-likelihood emulator distribution is Gaussian, so computes
    # Gaussian quantiles. Also supports rectified and truncated Gaussian quantiles.
    
    bounds <- .self$get_llik_bounds(input, ...)
    adjustment <- .self$emulator_model$get_dist_adjustment(adjustment, 
                                                           lower_bound=bounds$lower,
                                                           upper_bound=bounds$upper)
    if(!lower_tail) p <- 1 - p
    
    # First compute mean/variance of GP llik predictive distribution.
    llik_pred_list <- .self$predict(input, lik_par_val=lik_par_val, 
                                    return_mean=TRUE, return_var=TRUE,
                                    apply_shift=FALSE, adjustment=NULL,
                                    em_pred_list=em_pred_list)

    # Compute quantiles.
    if(is.null(adjustment)) { # Gaussian quantiles.
      q <- qnorm(p, llik_pred_list$mean, sqrt(llik_pred_list$var))
    } else if(isTRUE(adjustment == "truncated")) {
      q <- truncnorm::qtruncnorm(p, a=bounds$lower, b=bounds$upper, 
                                 mean=llik_pred_list$mean, sd=sqrt(llik_pred_list$var))
    } else if(isTRUE(adjustment == "rectified")) {
      q <- rect_norm_quantile(p, mean=llik_pred_list$mean, sd=sqrt(llik_pred_list$var), 
                              lower=bounds$lower, upper=bounds$upper)
    } else {
      stop("Quantile calculation unsupported with `adjustment` ", adjustment)
    }
    
    return(q)
  },
  
  .calc_lik_quantiles = function(input, p, lik_par_val=NULL, em_pred_list=NULL,
                                 log_scale=TRUE, conditional=default_conditional, 
                                 normalize=default_normalize, 
                                 lower_tail=TRUE, include_noise=TRUE,
                                 adjustment=NULL, ...) {
    # The log-likelihood emulator distribution is log-normal, so computes
    # log-normal quantiles. Also supports truncated log-normal quantiles,
    # and rectified (clipped) log-normal quantiles (the latter is estimated
    # via Monte Carlo). Note that quantiles will be shifted by the parent
    # function `calc_lik_quantiles()` so ensure the shifting is not done here.
    
    bounds <- .self$get_llik_bounds(input, ...)
    adjustment <- .self$emulator_model$get_dist_adjustment(adjustment, 
                                                           lower_bound=bounds$lower,
                                                           upper_bound=bounds$upper)
    if(!lower_tail) p <- 1 - p
    
    # First compute mean/variance of GP llik predictive distribution.
    llik_pred_list <- .self$predict(input, lik_par_val=lik_par_val, 
                                    return_mean=TRUE, return_var=TRUE,
                                    apply_shift=FALSE, adjustment=NULL,
                                    em_pred_list=em_pred_list, ...)
    
    # Compute log-normal or truncated log-normal quantile.
    if(is.null(adjustment)) {
      q <- qnorm(p)
      ln_quantiles <- llik_pred_list$mean + q * sqrt(llik_pred_list$var)
      if(!log_scale) ln_quantiles <- exp(ln_quantiles)
    } else if(adjustment == "truncated") {
      # Note that the bounds in EnvStats::qlnormTrunc are defined on the 
      # exponentiated scale.
      ln_quantiles <- EnvStats::qlnormTrunc(p, meanlog=llik_pred_list$mean,
                                            sdlog=sqrt(llik_pred_list$var),
                                            min=exp(bounds$lower),
                                            max=exp(bounds$upper))
      if(log_scale) ln_quantiles <- log(ln_quantiles)
    } else if(adjustment == "rectified") {
      # NOTE: this is computable but just using Monte Carlo estimate for now
      # for simplicity.
      samp <- .self$.sample(.self$get_input(input), N_samp=1e5, use_cov=FALSE,
                            adjustment="rectified", include_noise=include_noise,
                            conditional=conditional, normalize=normalize)
      ln_quantiles <- apply(exp(samp), 1, function(x) quantile(x, p))
      if(log_scale) ln_quantiles <- log(ln_quantiles)
    } else { 
      stop("Adjustment `", adjustment, 
           "` unsupported by llikEmulatorGP$calc_lik_quantile_approx().")
    }
    
    return(ln_quantiles)
  },
  
  calc_expected_lik_cond_var = function(input_eval, input_cond, include_noise=TRUE, 
                                        log_scale=TRUE, plugin=FALSE, ...) {
    # Since the emulator model directly emulates the log-likelihood in this 
    # class, the  options `plugin = FALSE` and `plugin = TRUE` are equivalent 
    # in this case. In either case, the analogous method for the GP 
    # `emulator_model` is dispatched. 
    
    .self$emulator_model$calc_expected_exp_cond_var(input_eval, input_cond, 
                                                    include_noise=include_noise, 
                                                    log_scale=log_scale, ...)
  }
)


# -----------------------------------------------------------------------------
# llikEmulatorExactGaussDiag class.
#    Implements an exact Gaussian likelihood with mean given by the output 
#    of a forward model, with a diagonal covariance matrix. Identical to 
#    llikEmulatorExactGauss except that `Sig` is assumed to be diagonal, so the 
#    `lik_par` for this class is defined to be the vector corresponding to the 
#    diagonal of this covariance matrix.`lik_par` may also be provided as a 
#    single value, which is interpreted as a homoskedastic variance. Concretely, 
#    implements a likelihood of the form, 
#        p(y1,...,yP|u,v1,...,vP) = prod_{p=1}^{P} prod_{n=1}^{N} N(yp_n|G(u)_p, vp),
#    where:
#    P = `N_output` 
#    N = `N_obs`
#    y1,...,yP are each vectors of length N.
#    v1,...,vP are variance parameters corresponding to each output, respectively. 
#
# In words, this likelihood assumes the forward model `G` has P outputs. A set of 
# observations yp is associated with each output p, and these observations are 
# assumed to be iid draws from N(G(u)_p, vp). The argument `y_obs` to the class 
# constructor is a (N,P) matrix storing each set of observations as columns.
# `fwd_model_vectorized()` is a vectorized version of `fwd_model()`, meaning 
# that it can accept a matrix of inputs (one per row). If not explicitly 
# passed, the method `vectorize_fwd_model()` will be utilized.
# 
# TODO: generalize to allow ragged arrays, where the yp need not be the same length.
# -----------------------------------------------------------------------------

llikEmulatorExactGaussDiag <- setRefClass(
  Class = "llikEmulatorExactGaussDiag", 
  contains = "llikEmulator",
  fields = list(fwd_model="ANY", fwd_model_vectorized="ANY", y="matrix", 
                N_output="integer", N_obs="integer")
)

llikEmulatorExactGaussDiag$methods(
  
  initialize = function(llik_lbl, y_obs, dim_par, fwd_model=NULL, 
                        fwd_model_vectorized=NULL, sig2=NULL, 
                        default_conditional=FALSE, default_normalize=FALSE, 
                        par_names=NULL, ...) {

    # Forward model must be provided, either vectorized or single-input version. 
    assert_that(is.null(fwd_model) || is.function(fwd_model))
    assert_that(is.null(fwd_model_vectorized) || is.function(fwd_model_vectorized))
    assert_that(is.function(fwd_model) || is.function(fwd_model_vectorized))
    assert_that(is.matrix(y_obs))
    initFields(fwd_model=fwd_model, fwd_model_vectorized=fwd_model_vectorized, 
               N_output=ncol(y_obs), N_obs=nrow(y_obs), y=y_obs)
    
    # Set parameter names. 
    if(is.null(par_names)) par_names <- paste0("input", 1:dim_par)
    assert_that(length(par_names) == dim_par)
    
    # Variance parameters for Gaussian likelihood. 
    if(!is.null(sig2)) {
      assert_that(is.numeric(sig2) && ((length(sig2)==N_output) || 
                                         (length(sig2)==1)) && all(sig2>0),
                  msg=paste0("`sig2` must be either vector of length `N_output`",
                             "or 1 and only contain positive numbers."))
    }
    
    callSuper(emulator_model=NULL, llik_label=llik_lbl, lik_par=sig2, 
              dim_input=dim_par, default_conditional=default_conditional, 
              input_names=par_names, default_normalize=default_normalize, 
              lik_description="Exact Gaussian likelihood, diagonal covariance structure.",
              emulator_description="No emulation.", exact_llik=TRUE, ...)
  },
  
  run_fwd_model = function(input, ...) {
    # `input` is an (M,D) matrix with input parameter values stacked in the rows. 
    # Returns (M, N_output) dimensional output. 
    
    input <- .self$get_input(input)
    if(is.null(.self$fwd_model_vectorized)) return(vectorize_fwd_model(input, ...))
    return(.self$fwd_model_vectorized(input, ...))
  },
  
  vectorize_fwd_model = function(input, ...) {
    model_output <- matrix(nrow=nrow(input), ncol=N_output)
    for(i in 1:nrow(input)) model_output[i,] <- .self$fwd_model(input[i,], ...)
    return(model_output)
  },
  
  get_lik_par = function(lik_par_val=NULL, ...) {
    # If the passed value is a scalar but the output dimension is larger than
    # one, then the value is interpreted as a homoscedastic variance parameter.
    # Otherwise, `lik_par_val` should be a vector representing the diagonal
    # of the observation covariance matrix.
    
    if(is.null(lik_par_val)) lik_par_val <- .self$lik_par
    
    assert_that(!is.null(lik_par_val), 
                msg="Missing likelihood parameter.")
    
    if((length(lik_par_val)==1) && (.self$N_output > 1)) {
      return(rep(lik_par_val, .self$N_output))
    }
    
    assert_that(length(lik_par_val)==.self$N_output, 
                msg="`lik_par_val` length not equal to 1 or `N_output`.")
    return(lik_par_val)
  },
  
  compute_log_det = function(lik_par_val=NULL) {
    # Computes the log of the normalizing constant of the Gaussian density. 
    # Specifically, for covariance matrix C = diag{sig_1*I, ..., sig_P*I}, 
    # where P = N_output and I is the N := N_obs dimensional identity,
    # computes:
    # 0.5 * log{det(2*pi*C)} = 0.5*N*P*log(2pi) + 0.5*N*sum_{p=1}^{P} log(sig2_p).
    # Note that `lik_par_val` is the diagonal of the covariance matrix C.
    # Note also that this log determinant term is also an upper bound for the 
    # Gaussian log-likelihood; i.e., log N(y|G(u),C) <= 0.5 * log{det(2*pi*C)}.
    
    sig2_val <- .self$get_lik_par(lik_par_val)
    -0.5*.self$N_obs * .self$N_output * log(2*pi) - 0.5*.self$N_obs * sum(log(sig2_val))
  },
  
  .assemble_llik = function(input, lik_par_val=NULL, 
                            conditional=default_conditional, 
                            normalize=default_normalize, ...) {
    # `input` should have dimension (N_input, D). Returns vector of length `N_input.` 
    # NOTE: currently ignoring conditional/normalize here, as I decide what 
    # to do with these arguments. Just always normalizing.

    # Fetch the variance parameters. 
    sig2_val <- lik_par_val
    
    # Run forward model. 
    fwd_model_vals <- .self$run_fwd_model(input, ...)
    
    # Construct log likelihood.
    llik <- vector(mode="numeric", length=nrow(input))
    for(i in seq_along(llik)) {
      llik[i] <- -0.5 * sum(mult_vec_with_mat_rows(1/sig2_val, 
                            add_vec_to_mat_rows(-fwd_model_vals[i,], .self$y)^2))
    }
    
    llik + .self$compute_log_det(sig2_val)
  }, 
  
  get_llik_bounds = function(lik_par_val=NULL, ...) {
    upper <- .self$compute_log_det(lik_par_val)
    c(-Inf,upper)
  },
  
  predict_emulator = function(input, em_pred_list=NULL, ...) {
    # No emulator to predict with, just returns the inputs as the mean 
    # prediction.

    list(mean=.self$get_input(input))
  },
  
  sample_emulator = function(input, em_pred_list=NULL, N_samp=1, ...) {
    # No emulator to sample from, simply return input. The argument 
    # `em_pred_list` is only present for consistency with other llikEmulator
    # classes. 
    
    .self$get_input(input)
  },
  
  update_emulator = function(...) {
    # There is no emulator, so this does nothing. 
  },
  
  .sample = function(input, lik_par_val=NULL, em_pred_list=NULL, N_samp=1L, 
                     conditional=default_conditional, 
                     normalize=default_normalize, ...) {
    
    # Compute unnormalized or normalized log-likelihood (exact, deterministic 
    # calculation - no sampling is actually performed). For consistency with 
    # other classes, duplicates the exact likelihood calculation when `N_samp`>1.

    matrix(.self$assemble_llik(input, lik_par_val, conditional, normalize), 
           nrow=nrow(input), ncol=N_samp)
  }, 
  
  .predict = function(input, lik_par_val=NULL, return_mean=TRUE, return_var=TRUE, 
                      return_cov=FALSE, return_cross_cov=FALSE, 
                      conditional=default_conditional, 
                      normalize=default_normalize, ...) {
    # Also just evaluates the exact log-likelihood. For consistency with other 
    # llikEmulator classes, these exact evaluations are stored as the "mean" 
    # element of a list. If requested, predictive variances and covariances 
    # are set to 0.

    return_list <- list()
    if(return_mean) {
      return_list$mean <- .self$assemble_llik(input, lik_par_val, 
                                              conditional, normalize)
    }
    if(return_var) return_list$var <- rep(0, nrow(input))
    if(return_cov) return_list$cov <- rep(0, nrow(input))
    
    return(return_list)
  },
  
  .predict_lik = function(input, lik_par_val=NULL, return_mean=TRUE, return_var=TRUE, 
                          return_cov=FALSE, return_cross_cov=FALSE, 
                          conditional=default_conditional, log_scale=TRUE,
                          normalize=default_normalize, ...) {
    # Same as `.predict()` but optionally exponentiates the values.
    
    l <- .self$predict(input, lik_par_val=lik_par_val, return_mean=return_mean,
                       return_var=return_var, ...)
    if(!log_scale) l$mean <- exp(l$mean)
    
    return(l)
  }

)

# -----------------------------------------------------------------------------
# llikEmulatorGPFwdGauss class.
#
# Implements a surrogate likelihood for the Gaussian inverse problem 
# y|u ~ N(G(u), Sig) where the forward model G is replaced by a 
# Gaussian process (GP) emulator G ~ GP(m,k). The forward model may have  
# multiple outputs, in which case the GP is assumed to be a multi-output GP 
# with the  same number of outputs as the observation space dimension.
# The attribute `N_output` is the number of outputs (i.e. the 
# dimension of the vector returned by G). This class does not assume any 
# special structure in the covariance `Sig`, so it is primarily intended for 
# dense observation covariances, or for multi-output GPs with dense output 
# covariances. Note that the expectation of the random likelihood (i.e., 
# the marginal approximation) is given by 
# E{N(y|G(u), Sig)} = N(y|m(u), Sig + k(u)). So the covariance matrix of the 
# marginal approximation can be dense when either Sig or k(u) are dense.
# The class llikEmulatorGPFwdGaussDiag treats the special case where both 
# Sig and k(u) are diagonal.
#
# At present, this class only supports emulators of the gpWrapperSum class.
# This is due to the fact that it is assuming the emulator is a multi-output
# GP, and the "include_output_cov" argument is only currently available in this
# class. This can be generalized to other gpWrapper classes by including 
# `include_output_cov`, which would just be a diagonal matrix for independent
# GPs.
#
# Likelihood bounds:
# Note that, unlike with `llikEmulatorGP`, there is no need to explicitly 
# constrain the log-likelihood emulator. Since the forward model emulator is
# being plugged into a Gaussian likelihood, then the Gaussian likelihood
# bounds will automatically be respected. Note that it is still possible
# to pass truncated/rectified adjustments via `...` to methods like 
# `sample_emulator()`. This is typically not recommended, since certain methods
# (e.g., `predict_lik()`) assume that the forward model emulator has a Gaussian
# predictive distribution.
# -----------------------------------------------------------------------------

llikEmulatorGPFwdGauss <- setRefClass(
  Class = "llikEmulatorGPFwdGauss", 
  contains = "llikEmulator",
  fields = list(y="ANY", N_output="integer", N_obs="integer", L_Cov="matrix")
                
)

llikEmulatorGPFwdGauss$methods(
  
  initialize = function(llik_lbl, gp_model, y_obs, Cov=NULL, 
                        default_conditional=FALSE, default_normalize=FALSE, 
                        par_names=NULL, ...) {
    # `Cov` can be a matrix of shape (N_output, N_output), a vector of length 
    # N_output (interpreted as diagonal of covariance), or a scalar 
    # (interpreted as homoscedastic variance parameter).
    
    assert_that(inherits(gp_model, "gpWrapper"), 
                msg="`gp_model` must inherit from `gpWrapper` class.")
    # assert_that(isTRUE(class(gp_model)=="gpWrapperSum"),
    #             msg="At present, only gpWrapperSum objects are allowed for the ",
    #                 "emulator model in llikEmulatorGPFwdGauss.")
    assert_that(is.numeric(y_obs) || is.matrix(y_obs))
    y_obs <- drop(y_obs)
    assert_that(length(y_obs) == gp_model$Y_dim)
    initFields(y=y_obs, N_output=length(y_obs))
    
    # Set parameter names. 
    dim_par <- gp_model$X_dim
    if(is.null(par_names)) par_names <- paste0("input", 1:dim_par)
    assert_that(length(par_names) == dim_par)
    
    # Covariance matrix of Gaussian likelihood. 
    if(!is.null(Cov)) {
      if(length(drop(Cov)) %in% c(1L, .self$N_output)) {
        Cov <- diag(Cov, nrow=.self$N_output, ncol=.self$N_output)
      }
      
      assert_that(is.matrix(Cov) && (nrow(Cov)==N_output) && (ncol(Cov)==N_output),
                  msg="`Cov` must be a positive definite matrix with dim `N_output` x `N_output`")
      initFields(L_Cov=t(chol(Cov))) 
    }
    
    # Upper bound for Gaussian log-likelihood.
    llik_bound_func <- function(input, lik_par_val=NULL, ...) {
      upper <- .self$compute_log_det(lik_par_val)
      list(lower=-Inf, upper=upper)
    }
    
    callSuper(emulator_model=gp_model, llik_label=llik_lbl, lik_par=Cov, 
              dim_input=dim_par, default_conditional=default_conditional, 
              input_names=par_names, default_normalize=default_normalize, 
              lik_description="Gaussian likelihood",
              emulator_description="Forward model GP emulator", 
              exact_llik=FALSE, llik_bounds=llik_bound_func, ...)
  },
  
  get_lik_par = function(lik_par_val=NULL, return_chol=FALSE, ...) {
    
    # If covariance matrix is passed, convert to Cholesky factor is requested.
    if(!is.null(lik_par_val)) {
      if(return_chol) lik_par_val <- t(chol(lik_par_val))
    } else {
      if(return_chol) lik_par_val <- .self$L_Cov
      else lik_par_val <- .self$lik_par
    }

    assert_that(!is.null(lik_par_val), 
                msg="Likelihood parameter (noise covariance) missing.")
  
    return(lik_par_val)
  },
  
  get_design_inputs = function(...) {
    .self$emulator_model$X
  },
  
  .get_design_llik = function(lik_par_val=NULL, conditional=default_conditional, 
                              normalize=default_normalize, ...) {
    .self$assemble_llik(.self$emulator_model$Y, lik_par_val=lik_par_val, 
                        conditional=conditional, normalize=normalize, ...)
  },
  
  compute_log_det = function(lik_par_val=NULL, lik_par_chol=NULL) {
    # Computes the log of the normalizing constant of the Gaussian density. 
    # Specifically, for covariance matrix C, computes -0.5 * log{det(2*pi*C)}.
    # Note that `lik_par_val` is the covariance matrix C. Optionally can pass 
    # `lik_par_chol`, the lower Cholesky factor of C. Note also that this
    # log-determinant term is an upper bound on the log-likelihood; i.e.,
    # log N(y|G(u),C) <= -0.5 * log{det(2*pi*C)}.
    
    if(is.null(lik_par_chol)) {
      lik_par_chol <- .self$get_lik_par(lik_par_val, return_chol=TRUE)
    }
    
    -0.5 * .self$N_output * log(2*pi) - sum(log(diag(lik_par_chol)))
  },
  
  .assemble_llik = function(fwd_model_vals, lik_par_val=NULL, lik_par_chol=NULL,
                            conditional=default_conditional, 
                            normalize=default_normalize, ...) {
    # `fwd_model_vals` should be of dimension (N_inputs, N_output). Returns 
    # vector of length N_inputs. `lik_par_val` is the covariance matrix of 
    # the Gaussian density. The density is computed using the Cholesky factor
    # of the covariance. If the Cholesky factor has already been computed, 
    # it can be passed via `lik_par_chol`.
    
    # Fetch the lower triangular Cholesky factor of the covariance matrix.
    if(is.null(lik_par_chol)) {
      L <- .self$get_lik_par(lik_par_val, return_chol=TRUE)
    } else {
      L <- lik_par_chol
    }
    
    # Construct log likelihood.
    llik <- -0.5 * colSums(solve(L, add_vec_to_mat_cols(.self$y, -t(fwd_model_vals)))^2)
    llik + .self$compute_log_det(lik_par_chol=L)
  },
  
  .sample = function(input, lik_par_val=NULL, em_pred_list=NULL, N_samp=1L, 
                     use_cov=TRUE, conditional=default_conditional, 
                     normalize=default_normalize, ...) {
    # Sample the log-likelihood emulator at specified inputs. `input` is M x D 
    # (M input vectors). Returns array of dimension (M, N_samp). `use_cov`
    # refers to covariance across inputs; inclusion of output covariance 
    # falls back on default behavior of .self$sample_emulator(), which should 
    # typically be to always include output covariance. See, e.g., 
    # gpWrapperSum$sample() for an example.
    
    fwd_model_samp <- .self$sample_emulator(input, N_samp=N_samp, 
                                            use_cov=use_cov, 
                                            pred_list=em_pred_list, ...)
    llik_samp <- matrix(nrow=nrow(input), ncol=N_samp)
    
    for(i in 1:N_samp) {
      llik_samp[,i] <- .self$assemble_llik(matrix(fwd_model_samp[,i,], 
                                                  nrow=nrow(input),
                                                  ncol=.self$N_output), 
                                           lik_par_val, conditional, 
                                           normalize, ...)
    }
    
    return(llik_samp)
  },
  
  sample_emulator = function(input, em_pred_list=NULL, N_samp=1L, 
                             use_cov=TRUE, ...) {
    # Sample the forward model emulator at specified inputs. `input` is M x D 
    # (M input vectors). Returns array of dimension (M, N_samp, N_output). 
    
    .self$emulator_model$sample(.self$get_input(input), use_cov=use_cov, 
                                include_noise=include_noise, 
                                N_samp=N_samp, pred_list=em_pred_list, ...)
  }, 
  
  .predict = function(input, lik_par_val=NULL, em_pred_list=NULL, return_mean=TRUE,  
                      return_var=TRUE, return_cov=FALSE, return_cross_cov=FALSE, 
                      input_cross=NULL, conditional=default_conditional,  
                      normalize=default_normalize, include_noise=TRUE, ...) {
    # Let M := nrow(input) (number of input points). Currently this method 
    # only supports computation of pointwise means and variances. Note that the 
    # both of these quantities require the underlying GP means and GP pointwise 
    # covariances in the output space.
    #
    # Return formatting:
    # "mean": numeric vector of length M.
    # "var": numeric vector of length M.
    # Covariance return options "cov" and "cross_cov" not yet supported.
    
    if(return_cross_cov || return_cov) {
      message("`return_cross_cov` and `return_cov` not yet supported for ",
              "`llikEmulatorGPFwdGauss$predict()`. Will not return these quantities.")
      return_cross_cov <- FALSE
      return_cov <- FALSE
    }
    
    L <- .self$get_lik_par(lik_par_val, return_chol=TRUE)
    llik_pred_list <- list()
    
    # Forward model emulator predictions. 
    if(is.null(em_pred_list)) {
      em_pred_list <- .self$predict_emulator(input, return_mean=TRUE, 
                                             return_var=FALSE, return_cov=FALSE, 
                                             return_output_cov=TRUE,
                                             return_cross_cov=return_cross_cov, 
                                             input_cross=input_cross, 
                                             include_noise=include_noise, ...)
    } else {
      assert_that(!is.null(em_pred_list$mean) && !is.null(em_pred_list$output_cov), 
                  msg=paste0("Log-likelihood predictive quantities require ", 
                             "forward model emulator mean and output covariance."))
    }
    
    # Log-likelihood mean can be shown to be equal to the plug-in mean plus 
    # a variance inflation term depending on the likelihood covariance and 
    # GP covariance. The inflation term is added below.
    if(return_mean) {
      llik_pred_list$mean <- vector(mode="numeric", length=nrow(input))
      llik_plug_in_mean <- .self$assemble_llik(em_pred_list$mean, 
                                               lik_par_val=lik_par_val, 
                                               conditional=conditional, 
                                               normalize=normalize, ...)
    }
    
    if(return_var) {
      llik_pred_list$var <- vector(mode="numeric", length=nrow(input))
    }
    
    if(return_mean || return_var) {
      for(i in 1:nrow(input)) {
        sig_inv_K <- backsolve(t(L), forwardsolve(L, em_pred_list$output_cov[i,,]))
          
        # Compute mean.
        if(return_mean) {
          var_inflation <- sum(diag(sig_inv_K))
          llik_pred_list$mean[i] <- llik_plug_in_mean[i] + var_inflation
        }
          
        # Compute variance.
        if(return_var) {
          resid_mean <- .self$y - em_pred_list$mean[i,]
          sig_inv_K_sig_inv <- t(backsolve(t(L), forwardsolve(L, t(sig_inv_K))))
          trace_term <- 0.5 * sum(sig_inv_K * t(sig_inv_K))
          quad_form_term <- 2 * sum(resid_mean * (sig_inv_K_sig_inv %*% resid_mean))
          llik_pred_list$var[i] <- trace_term + quad_form_term
        }
      }
    }
   
    return(llik_pred_list)
  }, 
  
  .predict_lik = function(input, lik_par_val=NULL, em_pred_list=NULL, 
                          return_mean=TRUE, return_var=TRUE, return_cov=FALSE, 
                          return_cross_cov=FALSE, input_cross=NULL, 
                          conditional=default_conditional,  
                          normalize=default_normalize, log_scale=FALSE, 
                          include_noise=TRUE, ...) {
    
    if(return_cross_cov || return_cov) {
      message("`return_cross_cov` and `return_cov` not yet supported for ",
              "`llikEmulatorGPFwdGauss$predict_lik()`. Will not return these quantities.")
      return_cross_cov <- FALSE
      return_cov <- FALSE
    }
    
    lik_pred_list <- list()

    # Forward model emulator predictions. 
    if(is.null(em_pred_list)) {
      em_pred_list <- .self$predict_emulator(input, 
                                             return_mean=TRUE, 
                                             return_var=TRUE, 
                                             return_cov=return_cov, 
                                             return_cross_cov=return_cross_cov, 
                                             input_cross=input_cross, 
                                             include_output_cov=TRUE,
                                             include_noise=include_noise, ...)
    } else {
      assert_that(!is.null(em_pred_list$mean), 
                  msg="Likelihood predictive quantities require emulator mean")
      assert_that(!(return_var && is.null(em_pred_list$output_cov)),
                  msg="Likelihood variance requires emulator output covariance.")
    }
    
    n_input <- nrow(em_pred_list$mean)
    Sig <- .self$get_lik_par(lik_par_val, return_chol=FALSE)
    
    if(return_mean) {
      lik_pred_list$log_mean <- vector(mode="numeric", length=n_input)
    }
    
    if(return_var) {
      lik_pred_list$log_var <- vector(mode="numeric", length=n_input)
    }
    
    for(i in 1:n_input) {
      K <- em_pred_list$output_cov[i,,]
      L_Sig_plus_K <- t(chol(Sig+K))
      
      if(return_mean) {
        lik_pred_list$log_mean[i] <- .self$assemble_llik(em_pred_list$mean[i,,drop=FALSE], 
                                                         lik_par_chol=L_Sig_plus_K,
                                                         conditional=conditional,
                                                         normalize=normalize)
      }
      
      if(return_var) {
        log_num1 <- .self$assemble_llik(em_pred_list$mean[i,,drop=FALSE],
                                        lik_par_val=0.5*Sig + K,
                                        conditional=conditional,
                                        normalize=normalize)
        log_num2 <- .self$assemble_llik(em_pred_list$mean[i,,drop=FALSE],
                                        lik_par_chol=L_Sig_plus_K/sqrt(2),
                                        conditional=conditional,
                                        normalize=normalize)
        log_denom1 <- .self$compute_log_det(lik_par_val=lik_par_val)
        log_denom2 <- .self$compute_log_det(lik_par_chol=L_Sig_plus_K)
        lik_pred_list$log_var[i] <- log_diff_exp(log_num1 - log_denom1, 
                                                 log_num2 - log_denom2)
      }
    }
    
    if(!log_scale) {
      if(return_mean) lik_pred_list$mean <- exp(lik_pred_list$log_mean)
      if(return_var) lik_pred_list$var <- exp(lik_pred_list$log_var)
    }

    return(lik_pred_list)
  }
)


# -----------------------------------------------------------------------------
# llikEmulatorGPFwdGaussDiag class.
#
# A special case of `llikEmulatorGPFwdGauss` where the covariance matrix is 
# assumed to be diagonal. See llikEmulatorGPFwdGauss` for details. Assumes 
# likelihood of the form: 
#
# p(y|u) = prod_{i=1}^{N_obs} N(y_i | G(u), diag{sig2_1, ..., sig2_P})
# 
# where P = N_output and G(u) is P-dimensional.
#
# The `lik_par` here is defined to be the vector of variances `sig2`. If `sig2` 
# is passed as a scalar value (but `N_output` is greater than 1), then the 
# same variance is used for all outputs.
# -----------------------------------------------------------------------------

llikEmulatorGPFwdGaussDiag <- setRefClass(
  Class = "llikEmulatorGPFwdGaussDiag", 
  contains = "llikEmulator",
  fields = list(y="matrix", N_output="integer", N_obs="integer")
  
)

llikEmulatorGPFwdGaussDiag$methods(
  
  initialize = function(llik_lbl, gp_model, y_obs, sig2=NULL, 
                        default_conditional=FALSE, 
                        default_normalize=FALSE, ...) {
    
    assert_that(inherits(gp_model, "gpWrapper"), 
                msg="`gp_model` must inherit from `gpWrapper` class.")
    assert_that(is.matrix(y_obs))
    assert_that(ncol(y_obs) == gp_model$Y_dim)
    initFields(y=y_obs, N_output=ncol(y_obs), N_obs=nrow(y_obs))
    
    # Set parameter names. 
    dim_par <- gp_model$X_dim
    par_names <- gp_model$X_names

    # Variance parameters for Gaussian likelihood. 
    if(!is.null(sig2)) {
      assert_that(is.numeric(sig2) && ((length(sig2)==N_output) ||
                                       (length(sig2)==1)) && all(sig2>0),
                  msg=paste0("`sig2` must be either vector of length `N_output`",
                             "or 1 and only contain positive numbers."))
    }
    
    # Upper bound for Gaussian log-likelihood.
    llik_bound_func <- function(input, lik_par_val=NULL, ...) {
      upper <- .self$compute_log_det(lik_par_val)
      list(lower=-Inf, upper=upper)
    }
    
    callSuper(emulator_model=gp_model, llik_label=llik_lbl, lik_par=sig2, 
              dim_input=dim_par, default_conditional=default_conditional, 
              input_names=par_names, default_normalize=default_normalize, 
              lik_description="Gaussian likelihood, diagonal covariance.",
              emulator_description="Forward model GP emulator", 
              exact_llik=FALSE, llik_bounds=llik_bound_func, ...)
  },
  
  get_lik_par = function(lik_par_val=NULL, ...) {
    # The likelihood parameter is the diagonal of the observation covariance.
    # If a scalar value is provided, then it is interpreted as a 
    # homoscedastic variance that will be replicated to fill out the diagonal.
    
    if(is.null(lik_par_val)) lik_par_val <- .self$lik_par
    
    assert_that(!is.null(lik_par_val), 
                msg="Likelihood parameter missing.")
    
    if((length(lik_par_val)==1) && (.self$N_output > 1)) {
      lik_par_val <- rep(lik_par_val, .self$N_output)
    }
    
    return(lik_par_val)
  },
  
  compute_log_det = function(lik_par_val=NULL) {
    # Computes the log of the normalizing constant of the Gaussian density. 
    # Specifically, for covariance matrix C = diag{sig_1*I, ..., sig_P*I}, 
    # where P = N_output and I is the N := N_obs dimensional identity,
    # computes:
    # 0.5 * log{det(2*pi*C)} = 0.5*N*P*log(2pi) + 0.5*N*sum_{p=1}^{P} log(sig2_p).
    # Note that `lik_par_val` is the diagonal of the covariance matrix C.
    # Note also that this log determinant term is also an upper bound for the 
    # Gaussian log-likelihood; i.e., log N(y|G(u),C) <= 0.5 * log{det(2*pi*C)}.
    
    sig2_val <- .self$get_lik_par(lik_par_val)
    -0.5*.self$N_obs * .self$N_output * log(2*pi) - 0.5*.self$N_obs * sum(log(sig2_val))
  },

  .assemble_llik = function(fwd_model_vals, lik_par_val=NULL, 
                            conditional=default_conditional, 
                            normalize=default_normalize, 
                            var_inflation_vals=NULL, ...) {
    # `fwd_model_vals` should be of dimension `M` x `N_output`, where `M` 
    # corresponds to the number of instances of forward model output; e.g., 
    # might be the forward model evaluated at `M` different inputs, or `M` 
    # samples from the forward model at the same input, or some combination. 
    # Returns numeric vector of length `M`. Note that `sample_emulator` returns 
    # an array of dimension `M` x `N_samp` x `N_output`. Therefore, this 
    # function also allows `fwd_model_vals` to have three dimensions, as long 
    # as the second dimension equals 1 (as is the case with a default call to 
    # `sample_emulator`). `var_inflation_vals` adds values to `sig2`
    # ("variance inflation"), and is used in computing the marginal likelihood 
    # approximation. If shape (`M`, `N_output`) then each row of 
    # `var_inflation_vals` is added to to `sig2` for the corresponding row of 
    # `fwd_model_vals`. If a numeric vector of length `N_output`, then the same 
    # variance inflation is applied to the likelihood evaluations at all `M` 
    # values of the forward model output. 

    # If `fwd_model_vals` has three dimensions, with second dimension equal to 
    # 1, then collapse  to a matrix along this dimension. 
    if(length(dim(fwd_model_vals)) == 3L) {
      if(dim(fwd_model_vals)[2] != 1L) {
        stop("`fwd_model_vals` does not have correct dimensions.")
      }
      fwd_model_vals <- matrix(fwd_model_vals[,1,], nrow=dim(fwd_model_vals)[1], 
                               ncol=dim(fwd_model_vals)[3])
    }
    
    assert_that(length(dim(fwd_model_vals)) == 2L)
    assert_that(ncol(fwd_model_vals) == .self$N_output)
    
    # Validate variance inflation argument. 
    inflate_var <- FALSE
    if(!is.null(var_inflation_vals)) {
      inflate_var <- TRUE
      if(is.null(dim(var_inflation_vals))) {
        assert_that(length(var_inflation_vals) == N_output)
        var_inflation_vals <- matrix(var_inflation_vals, 
                                     nrow=nrow(fwd_model_vals), 
                                     ncol=N_output, byrow=TRUE)
      }
      assert_that(all(dim(var_inflation_vals) == dim(fwd_model_vals)))
    }
    
    # Fetch the variance parameters. 
    sig2_val <- .self$get_lik_par(lik_par_val)
    
    # Construct log likelihood.
    llik <- vector(mode="numeric", length=nrow(fwd_model_vals))
    for(i in seq_along(llik)) {
      sig2_val_i <- ifelse(inflate_var, sig2_val + var_inflation_vals[i,], sig2_val)
      llik[i] <- -0.5 * sum(mult_vec_with_mat_rows(1/sig2_val_i,
                                               add_vec_to_mat_rows(-fwd_model_vals[i,], .self$y)^2))
      if(normalize || !conditional) llik[i] <- llik[i] - 0.5 * .self$N_obs * sum(log(sig2_val_i))
    }
    if(normalize) llik <- llik - 0.5*.self$N_obs * .self$N_output * log(2*pi)

    # Old approach.
    # llik <- rep(0, nrow(fwd_model_vals))
    # stdev <- sqrt(sig2_val)
    # for(i in seq_along(llik)) {
    #   for(j in seq_len(.self$N_obs)) {
    #     llik[i] <- llik[i] + sum(dnorm(.self$y[j,], fwd_model_vals[i,], stdev, log=TRUE))
    #   }
    # }
    
    return(llik)
  }, 
  
  get_design_inputs = function(...) {
    emulator_model$X
  },
  
  .get_design_llik = function(lik_par_val=NULL, conditional=default_conditional, 
                              normalize=default_normalize, ...) {
    .self$assemble_llik(.self$emulator_model$Y, lik_par_val=lik_par_val, 
                        conditional=conditional, normalize=normalize, ...)
  },
  
  sample_emulator = function(input, em_pred_list=NULL, N_samp=1L, use_cov=TRUE, 
                             include_noise=TRUE, ...) {
    # Sample the forward model emulator at specified inputs. `input` has  
    # dimension (M,D) (where M is the number of inputs). Returns array of 
    # dimension (M, N_samp, N_output). 
    
    input <- .self$get_input(input)
    
    .self$emulator_model$sample(input, use_cov=use_cov, 
                                include_noise=include_noise, 
                                N_samp=N_samp, pred_list=em_pred_list, ...)
  },
  
  .sample = function(input, lik_par_val=NULL, em_pred_list=NULL, N_samp=1L,
                    use_cov=TRUE, conditional=default_conditional, 
                    normalize=default_normalize, ...) {
    # Sample the log-likelihood emulator at specified inputs. `input` has 
    # dimension (M,D) (where M is the number of inputs). Returns matrix of 
    # dimension (M, N_samp). 
    
    fwd_model_samp <- .self$sample_emulator(input, em_pred_list, N_samp=N_samp, 
                                            use_cov=use_cov, ...)
    llik_samp <- matrix(nrow=nrow(input), ncol=N_samp)

    for(i in 1:N_samp) {
      llik_samp[,i] <- .self$assemble_llik(matrix(fwd_model_samp[,i,], nrow=nrow(input), ncol=.self$N_output), 
                                           lik_par_val, conditional, normalize)
    }
    
    return(llik_samp)
  }, 
  
  .predict = function(input, lik_par_val=NULL, em_pred_list=NULL, return_mean=TRUE,  
                      return_var=TRUE, return_cov=FALSE, return_cross_cov=FALSE, 
                      input_cross=NULL, conditional=default_conditional,  
                      normalize=default_normalize, include_noise=TRUE, ...) {
    
    if(return_cross_cov || return_cov) {
      stop("`return_cross_cov` and `return_cov` are not yet supported for `llikEmulatorGP$predict()`.")
    }
    
    # Forward model emulator predictions. 
    if(is.null(em_pred_list)) {
      em_pred_list <- .self$predict_emulator(input, 
                                             return_mean=TRUE, 
                                             return_var=TRUE, 
                                             return_cov=return_cov, 
                                             return_cross_cov=return_cross_cov, 
                                             input_cross=input_cross, 
                                             include_noise=include_noise, ...)
    } else {
      assert_that(!is.null(em_pred_list$mean) && !is.null(em_pred_list$var), 
                  msg=paste0("Log-likelihood predictive quantities require ", 
                             "forward model emulator mean and variance."))
    }
    
    # Compute induced likelihood emulator predictions. 
    llik_pred_list <- list()
    sig2_val <- .self$get_lik_par(lik_par_val)
    
    if(return_mean) {
      llik_plug_in_mean <- .self$assemble_llik(em_pred_list$mean, 
                                               lik_par_val=lik_par_val, 
                                               conditional=conditional, 
                                               normalize=normalize, ...)
      llik_var_inflation <- -0.5 * .self$N_obs * rowSums(exp(log(em_pred_list$var) - rep(log(sig2_val), each=nrow(input))))
      llik_pred_list$mean <- llik_plug_in_mean + llik_var_inflation
    }
    
    if(return_var) {
      vars <- vector(mode="numeric", length=nrow(input))
      for(i in seq_along(vars)) {
        mse <- colSums(add_vec_to_mat_rows(-em_pred_list$mean[i,], .self$y)^2)
        vars[i] <- .self$N_obs * sum((em_pred_list$var[i,] / sig2_val)^2) + 
                    sum(mse * em_pred_list$var[i,] / sig2_val^2)       
      }
      
      llik_pred_list$var <- vars
    }
    
    return(llik_pred_list)
  },
  
  # TODO: need to think about normalization here. For the marginal approximation, the 
  # determinant part of the Gaussian likelihood depends on u through the GP predictive 
  # variance. 
  # TODO: currently the variance calculations below implicitly assume a normalized 
  # likelihood; need to fix this. 
  .predict_lik = function(input, lik_par_val=NULL, em_pred_list=NULL, return_mean=TRUE,  
                          return_var=TRUE, return_cov=FALSE, return_cross_cov=FALSE, 
                          input_cross=NULL, conditional=default_conditional,  
                          normalize=default_normalize, log_scale=FALSE, 
                          include_noise=TRUE, ...) {
    # The mean/variance of the likelihood emulator in this case are available 
    # in closed-form, though note that the distribution is not known, so the 
    # mean and variance might not provide a full characterization of the 
    # emulator distribution. The mean returned here corresponds to the 
    # "marginal" likelihood approximation under a Gaussian error model, which 
    # implies that the forward model GP variance is simply added to the 
    # likelihood variance. 

    assert_that(normalize && !conditional, 
                msg=paste("`predict_lik()` for llikEmulatorGP currently",
                          "assumes normalized and not conditional likelihood."))
    
    if(return_cross_cov || return_cov) {
      stop(paste("`return_cross_cov` and `return_cov` are not yet supported", 
                 "for `llikEmulatorGP$predict_lik()`."))
    }
    
    # Forward model emulator predictions. 
    if(is.null(em_pred_list)) {
      em_pred_list <- .self$predict_emulator(input, return_mean=TRUE, 
                                             return_var=TRUE, return_cov=return_cov, 
                                             return_cross_cov=return_cross_cov, 
                                             input_cross=input_cross, 
                                             include_noise=include_noise, ...)
    } else {
      assert_that(!is.null(em_pred_list$mean) && !is.null(em_pred_list$var), 
                  msg=paste("Likelihood predictive quantities require forward",
                            "model emulator mean and variance."))
    }
    
    # Compute induced likelihood emulator predictions. 
    lik_pred_list <- list()
    sig2 <- .self$get_lik_par(lik_par_val)
    if(return_mean) {
      lik_pred_list$log_mean <- .self$assemble_llik(em_pred_list$mean, lik_par_val=sig2, 
                                                    conditional=conditional, 
                                                    normalize=normalize, 
                                                    var_inflation_vals=em_pred_list$var, ...)
      if(!log_scale) lik_pred_list$mean <- exp(lik_pred_list$log_mean)
    }
    
    if(return_var) {
      log_num1 <- .self$assemble_llik(em_pred_list$mean, lik_par_val=0.5 * sig2, 
                                      conditional=conditional, normalize=normalize, 
                                      var_inflation_vals=em_pred_list$var, ...)
      log_num2 <- .self$assemble_llik(em_pred_list$mean, lik_par_val=0.5 * sig2, 
                                      conditional=conditional, normalize=normalize, 
                                      var_inflation_vals=0.5 * em_pred_list$var, ...)
      log_2pi_term <- 0.5 * .self$N_obs * .self$N_output * log(2*pi)
      log_denom1 <-  log_2pi_term  + 0.5 * .self$N_obs * sum(log(sig2))
      log_denom2 <- log_2pi_term + 0.5 * .self$N_obs * rowSums(log(add_vec_to_mat_rows(sig2, em_pred_list$var)))
      lik_pred_list$log_var <- log_diff_exp(log_num1 - log_denom1, log_num2 - log_denom2)
      if(!log_scale) lik_pred_list$var <- exp(lik_pred_list$log_var)
    }
  
    return(lik_pred_list)
  }, 
  
  .calc_quantiles = function(input, p, lik_par_val=NULL, em_pred_list=NULL,
                             conditional=default_conditional, normalize=default_normalize,
                             lower_tail=TRUE, include_noise=TRUE, n_mc=1e3, ...) {
    # Monte Carlo approximation.
    
    if(!lower_tail) p <- 1 - p
    llik_samp <- .self$sample(input_grid, N_samp=n_mc, use_cov=FALSE, apply_shift=FALSE,
                              conditional=conditional, normalize=normalize, 
                              include_noise=include_noise, em_pred_list=em_pred_list, ...)
    q <- apply(llik_samp, 1L, function(x) quantile(x, p))
  },
  
  .calc_lik_quantiles = function(input, p, lik_par_val=NULL, em_pred_list=NULL,
                                 conditional=default_conditional, normalize=default_normalize,
                                 lower_tail=TRUE, include_noise=TRUE, log_scale=TRUE, n_mc=1e3, ...) {
    # Monte Carlo approximation. Note that this is the same as the exponential of the
    # log-likelihood quantiles, since the exponential is strictly increasing.
    
    log_q <- .self$calc_quantiles(input, p, lik_par_val=lik_par_val, em_pred_list=em_pred_list,
                                  conditional=conditional, normalize=normalize, lower_tail=lower_tail,
                                  include_noise=include_noise, n_mc=n_mc, ...)
    
    if(log_scale) return(log_q)
    return(exp(log_q))
  },
  
  calc_lik_marginal_approx = function(em_pred_list=NULL, input=NULL,
                                      lik_par_val=NULL, 
                                      conditional=default_conditional, 
                                      normalize=default_normalize,
                                      log_scale=TRUE, include_noise=TRUE, ...) {
    
    sel <- ifelse(log_scale, "log_mean", "mean")
    .self$predict_lik(input, lik_par_val=lik_par_val, em_pred_list=em_pred_list, 
                      return_mean=TRUE, return_var=FALSE, conditional=conditional,  
                      normalize=normalize, log_scale=log_scale, 
                      include_noise=include_noise, ...)[[sel]]
  },
  
  calc_expected_lik_cond_var = function(input_eval, input_cond, lik_par_val=NULL,
                                        include_noise=TRUE, log_scale=TRUE, plugin=FALSE, 
                                        conditional=default_conditional, 
                                        normalize=default_normalize, ...) {
    assert_that(plugin, 
                msg=paste("llikEmulatorGPFwdGaussDiag$calc_expected_lik_cond_var()",
                          "requires `plugin=TRUE`."))
    
    N_eval <- nrow(input_eval)
    sig2 <- .self$get_lik_par(lik_par_val)
    gp_copy <- .self$emulator_model$copy(shallow=FALSE)
    
    # GP emulator predictions at evaluation locations and at conditioning locations. 
    gp_pred <- gp_copy$predict(input_cond, return_mean=TRUE, return_var=FALSE, ...)
    gp_pred_eval <- gp_copy$predict(input_eval, return_mean=TRUE, return_var=TRUE, ...)
    
    # Update the GP model, treating the predictive mean as the observed response  
    # at the evaluation locations. This leaves the predictive mean unchanged, 
    # but updates the predictive variance. 
    gp_copy$update(input_cond, gp_pred$mean, update_hyperpar=FALSE, ...)
    
    # Predict with the conditional ("cond") GP (i.e., the updated GP) at the  
    # evaluation locations.
    gp_pred_cond <- gp_copy$predict(input_eval, return_mean=FALSE, 
                                    return_var=TRUE, ...)
    
    # Compute the log of the numerators of the first and second terms. 
    log_num1 <- assemble_llik(gp_pred_eval$mean, lik_par_val=0.5*sig2, 
                              conditional=conditional, normalize=normalize, 
                              var_inflation_vals=gp_pred_eval$var, ...)
    log_num2 <- assemble_llik(gp_pred_eval$mean, lik_par_val=0.5*sig2, 
                              conditional=conditional, normalize=normalize, 
                              var_inflation_vals=gp_pred_eval$var-0.5*gp_pred_cond$var, ...)
  
    # Compute the log of the denominators of the first and second terms. 
    log_2pi_term <- 0.5 * .self$N_obs * .self$N_output * log(2*pi)
    log_denom1 <-  log_2pi_term  + 0.5 * .self$N_obs * sum(log(sig2))
    log_denom2 <- log_2pi_term + 0.5 * .self$N_obs * rowSums(log(add_vec_to_mat_rows(sig2, gp_pred_cond$var)))
    log_evar <- log_diff_exp(log_num1-log_denom1, log_num2-log_denom2)
    
    if(log_scale) return(log_evar)
    return(exp(log_evar))
  }
  
)


# -----------------------------------------------------------------------------
# OLD/IN DEVELOPMENT: 
# These classes are not in a functional state. 
#  (1) llikSumEmulator: a sum of llikEmulator objects, representing a llik 
#                       that can be written as a sum of terms, in which each 
#                       term is emulated independently.
#  (2) llikEmulatorMultGausGP: Gaussian likelihood, GP emulating the quadratic 
#                              model-data misfit. 
#  (3) llikEmulatorExactGauss: Exact Gaussian likelihood with forward model 
#                              and general (dense) covariance matrix.
#  (4) llikSumEmulatorMultGausGP: A different approach `llikEmulatorMultGausGP`.
#                                 Will probably be deleted.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# llikSumEmulator Class 
# -----------------------------------------------------------------------------

#
# TODO: this needs lots of updates
#    - Add argument `em_pred_list` to the relevant functions. 
#    - Update quantile/CI functions. 
#    - Update plotting functions. 
#

llikSumEmulator <- setRefClass(
  Class = "llikSumEmulator", 
  contains = "llikEmulator",
  fields = list(llik_emulator_terms="list", N_terms="integer")
)

llikSumEmulator$methods(
  
  initialize = function(llik_emulator_list, default_conditional=FALSE, default_normalize=FALSE, 
                        lik_description="Sum of llik terms.", emulator_description="", ...) {
    
    assert_that(all(sapply(llik_emulator_list, function(obj) inherits(obj, "llikEmulator"))), 
                msg="`llikEmulator_list` must be a list of llikEmulator objects.")
    assert_that(length(unique(sapply(llik_emulator_list, function(obj) obj$dim_input)))==1, 
                msg="Inputs of llik terms have differing dimensions.")
    
    # Set labels for each llik emulator term, ensuring they are unique. The names attribute of 
    # the list in the field `llik_emulator_terms` is also assigned the vector of labels, ensuring
    # the order of the two is the same. 
    term_lbls <- sapply(llik_emulator_list, function(obj) obj$llik_label)
    assert_that(length(unique(term_lbls))==length(llik_emulator_list), 
                msg="Found duplicate `llik_label` attributes in `llik_emulator_list`.")
    names(llik_emulator_list) <- term_lbls
    initFields(llik_emulator_terms=llik_emulator_list, N_terms=length(llik_emulator_list))
    
    # Ensure that the `input_names` are the same for all llik terms. They are allowed 
    # to be ordered differently for each term, but viewed as sets they must be equal. 
    equal_set <- function(x,y) if(setequal(x,y)) x else FALSE
    input_names_list <- lapply(llik_emulator_list, function(obj) obj$input_names)
    assert_that(!isFALSE(Reduce(equal_set, input_names_list)), msg="Different `input_names` found for different llik terms.")
    
    # Since the llikEmulator terms are assumed to have independent predictive distributions, then the 
    # llikSumEmulator predictive distribution is Gaussian provided all of the terms are Gaussian. Otherwise, 
    # the distribution cannot be determined without further information. 
    llik_term_dists <- sapply(llik_emulator_list, function(x) x$llik_pred_dist)
    if(all(llik_term_dists == "Gaussian")) llik_sum_dist <- "Gaussian"
    else llik_sum_dist <- "unspecified"
    
    # The llikSumEmulator sum emulator is only exact if all of its component terms are exact. 
    llik_terms_exact <- sapply(llik_emulator_list, function(x) x$exact_llik)
    llik_sum_exact <- all(llik_terms_exact)
    
    callSuper(llik_label=term_lbls, lik_description=lik_description, emulator_description=emulator_description,
              emulator_model=NULL, default_conditional=default_conditional, 
              default_normalize=default_normalize, dim_input=llik_emulator_list[[1]]$dim_input,
              input_names=llik_emulator_list[[1]]$input_names,
              lik_par=NULL, llik_pred_dist=llik_sum_dist, exact_llik=llik_sum_exact)
    
  },
  
  get_llik_attr = function(attr_name, labels=llik_label) {
    
    attr_list <- list()
    for(i in seq_along(labels)) {
      lbl <- labels[i]
      attr_list[[lbl]] <- llik_emulator_terms[[lbl]]$field(attr_name)
    }
    
    return(attr_list)
    
  },
  
  get_lik_par = function(lik_par_val=NULL, labels=llik_label) {
    lik_par_list <- list()
    for(i in seq_along(labels)) {
      lbl <- labels[i]
      lik_par_list[[lbl]] <- llik_emulator_list[[lbl]]$get_lik_par(lik_par_val[[lbl]])
    }
    
    return(lik_par_list)
  },
  
  input_designs_equal = function(labels=llik_label, design_input_list=NULL, ...) {
    equal_design <- function(x,y) if(all.equal(x[,input_names,drop=FALSE], y[,input_names,drop=FALSE])) x else FALSE
    
    if(is.null(design_input_list)) {
      design_input_list <- list()
      for(lbl in labels) {
        design_input_list[[lbl]] <- llik_emulator_terms[[lbl]]$get_design_inputs(...)
      }
    }
    
    !isFALSE(Reduce(equal_design, design_input_list))
  },
  
  sample_emulator = function(input, N_samp=1, labels=llik_label, ...) {
    
    samp_list <- list()
    for(i in seq_along(labels)) {
      lbl <- labels[i]
      samp_list[[lbl]] <- llik_emulator_terms[[lbl]]$sample_emulator(input, N_samp=N_samp, ...)
    }   
    
    return(samp_list)
  },
  
  assemble_llik = function(emulator_vals_list, lik_par=NULL, conditional=default_conditional, 
                           normalize=default_normalize, sum_terms=TRUE, labels=names(emulator_vals_list), ...) {
    
    llik_list <- list()
    for(lbl in labels) {
      llik_list[[lbl]] <- llik_emulator_terms[[lbl]]$assemble_llik(emulator_vals_list[[lbl]], lik_par[[lbl]], 
                                                                   conditional=conditional, normalize=normalize)
    }
    
    # Combine into single array. Term label names assigned to third dimension. 
    llik_vals <- abind(llik_list, along=3)
    
    if(sum_terms) llik_vals <- rowSums(llik_vals, dims=2)
    return(llik_vals)
  },
  
  
  predict = function(input, lik_par_val=NULL, return_mean=TRUE, return_var=TRUE, 
                     return_cov=FALSE, return_cross_cov=FALSE, input_cross=NULL,
                     conditional=default_conditional, 
                     normalize=default_normalize, sum_terms=TRUE, 
                     labels=llik_label, ...) {
    
    predict_list <- list()
    for(lbl in labels) {
      predict_list[[lbl]] <- llik_emulator_terms[[lbl]]$predict(input, lik_par_val=lik_par_val[[lbl]],
                                                                return_mean=return_mean, return_var=return_var,
                                                                return_cov=return_cov, return_cross_cov=return_cross_cov,
                                                                input_cross=input_cross, conditional=conditional,
                                                                normalize=normalize, ...)
    }
    
    if(!sum_terms) return(predict_list)
    
    sum_list <- list()
    if(return_mean) sum_list$mean <- Reduce("+", lapply(predict_list, function(l) l$mean))
    if(return_cov) {
      sum_list$cov <- Reduce("+", lapply(predict_list, function(l) l$cov))
      sum_list$var <- diag(sum_list$cov)
    } else if(return_var) {
      sum_list$var <- Reduce("+", lapply(predict_list, function(l) l$var))
    }
    
    if(return_cross_cov) {
      sum_list$cross_cov <- Reduce("+", lapply(predict_list, function(l) l$cross_cov))
    }
    
    return(sum_list)
  },
  
  get_design_inputs = function(return_list=FALSE, labels=llik_label, ...) {
    design_input_list <- list()
    for(lbl in labels) {
      design_input_list[[lbl]] <- llik_emulator_terms[[lbl]]$get_design_inputs(...)
    }
    
    if(return_list) return(design_input_list)
    
    # If not returning list, then requires the design points to be the same across 
    # all emulator terms. 
    assert_that(input_designs_equal(design_input_list=design_input_list), 
                msg="`return_list==FALSE` requires that the design input points are equal for each llik emulator term.")
    
    return(design_input_list[[1]][,input_names, drop=FALSE])
    
  },
  
  get_design_llik = function(lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, 
                             return_list=FALSE, labels=llik_label, ...) {
    design_llik_list <- list()
    for(lbl in labels) {
      design_llik_list[[lbl]] <- llik_emulator_terms[[lbl]]$get_design_llik(lik_par_val=lik_par_val[[lbl]], 
                                                                            conditional=conditional,
                                                                            normalize=normalize, ...)
    }
    
    if(return_list) return(design_llik_list)
    
    # If not returning list, then requires the design input points to be the same across 
    # all emulator terms.
    assert_that(input_designs_equal(labels=labels), 
                msg="`return_list==FALSE` requires that the design input points are equal for each llik emulator term.")    
    return(Reduce("+", design_llik_list))
  },
  
  sample = function(input, lik_par_val=NULL, N_samp=1, conditional=default_conditional, 
                    normalize=default_normalize, sum_terms=TRUE, labels=llik_label, ...) {
    
    emulator_samp_list <- sample_emulator(input, N_samp=N_samp, labels=labels, ...)
    assemble_llik(emulator_samp_list, lik_par_val, conditional=conditional, normalize=normalize,
                  sum_terms=sum_terms, labels=labels)
  }, 
  
  # TODO: need to update this to reflect changes in `llikEmulator$plot_llik_samp_1d`; i.e. the 
  # addition of `plot_type` argument and the changes due to this. 
  plot_llik_samp_1d = function(input, lik_par_val=NULL, N_samp=1, conditional=default_conditional, 
                               normalize=default_normalize, true_llik=NULL, include_design=FALSE, 
                               labels=llik_label, sum_terms=TRUE, ...) {
    # `true_llik` should be a vector (or 1 col matrix) if `sum_terms` is TRUE. Otherwise should
    # be a matrix with `N_terms` cols, with colnames set to llik labels. 
    
    assert_that(dim_input==1, msg=paste0("plot_llik_samp_1d() requires 1d input space. input_dim = ", dim_input))
    input <- get_input(input)
    
    if(sum_terms) {
      llik_samp <- .self$sample(input, lik_par=lik_par_val, N_samp=N_samp, conditional=conditional,
                                normalize=normalize, sum_terms=TRUE, labels=labels, ...)
      plt <- ggmatplot(input, llik_samp, plot_type="line", color="gray") +
        theme(legend.position = "none") +
        ggtitle("Log Likelihood Samples") +
        xlab(input_names) + ylab("Log Likelihood")
      
      if(!is.null(true_llik)) {
        true_llik <- drop(true_llik)
        assert_that(is.numeric(true_llik) && length(true_llik)==nrow(input), 
                    msg="With `sum_terms==TRUE`, `true_llik` must be vector of length `nrow(input)`.")
        df <- data.frame(x=input[,1], y=true_llik)
        plt <- plt + geom_line(aes(x=x, y=y), df, inherit.aes=FALSE, color="red")
      }
      
      if(include_design) {
        design_df <- data.frame(x=drop(get_design_inputs(return_list=FALSE, labels=labels)),
                                y=drop(get_design_llik(lik_par_val, conditional, normalize, return_list=FALSE, labels=labels)))
        plt <- plt + geom_point(aes(x=x, y=y), design_df, inherit.aes=FALSE, color="red")
      }
      
      return(plt)
    }
    
    # If not summing terms, produce one plot per llik term. Simply call the plot function 
    # for each llik emulator term. 
    plts <- list()
    for(lbl in labels) {
      plts[[lbl]] <- llik_emulator_terms[[lbl]]$plot_llik_samp_1d(input, lik_par_val[[lbl]], N_samp,
                                                                  conditional, normalize, true_llik[[lbl]],
                                                                  include_design, ...)
    }
    
    return(plts)
  },
  
  plot_llik_pred_1d = function(input, lik_par_val=NULL, conditional=default_conditional,
                               normalize=default_normalize, include_CI=FALSE, CI_prob=0.9, 
                               llik_pred_list=NULL, true_llik=NULL, sum_terms=TRUE, labels=llik_label, 
                               plot_title=NULL, ...) {
    
    assert_that(dim_input==1, msg=paste0("plot_llik_pred_1d() requires 1d input space. dim_input = ", dim_input))
    
    # If `sum_terms==TRUE` produce plot based on the llikSumEmulator predictions directly. 
    if(sum_terms) {
      if(is.null(llik_pred_list)) {
        llik_pred_list <- .self$predict(input, lik_par_val=lik_par_val, return_mean=TRUE, return_var=TRUE, 
                                        conditional=conditional, normalize=normalize, sum_terms=TRUE, labels=labels, ...)
      }
      
      # Plot title and labels.
      if(is.null(plot_title)) {
        plot_title <- "Log Likelihood Sum Emulator Predictions"
        if(include_CI) plot_title <- paste0(plot_title, ", ", 100*CI_prob, "% CI")
      }
      xlab <- input_names
      ylab <- "llik sum"
      if(!normalize) {
        ylab <- paste0(ylab, ", ", ifelse(conditional, "unnormalized/conditional", "unnormalized"))
      }
      
      # Compute confidence interval. 
      if(include_CI) CI_list <- .self$calc_confidence_interval(llik_pred_list=llik_pred_list, CI_prob=CI_prob, ...)
      
      # Produce plot. 
      plt <- plot_pred_1d_helper(X_new=drop(input), pred_mean=drop(llik_pred_list$mean), 
                                 include_CI=include_CI, CI_lower=CI_list$lower, CI_upper=CI_list$upper, 
                                 y_new=drop(true_llik), X_design=drop(.self$get_design_inputs(return_list=FALSE, labels=labels, ...)), 
                                 y_design=drop(get_design_llik(lik_par_val=lik_par_val, conditional=conditional,  
                                                               normalize=normalize, return_list=FALSE, labels=labels, ...)), 
                                 plot_title=plot_title, xlab=xlab, ylab=ylab)
      return(plt)
    }
    
    # If not summing terms, produce one plot per llik term. Simply call the plot function 
    # for each llik emulator term. 
    plts <- list()
    for(lbl in labels) {
      plts[[lbl]] <- llik_emulator_terms[[lbl]]$plot_llik_pred_1d(input, lik_par_val=lik_par_val[[lbl]],
                                                                  conditional=conditional, normalize=normalize, 
                                                                  true_llik=true_llik[[lbl]], include_CI=include_CI,
                                                                  CI_prob=CI_prob, llik_pred_list=llik_pred_list[[lbl]], 
                                                                  ylab=paste0("llik term: ", lbl), ...)
    }
    
    return(plts)
  },
  
  calc_quantiles = function(p, input=NULL, lik_par_val=NULL, conditional=default_conditional, 
                            normalize=default_normalize, llik_pred_list=NULL,
                            lower_tail=TRUE, labels=llik_label, sum_terms=TRUE, ...) {
    # If `sum_terms==TRUE` then `llik_pred_list` must be summed across terms. If `sum_terms==FALSE` then 
    # it must be a list of the llik pred lists for each term. 
    
    # If `sum_terms==TRUE` then only proceed if the predictive distribution of the llik sum is 
    # Gaussian. Otherwise, the llik distribution cannot be determined from the llik term distributions. 
    if(sum_terms) {
      assert_that(llik_pred_dist=="Gaussian", 
                  msg="`calc_quantiles()` method of `llikSumEmulator` currently only supports `sum_terms==TRUE` when `llik_pred_dist=='Gaussian'`")
      
      if(is.null(llik_pred_list)) {
        llik_pred_list <- .self$predict(input, lik_par_val=lik_par_val, return_mean=TRUE, return_var=TRUE, 
                                        conditional=conditional, normalize=normalize, sum_terms=TRUE, labels=labels, ...)
      }
      
      return(qnorm(p, drop(llik_pred_list$mean), sqrt(drop(llik_pred_list$var)), lower.tail=lower_tail))
    }
    
    # Otherwise call the quantile function for each term separately and return list of the results. 
    quantiles_list <- list()
    for(lbl in labels) {
      quantiles_list[[lbl]] <- llik_emulator_terms[[lbl]]$calc_quantiles(p, input=input, lik_par_val=lik_par_val[[lbl]], 
                                                                         conditional=conditional, normalize=normalize, 
                                                                         lower_tail=lower_tail, 
                                                                         llik_pred_list=llik_pred_list[[lbl]], ...)
    }
    
    return(quantiles_list)
  } 
  
)

# -----------------------------------------------------------------------------
# llikEmulatorMultGausGP class
#    Surrogate for a multiplicative Gaussian likelihood with a single variance
#    parameter. A GP emulator models the map u -> ||y - G(u)||^2.
# -----------------------------------------------------------------------------

llikEmulatorMultGausGP <- setRefClass(
  Class = "llikEmulatorMultGausGP", 
  contains = "llikEmulator",
  fields = list(N_obs="integer")
)

llikEmulatorMultGausGP$methods(
  
  initialize = function(gp_model, llik_lbl, N_obs, sig2=NULL, default_conditional=FALSE, 
                        default_normalize=FALSE, ...) {
    assert_that(inherits(gp_model, "gpWrapper"), msg="`gp_model` must inherit from `gpWrapper` class.")
    assert_that(is.integer(N_obs) && (length(N_obs) == 1) && (N_obs > 0),
                msg="`N_obs` must be an integer greater than 0.")
    assert_that(gp_model$Y_dim==1, msg="`llikEmulatorMultGausGP` only supports single-output GP emulator.")
    assert_that(!is.null(gp_model$X_names) && noNA(gp_model$X_names), 
                msg="`llikEmulatorMultGausGP` requires that `gp_model` has `X_names` field set.")
    if(!is.null(sig2)) assert_that(is.numeric(sig2) && (sig2>0), msg="`sig2` must be NULL or numeric positive value.")
    
    initFields(N_obs=N_obs)
    callSuper(emulator_model=gp_model, llik_label=llik_lbl, lik_par=sig2, input_names=gp_model$X_names,
              dim_input=gp_model$X_dim, default_conditional=default_conditional, 
              default_normalize=default_normalize, 
              lik_description="Multiplicative Gaussian.",
              emulator_description="GP emulating sum of squared error function.", 
              llik_pred_dist="Gaussian", exact_llik=FALSE, ...)
  }, 
  
  assemble_llik = function(SSR, lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize) {
    # SSR should be N_input x N_samp. 
    
    # Fetch the variance parameter. 
    sig2 <- get_lik_par(lik_par_val)
    
    # Construct likelihood using SSR.  
    llik <- -0.5 * SSR / sig2
    if(normalize || !conditional) llik <- llik - 0.5*N_obs*log(sig2)
    if(normalize) llik <- llik - 0.5*N_obs*log(2*pi)
    
    return(drop(llik))
    
  },
  
  get_design_inputs = function(...) {
    emulator_model$X
  },
  
  get_design_llik = function(lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, ...) {
    assemble_llik(emulator_model$Y, lik_par_val=lik_par_val, conditional=conditional, normalize=normalize)
  },
  
  sample_emulator = function(input, em_pred_list=NULL, N_samp=1, use_cov=FALSE, 
                             include_noise=TRUE, adjustment="rectified", ...) {
    
    emulator_model$sample(get_input(input), use_cov=use_cov, 
                          include_noise=include_noise, 
                          N_samp=N_samp, adjustment=adjustment, 
                          pred_list=em_pred_list)[,,1,drop=FALSE]             
  },
  
  sample = function(input, lik_par_val=NULL, em_pred_list=NULL, N_samp=1, 
                    use_cov=FALSE, include_noise=TRUE, 
                    conditional=default_conditional, normalize=default_normalize, ...) {
    # Sample SSR. 
    samp <- sample_emulator(input, em_pred_list, N_samp, use_cov, 
                            include_noise, ...)
    
    # Compute unnormalized or normalized log-likelihood. 
    assemble_llik(samp, lik_par_val, conditional, normalize)
  }, 
  
  predict = function(input, lik_par_val=NULL, em_pred_list=NULL, return_mean=TRUE, 
                     return_var=TRUE,  return_cov=FALSE, return_cross_cov=FALSE, 
                     input_cross=NULL, conditional=default_conditional, 
                     normalize=default_normalize, include_noise=TRUE, ...) {
    
    if(is.null(em_pred_list)) {
      pred_list <- .self$emulator_model$predict(get_input(input), return_mean=return_mean, return_var=return_var,
                                                return_cov=return_cov, return_cross_cov=return_cross_cov,
                                                X_cross=input_cross, include_noise=include_noise, ...)
    } else {
      pred_list <- em_pred_list
    }
    
    if(return_mean) {
      pred_list$mean <- .self$assemble_llik(pred_list$mean, lik_par=lik_par_val, conditional, normalize)
    }
    
    if(return_cov) {
      pred_list$cov <- pred_list$cov[,,1] / (4*lik_par_val^2)
      pred_list$var <- diag(pred_list$cov)
    } else if(return_var) {
      pred_list$var <- pred_list$var / (4*lik_par_val^2)
    }
    
    if(return_cross_cov) {
      pred_list$cross_cov <- pred_list$cov / (4*lik_par_val^2)
    }
    
    return(pred_list)
  },
  
  
  calc_quantiles = function(p, input=NULL, lik_par_val=NULL, em_pred_list=NULL, llik_pred_list=NULL,
                            target="llik", conditional=default_conditional, normalize=default_normalize, 
                            lower_tail=TRUE, include_noise=TRUE, ...) {
    # The log-likelihood emulator distribution is Gaussian, and the likelihood emulator is log-normal. 
    # Both the Gaussian and Log-normal quantile functions are parameterized in terms of the underlying 
    # Gaussian mean/variance, so only `llik_pred_list` is required, regardless of whether `target`
    # is "llik" or "lik". 
    
    assert_that(target %in% c("llik", "lik"))
    
    # Log-likelihood or likelihood predictions. 
    if(is.null(llik_pred_list)) {
      llik_pred_list <- .self$predict(input, lik_par_val=lik_par_val, em_pred_list=em_pred_list, 
                                      return_mean=TRUE, return_var=TRUE, conditional=conditional,
                                      normalize=normalize, log_scale=FALSE, ...)
    } else {
      assert_that(!is.null(llik_pred_list$mean) && !is.null(llik_pred_list$var))
    }
    
    # Log-likelihood emulator distribution is Gaussian. Likelihood emulator is log-normal. 
    if(target == "llik") q <- qnorm(p, drop(llik_pred_list$mean), sqrt(drop(llik_pred_list$var)), lower.tail=lower_tail)
    else q <- qlnorm(p, drop(llik_pred_list$mean), sqrt(drop(llik_pred_list$var)), lower.tail=lower_tail)
    
    return(q)
  }
  
)


# -----------------------------------------------------------------------------
# llikEmulatorExactGauss class
#
# This simply implements the exact likelihood corresponding to a  
# Gaussian inverse problem: y|u ~ N(G(u), Sig), where G may be nonlinear.
# The forward map G is defined by an attribute called `fwd_model` in the class, 
# which is a function that users pass in when instantiating the class. 
# This is exact in the sense that there is no emulation here - the reason for 
# implementing this as a llikEmulator class is to use it for algorithm testing; 
# e.g. ensuring the correctness of an MCMC implementation. The `lik_par` here 
# is defined to be the covariance matrix `Sig`. 
#
# In the case that `G` is linear and `u` is assigned a Gaussian prior, then 
# this yields a linear Gaussian inverse problem, implying that the posterior 
# u|y is Gaussian (when `Sig` is fixed). When `Sig` is not fixed and instead 
# assigned an inverse Wishart prior then the joint posterior u, Sig|y is 
# Normal Inverse Wishart.
#
# Note that the class llikEmulatorExactGaussDiag should be used in the 
# special case where `Sig` is constrained to be diagonal. 
# -----------------------------------------------------------------------------

llikEmulatorExactGauss <- setRefClass(
  Class = "llikEmulatorExactGauss", 
  contains = "llikEmulator",
  fields = list(fwd_model="ANY", fwd_model_vectorized="ANY", y="numeric", 
                N_obs="integer", L_Cov="matrix")
)

llikEmulatorExactGauss$methods(
  
  initialize = function(llik_lbl, y_obs, dim_par, fwd_model=NULL, fwd_model_vectorized=NULL, Cov=NULL, 
                        default_conditional=FALSE, default_normalize=FALSE, 
                        par_names=NULL, ...) {
    
    # Forward model must be provided, either vectorized or single-input version. 
    assert_that(is.null(fwd_model) || is.function(fwd_model))
    assert_that(is.null(fwd_model_vectorized) || is.function(fwd_model_vectorized))
    assert_that(is.function(fwd_model) || is.function(fwd_model_vectorized))
    initFields(fwd_model=fwd_model, fwd_model_vectorized=fwd_model_vectorized, 
               N_obs=length(drop(y_obs)), y=drop(y_obs))
    
    # Set parameter names. 
    if(is.null(par_names)) par_names <- paste0("input", 1:dim_par)
    assert_that(length(par_names) == dim_par)
    
    # Covariance matrix of Gaussian likelihood. 
    if(!is.null(Cov)) {
      assert_that(is.matrix(Cov) && (nrow(Cov)==N_obs) && (ncol(Cov)==N_obs),
                  msg="`Cov` must be a positive definite matrix with dim `N_obs` x `N_obs`")
      initFields(L_Cov=t(chol(Cov))) 
    }
    
    callSuper(emulator_model=NULL, llik_label=llik_lbl, lik_par=Cov, dim_input=dim_par,
              default_conditional=default_conditional, input_names=par_names,
              default_normalize=default_normalize, 
              lik_description="Exact Gaussian likelihood.",
              emulator_description="No emulation.", exact_llik=TRUE, ...)
  },
  
  run_fwd_model = function(input, ...) {
    # `input` is an M x D matrix with input parameter values stacked in the rows. 
    # `fwd_model_vectorized(input)` returns N_obs x M dimensional output. 
    
    if(is.null(.self$fwd_model_vectorized)) return(vectorize_fwd_model(input, ...))
    return(.self$fwd_model_vectorized(input, ...))
  },
  
  vectorize_fwd_model = function(input, ...) {
    model_output <- matrix(nrow=N_obs, ncol=nrow(input))
    for(i in 1:nrow(input)) model_output[,i] <- .self$fwd_model(input[i,], ...)
    return(model_output)
  },
  
  get_lik_par = function(lik_par_val=NULL, return_chol=FALSE, ...) {
    # If covariance matrix is passed, convert to Cholesky factor is requested.
    if(!is.null(lik_par_val)) {
      if(return_chol) lik_par_val <- t(chol(lik_par_val))
    } else {
      if(return_chol) lik_par_val <- .self$L_Cov
      else lik_par_val <- .self$lik_par
    }
    
    assert_that(!is.null(lik_par_val), 
                msg="Likelihood parameter (noise covariance) missing.")
    
    return(lik_par_val)
  },
  
  assemble_llik = function(input, lik_par_val=NULL, 
                           conditional=default_conditional, 
                           normalize=default_normalize, ...) {
    # `input` should be N_input x N_samp. 
    
    # Fetch the lower triangular Cholesky factor of the covariance matrix.
    L <- .self$get_lik_par(lik_par_val, return_chol=TRUE)
    
    # Construct log likelihood. 
    llik <- -0.5 * colSums(solve(L, y - .self$run_fwd_model(input, ...))^2, na.rm=TRUE)
    if(normalize || !conditional) llik <- llik - sum(log(diag(L)))
    if(normalize) llik <- llik - 0.5*N_obs*log(2*pi)
    
    return(drop(llik))
  }, 
  
  sample_emulator = function(input, em_pred_list=NULL, N_samp=1, ...) {
    # No emulator to sample from, simply return input. The argument 
    # `em_pred_list` is only present for consistency with other llikEmulator
    # classes. 
    
    input
  },
  
  update_emulator = function(...) {
    # There is no emulator, so this does nothing. 
  },
  
  sample = function(input, lik_par=NULL, em_pred_list=NULL, N_samp=1, 
                    conditional=default_conditional, normalize=default_normalize, ...) {
    # Compute unnormalized or normalized log-likelihood (exact, deterministic 
    # calculation - no sampling is actually performed). For consistency with 
    # other classes, duplicates the exact likelihood calculation when `N_samp`>1.
    matrix(assemble_llik(get_input(input), lik_par, conditional, normalize), 
           nrow=nrow(input), ncol=N_samp)
  }
  
)


# -----------------------------------------------------------------------------
# llikSumEmulatorMultGausGP: Encapsulates an approximation of a multiplicative 
# Gaussian likelihood where the "sum of squared error" functions have been 
# approximated by independent GPs. In particular, considers log-likelihood 
# of the form: 
#   log p(Y|u, Sig) = C - 0.5 * sum_{p=1}^{P} sum_{t=1}^{T_p} (Y_{tp} - G(u)_{tp})^2 / sig2_p
# where C is a constant and G a function. The mappings 
#   Phi(u) := sum_{t=1}^{T_p} (Y_{tp} - G(u)_{tp})^2
# are assumed to be emulated by independent GPs. Hence the `model` field of this 
# class is required to inherit from the `gpWrapper` class. This class interprets
# the `input` argument as `u` (typically the forward model parameters, which 
# are of primary interest). This class also requires a second argument `sig2` to its 
# core functions, which is the vector of variance parameters. 
#
# By default, the class represents a "conditional" llik (meaning conditional 
# on some current value of `sig2`), which practically means that setting 
# `normalize=FALSE` will drop additive constants that depend on `sig2`. 
# By setting `default_conditional=FALSE`, this behavior is reversed and 
# terms involving `sig2` are no longer interpreted as normalizing constants. 
# -----------------------------------------------------------------------------

#
# TODO: should probably also have method for computing the exact likelihood. i.e., 
#       pass in either the sufficient statistic or model output and it will compute. 
#

llikSumEmulatorMultGausGP <- setRefClass(
  Class = "llikSumEmulatorMultGausGP", 
  contains = "llikEmulator",
  fields = list(N_obs="integer", N_output="integer", sum_obs="integer")
)

llikSumEmulatorMultGausGP$methods(
  
  initialize = function(gp_model, N_obs, N_output, sig2=NULL, default_conditional=TRUE, 
                        default_normalize=FALSE, ...) {
    assert_that(inherits(gp_model, "gpWrapper"), msg="`gp_model` must inherit from `gpWrapper` class.")
    assert_that(is.integer(N_output) && N_output>0, msg="`N_output` must be an integer greater than 0.")
    assert_that(is.integer(N_obs) && (length(N_obs)==N_output), 
                msg="`N_obs` must be integer vector of length equal to `N_output`.")
    assert_that(gp_model$Y_dim==N_output, msg="Number of independent GP emulators must equal `N_output`.")
    if(!is.null(sig2)) {
      assert_that(is.numeric(sig2) && (length(sig2)==N_output) && all(sig2>0), 
                  msg="`sig2` must be NULL or numeric vector of length `N_output`.")
    }
    
    initFields(N_obs=N_obs, N_output=N_output, default_conditional=default_conditional, 
               default_normalize=default_normalize, sum_obs=sum(N_obs))
    callSuper(lik_description="Multiplicative Gaussian.",
              emulator_description="Independent GPs emulating sum of squared error functions.",
              emulator_model=gp_model, ...)
  },
  
  
  sample = function(input, sig2, N_samp=1, use_cov=FALSE, include_noise=TRUE, sum_output_llik=TRUE,
                    conditional=default_conditional, normalize=FALSE, ...) {
    samp <- emulator_model$sample(input, use_cov=use_cov, include_noise=include_noise, N_samp=N_samp)
    
    for(j in 1:N_output) {
      samp[,,j] <- -0.5 * samp[,,j] / sig2[j]
      if(normalize || !conditional) samp[,,j] <- samp[,,j] - 0.5*N_obs[j]*log(sig2[j])
    }
    
    if(normalize) samp <- samp - 0.5*sum_obs*log(2*pi)
    if(sum_output_llik) return(rowSums(samp, dims=2))
    return(samp)
  }
  
)



