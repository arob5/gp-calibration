#
# llikEmulator.r
# Class definitions for the `llikEmulator`, `llikSumEmulator`, and classes which 
# inherit from these classes. All classes defined using reference classes. 
#
# Andrew Roberts
#
# Depends on: general_helper_functions.r
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
# surrogate/emulator model, which is stored in the field 
# `emulator_model` of the class. Currently, this model is almost 
# always a Gaussian process (GP), in which case `emulator_model`
# is a `gpWrapper` object. 
# `lik_par` is stores the likelihood parameters. These 
# can be stored in advance if the likelihood parameters will be fixed
# throughout an analysis; otherwise, this field can be left NULL. Note 
# that the likelihood approximation may be a function of the likelihood 
# parameters (e.g. the likelihood parameters are included as inputs to 
# a GP emulator). 
#
# This class organizes the pieces comprising the likelihood into 
# (1) "inputs" or "par": the primary parameters of interest; e.g. the parameters 
# characterizing the forward model in an inverse problem; and (2) "lik_par", 
# other likelihood parameters (e.g. variance parameters) that define 
# the observation model, may be unknown, but are not of primary interest. 
#
# The core methods for the `llikEmulator` class are `get_lik_par()`, 
# `get_input()`, `assemble_llik()`, `sample()`, and `predict()`. 
# 
# Implementing a deterministic approximation: 
#    TODO
# Implementing an exact llik:
#    TODO
# -----------------------------------------------------------------------------

llikEmulator <- setRefClass(
  Class = "llikEmulator", 
  fields = list(emulator_model="ANY", lik_description="character", llik_label="character",
                emulator_description="character", default_conditional="logical",
                default_normalize="logical", use_fixed_lik_par="logical", lik_par="ANY",
                input_names="character", dim_input="integer", llik_pred_dist="character", 
                exact_llik="logical")
)

llikEmulator$lock("llik_label")
llikEmulator$lock("llik_pred_dist")
llikEmulator$lock("exact_llik")

llikEmulator$methods(
  
  initialize = function(llik_label, input_names, lik_description, emulator_description, dim_input,  
                        emulator_model=NULL, default_conditional=FALSE, default_normalize=FALSE,
                        use_fixed_lik_par=FALSE, lik_par=NULL, llik_pred_dist="unspecified", 
                        exact_llik=FALSE, ...) {

    initFields(llik_label=llik_label, input_names=input_names, lik_description=lik_description, 
               dim_input=dim_input, emulator_description=emulator_description,
               emulator_model=emulator_model, default_conditional=default_conditional,
               default_normalize=default_normalize, use_fixed_lik_par=use_fixed_lik_par, 
               lik_par=lik_par, llik_pred_dist=llik_pred_dist, exact_llik=exact_llik)  
        
    if(use_fixed_lik_par && is.null(lik_par)) stop("Fixed `lik_par` not passed but `use_fixed_lik_par` is TRUE.")
  }, 
  
  get_lik_par = function(lik_par_val=NULL, ...) {
    if(use_fixed_lik_par) return(lik_par)
    assert_that(!is.null(lik_par_val), 
                msg="`lik_par_val` arg must be non-NULL if `use_fixed_lik_par` is FALSE.")
    return(lik_par_val)
  },
  
  get_llik_term_attr = function(attr_name, ...) {
    setNames(.self$field(attr_name), llik_label)
  },
  
  get_input = function(input, ...) {
    assert_that(is.matrix(input) && (ncol(input)==dim_input), 
                msg="`input` must be matrix with ncol equal to `dim_input`.")
    assert_that(!is.null(colnames(input)) && all(is.element(colnames(input), input_names)),
                msg="`input` must have colnames set to subset of `input_names`.")

    return(input[,input_names, drop=FALSE])
  },
  
  sample_emulator = function(input, N_samp=1, ...) {
    .NotYetImplemented()
  },
  
  predict_emulator = function(input, lik_par_val=NULL, emulator_pred_list=NULL, return_mean=TRUE, 
                              return_var=TRUE, return_cov=FALSE, return_cross_cov=FALSE, input_cross=NULL, ...) {
    # This default method assumes the emulator inherits from `gpWrapper`, but different llikEmulator 
    # classes may want to override this default method. 
    
    assert_that(inherits(.self$emulator_model, "gpWrapper"))
    if(!is.null(emulator_pred_list)) return(emulator_pred_list)
    
    .self$emulator_model$predict(.self$get_input(input), return_mean=return_mean, 
                                 return_var=return_var, return_cov=return_cov, 
                                 return_cross_cov=return_cross_cov, input_cross=input_cross)
  },
  
  assemble_llik = function(input, lik_par_val=NULL, conditional=default_conditional, 
                           normalize=default_normalize, ...) {
    # Since a log-likelihood is scalar-valued, this should always return a numeric vector
    # of length equal to the number of rows in `input`. 
    .NotYetImplemented()
  },
  
  sample = function(input, lik_par_val=NULL, N_samp=1, conditional=default_conditional, 
                    normalize=default_normalize, ...) {
    .NotYetImplemented()
  },
  
  predict = function(input, lik_par_val=NULL, return_mean=TRUE, return_var=TRUE, 
                     return_cov=FALSE, return_cross_cov=FALSE, input_cross=NULL,
                     conditional=default_conditional, normalize=default_normalize, ...) {
    .NotYetImplemented()
  },
  
  predict_lik = function(input, lik_par_val=NULL, emulator_pred_list=NULL, return_mean=TRUE,  
                         return_var=TRUE, return_cov=FALSE, return_cross_cov=FALSE, 
                         input_cross=NULL, conditional=default_conditional,  
                         normalize=default_normalize, log_scale=FALSE, ...) {
    .NotYetImplemented()
  },
  
  calc_quantiles = function(p, input=NULL, lik_par_val=NULL, conditional=default_conditional, 
                            normalize=default_normalize, llik_pred_list=NULL,
                            lower_tail=TRUE, ...) {
    .NotYetImplemented()
  },
  
  calc_lik_approx = function(input, approx_type, lik_par_val=NULL, emulator_pred_list=NULL,
                             conditional=default_conditional, normalize=default_normalize,
                             log_scale=FALSE, ...) {
    # A helper method that will call the method `calc_lik_<approx_type>_approx()`, which 
    # returns an approximation of the likelihood function evaluated at the inputs `input`, 
    # using the likelihood approximation type `approx_type`. 
    # See description on usingMethods here: 
    # https://docs.tibco.com/pub/enterprise-runtime-for-R/6.0.0/doc/html/Language_Reference/methods/setRefClass.html
    
    # Note that for some reason this doesn't work if called like `.self$usingMethods()`. 
    usingMethods(calc_lik_marginal_approx, calc_lik_mean_approx, calc_lik_sample_approx)
    
    if(is.null(emulator_pred_list)) {
      emulator_pred_list <- .self$predict_emulator(input, lik_par_val=lik_par_val, ...)
    }

    lik_approx_method_name <- paste0("calc_lik_", approx_type, "_approx")
    get(lik_approx_method_name)(input=input, lik_par_val=lik_par_val, emulator_pred_list=emulator_pred_list, 
                                conditional=conditional, normalize=normalize, log_scale=log_scale, ...)
  },
  
  calc_lik_approx_comparison = function(input, approx_types, lik_par_val=NULL, emulator_pred_list=NULL,
                                        conditional=default_conditional, normalize=default_normalize,
                                        include_nugget=TRUE, log_scale=FALSE, ...) {
    # A wrapper around `calc_lik_approx` that allows multiple likelihood approximation types 
    # to be computed. Returns a matrix of dimension `(nrow(input), length(approx_types))`, 
    # each column containing a different likelihood approximation type. 
    
    if(is.null(emulator_pred_list)) {
      emulator_pred_list <- .self$predict_emulator(input, lik_par_val=lik_par_val, ...)
    }
    
    lik_approx_vals <- matrix(nrow=nrow(input), ncol=length(approx_types),
                              dimnames=list(NULL, approx_types))
    for(i in seq_along(approx_types)) {
      lik_approx_vals[,i] <-   calc_lik_approx(input=input, approx_type=approx_types[i], lik_par_val=lik_par_val, 
                                               emulator_pred_list=emulator_pred_list, conditional=conditional, 
                                               normalize=normalize, include_nugget=TRUE, log_scale=log_scale, ...)
    }
    
    return(lik_approx_vals)
  },  
  
  calc_lik_marginal_approx = function(input, lik_par_val=NULL, emulator_pred_list=NULL, conditional=default_conditional,
                                      normalize=default_normalize, include_nugget=TRUE, log_scale=FALSE, ...) {
    .NotYetImplemented()  
  },
  
  calc_lik_sample_approx = function(input, lik_par_val=NULL, emulator_pred_list=NULL, conditional=default_conditional,
                                    normalize=default_normalize, include_nugget=TRUE, log_scale=FALSE, ...) {
    .NotYetImplemented()  
  },
  
  calc_lik_quantile_approx = function(input, alpha=0.9, lik_par_val=NULL, emulator_pred_list=NULL, 
                                      conditional=default_conditional, normalize=default_normalize,
                                      include_nugget=TRUE, log_scale=FALSE, ...) {
    .NotYetImplemented()
  },
  
  calc_lik_mean_approx = function(input, lik_par_val=NULL, emulator_pred_list=NULL, conditional=default_conditional,
                                  normalize=default_normalize, include_nugget=TRUE, log_scale=FALSE, ...) {
    # `Mean approx` means that the `emulator_model` predictive mean is computed at inputs `input`, then 
    # the predictive mean is passed to `assemble_llik()` and the result is exponentiated (if 
    # `log_scale` is FALSE). This implements the "plug-in emulator mean" approximation. This default method 
    # can be written for llikEmulator classes with special structure that differs from this. 
    
    input <- .self$get_input(input)
    
    if(is.null(emulator_pred_list)) {
      emulator_pred_list <- .self$predict_emulator(input, lik_par_val=lik_par_val, return_var=FALSE, 
                                                   include_nugget=include_nugget, ...)
    } else {
      assert_that(!is.null(emulator_pred_list$mean))
    }
      
    llik_pred <- .self$assemble_llik(emulator_pred_list$mean, lik_par_val=lik_par_val, 
                                     conditional=conditional, normalize=normalize, ...)
    if(log_scale) return(llik_pred)
    return(exp(llik_pred))
  },
  
  calc_expected_lik_cond_var = function(input_eval, input_cond, include_nugget=TRUE, log_scale=TRUE, plugin=FALSE, ...) {
    # Under the default setting `plugin = FALSE`, computes the log of 
    # E_l Var[exp(L(input_eval))|L(input_cond)=l]
    # where L is the log-likelihood emulator and the expectation is with respect to the current 
    # log-likelihood emulator predictive distribution (i.e., this is a posterior predictive quantity). 
    # If `plugin = TRUE`, then the function will instead compute 
    # E_g Var[exp(L(input_eval; G|{input_cond,g}))|g]
    # where `G` denotes the underlying emulator `emulator_model`. In words, the underlying emulator 
    # will first be conditioned at the new inputs `input_cond`- G|{input_cond,g} - where `g` is assumed to 
    # be distributed according to the current emulator model predictive distribution. This conditioned
    # emulator is plugged into the log-likelihood and then the expectation of the exponentiated 
    # log-likelihood is computed with respect to `g`. Certain classes may implement this method 
    # with only one of the two options for `plugin`, in which case an error should be thrown if 
    # the other option is passed. Other classes may not have this method at all, or return an 
    # approximation to this expected conditional variance.
    #
    # Returns vector of length `length(X_eval)` containing the expected conditional variance 
    # at the evaluation points `input_eval`. 
    .NotYetImplemented()
  },
  
  get_pred_interval = function(input, lik_par_val=NULL, emulator_pred_list=NULL, target_pred_list=NULL, target="llik",  
                               method="pm_std_dev", N_std_dev=1, CI_prob=0.9, 
                               conditional=default_conditional, normalize=default_normalize, include_nugget=TRUE, ...) {
    # Options for `target`: "llik" and "lik". 
    # Options for `method`: "pm_std_dev" (pm = "plus-minus"), "CI". The former computes the bounds by 
    # adding and subtracting `N_std_dev` standard deviations from the predictive mean. The latter 
    # computes a 100*`CI_prob`% confidence interval.
    # `target_pred_list` is either the pred list returned by `predict` or `predict_lik`; which one it is 
    # should align with the value of `target`. 
    
    assert_that(target %in% c("llik", "lik"))
    assert_that(method %in% c("pm_std_dev", "CI"))
    interval_list <- list()
    
    # Log likelihood or likelihood predictions. 
    if(is.null(target_pred_list)) {
      if(target == "llik") {
        target_pred_list <- .self$predict(input, lik_par_val=lik_par_val, emulator_pred_list=emulator_pred_list, 
                                          return_mean=TRUE, return_var=TRUE, conditional=conditional,
                                          normalize=normalize, log_scale=FALSE, ...)
      } else {
        target_pred_list <- .self$predict_lik(input, lik_par_val=lik_par_val, emulator_pred_list=emulator_pred_list, 
                                              return_mean=TRUE, return_var=TRUE, conditional=conditional,
                                              normalize=normalize, ...)
      }
    } else {
      assert_that(!is.null(target_pred_list$mean) && !is.null(target_pred_list$var))
    }
    
    # Plus/minus standard deviation method 
    if(method == "pm_std_dev") {
      interval_list$lower <- target_pred_list$mean - N_std_dev * sqrt(target_pred_list$var)
      interval_list$upper <- target_pred_list$mean + N_std_dev * sqrt(target_pred_list$var)
    } else {
      CI_tail_prob <- 0.5 * (1-CI_prob)
      interval_list$lower <- .self$calc_quantiles(p=CI_tail_prob, input=input, lik_par_val=lik_par_val, 
                                                  emulator_pred_list=emulator_pred_list, 
                                                  target=target, conditional=conditional, normalize=normalize, 
                                                  lower_tail=FALSE, include_nugget=include_nugget, ...)
      interval_list$upper <- .self$calc_quantiles(p=CI_tail_prob, input=input, lik_par_val=lik_par_val, 
                                                  emulator_pred_list=emulator_pred_list, 
                                                  target=target, conditional=conditional, normalize=normalize, 
                                                  lower_tail=TRUE, include_nugget=include_nugget, ...)
    }
    
    return(interval_list)
    
  },
  
  get_design_inputs = function(...) {
    .NotYetImplemented()
  },
  
  
  get_design_llik = function(lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, ...) {
    .NotYetImplemented()
  },
  
  
  plot_samp_1d = function(input, lik_par_val=NULL, emulator_pred_list=NULL, N_samp=1, plot_type="llik",  
                          conditional=default_conditional, normalize=default_normalize, 
                          true_llik=NULL, include_design=TRUE, use_cov=TRUE, ground_truth_col="black",
                          design_col="black", ...) {
    # `plot_type` options: "llik", "lik". 
    
    assert_that(dim_input==1, msg=paste0("plot_llik_samp_1d() requires 1d input space. input_dim = ", dim_input))
    assert_that(plot_type %in% c("llik", "lik"))
    
    input <- get_input(input)
    samp <- .self$sample(input, lik_par=lik_par_val, emulator_pred_list=emulator_pred_list,
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
      plt <- plt + geom_line(aes(x=x, y=y), df, inherit.aes=FALSE, color=ground_truth_col)
    }
    
    if(include_design) {
      design_response_vals <- drop(get_design_llik(lik_par_val, conditional, normalize))
      if(plot_type == "lik") design_response_vals <- exp(design_response_vals)
      
      design_df <- data.frame(x=drop(get_design_inputs()), y=design_response_vals)
      plt <- plt + geom_point(aes(x=x, y=y), design_df, inherit.aes=FALSE, color=design_col)
    }

    return(plt)
    
  },
  
  plot_pred_1d = function(input, lik_par_val=NULL, emulator_pred_list=NULL, plot_type="llik", 
                          conditional=default_conditional, normalize=default_normalize,
                          include_interval=TRUE, interval_method="pm_std_dev",
                          N_std_dev=1, CI_prob=0.9, true_llik=NULL, 
                          include_design=TRUE, xlab=NULL, ylab=NULL, plot_title=NULL, ...) {
    
    assert_that(dim_input==1, msg=paste0("plot_llik_pred_1d() requires 1d input space. dim_input = ", dim_input))
    assert_that(plot_type %in% c("llik", "lik"))
    
    input <- .self$get_input(input)
    
    # Default axis labels. 
    if(is.null(xlab)) xlab <- input_names
    if(is.null(ylab)) ylab <- plot_type
    
    # Compute required predictive quantities if not already provided. 
    true_vals <- true_llik
    if(plot_type == "llik") {
      pred_list <- .self$predict(input, lik_par_val=lik_par_val, emulator_pred_list=emulator_pred_list,  
                                 return_mean=TRUE, return_var=include_interval, 
                                 conditional=conditional, normalize=normalize, ...)
    } else {
      pred_list <- .self$predict_lik(input, lik_par_val=lik_par_val, emulator_pred_list=emulator_pred_list,  
                                     return_mean=TRUE, return_var=include_interval, conditional=conditional, 
                                     normalize=normalize, log_scale=FALSE, ...)
      if(!is.null(true_vals)) true_vals <- exp(true_vals)
    }
    
    # Plot title and labels.
    if(is.null(plot_title)) {
      plot_title <- "Likelihood Emulator Predictions"
      if(plot_type == "llik") plot_title <- paste0("Log ", plot_title)
      if(include_interval && (interval_method=="CI")) plot_title <- paste0(plot_title, ", ", 100*CI_prob, "% CI")
      if(include_interval && (interval_method=="pm_std_dev")) plot_title <- paste0(plot_title, ", +/- ", N_std_dev, " std dev")
    }
    
    if(!normalize) {
      ylab <- paste0(ylab, ", ", ifelse(conditional, "unnormalized/conditional", "unnormalized"))
    }

    # Compute prediction interval.  
    if(include_interval) interval_list <- .self$get_pred_interval(input, lik_par_val=lik_par_val, target_pred_list=pred_list, 
                                                                  target=plot_type, method=interval_method, N_std_dev=N_std_dev,
                                                                  CI_prob=CI_prob, conditional=conditional, normalize=normalize, ...)
    else interval_list <- NULL

    # Design points. 
    design_inputs <- NULL
    design_response_vals <- NULL 
    if(include_design) {
      design_inputs <- drop(.self$get_design_inputs(...))
      design_response_vals <- drop(get_design_llik(lik_par_val, conditional, normalize, ...))
      if(plot_type == "lik") design_response_vals <- exp(design_response_vals)
    } 
    
    # Produce plot. 
    plt <- plot_pred_1d_helper(X_new=drop(input), pred_mean=drop(pred_list$mean), 
                               include_CI=include_interval, CI_lower=interval_list$lower,  
                               CI_upper=interval_list$upper, y_new=drop(true_vals), 
                               X_design=design_inputs, y_design=design_response_vals, 
                               plot_title=plot_title, xlab=xlab, ylab=ylab, ...)
    
    return(plt)                     
  }, 
  
  
  plot_pred_validation = function(input, true_llik, lik_par_val=NULL, emulator_pred_list=NULL, 
                                  plot_type="llik", conditional=default_conditional, 
                                  normalize=default_normalize, include_interval=TRUE,
                                  interval_method="pm_std_dev", N_std_dev=1, CI_prob=0.9, 
                                  include_design=TRUE, xlab=NULL, ylab=NULL, plot_title=NULL, ...) {
    
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
      pred_list <- .self$predict(input, lik_par_val=lik_par_val, emulator_pred_list=emulator_pred_list,  
                                 return_mean=TRUE, return_var=include_interval, 
                                 conditional=conditional, normalize=normalize, ...)
    } else {
      pred_list <- .self$predict_lik(input, lik_par_val=lik_par_val, emulator_pred_list=emulator_pred_list,  
                                     return_mean=TRUE, return_var=include_interval, conditional=conditional, 
                                     normalize=normalize, log_scale=FALSE, ...)
      if(!is.null(true_vals)) true_vals <- exp(true_vals)
    }
    
    # Compute prediction interval.  
    if(include_interval) interval_list <- .self$get_pred_interval(input, lik_par_val=lik_par_val, target_pred_list=pred_list, 
                                                                  target=plot_type, method=interval_method, N_std_dev=N_std_dev,
                                                                  CI_prob=CI_prob, conditional=conditional, normalize=normalize, ...)
    else interval_list <- NULL
    
    # Design points. 
    design_inputs <- NULL
    design_response_vals <- NULL 
    if(include_design) {
      design_response_vals <- drop(get_design_llik(lik_par_val, conditional, normalize, ...))
      if(plot_type == "lik") design_response_vals <- exp(design_response_vals)
    } 
    
    # Produce plot. 
    plt <- plot_true_pred_scatter(y_pred=drop(pred_list$mean), y_true=drop(true_vals),
                                  include_CI=include_interval, CI_lower=interval_list$lower,  
                                  CI_upper=interval_list$upper, y_design=design_response_vals, 
                                  plot_title=plot_title, xlab=xlab, ylab=ylab, ...)
    
    return(plt)
  }
  
)


# -----------------------------------------------------------------------------
# llikSumEmulator Class 
# -----------------------------------------------------------------------------

#
# TODO: this needs lots of updates
#    - Add argument `emulator_pred_list` to the relevant functions. 
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
              use_fixed_lik_par=all(sapply(llik_emulator_list, function(x) x$use_fixed_lik_par)), 
              lik_par=NULL, llik_pred_dist=llik_sum_dist, exact_llik=llik_sum_exact)
              
  },
  
  get_llik_term_attr = function(attr_name, labels=llik_label) {
    
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
# for reference and validation purposes. Similar to 
# `normalize` and `conditional`, an error is thrown 
# if the user tries to pass a likelihood parameter that differs from the 
# fixed value set when initializing the class.
# -----------------------------------------------------------------------------

llikEmulatorGP <- setRefClass(
  Class = "llikEmulatorGP", 
  contains = "llikEmulator",
  fields = list(N_obs="integer")
)

llikEmulatorGP$methods(
  
  initialize = function(llik_lbl, gp_model, default_conditional, default_normalize,
                        use_fixed_lik_par=TRUE, lik_par=NULL, ...) {
                        
    assert_that(use_fixed_lik_par, msg="llikEmulatorGP does not yet support emulation as a function of `lik_par`.")
    assert_that(inherits(gp_model, "gpWrapper"), msg="`gp_model` must inherit from `gpWrapper` class.")
    assert_that(gp_model$Y_dim==1, msg="`llikEmulatorMultGausGP` only supports single-output GP emulator.")
    assert_that(!is.null(gp_model$X_names) && noNA(gp_model$X_names), 
                msg="`llikEmulatorMultGausGP` requires that `gp_model` has `X_names` field set.")
    
    callSuper(emulator_model=gp_model, llik_label=llik_lbl, lik_par=lik_par, input_names=gp_model$X_names,
              dim_input=gp_model$X_dim, default_conditional=default_conditional, 
              default_normalize=default_normalize, use_fixed_lik_par=use_fixed_lik_par, 
              lik_description="Generic log likelihood",
              emulator_description="GP directly emulating the log-likelihood.", 
              llik_pred_dist="Gaussian", exact_llik=FALSE, ...)
  }, 
  
  check_fixed_quantities = function(conditional=NULL, normalize=NULL, lik_par_val=NULL) {
    if(!is.null(conditional)) {
      assert_that(conditional==default_conditional, 
                  msg="`llikEmulatorGP` class requires `conditional` to agree with `default_conditional`.")
    }
    
    if(!is.null(normalize)) {
      assert_that(normalize==default_normalize, 
                  msg="`llikEmulatorGP` class requires `normalize` to agree with `default_normalize`.")
    }
    
    if(!is.null(lik_par_val)) {
      assert_that(lik_par_val==get_lik_par(lik_par), 
                  msg="`llikEmulatorGP` class requires `lik_par_val` to agree with `lik_par`.")
    }
  },
  
  assemble_llik = function(llik, lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, ...) {
    # llik should be N_input x 1. Since the llik is emulated directly, then this function simply returns 
    # the argument `llik` after performing argument validation. 
    .self$check_fixed_quantities(conditional, normalize, lik_par_val)
    return(drop(llik))
  },
  
  get_design_inputs = function(...) {
    emulator_model$X
  }, 
  
  get_design_llik = function(lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, ...) {
    # Returns the response values in the design of the GP, since the response is the llik in this case. 
    .self$check_fixed_quantities(conditional, normalize, lik_par_val)
    return(emulator_model$Y)
  },
  
  sample_emulator = function(input, emulator_pred_list=NULL, N_samp=1, use_cov=FALSE, include_nugget=TRUE, ...) {
    emulator_model$sample(get_input(input), use_cov=use_cov, include_nugget=include_nugget, 
                          N_samp=N_samp, pred_list=emulator_pred_list, ...)[,,1,drop=FALSE]             
  },
  
  sample = function(input, lik_par_val=NULL, emulator_pred_list=NULL, N_samp=1, use_cov=FALSE, 
                    include_nugget=TRUE, conditional=default_conditional, normalize=default_normalize, ...) {
    # Directly returns the emulator samples, since these are llik samples. 
    .self$check_fixed_quantities(conditional, normalize, lik_par_val)
    sample_emulator(input, emulator_pred_list, N_samp, use_cov, include_nugget, ...)[,,1]
  }, 
  
  predict = function(input, lik_par_val=NULL, emulator_pred_list=NULL, return_mean=TRUE, 
                     return_var=TRUE, return_cov=FALSE, return_cross_cov=FALSE, 
                     input_cross=NULL, conditional=default_conditional, 
                     normalize=default_normalize, include_nugget=TRUE, ...) {
    # Log-likelihood emulator mean/var/cov predictions. Since the GP directly emulates 
    # the log-likelihood, then simply return the GP predictions directly. 
    
    .self$check_fixed_quantities(conditional, normalize, lik_par_val)
    
    if(!is.null(emulator_pred_list)) return(emulator_pred_list)
    
    .self$emulator_model$predict(get_input(input), return_mean=return_mean, return_var=return_var,
                                 return_cov=return_cov, return_cross_cov=return_cross_cov,
                                 X_cross=input_cross, include_nugget=include_nugget, ...)
  }, 
  
  predict_lik = function(input, lik_par_val=NULL, emulator_pred_list=NULL, return_mean=TRUE,  
                         return_var=TRUE, return_cov=FALSE, return_cross_cov=FALSE, 
                         input_cross=NULL, conditional=default_conditional, 
                         normalize=default_normalize, include_nugget=TRUE, log_scale=FALSE, ...) {
    # Likelihood emulator mean/var/cov predictions. For this class (under direct GP emulation
    # of the log-likelihood), the likelihood emulator is a log-normal process. Thus, 
    # the GP log-likelihood predictions can simply be transformed to obtain log-normal 
    # likelihood predictions. Currently `return_cross_cov` is not supported. Note that both the 
    # GP mean and variance are required to compute the log-normal mean and/or variance. 
    
    if(return_cross_cov) {
      stop("`return_cross_cov` is not yet supported for `llikEmulatorGP$predict_lik()`.")
    }
    
    llik_pred <- .self$predict(input=input, lik_par_val=lik_par_val, emulator_pred_list=emulator_pred_list, 
                               return_mean=TRUE, return_var=TRUE, return_cov=return_cov,  
                               return_cross_cov=return_cross_cov, input_cross=input_cross, 
                               conditional=conditional, normalize=normalize, include_nugget=include_nugget, ...)
                               
    convert_Gaussian_to_LN(mean_Gaussian=llik_pred$mean, var_Gaussian=llik_pred$var, 
                           cov_Gaussian=llik_pred$cov, return_mean=return_mean,
                           return_var=return_var, return_cov=return_cov, log_scale=log_scale)
  }, 
  

  calc_lik_marginal_approx = function(input, lik_par_val=NULL, emulator_pred_list=NULL, 
                                      conditional=default_conditional, normalize=default_normalize,
                                      include_nugget=TRUE, log_scale=FALSE, ...) {
    
    selector <- ifelse(log_scale, "log_mean", "mean")
    
    predict_lik(input, lik_par_val=lik_par_val, emulator_pred_list=emulator_pred_list, 
                return_mean=TRUE, return_var=FALSE, conditional=conditional, normalize=normalize,
                log_scale=log_scale, include_nugget=include_nugget, ...)[[selector]]
    
  },
  
  
  calc_lik_quantile_approx = function(input, alpha=0.9, lik_par_val=NULL, emulator_pred_list=NULL, 
                                      conditional=default_conditional, normalize=default_normalize,
                                      include_nugget=TRUE, log_scale=FALSE, ...) {
    
    input <- .self$get_input(input)
    
    if(is.null(emulator_pred_list)) {
      emulator_pred_list <- .self$predict_emulator(input, lik_par_val=lik_par_val, return_var=TRUE, ...)
    } else {
      assert_that(!is.null(emulator_pred_list$mean))
      assert_that(!is.null(emulator_pred_list$var))
    }
    
    q <- qnorm(alpha)
    emulator_pred_list$mean + q * sqrt(emulator_pred_list$var)
    
  },
  
  
  calc_quantiles = function(p, input=NULL, lik_par_val=NULL, emulator_pred_list=NULL, llik_pred_list=NULL,
                            target="llik", conditional=default_conditional, normalize=default_normalize, 
                            lower_tail=TRUE, include_nugget=TRUE, ...) {
    # The log-likelihood emulator distribution is Gaussian, and the likelihood emulator is log-normal.
    # Both the Gaussian and Log-normal quantile functions are parameterized in terms of the underlying 
    # Gaussian mean/variance, so only `llik_pred_list` is required, regardless of whether `target`
    # is "llik" or "lik". 
    
    assert_that(target %in% c("llik", "lik"))
    
    # Log-likelihood or likelihood predictions. 
    if(is.null(llik_pred_list)) {
      llik_pred_list <- .self$predict(input, lik_par_val=lik_par_val, emulator_pred_list=emulator_pred_list, 
                                      return_mean=TRUE, return_var=TRUE, conditional=conditional,
                                      normalize=normalize, log_scale=FALSE, ...)
    } else {
      assert_that(!is.null(llik_pred_list$mean) && !is.null(llik_pred_list$var))
    }
    
    # Log-likelihood emulator distribution is Gaussian. Likelihood emulator is log-normal. 
    if(target == "llik") q <- qnorm(p, drop(llik_pred_list$mean), sqrt(drop(llik_pred_list$var)), lower.tail=lower_tail)
    else q <- qlnorm(p, drop(llik_pred_list$mean), sqrt(drop(llik_pred_list$var)), lower.tail=lower_tail)
    
    return(q)
  }, 
  
  calc_expected_lik_cond_var = function(input_eval, input_cond, include_nugget=TRUE, 
                                        log_scale=TRUE, plugin=FALSE, ...) {
    # Since the emulator model directly emulates the log-likelihood in this class, the 
    # options `plugin = FALSE` and `plugin = TRUE` are equivalent in this case. In either
    # case, the analogous method for the GP `emulator_model` is dispatched. 
    
    .self$emulator_model$calc_expected_exp_cond_var(input_eval, input_cond, include_nugget=include_nugget, 
                                                    log_scale=log_scale, ...)
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
                        default_normalize=FALSE, use_fixed_lik_par=FALSE, ...) {
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
              default_normalize=default_normalize, use_fixed_lik_par=use_fixed_lik_par, 
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
  
  sample_emulator = function(input, emulator_pred_list=NULL, N_samp=1, use_cov=FALSE, 
                             include_nugget=TRUE, adjustment="rectified", ...) {
    emulator_model$sample(get_input(input), use_cov=use_cov, include_nugget=include_nugget, 
                          N_samp=N_samp, adjustment=adjustment, pred_list=emulator_pred_list)[,,1,drop=FALSE]             
  },
  
  sample = function(input, lik_par_val=NULL, emulator_pred_list=NULL, N_samp=1, use_cov=FALSE, 
                    include_nugget=TRUE, conditional=default_conditional, normalize=default_normalize, ...) {
    # Sample SSR. 
    samp <- sample_emulator(input, emulator_pred_list, N_samp, use_cov, include_nugget, ...)
    
    # Compute unnormalized or normalized log-likelihood. 
    assemble_llik(samp, lik_par_val, conditional, normalize)
  }, 
  
  predict = function(input, lik_par_val=NULL, emulator_pred_list=NULL, return_mean=TRUE, 
                     return_var=TRUE,  return_cov=FALSE, return_cross_cov=FALSE, 
                     input_cross=NULL, conditional=default_conditional, 
                     normalize=default_normalize, include_nugget=TRUE, ...) {
    
    if(is.null(emulator_pred_list)) {
      pred_list <- .self$emulator_model$predict(get_input(input), return_mean=return_mean, return_var=return_var,
                                                return_cov=return_cov, return_cross_cov=return_cross_cov,
                                                X_cross=input_cross, include_nugget=include_nugget, ...)
    } else {
      pred_list <- emulator_pred_list
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

  
  calc_quantiles = function(p, input=NULL, lik_par_val=NULL, emulator_pred_list=NULL, llik_pred_list=NULL,
                             target="llik", conditional=default_conditional, normalize=default_normalize, 
                             lower_tail=TRUE, include_nugget=TRUE, ...) {
    # The log-likelihood emulator distribution is Gaussian, and the likelihood emulator is log-normal. 
    # Both the Gaussian and Log-normal quantile functions are parameterized in terms of the underlying 
    # Gaussian mean/variance, so only `llik_pred_list` is required, regardless of whether `target`
    # is "llik" or "lik". 
    
    assert_that(target %in% c("llik", "lik"))
    
    # Log-likelihood or likelihood predictions. 
    if(is.null(llik_pred_list)) {
      llik_pred_list <- .self$predict(input, lik_par_val=lik_par_val, emulator_pred_list=emulator_pred_list, 
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
# This simply implements the exact likelihood corresponding to a  
# Gaussian inverse problem: y|u ~ N(G(u), Sig), where G may be nonlinear.
# The forward map G is defined by an attribute called `fwd_model` in the class, 
# which is a function that users pass in when instantiating the class. 
# This is exact in the sense that there is no emulation here -
# the reason for implementing this as a 
# llikEmulator class is to use it for algorithm testing; e.g. ensuring the 
# correctness of an MCMC implementation. The `lik_par` here is defined to 
# be the covariance matrix `Sig`. 
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
  fields = list(fwd_model="ANY", fwd_model_vectorized="ANY", y="numeric", N_obs="integer", L_Cov="matrix")
)

llikEmulatorExactGauss$methods(
  
  initialize = function(llik_lbl, y_obs, dim_par, fwd_model=NULL, fwd_model_vectorized=NULL, Cov=NULL, 
                        default_conditional=FALSE, default_normalize=FALSE, use_fixed_lik_par=FALSE, 
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
              default_normalize=default_normalize, use_fixed_lik_par=use_fixed_lik_par, 
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
    if(use_fixed_lik_par) {
      if(return_chol) return(L_Cov)
      else return(lik_par)
    }
      
    assert_that(!is.null(lik_par_val), 
                msg="`lik_par_val` arg must be non-NULL if `use_fixed_lik_par` is FALSE.")
    if(return_chol) return(t(chol(lik_par_val)))
    else return(lik_par_val)
  },
  
  assemble_llik = function(input, lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, ...) {
    # `input` should be N_input x N_samp. 
    
    # Fetch the lower triangular Cholesky factor of the covariance matrix.
    L <- get_lik_par(lik_par_val, return_chol=TRUE)
    
    # Construct log likelihood. 
    llik <- -0.5 * colSums(solve(L, y - run_fwd_model(input, ...))^2, na.rm=TRUE)
    if(normalize || !conditional) llik <- llik - sum(log(diag(L)))
    if(normalize) llik <- llik - 0.5*N_obs*log(2*pi)

    return(drop(llik))
  }, 
  
  sample_emulator = function(input, emulator_pred_list=NULL, N_samp=1, ...) {
    # No emulator to sample from, simply return input. The argument 
    # `emulator_pred_list` is only present for consistency with other llikEmulator
    # classes. 
    
    input
  },
  
  sample = function(input, lik_par=NULL, emulator_pred_list=NULL, N_samp=1, 
                    conditional=default_conditional, normalize=default_normalize, ...) {
    # Compute unnormalized or normalized log-likelihood (exact, deterministic 
    # calculation - no sampling is actually performed). For consistency with 
    # other classes, duplicates the exact likelihood calculation when `N_samp`>1.
    matrix(assemble_llik(get_input(input), lik_par, conditional, normalize), 
           nrow=nrow(input), ncol=N_samp)
  }
  
)


# -----------------------------------------------------------------------------
# llikEmulatorExactGaussDiag class.
#    Identical to llikEmulatorExactGauss except that `Sig` is assumed to 
#    be diagonal, so the `lik_par` for this class is defined to be the vector 
#    corresponding to the diagonal of this covariance matrix.`lik_par` may 
#    also be provided as a single value, which is interpreted as a 
#    homoskedastic variance. Concretely, implements a likelihood of the form, 
#    likelihood of the form
#        p(y1,...,yP|u,v1,...,vP) = prod_{p=1}^{P} prod_{n=1}^{N} N(yp_n|G(u)_p, vp),
#    where:
#    P = `N_output` and N = `N_obs`
#    y1,...,yP are each vectors of length N.
#    v1,...,vP are variance parameters corresponding to each output, respectively. 
#
# In words, this likelihood assumes the forward model `G` has P outputs. A set of 
# observations yp is associated with each output p, and these observations are 
# assumed to be iid draws from N(G(u)_p, vp). The argument `y_obs` to the class 
# constructor is a (N,P) matrix storing each set of observations as columns.  
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
  
  initialize = function(llik_lbl, y_obs, dim_par, fwd_model=NULL, fwd_model_vectorized=NULL, 
                        sig2=NULL, default_conditional=FALSE, default_normalize=FALSE, 
                        use_fixed_lik_par=FALSE, par_names=NULL, ...) {
                        
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
      assert_that(is.numeric(sig2) && ((length(sig2)==N_output) || (length(sig2)==1)) && all(sig2>0),
                  msg="`sig2` must be either vector of length `N_output` or 1 and only contain positive numbers.")
    }
    
    callSuper(emulator_model=NULL, llik_label=llik_lbl, lik_par=sig2, dim_input=dim_par,
              default_conditional=default_conditional, input_names=par_names,
              default_normalize=default_normalize, use_fixed_lik_par=use_fixed_lik_par, 
              lik_description="Exact Gaussian likelihood, diagonal covariance structure.",
              emulator_description="No emulation.", exact_llik=TRUE, ...)
  },
  
  run_fwd_model = function(input, ...) {
    # `input` is an (M,D) matrix with input parameter values stacked in the rows. 
    # Returns (M, N_output) dimensional output. 
    
    if(is.null(.self$fwd_model_vectorized)) return(vectorize_fwd_model(input, ...))
    return(.self$fwd_model_vectorized(input, ...))
  },
  
  vectorize_fwd_model = function(input, ...) {
    model_output <- matrix(nrow=nrow(input), ncol=N_output)
    for(i in 1:nrow(input)) model_output[i,] <- .self$fwd_model(input[i,], ...)
    return(model_output)
  },
  
  get_lik_par = function(lik_par_val=NULL, ...) {
    if(use_fixed_lik_par) lik_par_val <- lik_par
    
    assert_that(!is.null(lik_par_val), 
                msg="`lik_par_val` arg must be non-NULL if `use_fixed_lik_par` is FALSE.")
    
    if((length(lik_par_val)==1) && (.self$N_output > 1)) return(rep(lik_par_val, N_output))
    
    assert_that(length(lik_par_val)==N_output, msg="`lik_par_val` length not equal to 1 or `N_output`.")
    return(lik_par_val)
  },
  
  assemble_llik = function(input, lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, ...) {
    # `input` should have dimension (N_input, D). Returns vector of length `N_input.` 

    # Fetch the variance parameters. 
    sig2_val <- get_lik_par(lik_par_val)
    
    # Run forward model. 
    fwd_model_vals <- .self$run_fwd_model(input, ...)
    
    # Construct log likelihood.
    llik <- vector(mode="numeric", length=nrow(input))
    for(i in seq_along(llik)) {
      llik[i] <- -0.5 * sum(mult_vec_with_mat_rows(1/sig2_val, 
                            add_vec_to_mat_rows(-fwd_model_vals[i,], .self$y)^2))
    }
    
    if(normalize || !conditional) llik <- llik - 0.5 * .self$N_obs * sum(log(sig2_val))
    if(normalize) llik <- llik - 0.5*.self$N_obs * .self$N_output * log(2*pi)
    
    return(llik)
  }, 
  
  sample_emulator = function(input, emulator_pred_list=NULL, N_samp=1, ...) {
    # No emulator to sample from, simply return input. The argument 
    # `emulator_pred_list` is only present for consistency with other llikEmulator
    # classes. 
    
    input
  },
  
  sample = function(input, lik_par=NULL, emulator_pred_list=NULL, N_samp=1, 
                    conditional=default_conditional, normalize=default_normalize, ...) {
    
    # Compute unnormalized or normalized log-likelihood (exact, deterministic 
    # calculation - no sampling is actually performed). For consistency with 
    # other classes, duplicates the exact likelihood calculation when `N_samp`>1.
    matrix(assemble_llik(get_input(input), lik_par, conditional, normalize), 
           nrow=nrow(input), ncol=N_samp)
  }
  
)


# -----------------------------------------------------------------------------
# llikEmulatorFwdGauss class.
#
# Implements a surrogate likelihood for the Gaussian inverse problem 
# y|u ~ N(G(u), Sig) where the forward model G is replaced by a 
# Gaussian process (GP) emulator. The forward model may have multiple 
# outputs, in which case independent GP emulators are fit for 
# each output. The attribute `N_output` is the number of outputs (i.e. the 
# dimension of the vector returned by G). This class also allows for 
# multiple independent observations, in which case the likelihood takes the form 
# prod_{n=1}^{N_obs} N(y_n | G(u), Sig), where `N_obs` is the number of 
# independent observations. Note that G(u) and Sig are assumed constant over 
# all observations. The `lik_par` for this class is defined to be the 
# covariance matrix `Sig`. The attribute `y` is the `N_obs` x `N_output`
# data matrix. If `y` is a numeric vector, then it is assumed that 
# `N_obs = 1` and the vector is reshaped into a 1 x `N_output` matrix. 
# -----------------------------------------------------------------------------

llikEmulatorFwdGauss <- setRefClass(
  Class = "llikEmulatorFwdGauss", 
  contains = "llikEmulator",
  fields = list(y="ANY", N_output="integer", N_obs="integer", L_cov="matrix")
                
)

llikEmulatorFwdGauss$methods(
  
  initialize = function(llik_lbl, gp_model, y_obs, Cov=NULL, default_conditional=FALSE, 
                        default_normalize=FALSE, use_fixed_lik_par=FALSE, par_names=NULL, ...) {
    
    assert_that(inherits(gp_model, "gpWrapper"), msg="`gp_model` must inherit from `gpWrapper` class.")
    assert_that(is.numeric(y_obs) || is.matrix(y_obs))
    if(is.numeric(y_obs)) y_obs <- matrix(y_obs, nrow=1)
    assert_that(ncol(y_obs) == gp_model$Y_dim)
    initFields(y=y_obs, N_output=ncol(y_obs), N_obs=nrow(y_obs))
    
    # Set parameter names. 
    dim_par <- gp_model$X_dim
    if(is.null(par_names)) par_names <- paste0("input", 1:dim_par)
    assert_that(length(par_names) == dim_par)
    
    # Covariance matrix of Gaussian likelihood. 
    if(!is.null(Cov)) {
      assert_that(is.matrix(Cov) && (nrow(Cov)==N_output) && (ncol(Cov)==N_output),
                  msg="`Cov` must be a positive definite matrix with dim `N_output` x `N_output`")
      initFields(L_Cov=t(chol(Cov))) 
    }
    
    callSuper(emulator_model=gp_model, llik_label=llik_lbl, lik_par=Cov, dim_input=dim_par,
              default_conditional=default_conditional, input_names=par_names,
              default_normalize=default_normalize, use_fixed_lik_par=use_fixed_lik_par, 
              lik_description="Gaussian likelihood",
              emulator_description="Forward model GP emulator", exact_llik=FALSE, ...)
              
  },
  
  get_lik_par = function(lik_par_val=NULL, return_chol=FALSE, ...) {
    if(use_fixed_lik_par) {
      if(return_chol) return(L_Cov)
      else return(lik_par)
    }
    
    assert_that(!is.null(lik_par_val), 
                msg="`lik_par_val` arg must be non-NULL if `use_fixed_lik_par` is FALSE.")
    if(return_chol) return(t(chol(lik_par_val)))
    else return(lik_par_val)
  },
  
  assemble_llik = function(fwd_model_vals, lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, ...) {
    # `fwd_model_vals` should be of dimension `N_inputs` x `N_output`.
    
    # Fetch the lower triangular Cholesky factor of the covariance matrix.
    L <- get_lik_par(lik_par_val, return_chol=TRUE)
    
    # Construct log likelihood.
    llik <- rep(0, nrow(fwd_model_vals))
    for(n in seq_along(.self$N_obs)) {
      llik <- llik + colSums(solve(L, add_vec_to_mat_cols(.self$y[n,], -t(fwd_model_vals)))^2)
    }
    if(normalize || !conditional) llik <- llik - .self$N_obs * sum(log(diag(L)))
    if(normalize) llik <- llik - 0.5 * .self$N_obs * .self$N_output * log(2*pi)

    return(llik)
  }, 
  
  sample_emulator = function(input, emulator_pred_list=NULL, N_samp=1, use_cov=FALSE, ...) {
    # Sample the forward model emulator at specified inputs. `input` is M x D 
    # (M input vectors). Returns array of dimension (M, N_samp, N_output). 
    
    emulator_model$sample(get_input(input), use_cov=use_cov, include_nugget=include_nugget, 
                          N_samp=N_samp, pred_list=emulator_pred_list, ...)
  },
  
  sample = function(input, lik_par_val=NULL, emulator_pred_list=NULL, N_samp=1, use_cov=FALSE, 
                    conditional=default_conditional, normalize=default_normalize, ...) {
    # Sample the log-likelihood emulator at specified inputs. `input` is M x D 
    # (M input vectors). Returns array of dimension (M, N_samp). 
    
    input <- get_input(input)
    fwd_model_samp <- .self$sample_emulator(input, emulator_pred_list, N_samp=N_samp, 
                                            use_cov=use_cov, ...)
    llik_samp <- matrix(nrow=nrow(input), ncol=N_samp)
    
    for(i in 1:N_samp) {
      llik_samp[,i] <- assemble_llik(matrix(fwd_model_samp[,i,], nrow=nrow(input), ncol=N_output), 
                                     lik_par_val, conditional, normalize, ...)
    }
    
    return(llik_samp)
  }
  
)


# -----------------------------------------------------------------------------
# llikEmulatorFwdGaussDiag class.
#
# A special case of `llikEmulatorFwdGauss` where the covariance matrix is 
# assumed to be diagonal. The `lik_par` here is defined to be the 
# vector of variances `sig2`. If `sig2` is passed as a scalar value (but 
# `N_output` is greater than 1), then the same variance is used for all outputs. 
# -----------------------------------------------------------------------------

llikEmulatorFwdGaussDiag <- setRefClass(
  Class = "llikEmulatorFwdGaussDiag", 
  contains = "llikEmulator",
  fields = list(y="matrix", N_output="integer", N_obs="integer")
  
)

llikEmulatorFwdGaussDiag$methods(
  
  initialize = function(llik_lbl, gp_model, y_obs, sig2=NULL, default_conditional=FALSE, 
                        default_normalize=FALSE, use_fixed_lik_par=FALSE, par_names=NULL, ...) {
    
    assert_that(inherits(gp_model, "gpWrapper"), msg="`gp_model` must inherit from `gpWrapper` class.")
    assert_that(is.matrix(y_obs))
    assert_that(ncol(y_obs) == gp_model$Y_dim)
    initFields(y=y_obs, N_output=ncol(y_obs), N_obs=nrow(y_obs))
    
    # Set parameter names. 
    dim_par <- gp_model$X_dim
    if(is.null(par_names)) par_names <- paste0("input", 1:dim_par)
    assert_that(length(par_names) == dim_par)
    
    # Variance parameters for Gaussian likelihood. 
    if(!is.null(sig2)) {
      assert_that(is.numeric(sig2) && ((length(sig2)==N_output) || (length(sig2)==1)) && all(sig2>0),
                  msg="`sig2` must be either vector of length `N_output` or 1 and only contain positive numbers.")
    }
    
    callSuper(emulator_model=gp_model, llik_label=llik_lbl, lik_par=sig2, dim_input=dim_par,
              default_conditional=default_conditional, input_names=par_names,
              default_normalize=default_normalize, use_fixed_lik_par=use_fixed_lik_par, 
              lik_description="Gaussian likelihood, diagonal covariance.",
              emulator_description="Forward model GP emulator", exact_llik=FALSE, ...)
  },
  
  get_lik_par = function(lik_par_val=NULL, ...) {
    if(use_fixed_lik_par) lik_par_val <- lik_par
    
    assert_that(!is.null(lik_par_val), 
                msg="`lik_par_val` arg must be non-NULL if `use_fixed_lik_par` is FALSE.")
    
    if((length(lik_par_val)==1) && (.self$N_output > 1)) return(rep(lik_par_val, N_output))
    
    assert_that(length(lik_par_val)==N_output, msg="`lik_par_val` length not equal to 1 or `N_output`.")
    return(lik_par_val)
  },
  
  assemble_llik = function(fwd_model_vals, lik_par_val=NULL, conditional=default_conditional, 
                           normalize=default_normalize, var_inflation_vals=NULL, ...) {
    # `fwd_model_vals` should be of dimension `M` x `N_output`, where `M` corresponds to the number of 
    # instances of forward model output; e.g., might be the forward model evaluated at `M` different inputs, 
    # or `M` samples from the forward model at the same input, or some combination. Returns numeric vector of length `M`. 
    # Note that `sample_emulator` returns an array of dimension `M` x `N_samp` x `N_output`. Therefore, this 
    # function also allows `fwd_model_vals` to have three dimensions, as long as the second dimension 
    # equals 1 (as is the case with a default call to `sample_emulator`). 
    # `var_inflation_vals` allows adds values to `sig2`, and is used in computing the marginal likelihood 
    # approximation. If shape (`M`, `N_output`) then each row of `var_inflation_vals` to `sig2` for the 
    # corresponding row of `fwd_model_vals`. If a numeric vector of length `N_output`, then the same variance 
    # inflation is applied to the likelihood evaluations at all `M` values of the forward model output. 

    # If `fwd_model_vals` has three dimensions, with second dimension equal to 1, then collapse 
    # to a matrix along this dimension. 
    if(length(dim(fwd_model_vals)) == 3L) {
      if(dim(fwd_model_vals)[2] != 1L) stop("`fwd_model_vals` does not have correct dimensions.")
      fwd_model_vals <- matrix(fwd_model_vals[,1,], nrow=dim(fwd_model_vals)[1], ncol=dim(fwd_model_vals)[3])
    }
    
    assert_that(length(dim(fwd_model_vals)) == 2L)
    assert_that(ncol(fwd_model_vals) == .self$N_output)
    
    # Validate variance inflation argument. 
    inflate_var <- FALSE
    if(!is.null(var_inflation_vals)) {
      inflate_var <- TRUE
      if(is.null(dim(var_inflation_vals))) {
        assert_that(length(var_inflation_vals) == N_output)
        var_inflation_vals <- matrix(var_inflation_vals, nrow=nrow(fwd_model_vals), 
                                     ncol=N_output, byrow=TRUE)
      }
      assert_that(all(dim(var_inflation_vals) == dim(fwd_model_vals)))
    }
    
    # Fetch the variance parameters. 
    sig2_val <- get_lik_par(lik_par_val)
    
    # Construct log likelihood.
    llik <- vector(mode="numeric", length=nrow(fwd_model_vals))
    for(i in seq_along(llik)) {
      sig2_val_i <- ifelse(inflate_var, sig2_val + var_inflation_vals[i,], sig2_val)
      llik[i] <- -0.5 * sum(mult_vec_with_mat_rows(1/sig2_val_i, 
                                               add_vec_to_mat_rows(-fwd_model_vals[i,], .self$y)^2))  
      if(normalize || !conditional) llik[i] <- llik[i] - 0.5 * .self$N_obs * sum(log(sig2_val_i))
    }
    if(normalize) llik <- llik - 0.5*.self$N_obs * .self$N_output * log(2*pi)
    
    return(llik)
  }, 
  
  get_design_inputs = function(...) {
    emulator_model$X
  },
  
  get_design_llik = function(lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, ...) {
    assemble_llik(emulator_model$Y, lik_par_val=lik_par_val, conditional=conditional, normalize=normalize, ...)
  },
  
  sample_emulator = function(input, emulator_pred_list=NULL, N_samp=1, use_cov=FALSE, include_nugget=TRUE, ...) {
    # Sample the forward model emulator at specified inputs. `input` has dimension (M,D) (where M is the 
    # number of inputs). Returns array of dimension (M, N_samp, N_output). 
    
    emulator_model$sample(get_input(input), use_cov=use_cov, include_nugget=include_nugget, 
                          N_samp=N_samp, pred_list=emulator_pred_list, ...)
  },
  
  sample = function(input, lik_par_val=NULL, emulator_pred_list=NULL, N_samp=1, use_cov=TRUE, 
                    conditional=default_conditional, normalize=default_normalize, ...) {
    # Sample the log-likelihood emulator at specified inputs. `input` has dimension (M,D) (where M is the 
    # number of inputs). Returns matrix of dimension (M, N_samp). 
    
    input <- get_input(input)
    fwd_model_samp <- .self$sample_emulator(input, emulator_pred_list, N_samp=N_samp, use_cov=use_cov, ...)
    llik_samp <- matrix(nrow=nrow(input), ncol=N_samp)

    for(i in 1:N_samp) {
      llik_samp[,i] <- assemble_llik(matrix(fwd_model_samp[,i,], nrow=nrow(input), ncol=N_output), 
                                     lik_par_val, conditional, normalize)
    }
    
    return(llik_samp)
  }, 
  
  predict = function(input, lik_par_val=NULL, emulator_pred_list=NULL, return_mean=TRUE,  
                     return_var=TRUE, return_cov=FALSE, return_cross_cov=FALSE, 
                     input_cross=NULL, conditional=default_conditional,  
                     normalize=default_normalize, log_scale=FALSE, include_nugget=TRUE, ...) {
    
    if(return_cross_cov || return_cov) {
      stop("`return_cross_cov` and `return_cov` are not yet supported for `llikEmulatorGP$predict()`.")
    }
    
    # Forward model emulator predictions. 
    if(is.null(emulator_pred_list)) {
      emulator_pred_list <- .self$predict_emulator(get_input(input), return_mean=TRUE, return_var=TRUE,
                                                   return_cov=return_cov, return_cross_cov=return_cross_cov, 
                                                   input_cross=input_cross, include_nugget=include_nugget, ...)
    } else {
      assert_that(!is.null(emulator_pred_list$mean) && !is.null(emulator_pred_list$var), 
                  msg="Log-likelihood predictive quantities require forward model emulator mean and variance.")
    }
    
    # Compute induced likelihood emulator predictions. 
    llik_pred_list <- list()
    sig2_val <- get_lik_par(lik_par_val)
    
    if(return_mean) {
      llik_plug_in_mean <- assemble_llik(emulator_pred_list$mean, lik_par_val=lik_par_val, 
                                         conditional=conditional, normalize=normalize, ...)
      llik_var_inflation <- -0.5 * .self$N_obs * rowSums(exp(log(emulator_pred_list$var) - rep(log(sig2_val), each=nrow(input))))
      llik_pred_list$mean <- llik_plug_in_mean + llik_var_inflation
    }
    
    if(return_var) {
      vars <- vector(mode="numeric", length=nrow(input))
      for(i in seq_along(vars)) {
        mse <- colSums(add_vec_to_mat_rows(-emulator_pred_list$mean[i,], .self$y)^2)
        vars[i] <- .self$N_obs * sum((emulator_pred_list$var[i,] / sig2_val)^2) + 
                    sum(mse * emulator_pred_list$var[i,] / sig2_val^2)       
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
  predict_lik = function(input, lik_par_val=NULL, emulator_pred_list=NULL, return_mean=TRUE,  
                         return_var=TRUE, return_cov=FALSE, return_cross_cov=FALSE, 
                         input_cross=NULL, conditional=default_conditional,  
                         normalize=default_normalize, log_scale=FALSE, include_nugget=TRUE, ...) {
    # The mean/variance of the likelihood emulator in this case are available in closed-form, though 
    # note that the distribution is not known, so the mean and variance might not provide a full 
    # characterization of the emulator distribution. The mean returned here corresponds to the 
    # "marginal" likelihood approximation under a Gaussian error model, which implies that the 
    # forward model GP variance is simply added to the likelihood variance. 
    
    assert_that(normalize && !conditional, 
                msg="`predict_lik()` for llikEmulatorGP currently assumes normalized and not conditional likelihood.")
    
    if(return_cross_cov || return_cov) {
      stop("`return_cross_cov` and `return_cov` are not yet supported for `llikEmulatorGP$predict_lik()`.")
    }
    
    # Forward model emulator predictions. 
    if(is.null(emulator_pred_list)) {
      emulator_pred_list <- .self$predict_emulator(get_input(input), return_mean=TRUE, return_var=TRUE,
                                                   return_cov=return_cov, return_cross_cov=return_cross_cov, 
                                                   input_cross=input_cross, include_nugget=include_nugget, ...)
    } else {
      assert_that(!is.null(emulator_pred_list$mean) && !is.null(emulator_pred_list$var), 
                  msg="Likelihood predictive quantities require forward model emulator mean and variance.")
    }
    
    # Compute induced likelihood emulator predictions. 
    lik_pred_list <- list()
    sig2 <- .self$get_lik_par(lik_par_val)
    if(return_mean) {
      lik_pred_list$log_mean <- assemble_llik(emulator_pred_list$mean, lik_par_val=sig2, 
                                              conditional=conditional, normalize=normalize, 
                                              var_inflation_vals=emulator_pred_list$var, ...)
      if(!log_scale) lik_pred_list$mean <- exp(lik_pred_list$log_mean)
    }
    
    if(return_var) {
      log_num1 <- assemble_llik(emulator_pred_list$mean, lik_par_val=0.5 * sig2, 
                                conditional=conditional, normalize=normalize, 
                                var_inflation_vals=emulator_pred_list$var, ...)
      log_num2 <- assemble_llik(emulator_pred_list$mean, lik_par_val=0.5 * sig2, 
                                conditional=conditional, normalize=normalize, 
                                var_inflation_vals=0.5 * emulator_pred_list$var, ...)
      log_2pi_term <- 0.5 * .self$N_obs * .self$N_output * log(2*pi)
      log_denom1 <-  log_2pi_term  + 0.5 * .self$N_obs * sum(log(sig2))
      log_denom2 <- log_2pi_term + 0.5 * .self$N_obs * rowSums(log(add_vec_to_mat_rows(sig2, emulator_pred_list$var)))
      lik_pred_list$log_var <- log_diff_exp(log_num1 - log_denom1, log_num2 - log_denom2)
      if(!log_scale) lik_pred_list$var <- exp(lik_pred_list$log_var)
    }
  
    return(lik_pred_list)
      
  }, 
  
  
  calc_lik_marginal_approx = function(input, lik_par_val=NULL, emulator_pred_list=NULL, 
                                      conditional=default_conditional, normalize=default_normalize,
                                      log_scale=FALSE, include_nugget=TRUE, ...) {
    
    predict_lik(input, lik_par_val=lik_par_val, emulator_pred_list=emulator_pred_list, 
                return_mean=TRUE, return_var=FALSE, conditional=conditional,  
                normalize=normalize, log_scale=FALSE, include_nugget=include_nugget, ...)$mean
    
  }, 
  
  calc_expected_lik_cond_var = function(input_eval, input_cond, lik_par_val=NULL,
                                        include_nugget=TRUE, log_scale=TRUE, plugin=FALSE, 
                                        conditional=default_conditional, normalize=default_normalize, ...) {
    assert_that(plugin, msg="llikEmulatorFwdGaussDiag$calc_expected_lik_cond_var() requires `plugin=TRUE`.")
    
    N_eval <- nrow(input_eval)
    sig2 <- .self$get_lik_par(lik_par_val)
    gp_copy <- .self$emulator_model$copy(shallow=FALSE)
    
    # GP emulator predictions at evaluation locations and at conditioning locations. 
    gp_pred <- gp_copy$predict(input_cond, return_mean=TRUE, return_var=FALSE, ...)
    gp_pred_eval <- gp_copy$predict(input_eval, return_mean=TRUE, return_var=TRUE, ...)
    
    # Update the GP model, treating the predictive mean as the observed response at 
    # the evaluation locations. This leaves the predictive mean unchanged, but updates
    # the predictive variance. 
    gp_copy$update(input_cond, gp_pred$mean, update_hyperpar=FALSE, ...)
    
    # Predict with the conditional ("cond") GP (i.e., the updated GP) at the  
    # evaluation locations.
    gp_pred_cond <- gp_copy$predict(input_eval, return_mean=FALSE, return_var=TRUE, ...)
    
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
  
  
  sample = function(input, sig2, N_samp=1, use_cov=FALSE, include_nugget=TRUE, sum_output_llik=TRUE,
                    conditional=default_conditional, normalize=FALSE, ...) {
    samp <- emulator_model$sample(input, use_cov=use_cov, include_nugget=include_nugget, N_samp=N_samp)
    
    for(j in 1:N_output) {
      samp[,,j] <- -0.5 * samp[,,j] / sig2[j]
      if(normalize || !conditional) samp[,,j] <- samp[,,j] - 0.5*N_obs[j]*log(sig2[j])
    }
    
    if(normalize) samp <- samp - 0.5*sum_obs*log(2*pi)
    if(sum_output_llik) return(rowSums(samp, dims=2))
    return(samp)
  }
  
)

