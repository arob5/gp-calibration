#
# gp.r
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
# just the exact likelihood, which is useful for testing code). Given the 
# generality, not all of the core methods may apply to all approximate 
# likelihoods. For example, the `mean_lik()` method is intended to return 
# the expectation of the likelihood (i.e. the expectation of the 
# exponentiated log-likelihood). In certain cases, there may not be a closed
# form for this expression and the user may just implement this method 
# to return an error if it is called. Alternatively, the user may want to 
# implement this function to return a numerical approximation of this 
# quantity. The field `model` is intended to be a model object of some 
# sort (e.g. a gpWrapper instance) which is used to define the llik 
# approximation. `lik_par` is stores the likelihood parameters. These 
# can be stored in advance if the likelihood parameters will be fixed
# throughout an analysis; otherwise, this field can be left NULL. Note 
# that the likelihood approximation may be a function of the likelihood 
# parameters (e.g. the likelihood parameters are included as inputs to 
# a GP emulator). 
#
# This class organizes the pieces comprising the likelihood into 
# (1) "inputs": the primary parameters of interest; e.g. the parameters 
# characterizing the forward model in an inverse problem; (2) "lik_par", 
# other likelihood parameters (e.g. variance parameters) that define 
# the observation model, may be unknown, but are not of primary interest; 
# and (3) "lik_components", which are fixed quantities such as a 
# sample size. 
# 
# Implementing a deterministic approximation: 
# There are different ways to approach this. One is to just have the 
# methods `mean_log()`, `mean_lik()`, `sample()` all return the 
# determistic approximation. It may make sense to implement 
# `var_log()` and `var_lik()` to return 0 or throw an error, depending 
# on the specific needs. Note that this same method could be used 
# to implement the exact (not approximate) likelihood as well, which 
# could then be used for testing purposes. 
# -----------------------------------------------------------------------------

#
# TODO: update this so that it is defined by default to be a CONDITIONAL likelihood, 
# viewed as a function of `input` with `lik_par` being the other parameters being 
# conditioned on, and `lik_component` being other fixed stuctural information (such 
# as sample size, number of outputs, etc.). Thus, by default all functions should 
# output unnormalized values, dropping terms that only depend on `lik_par`. Overriding 
# the default would then include these terms. An alternative is to have a class 
# attribute that stores the default behavior.

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
  
  assemble_llik = function(input, lik_par_val=NULL, N_samp=1, conditional=default_conditional, 
                           normalize=default_normalize, ...) {
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
  
  calc_quantiles = function(p, input=NULL, lik_par_val=NULL, conditional=default_conditional, 
                             normalize=default_normalize, llik_pred_list=NULL,
                             lower_tail=TRUE, ...) {
    .NotYetImplemented()
  },
  
  calc_confidence_interval = function(input=NULL, lik_par_val=NULL, conditional=default_conditional, 
                                      normalize=default_normalize, llik_pred_list=NULL, CI_prob=0.9, ...) {

    if(is.null(llik_pred_list)) {
      llik_pred_list <- .self$predict(input, lik_par_val=lik_par_val, conditional=conditional, normalize=normalize, ...)
    }
    
    CI_list <- list()
    CI_tail_prob <- 0.5 * (1-CI_prob)
    CI_list$upper <- .self$calc_quantiles(p=CI_tail_prob, input=input, lik_par_val=lik_par_val,  
                                          conditional=conditional, normalize=normalize, llik_pred_list=llik_pred_list,
                                          lower_tail=TRUE, ...)
    CI_list$lower <- .self$calc_quantiles(p=CI_tail_prob, input=input, lik_par_val=lik_par_val,  
                                          conditional=conditional, normalize=normalize, llik_pred_list=llik_pred_list,
                                          lower_tail=FALSE, ...)
    
    return(CI_list)
  },

  get_design_inputs = function(...) {
    .NotYetImplemented()
  },
  
  get_design_llik = function(lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, ...) {
    .NotYetImplemented()
  },
  
  plot_llik_samp_1d = function(input, lik_par_val=NULL, N_samp=1, conditional=default_conditional, 
                               normalize=default_normalize, true_llik=NULL, include_design=FALSE, ...) {
    
    assert_that(dim_input==1, msg=paste0("plot_llik_samp_1d() requires 1d input space. input_dim = ", dim_input))
    
    input <- get_input(input)
    llik_samp <- .self$sample(input, lik_par=lik_par_val, N_samp=N_samp, ...)
    
    plt <- ggmatplot(input, llik_samp, plot_type="line", color="gray") + 
            theme(legend.position = "none") + 
            ggtitle("Log Likelihood Samples") + 
            xlab(input_names) + ylab(paste0("Log Likelihood: ", llik_label))

    if(!is.null(true_llik)) {
      df <- data.frame(x=input[,1], y=drop(true_llik))
      plt <- plt + geom_line(aes(x=x, y=y), df, inherit.aes=FALSE, color="red")
    }
    
    if(include_design) {
      design_df <- data.frame(x=drop(get_design_inputs()), 
                              y=drop(get_design_llik(lik_par_val, conditional, normalize)))
      plt <- plt + geom_point(aes(x=x, y=y), design_df, inherit.aes=FALSE, color="red")
    }

    return(plt)
    
  },
  
  plot_llik_pred_1d = function(input, lik_par_val=NULL, conditional=default_conditional,
                               normalize=default_normalize, include_CI=FALSE, 
                               CI_prob=0.9, llik_pred_list=NULL, true_llik=NULL, 
                               xlab=input_names, ylab="llik", plot_title=NULL, ...) {
    assert_that(dim_input==1, msg=paste0("plot_llik_pred_1d() requires 1d input space. dim_input = ", dim_input))
    
    # Compute required predictive quantities if not already provided. 
    if(is.null(llik_pred_list)) {
      llik_pred_list <- .self$predict(input, lik_par_val=lik_par_val, return_mean=TRUE, return_var=TRUE, 
                                      conditional=conditional, normalize=normalize, ...)
    }
    
    # Plot title and labels.
    if(is.null(plot_title)) {
      plot_title <- "Log Likelihood Emulator Predictions"
      if(include_CI) plot_title <- paste0(plot_title, ", ", 100*CI_prob, "% CI")
    }
    if(!normalize) {
      ylab <- paste0(ylab, ", ", ifelse(conditional, "unnormalized/conditional", "unnormalized"))
    }
    
    # Compute confidence interval. 
    if(include_CI) CI_list <- .self$calc_confidence_interval(llik_pred_list=llik_pred_list, CI_prob=CI_prob)
      
    # Produce plot. 
    plt <- plot_pred_1d_helper(X_new=drop(input), pred_mean=drop(llik_pred_list$mean), 
                               include_CI=include_CI, CI_lower=CI_list$lower, CI_upper=CI_list$upper, 
                               y_new=drop(true_llik), X_design=drop(.self$get_design_inputs(...)), 
                               y_design=drop(get_design_llik(lik_par_val=lik_par_val, conditional=conditional, normalize=normalize, ...)), 
                               plot_title=plot_title, xlab=xlab, ylab=ylab)
                        
    return(plt)                     
  }
  
)


# -----------------------------------------------------------------------------
# llikSumEmulator Class 
# -----------------------------------------------------------------------------

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
  }, 
  
  predict_exp = function(input, lik_par_val=NULL, return_mean=TRUE, return_var=TRUE, 
                         return_cov=FALSE, return_cross_cov=FALSE, input_cross=NULL,
                         conditional=default_conditional, normalize=default_normalize,
                         log_moments=FALSE, llik_pred_list=NULL, sum_terms=FALSE,
                         labels=llik_label, ...) {
    # Note that the mean of the likelihood is required in to compute the var/cov of the 
    # likelihood, so it is always returned regardless of the value `return_mean`. Also 
    # the log moments are always returned. The logic here is that numerical overflow is 
    # common on the exponential scale; a costly prediction computation could be wasted 
    # if only the exponential is returned. Always returning the moments on the log scale
    # thus acts as a backup. 
    
    # In this case, simply call `predict_exp()` method for each llikEmulator term, and then return 
    # the results as a list. 
    if(!sum_terms) {
      pred_exp_list <- list()
      for(lbl in labels) {
        pred_exp_list[[lbl]] <- llik_emulator_terms[[lbl]]$predict_exp(input, lik_par_val=lik_par_val[[lbl]], 
                                                                       return_mean=return_mean, return_var=return_var,
                                                                       return_cov=return_cov, return_cross_cov=return_cross_cov,
                                                                       input_cross=input_cross, conditional=conditional, 
                                                                       normalize=normalize, log_moments=log_moments, 
                                                                       llik_pred_list=llik_pred_list[[lbl]])
      }
      return(pred_exp_list)
    }
    
    # If considering the llik sum, the `predict_exp` method is only valid when the distribution of 
    # the exponential of the sum is known. Currently, this means that the llikSumEmulator must have 
    # a Gaussian predictive distribution, so that its exponential is log normal. 
    assert_that(llik_pred_dist=="Gaussian", 
                msg="`llikSumEmulator$predict_exp(sum_terms=TRUE, ...)` is only valid when `llik_pred_dist=='Gaussian'`.")
    
    # Compute Gaussian llikSumEmulator predictions. 
    if(is.null(llik_pred_list)) {
      # Computing the log-normal moments requires both the Gaussian mean and var (or cov).
      llik_pred_list <- .self$predict(input, lik_par_val=lik_par_val, return_mean=TRUE, 
                                      return_var=!return_cov, return_cov=return_cov, 
                                      return_cross_cov=return_cross_cov, input_cross=input_cross,
                                      conditional=conditional, normalize=normalize, labels=labels, ...)
    }
    
    # The log of the likelihood expectation is always returned. 
    return_list <- list()
    return_list$log_mean <- llik_pred_list$mean + 0.5 * llik_pred_list$var
    
    if(return_cov) {
      return_list$log_cov <- log_exp_minus_1(llik_pred_list$cov) + 
                             outer(return_list$log_mean, return_list$log_mean, FUN="+")
      return_list$log_var <- diag(return_list$log_cov)
    } else if(return_var) {
      return_list$log_var <- log_exp_minus_1(llik_pred_list$var) + 2*return_list$log_mean
    }
    
    # Cross covariance also requires the predictive mean and variance at inputs `input_cross`. 
    if(return_cross_cov) {
      llik_pred_list_cross <- .self$predict(input_cross, lik_par_val, return_mean=TRUE, return_var=TRUE, 
                                            conditional=conditional, normalize=normalize, ...)
      return_list$log_mean_cross <- llik_pred_list_cross$mean + 0.5 * llik_pred_list_cross$var                                   
      
      return_list$log_cross_cov <- log_exp_minus_1(llik_pred_list$cross_cov) + 
                                   outer(return_list$log_mean, return_list$log_mean_cross, FUN="+")
    }
    
    return(return_list)
  }
  
)


# -----------------------------------------------------------------------------
# llikEmulatorMultGausGP class
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
  
  assemble_llik = function(SSR, lik_par=NULL, conditional=default_conditional, normalize=default_normalize) {
    # SSR should be N_input x N_samp. 
    
    # Fetch the variance parameter. 
    sig2 <- get_lik_par(lik_par)
    
    # Construct likelihood using SSR.  
    llik <- -0.5 * SSR / sig2
    if(normalize || !conditional) llik <- llik - 0.5*N_obs*log(sig2)
    if(normalize) llik <- llik - 0.5*N_obs*log(2*pi)
    
    return(llik)
    
  },
  
  get_design_inputs = function(...) {
    emulator_model$X
  },
  
  get_design_llik = function(lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, ...) {
    assemble_llik(emulator_model$Y, lik_par=lik_par_val, conditional=conditional, normalize=normalize)
  },
  
  sample_emulator = function(input, N_samp=1, use_cov=FALSE, include_nugget=TRUE, adjustment="rectified", ...) {
    emulator_model$sample(get_input(input), use_cov=use_cov, include_nugget=include_nugget, 
                          N_samp=N_samp, adjustment=adjustment)[,,1,drop=FALSE]             
  },
  
  sample = function(input, lik_par=NULL, N_samp=1, use_cov=FALSE, include_nugget=TRUE,
                    conditional=default_conditional, normalize=default_normalize, ...) {
    # Sample SSR. 
    samp <- sample_emulator(get_input(input), N_samp, use_cov, include_nugget)
    
    # Compute unnormalized or normalized log-likelihood. 
    assemble_llik(samp, lik_par, conditional, normalize)
  }, 
  
  predict = function(input, lik_par_val=NULL, return_mean=TRUE, return_var=TRUE, 
                     return_cov=FALSE, return_cross_cov=FALSE, input_cross=NULL,
                     conditional=default_conditional, 
                     normalize=default_normalize, include_nugget=TRUE, ...) {
    
    pred_list <- .self$emulator_model$predict(input, return_mean=return_mean, return_var=return_var,
                                              return_cov=return_cov, return_cross_cov=return_cross_cov,
                                              X_cross=input_cross, include_nugget=include_nugget)
    
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
  
  predict_exp = function(input, lik_par_val=NULL, return_mean=TRUE, return_var=TRUE, 
                         return_cov=FALSE, return_cross_cov=FALSE, input_cross=NULL,
                         conditional=default_conditional, 
                         normalize=default_normalize, include_nugget=TRUE, 
                         only_log_moments=FALSE, llik_pred_list=NULL, ...) {
    # Note that the mean of the likelihood is required in to compute the var/cov of the 
    # likelihood, so it is always returned regardless of the value `return_mean`. Also 
    # the log moments are always returned. The logic here is that numerical overflow is 
    # common on the exponential scale; a costly prediction computation could be wasted 
    # if only the exponential is returned. Always returning the moments on the log scale
    # thus acts as a backup. 
    
    if(is.null(llik_pred_list)) {
      # Computing the log-normal moments requires both the Gaussian mean and var (or cov).
      llik_pred_list <- .self$predict(input, lik_par_val=lik_par_val, return_mean=TRUE, return_var=!return_cov, 
                                      return_cov=return_cov, return_cross_cov=return_cross_cov, 
                                      input_cross=input_cross, conditional=conditional, 
                                      normalize=normalize, include_nugget=include_nugget, ...)
    }
    
    # The log of the likelihood expectation always returned. 
    return_list <- list()
    return_list$log_mean <- drop(llik_pred_list$mean) + 0.5 * drop(llik_pred_list$var)
    
    # The log of the likelihood variance is always returned provided `return_cov` or 
    # `return_var` is TRUE. Note that log covariance matrices are not computed since 
    # covariance matrices may contain negative values. 
    if(return_var || return_cov) {
      return_list$log_var <- log_exp_minus_1(llik_pred_list$var) + 2*return_list$log_mean
    }
    
    # If only values of the log scale are requested, then return them now. 
    if(only_log_moments) {
      if(return_cov || return_cross_cov) message("Not including `cov` or `cross_cov` on the log-scale since these matrices may contain negative values.")
      return(return_list)
    }
    
    # Otherwise, return the values on the exponential scale as well. Note that these computations may risk 
    # numerical overflow. 
    if(return_cov) {
      return_list$cov <- exp(outer(return_list$log_mean, return_list$log_mean, FUN="+")) * (exp(llik_pred_list$cov)-1)
      return_list$var <- diag(return_list$cov)
    } else if(return_var) {
      return_list_var <- exp(return_list$log_var)
    }
    
    if(return_mean) return_list$mean <- exp(return_list$log_mean)
    
    # Cross covariance also requires the predictive mean and variance at inputs `input_cross`. 
    if(return_cross_cov) {
      llik_pred_list_cross <- .self$predict(input_cross, lik_par_val=lik_par_val, return_mean=TRUE,  
                                            return_var=TRUE, conditional=conditional, normalize=normalize, ...)
      return_list$log_mean_cross <- llik_pred_list_cross$mean + 0.5 * llik_pred_list_cross$var                                   
      return_list$cross_cov <- exp(outer(return_list$log_mean, return_list$log_mean_cross, FUN="+")) * (exp(llik_pred_list$cross_cov)-1)
    }
    
    return(return_list)
  }, 
  
  
  calc_quantiles = function(p, input=NULL, lik_par_val=NULL, conditional=default_conditional, 
                             normalize=default_normalize, llik_pred_list=NULL,
                             lower_tail=TRUE, include_nugget=TRUE, ...) {
    if(is.null(llik_pred_list)) {
      llik_pred_list <- .self$predict(input, lik_par_val=lik_par_val, return_mean=TRUE, return_var=TRUE, 
                                      conditional=conditional, normalize=normalize, 
                                      include_nugget=include_nugget, ...)
    }
    
    qnorm(p, drop(llik_pred_list$mean), sqrt(drop(llik_pred_list$var)), lower.tail=lower_tail)
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

llikEmulatorExactLinGauss <- setRefClass(
  Class = "llikEmulatorExactLinGauss", 
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
              lik_description="Exact linear Gaussian likelihood.",
              emulator_description="No emulation.", exact_llik=TRUE, ...)
  },
  
  run_fwd_model = function(input, ...) {
    if(is.null(.self$fwd_model_vectorized)) return(vectorize_fwd_model(input, ...))
    return(.self$fwd_model_vectorized(input, ...))
  },
  
  vectorize_fwd_model = function(input, ...) {
    model_output <- vector(mode="numeric", length=nrow(input))
    for(i in 1:nrow(input)) model_output[i] <- .self$fwd_model(input[i,], ...)
    return(matrix(model_output, ncol=1))
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

    return(matrix(llik, ncol=1))
  }, 
  
  sample_emulator = function(input, N_samp=1, ...) {
    # No emulator to sample from, simply return input. 
    input
  },
  
  sample = function(input, lik_par=NULL, N_samp=1, conditional=default_conditional,
                    normalize=default_normalize, ...) {
    
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
#    homoskedastic variance. 
# -----------------------------------------------------------------------------

llikEmulatorExactGaussDiag <- setRefClass(
  Class = "llikEmulatorExactGaussDiag", 
  contains = "llikEmulator",
  fields = list(fwd_model="matrix", y="numeric", N_obs="integer")
)

llikEmulatorExactGaussDiag$methods(
  
  initialize = function(llik_lbl, fwd_model, y_obs, sig2=NULL, default_conditional=FALSE, 
                        default_normalize=FALSE, use_fixed_lik_par=FALSE, ...) {
    
    assert_that(is.matrix(fwd_model), msg="`fwd_model` must be a matrix.")
    initFields(fwd_model=fwd_model, N_obs=nrow(fwd_model))
    d <- ncol(fwd_model)
    
    if(!is.null(colnames(fwd_model))) input_names_val <- colnames(fwd_model)
    else input_names_val <- paste0("input", 1:d)
    
    if(!is.null(sig2)) {
      assert_that(is.numeric(sig2) && ((length(sig2)==N_obs) || (length(sig2)==1)) && all(sig2>0),
                  msg="`sig2` must be either vector of length `N_obs` or 1 containing positive numbers.")
    }
    
    y_obs <- drop(y_obs)
    assert_that(length(y_obs)==N_obs, 
                msg="Number of observations implied by `y_obs` and `fwd_model` disagree.")
    initFields(y=y_obs)
    
    callSuper(emulator_model=NULL, llik_label=llik_lbl, lik_par=sig2, dim_input=d,
              default_conditional=default_conditional, input_names=input_names_val,
              default_normalize=default_normalize, use_fixed_lik_par=use_fixed_lik_par, 
              lik_description="Exact linear Gaussian likelihood, diagonal covariance structure.",
              emulator_description="No emulation.", ...)
  },
  
  get_lik_par = function(lik_par_val=NULL, ...) {
    if(use_fixed_lik_par) lik_par_val <- lik_par
    
    assert_that(!is.null(lik_par_val), 
                msg="`lik_par_val` arg must be non-NULL if `use_fixed_lik_par` is FALSE.")
    
    if(length(lik_par_val)==1) return(rep(lik_par_val, N_obs))
    
    assert_that(length(lik_par_val)==N_obs, msg="`lik_par_val` length not equal to 1 or `N_obs`.")
    return(lik_par_val)
  },
  
  assemble_llik = function(input, lik_par_val=NULL, conditional=default_conditional, normalize=default_normalize, ...) {
    # `input` should be N_input x N_samp. 
    
    # Fetch the variance parameters. 
    sig2_val <- get_lik_par(lik_par_val)
    
    # Construct log likelihood. 
    llik <- -0.5 * colSums((y - fwd_model %*% t(input))^2 / sig2_val, na.rm=TRUE)
    if(normalize || !conditional) llik <- llik - 0.5 * sum(log(sig2_val))
    if(normalize) llik <- llik - 0.5*N_obs*log(2*pi)

    return(matrix(llik, ncol=1))
  }, 
  
  sample_emulator = function(input, N_samp=1, ...) {
    # No emulator to sample from, simply return input. 
    input
  },
  
  sample = function(input, lik_par=NULL, N_samp=1, conditional=default_conditional,
                    normalize=default_normalize, ...) {
    
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

