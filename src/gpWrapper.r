#
# gp.r
# Class definitions for the gp, gpInd, and gpLik classes. All classes defined 
# using reference classes. 
#
# Andrew Roberts
# 

library(assertthat)
library(abind)

# -----------------------------------------------------------------------------
# gpWrapper: the base Gaussian Process (GP) wrapper class.
#
# This class represents a "wrapper" over around an GP implementation from a
# different R package. It serves as an interface between these other GP 
# packages and PEcAn code. The idea is to define a modular object with 
# a standardized interface providing access to common GP functionality 
# (fit, predict, sample, plot, etc.) that can then be plugged into downstream 
# tasks. This allows different R packages to be swapped in without modifying 
# the downstream code. 
# -----------------------------------------------------------------------------

gpWrapper <- setRefClass(
   Class = "gpWrapper", 
   fields = list(gp_model="ANY", lib="character", X="matrix", Y="matrix", 
                 X_dim="integer", Y_dim="integer", non_na_idx="matrix",
                 scale_input="logical", normalize_output="logical",
                 X_bounds="matrix", Y_mean="numeric", Y_std="numeric",
                 X_names="character", Y_names="character",
                 X_train="matrix", Y_train="matrix", nugget="numeric", 
                 valid_kernels="character", valid_mean_funcs="character", 
                 kernel_name_map="list", mean_func_name_map="list")
)

gpWrapper$methods(
  
  initialize = function(X, Y, scale_input=FALSE, normalize_output=FALSE, 
                        x_names=NULL, y_names=NULL, 
                        nugget=sqrt(.Machine$double.eps), ...) {
    
    # Handle missing values. 
    assert_that(!any(is.na(X)), msg="NAs found in X.")
    if(any(is.na(Y))) message("NAs found in Y will be removed (for each output variable separately).")
    
    initFields(X=X, Y=Y, X_dim=ncol(X), Y_dim=ncol(Y), scale_input=scale_input, 
               normalize_output=normalize_output, X_bounds=apply(X, 2, range),
               non_na_idx=!is.na(Y), nugget=nugget)
    
    if(normalize_output) {
      sd_rm_na <- function(x) sd(x, na.rm=TRUE)
      initFields(Y_mean=colMeans(Y, na.rm=TRUE), Y_std=apply(Y, 2, sd_rm_na))
      initFields(Y_train=.self$normalize(Y))
    } else {
      initFields(Y_train=Y)
    }
    
    if(scale_input) {
      initFields(X_train=.self$scale(X))
    } else {
      initFields(X_train=X)
    }
    
    initFields(valid_kernels=c("Gaussian", "Matern5_2", "Matern3_2"), 
               valid_mean_funcs=c("constant", "linear", "quadratic"))
    
    if(is.null(x_names)) x_names <- paste0("x", 1:X_dim)
    if(is.null(y_names)) y_names <- paste0("y", 1:Y_dim)
    initFields(X_names=x_names, Y_names=y_names, ...)
  },
  
  scale = function(Xnew, inverse=FALSE) {
    if(inverse) {
      Xnew <- Xnew %*% diag(X_bounds[2,] - X_bounds[1,], X_dim) + 
              matrix(X_bounds[1,], nrow=nrow(Xnew), ncol=X_dim, byrow=TRUE)
    } else {
      Xnew <- (Xnew - matrix(X_bounds[1,], nrow=nrow(Xnew), ncol=X_dim, byrow=TRUE)) %*% 
              diag(1/(X_bounds[2,] - X_bounds[1,]), X_dim)
    }
    
    return(Xnew)
  }, 
  
  normalize = function(Ynew, inverse=FALSE) {
    if(inverse) {
      Ynew <- Ynew * matrix(Y_std, nrow=nrow(Ynew), ncol=Y_dim, byrow=TRUE) + 
              matrix(Y_mean, nrow=nrow(Ynew), ncol=Y_dim, byrow=TRUE)
    } else {
      Ynew <- (Ynew - matrix(Y_mean, nrow=nrow(Ynew), ncol=Y_dim, byrow=TRUE)) / 
              matrix(Y_std, nrow=nrow(Ynew), ncol=Y_dim, byrow=TRUE)
    }
    
    return(Ynew)
  },
  
  fit_package = function(X_fit, y_fit, kernel_name="Gaussian", mean_func_name="constant", fixed_pars=list(), ...) {
    err_msg <- "A fit_package() method must be implemented for each class inheriting from 
                gpWrapper. gpWrapper does not implement its own fit_package() method. The
                fit() method should ultimately call fit_package()."
    stop(err_msg)
  },
  
  
  fit = function(kernel_name="Gaussian", mean_func_name="constant", estimate_nugget=TRUE, fixed_pars=list(), ...) {
    fits_list <- vector(mode="list", length=Y_dim)
    for(i in seq_along(fits_list)) {
      idx_sel <- non_na_idx[,i]
      fits_list[[i]] <- fit_package(X_train[idx_sel,,drop=FALSE], Y_train[idx_sel, i], 
                                    kernel_name, mean_func_name, estimate_nugget, fixed_pars, ...)
    }

    gp_model <<- fits_list
  }, 

  fit_parallel = function(kernel_name="Gaussian", mean_func_name="constant", fixed_pars=list(), ...) {
    .NotYetImplemented()
  },
  
  
  predict_package = function(X_new, output_idx, return_mean=TRUE, return_var=TRUE, return_cov=FALSE, X_cov=NULL, include_nugget=TRUE) {
    err_msg <- "A predict_package() method must be implemented for each class inheriting from 
                gpWrapper. gpWrapper does not implement its own predict_package() method. The
                predict() method should ultimately call predict_package()."
    stop(err_msg)
  }, 
  
  predict = function(X_new, return_mean=TRUE, return_var=TRUE, return_cov=FALSE, 
                     X_cov=NULL, include_nugget=TRUE) {
    
    # Scale inputs, if required. 
    if(scale_input) {
      X_new <- .self$scale(X_new)
      if(return_cov && !is.null(X_cov)) X_cov <- .self$scale(X_cov)
    }
    
    # If covariance is requested, default to computing cov at inputs `X_new`. 
    if(return_cov && is.null(X_cov)) X_cov <- X_new
    
    # Predict for each independent GP. 
    pred_list <- vector(mode="list", length=Y_dim)
    for(i in seq_along(pred_list)) {
      pred_list[[i]] <- predict_package(X_new, i, return_mean, return_var, return_cov, X_cov, include_nugget)
    }
    
    names(pred_list) <- Y_names
    mean_pred <- do.call(cbind, lapply(pred_list, function(l) l$mean))
    var_pred <- do.call(cbind, lapply(pred_list, function(l) l$var))
    cov_pred <- abind(lapply(pred_list, function(l) l$cov), along=3)
    
    # Set negative variances to 0. 
    neg_var_idx <- (var_pred < 0)
    if(any(neg_var_idx)) {
      message("Thresholding negative variance predictions at 0.")
      var_pred[neg_var_idx] <- 0
    }
    
    # Return outputs to unnormalized scale. 
    if(normalize_output) {
      if(return_mean) mean_pred <- .self$normalize(mean_pred, inverse=TRUE)
      if(return_var) var_pred <- var_pred %*% diag(Y_std^2)
      if(return_cov) for(i in 1:Y_dim) cov_pred[,,i] <- Y_std[i]^2 * cov_pred[,,i]
    }
      
    return(list(mean=mean_pred, var=var_pred, cov=cov_pred))
  }, 
  
  
  predict_parallel = function(X_new, return_mean=TRUE, return_var=TRUE, 
                              return_cov=FALSE, X_cov=NULL, include_nugget=TRUE) {
    .NotYetImplemented()
  },
  
  
  sample = function(X_new, use_cov=FALSE, include_nugget=TRUE, N_samp=1, pred_list=NULL) {
    # Compute required predictive quantities if not already provided. 
    if(is.null(pred_list)) {
      pred_list <- predict(X_new, return_mean=TRUE, return_var=!use_cov, return_cov=use_cov, 
                           X_cov=X_new, include_nugget=include_nugget)
    }
    
    # Compute lower Cholesky factors of the predictive covariance matrices. 
    # If not using predictive cov, Cholesky factors are diagonal with standard devs on diag. 
    if(use_cov) {
      if(is.null(pred_list$chol_cov)) pred_list$chol_cov <- abind(lapply(1:Y_dim, 
                                                                         function(i) t(chol(pred_list$cov[,,i]))), along=3)
    } else {
      pred_list$chol_cov <- abind(lapply(1:Y_dim, function(i) diag(sqrt(pred_list$var[,i]))), along=3)
    }
    
    samp <- array(dim=c(nrow(X_new), N_samp, Y_dim))
    for(i in 1:Y_dim) {
      samp[,,i] <- sample_Gaussian_chol(pred_list$mean[,i], pred_list$chol_cov[,,i], N_samp)
    }
    
    return(samp)
    
  },
  
  plot_pred_1d = function(X_new, include_nugget=TRUE, CI_prob=0.9, pred_list=NULL, Y_new=NULL) {
  
    assert_that(X_dim==1, msg=paste0("plot_pred_1d() requires 1d input space. X_dim = ", X_dim))
    
    # Compute required predictive quantities if not already provided. 
    if(is.null(pred_list)) {
      pred_list <- predict(X_new, return_mean=TRUE, return_var=TRUE, include_nugget=include_nugget)
    }
    
    CI_tail_prob <- 0.5 * (1-CI_prob)
    plts <- vector(mode="list", length=Y_dim)
    for(i in 1:Y_dim) {
      df_train <- data.frame(x=X[non_na_idx[,i],1], y=Y[non_na_idx[,i],i])
      df_pred <- data.frame(x=X_new[,1], y_mean=pred_list$mean[,i], y_sd=sqrt(pred_list$var[,i]))
      df_pred$CI_upper <- qnorm(CI_tail_prob, df_pred$y_mean, df_pred$y_sd)
      df_pred$CI_lower <- qnorm(CI_tail_prob, df_pred$y_mean, df_pred$y_sd, lower.tail=FALSE)
      if(!is.null(Y_new)) df_pred$y <- Y_new[,i]
      
      plts[[i]] <- ggplot(df_pred) + 
                   geom_line(aes(x, y_mean), color="blue") + 
                   geom_line(aes(x, CI_upper), color="gray") + 
                   geom_line(aes(x, CI_lower), color="gray") +
                   geom_point(aes(x,y), df_train, color="red")
      if(!is.null(Y_new)) plts[[i]] <- plts[[i]] + geom_line(aes(x,y), linetype="dotted")
    }
    
    return(plts)
    
  },
  
  plot_pred = function(X_new, Y_new, include_CI=FALSE, include_nugget=TRUE, CI_prob=0.9, pred_list=NULL) {
    # Compute required predictive quantities if not already provided. 
    if(is.null(pred_list)) {
      pred_list <- predict(X_new, return_mean=TRUE, return_var=include_CI, include_nugget=include_nugget)
    }
    
    CI_tail_prob <- 0.5 * (1-CI_prob)
    plts <- vector(mode="list", length=Y_dim)
    
    for(i in 1:Y_dim) {
      # df_train <- data.frame(x=X[non_na_idx[,i],1], y=Y[non_na_idx[,i],i])
      df_pred <- data.frame(y=Y_new[,i], y_mean=pred_list$mean[,i])
      if(include_CI) {
        df_pred$y_sd <- sqrt(pred_list$var[,i])
        df_pred$CI_upper <- qnorm(CI_tail_prob, df_pred$y_mean, df_pred$y_sd)
        df_pred$CI_lower <- qnorm(CI_tail_prob, df_pred$y_mean, df_pred$y_sd, lower.tail=FALSE)
      }
      
      plts[[i]] <- ggplot(df_pred) + 
                   geom_point(aes(y, y_mean), color="black") + 
                   geom_abline(slope=1, intercept=0, color="red") + 
                   xlab("Observed") + ylab("Predicted")
    }
    
    return(plts)
  },
  
  plot_loocv = function() {
    .NotYetImplemented()
  },
  
  sample_Gaussian_chol = function(mu, C_chol_lower, N_samp=1) {
    drop(mu) + C_chol_lower %*% matrix(rnorm(nrow(C_chol_lower)*N_samp), ncol=N_samp)
  },
  
  map_kernel_name = function(kernel_name) {
    return(kernel_name_map[[kernel_name]])
  },
  
  map_mean_func_name = function(mean_func_name) {
    return(mean_func_name_map[[mean_func_name]])
  }

)


#
# Package-Specific Classes that Inherit from gpWrapper:
# Each time a new R GP library is to be added to the framework, a new class 
# must be defined for the package which inherits from the gpWrapper base 
# class. 
#
# Attributes defined by gpWrapper but that must be specified in the package-specific
# class:
#    kernel_name_map: list, to be interpreted as a dictionary with the names
#                     attribute equal to the dictionary keys, which should be 
#                     valid kernel names (the valid kernel names are defined 
#                     by the `valid_kernels` attribute of gpWrapper). The 
#                     values of the dictionary are the corresponding name 
#                     used by the specific GP package. 
#    mean_func_name_map: list, the analog of `kernel_name_map` for the mean 
#                        functions. Note that `valid_mean_funcs` in 
#                        gpWrapper provides the set of valid names. 
#    
#
# Required methods: 
#    fit_package
#    predict_package: The analog to `fit_package` but for prediction. 
#
# fit_package():
#    Handles the translation from the GP fit interface defined 
#    in gpWrapper to the format required by the GP package. 
#    fit_package`.
#
#    Args:
#        X_fit, y_fit: Matrix and numeric vector, respectively. These are the exact
#                      design matrix and response vector that will be used to fit the 
#                      GP model. Note that they may differ from X_train, Y_train in 
#                      that missing values may have been dropped (which is handled 
#                      by gpWrapper). 

# Methods to add: 
#    get_input_dim()
#    get_num_design()
#    get_output_dim() (for the indGPWrapper class)
#    get_pred_quantiles()
#    predict_parallel() (for indGPWRapper)
#    fit() (maybe standardize some inputs to fit method but allow it to be pretty flexible)

# -----------------------------------------------------------------------------
# gpHet: Class providing interface between calibration code and hetGP package.  
# -----------------------------------------------------------------------------

gpWrapperHet <- setRefClass(
  Class = "gpWrapperHet",
  contains = "gpWrapper"
)

gpWrapperHet$methods(
  
  initialize = function(X, Y, ...) {
    library(hetGP)
    
    initFields(kernel_name_map=list(Gaussian="Gaussian", Matern5_2="Matern5_2", Matern3_2="Matern3_2"), 
               mean_func_name_map=list(constant="beta0"))
    callSuper(X=X, Y=Y, lib="hetGP", ...)
    
    assert_that(all(names(kernel_name_map) %in% valid_kernels), 
                msg=paste0("Some kernel names not found in base gpWrapper class: ", 
                           paste(setdiff(names(kernel_name_map), valid_kernels), collapse=", ")))
    assert_that(all(names(mean_func_name_map) %in% valid_mean_funcs), 
                msg=paste0("Some mean function names not found in base gpWrapper class: ", 
                           paste(setdiff(names(mean_func_name_map), valid_mean_funcs), collapse=", ")))
  }, 
  
  fit_package = function(X_fit, y_fit, kernel_name="Gaussian", mean_func_name="constant", estimate_nugget=TRUE, 
                         fixed_pars=list(), ...) {
    
    # Deal with noiseless case. 
    if(estimate_nugget) {
      assert_that(!("g" %in% names(fixed_pars)), 
                  msg="`estimate_nugget` is TRUE but the nugget `g` is in `fixed_pars`.")
    } else {
      if(!("g" %in% names(fixed_pars))) fixed_pars[["g"]] <- nugget
    }
    
    if(length(fixed_pars) == 0) fixed_pars <- NULL
    
    gp_fits <- hetGP:::mleHomGP(X_fit, y_fit, covtype=map_kernel_name(kernel_name), known=fixed_pars, ...)
  },
  
  predict_package = function(X_new, output_idx, return_mean=TRUE, return_var=TRUE, 
                             return_cov=FALSE, X_cov=NULL, include_nugget=TRUE) {
    
    if(!return_cov) X_cov <- NULL
    
    pred <- hetGP:::predict(gp_model[[output_idx]], X_new, xprime=X_cov)
    return_list <- list()
    if(return_mean) return_list$mean <- pred$mean
    if(return_var) {
      return_list$var <- pred$sd2
      if(include_nugget) return_list$var <- return_list$var + pred$nugs
    }
    
    if(return_cov) {
      return_list$cov <- pred$cov
      # TODO: need to check that this is right. 
      if(include_nugget) diag(return_list$cov) <- diag(return_list$cov) + pred$nugs 
    }
    
    return(return_list)
  }
  
)


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
# Implementing a deterministic approximation: 
# There are different ways to approach this. One is to just have the 
# methods `mean_log()`, `mean_lik()`, `sample()` all return the 
# determistic approximation. It may make sense to implement 
# `var_log()` and `var_lik()` to return 0 or throw an error, depending 
# on the specific needs. Note that this same method could be used 
# to implement the exact (not approximate) likelihood as well, which 
# could then be used for testing purposes. 
# -----------------------------------------------------------------------------

llikEmulator <- setRefClass(
  Class = "llikEmulator", 
  fields = list(model="ANY", lik_par="ANY", lik_description="character", 
                emulator_description="character")
)

llikEmulator$methods(
  
  initialize = function(lik_description, emulator_description, model=NULL, base_lik_par=NULL, ...) {
    initFields(lik_description=lik_description, emulator_description=emulator_description,
               model=model, base_lik_par=base_lik_par)
  }, 
  
  sample = function(input, lik_par=NULL) {
    stop("`sample()` method not implemented.")
  }, 
  
  mean_log = function(input, lik_par=NULL) {
    stop("`mean_log()` method not implemented.")
  }, 
  
  var_log = function(input, lik_par=NULL) {
    stop("`var_log()` method not implemented.")
  }, 
  
  mean_lik = function(input, lik_par=NULL) {
    stop("`mean_lik()` method not implemented.")
  }, 
  
  var_lik = function(input, lik_par=NULL) {
    stop("`var_lik()` method not implemented.")
  }

)


























#
# TEST
#

# Testing base class. 
# X <- matrix(seq(10,20,length.out=5), ncol=1)
# Y <- X^2 + 0.2*matrix(rnorm(nrow(X)), ncol=1)
# 
# gp <- gpWrapper(X, Y, normalize_output=TRUE, scale_input=TRUE)
# gp$field("lib")
# gp$field("normalize_output")
# gp$field("scale_input")
# gp$field("X_names")
# gp$field("Y_names")
# 
# gp$field("Y")
# gp$field("Y_mean")
# gp$field("Y_std")
# gp$field("Y_train")
# gp$field("X")
# gp$field("X_bounds")
# gp$field("X_train")
# 
# gp$fit()
# 
# 
# # Testing hetGP wrapper.
# gpHet <- gpWrapperHet(X, Y, normalize_output=TRUE, scale_input=TRUE)
# 
# gpHet$field("lib")
# gpHet$field("normalize_output")
# gpHet$field("scale_input")
# gpHet$field("X_names")
# gpHet$field("Y_names")
# gpHet$field("Y")
# gpHet$field("Y_mean")
# gpHet$field("Y_std")
# gpHet$field("Y_train")
# gpHet$field("X")
# gpHet$field("X_bounds")
# gpHet$field("X_train")
# gpHet$mean_func_name_map
# gpHet$kernel_name_map
# gpHet$field("gp_model")




