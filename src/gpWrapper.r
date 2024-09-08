#
# gp.r
# Class definitions for the `gpEmulator` class, and classes which inherit from 
# this class. All classes defined using reference classes. 
#
# Andrew Roberts
# 
# Depends: statistical_helper_functions.r
#

library(assertthat)
library(tmvtnorm)
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
                 X_dim="integer", Y_dim="integer",
                 scale_input="logical", normalize_output="logical",
                 X_bounds="matrix", Y_mean="numeric", Y_std="numeric",
                 X_names="character", Y_names="character",
                 X_train="matrix", Y_train="matrix", default_nugget="numeric", 
                 valid_kernels="character", valid_mean_funcs="character", 
                 kernel_name_map="function", mean_func_name_map="function")
)

gpWrapper$methods(
  
  initialize = function(X, Y, scale_input=FALSE, normalize_output=FALSE, 
                        x_names=NULL, y_names=NULL, 
                        default_nugget=sqrt(.Machine$double.eps), ...) {
    
    # X and Y should only be missing when the generator is being called due to 
    # a call to `$copy()`. 
    # TODO: think of better long-term solution. 
    if(missing(X) || missing(Y)) return(NULL)
    
    # Handle missing values. 
    assert_that(is.matrix(X) && is.matrix(Y), msg="X and Y must be matrices.")
    assert_that(!any(is.na(X)), msg="NAs found in X.")
    assert_that(!any(is.na(Y)), msg="NAs found in Y.")
    initFields(X=X, Y=Y, X_dim=ncol(X), Y_dim=ncol(Y))
    
    if(is.null(x_names)) {
      x_names <- colnames(X)
      if(is.null(x_names)) x_names <- paste0("x", 1:X_dim)
    }

    if(is.null(y_names)) {
      y_names <- colnames(Y)
      if(is.null(y_names)) y_names <- paste0("y", 1:Y_dim)
    }
    
    initFields(X_names=x_names, Y_names=y_names, scale_input=scale_input,
               normalize_output=normalize_output, X_bounds=apply(X, 2, range),
               default_nugget=default_nugget)

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
    
    # Set input and output names. 
    colnames(X) <<- X_names
    colnames(Y) <<- Y_names
    colnames(X_train) <<- X_names
    colnames(Y_train) <<- Y_names
    
    initFields(valid_kernels=c("Gaussian", "Matern5_2", "Matern3_2", "Gaussian_plus_Quadratic"), 
               valid_mean_funcs=c("constant", "linear", "quadratic"), ...)
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
      fits_list[[i]] <- fit_package(X_fit=X_train, y_fit=Y_train[,i], kernel_name=kernel_name, 
                                    mean_func_name=mean_func_name, estimate_nugget=estimate_nugget, 
                                    fixed_pars=fixed_pars, ...)
    }

    gp_model <<- fits_list
  }, 

  fit_parallel = function(kernel_name="Gaussian", mean_func_name="constant", fixed_pars=list(), ...) {
    .NotYetImplemented()
  },
  
  predict_package = function(X_new, output_idx, return_mean=TRUE, return_var=TRUE, return_cov=FALSE, 
                             return_cross_cov=NULL, X_cross=NULL, include_nugget=TRUE) {
    err_msg <- "A predict_package() method must be implemented for each class inheriting from 
                gpWrapper. gpWrapper does not implement its own predict_package() method. The
                predict() method should ultimately call predict_package()."
    stop(err_msg)
  }, 
  
  predict = function(X_new, return_mean=TRUE, return_var=TRUE, return_cov=FALSE, 
                     return_cross_cov=FALSE, X_cross=NULL, include_nugget=TRUE, 
                     return_trend=TRUE, ...) {
    # Logic for all predict() functions:
    #   - `return_cov` refers to k(X_new, X_new) while `return_cross_cov` refers to k(X_new, X_cross).
    #   - The former is always diagonal, and `include_nugget==TRUE` will cause nugget variances to be 
    #     added to the diagonal. The latter may or may not be diagonal, and `include_nugget` has no
    #     effect on the cross cov matrix (nugget variances will NOT be added in either case). 
    #   - If `return_cov==TRUE`, then the variances will always be returned as well, by simply taking 
    #     the diagonal of the covariance matrix. 
    #   - `return_cross_cov` has no effect on the other settings; can be thought of as an optional add-on. 

    # Scale inputs, if required. 
    if(scale_input) {
      X_new <- .self$scale(X_new)
      if(return_cross_cov && !is.null(X_cross)) X_cross <- .self$scale(X_cross)
    }

    # If covariance is requested, default to computing scov at inputs `X_new`.
    if(return_cross_cov && is.null(X_cross)) stop("`return_cross_cov` is TRUE but `X_cross` is NULL.") 
      
    # Predict for each independent GP. 
    pred_list <- vector(mode="list", length=Y_dim)
    for(i in seq_along(pred_list)) {
      pred_list[[i]] <- predict_package(X_new=X_new, output_idx=i, return_mean=return_mean, 
                                        return_var=return_var, return_cov=return_cov, 
                                        return_cross_cov=return_cross_cov, X_cross=X_cross, 
                                        include_nugget=include_nugget, return_trend=return_trend, ...)
    }

    names(pred_list) <- Y_names
    mean_pred <- do.call(cbind, lapply(pred_list, function(l) l$mean))
    trend_pred <- do.call(cbind, lapply(pred_list, function(l) l$trend))
    var_pred <- do.call(cbind, lapply(pred_list, function(l) l$var))
    cov_pred <- abind(lapply(pred_list, function(l) l$cov), along=3)
    cross_cov_pred <- abind(lapply(pred_list, function(l) l$cross_cov), along=3)
    
    # Set negative variances to 0. 
    if(return_var) {
      neg_var_idx <- (var_pred < 0)
      if(any(neg_var_idx)) {
        message("Thresholding negative variance predictions at 0.")
        var_pred[neg_var_idx] <- 0
      }
    }
    
    if(return_cov) {
      for(i in 1:Y_dim) {
        neg_var_idx <- (diag(cov_pred[,,i]) < 0)
        if(any(neg_var_idx)) {
          eps <- .Machine$double.eps
          message("Thresholding negative diagonal values of predicted covariance at ",
                  eps, "; output ", i)
          diag(cov_pred[,,i])[neg_var_idx] <- eps
        }
      }
    }
    
    # Return outputs to unnormalized scale. 
    if(normalize_output) {
      if(return_mean) mean_pred <- .self$normalize(mean_pred, inverse=TRUE)
      if(return_trend && !is.null(trend_pred)) trend_pred <- .self$normalize(trend_pred, inverse=TRUE)
      if(return_var) var_pred <- var_pred %*% diag(Y_std^2, nrow=Y_dim)
      if(return_cov) for(i in 1:Y_dim) cov_pred[,,i] <- Y_std[i]^2 * cov_pred[,,i]
      if(return_cross_cov) for(i in 1:Y_dim) cross_cov_pred[,,i] <- Y_std[i]^2 * cross_cov_pred[,,i]
    }
      
    return(list(mean=mean_pred, var=var_pred, cov=cov_pred, 
                cross_cov=cross_cov_pred, trend=trend_pred))
  }, 
  
  
  predict_parallel = function(X_new, return_mean=TRUE, return_var=TRUE, return_cov=FALSE,
                              return_cross_cov=FALSE, X_cross=NULL, include_nugget=TRUE,
                              return_trend=TRUE) {
    .NotYetImplemented()
  },
  
  
  sample = function(X_new, use_cov=FALSE, include_nugget=TRUE, N_samp=1, pred_list=NULL, adjustment="none", ...) {
    # If `pred_list` is passed, it should have all the required components. 
    # Returns array of dimension (num input, num samp, Y_dim). 

    # Compute required predictive quantities if not already provided. 
    if(is.null(pred_list)) {
      pred_list <- predict(X_new, return_mean=TRUE, return_var=!use_cov, return_cov=use_cov, 
                           include_nugget=include_nugget)
    } else {
      if(use_cov) assert_that(!is.null(pred_list$cov))
    }
    
    # Compute lower Cholesky factors of the predictive covariance matrices. 
    # If not using predictive cov, Cholesky factors are diagonal with standard devs on diag. 
    # For truncated normal, `tmvnorm` doesn't accept Cholesky factor, so just pass covariance matrix. 
    if((adjustment=="truncated") && !use_cov) {
      pred_list$cov <- abind(lapply(1:Y_dim, function(i) diag(pred_list$var[,i], nrow=nrow(X_new))), along=3)
    } else if((adjustment != "truncated") && use_cov) {
      if(is.null(pred_list$chol_cov)) pred_list$chol_cov <- abind(lapply(1:Y_dim, function(i) t(chol(pred_list$cov[,,i]))), along=3)
    } else if(adjustment != "truncated") {
      pred_list$chol_cov <- abind(lapply(1:Y_dim, function(i) diag(sqrt(pred_list$var[,i]))), along=3)
    }
    
    samp <- array(dim=c(nrow(X_new), N_samp, Y_dim))
    for(i in 1:Y_dim) {
      if(adjustment=="truncated") { # Zero left-truncated Gaussian. 
        samp[,,i] <- t(rtmvnorm(N_samp, mean=pred_list$mean[,i], sigma=pred_list$cov[,,i], lower=rep(0, nrow(X_new))))
      } else {
        samp[,,i] <- sample_Gaussian_chol(pred_list$mean[,i], pred_list$chol_cov[,,i], N_samp)
        if(adjustment=="rectified") samp[,,i] <- pmax(samp[,,i], 0)
      }
    }
    
    return(samp)
    
  },
  
  update = function(X_new, Y_new, update_hyperpar=FALSE, ...) {
    N_design_curr <- nrow(X_train)
    .self$augment_design(X_new, Y_new)
    N_design_new <- nrow(X_train)
    
    for(i in 1:Y_dim) {
      gp_model[[i]] <<- update_package(X_train[(N_design_curr+1):N_design_new,,drop=FALSE], 
                                       Y_train[(N_design_curr+1):N_design_new, i], 
                                       output_idx=i, update_hyperpar=update_hyperpar, ...)
    }
  },
  
  update_package = function(X_new, Y_new, update_hyperpar=FALSE, ...) {
    err_msg <- "An update_package() method must be implemented for each class inheriting from 
                gpWrapper. gpWrapper does not implement its own update_package() method. The
                update() method should ultimately call update_package()."
    stop(err_msg)
  },
  
  plot_pred_1d = function(X_new, include_nugget=TRUE, include_interval=TRUE, interval_method="pm_std_dev",
                          N_std_dev=1, CI_prob=0.9, pred_list=NULL, Y_new=NULL, plot_title=NULL, 
                          xlab=NULL, ylab=NULL, ...) {
    
    assert_that(X_dim==1, msg=paste0("plot_pred_1d() requires 1d input space. X_dim = ", X_dim))
    
    # Compute required predictive quantities if not already provided. 
    if(is.null(pred_list)) {
      pred_list <- predict(X_new, return_mean=TRUE, return_var=TRUE, include_nugget=include_nugget)
    }
    
    # Default x-axis label.
    if(is.null(xlab)) xlab <- X_names

    plts <- vector(mode="list", length=Y_dim)
    for(i in 1:Y_dim) {
      
      # Default y-axis label. 
      if(is.null(ylab)) ylab <- Y_names[i]
      
      plts[[i]] <- plot_Gaussian_pred_1d(X_new=X_new[,1], pred_mean=pred_list$mean[,i], pred_var=pred_list$var[,i], 
                                         include_interval=include_interval, interval_method=interval_method, 
                                         N_std_dev=N_std_dev, CI_prob=CI_prob, y_new=Y_new[,i], X_design=X[,1], 
                                         y_design=Y[,i], xlab=xlab, ylab=ylab, plot_title=plot_title, ...)
    }
    
    return(plts)
    
  },
  
  
  plot_pred_2d = function(X_new, include_nugget=TRUE, pred_list=NULL, Y_new=NULL) {
    # Returns a list of heatmaps for the predictive mean and standard deviation. 
    
    assert_that(X_dim==2, msg=paste0("plot_pred_2d() requires 1d input space. X_dim = ", X_dim))
    stop("gpWrapper$plot_pred_2d not yet implemented.")
  },
  
  
  plot_pred = function(X_new, Y_new, include_CI=FALSE, include_nugget=TRUE, CI_prob=0.9, pred_list=NULL) {
    # Compute required predictive quantities if not already provided. 
    if(is.null(pred_list)) {
      pred_list <- predict(X_new, return_mean=TRUE, return_var=include_CI, include_nugget=include_nugget)
    }
    
    CI_tail_prob <- 0.5 * (1-CI_prob)
    plts <- vector(mode="list", length=Y_dim)
    
    for(i in 1:Y_dim) {
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
  
  map_kernel_name = function(kernel_name, ...) {
    # Additional arguments `...` sometimes include the design data to set empirical 
    # bounds/priors on the kernel hyperparameters. 
    
    assert_that(kernel_name %in% valid_kernels, 
                msg=paste0("Kernel ", kernel_name, " not found in gpWrapper$valid_kernels."))
    
    return(kernel_name_map(kernel_name, ...))
  },
  
  map_mean_func_name = function(mean_func_name, ...) {
    
    assert_that(mean_func_name %in% valid_mean_funcs, 
                msg=paste0("Mean function ", mean_func_name, " not found in gpWrapper$valid_mean_funcs"))
    
    return(mean_func_name_map(mean_func_name, ...))
  }, 
  
  augment_design = function(X_new, Y_new) {
    assert_that(is.matrix(X_new) && is.matrix(Y_new))
    assert_that(!any(is.na(X_new)) && !any(is.na(Y_new)))
    assert_that(nrow(X_new) == nrow(Y_new))
    assert_that(ncol(X_new) == X_dim)
    assert_that(ncol(Y_new) == Y_dim)
    
    X <<- rbind(X, X_new)
    Y <<- rbind(Y, Y_new)
    
    if(scale_input) X_train <<- rbind(X_train, .self$scale(X_new))
    else X_train <<- rbind(X_train, X_new)
    
    if(normalize_output) Y_train <<- rbind(Y_train, .self$normalize(Y_new))
    else Y_train <<- rbind(Y_train, Y_new)
  }, 
  
  calc_expected_cond_var = function(X_eval, X_cond, include_nugget=TRUE, log_scale=TRUE, ...) {
    # Computes the log of E_l Var[exp(f(X_eval))|f(X_cond)=l], the expected conditional variance 
    # due to conditioning on `X_cond`. The expected conditional variance is evaluated at inputs 
    # `X_eval`. The expectation is with respect to the current GP predictive distribution, 
    # which can be viewed as a posterior predictive expectation. When `gp_model` contains 
    # more than one GP, then this expectation is computed independently for each GP. 
    # Returns matrix of shape `(nrow(X_eval), Y_dim)`, where `Y_dim` is the number of 
    # independent GPs. 
    
    N_eval <- nrow(X_eval)
    gp_copy <- .self$copy(shallow=FALSE)
    
    # Emulator predictions at evaluation locations and at conditioning locations. 
    pred <- gp_copy$predict(X_cond, return_mean=TRUE, return_var=TRUE, ...)
    pred_eval <- gp_copy$predict(X_eval, return_mean=FALSE, return_var=TRUE, ...)
    
    # Update the GP model, treating the predictive mean as the observed 
    # response at the evaluation locations. 
    gp_copy$update(X_cond, pred$mean, update_hyperpar=FALSE, ...)
    
    # Predict with the conditional ("cond") GP (i.e., the updated GP) at the  
    # evaluation locations. Convert to the exponentiated scale to obtain 
    # log-normal predictive quantities.
    pred_cond <- gp_copy$predict(X_eval, return_mean=TRUE, return_var=TRUE, ...)
    lnp_cond_log_var <- matrix(nrow=N_eval, ncol=.self$Y_dim)
    for(j in 1:.self$Y_dim) {
      lnp_cond_log_var[,j] <- convert_Gaussian_to_LN(mean_Gaussian=pred_cond$mean[,j], 
                                                     var_Gaussian=pred_cond$var[,j],
                                                     return_mean=FALSE, return_var=TRUE, 
                                                     log_scale=TRUE)$log_var
    }
    
    # Compute the variance inflation term on the log scale. 
    log_var_inflation <- 2 * (pred_eval$var - pred_cond$var)
    
    # Compute log expected conditional variance.  
    log_evar <- lnp_cond_log_var + log_var_inflation
    
    if(log_scale) return(log_evar)
    return(exp(log_evar))
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
#    kernel_name_map: function, maps kernel name to either the associated kernel 
#                     object as required by the function. For some packages (e.g., 
#                     kergp) the output will be an actual object from a kernel class. 
#                     For others (e.g., hetGP), the output will be a string 
#                     specifying the kernel names used by hetGP. The `valid_kernels`
#                     attribute of `gpWrapper` provides the valid inputs to the 
#                     `kernel_name_map` function. 
#    mean_func_name_map: function, the analog of `kernel_name_map` for the mean 
#                        functions. The `valid_mean_funcs` attribute 
#                        is the analog of `valid_kernels`. 
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
    # X and Y should only be missing when the generator is being called due to 
    # a call to `$copy()`. 
    # TODO: think of better long-term solution. 
    if(missing(X) || missing(Y)) return(NULL)
    
    
    library(hetGP)
    
    kernel_map <- function(kernel_name) {
      if(kernel_name == "Gaussian") return("Gaussian")
      else if(kernel_name == "Materm5_2") return("Matern5_2")
      else if(kernel_name == "Matern3_2") return("Matern3_2")
      else stop("gpWrapperHet does not support kernel: ", kernel_name)
    }
    
    mean_func_map <- function(mean_func_name) {
      if(mean_func_name == "constant") return("beta0")
      else stop("gpWrapperHet does not support mean function: ", mean_func_name)
    }
    
    initFields(kernel_name_map=kernel_map, mean_func_name_map=mean_func_map)
    callSuper(X=X, Y=Y, lib="hetGP", ...)
  }, 
  
  fit_package = function(X_fit, y_fit, kernel_name="Gaussian", mean_func_name="constant", estimate_nugget=TRUE, 
                         fixed_pars=list(), ...) {
    # Deal with noiseless case. 
    if(estimate_nugget) {
      assert_that(!("g" %in% names(fixed_pars)), 
                  msg="`estimate_nugget` is TRUE but the nugget `g` is in `fixed_pars`.")
    } else {
      if(!("g" %in% names(fixed_pars))) fixed_pars[["g"]] <- default_nugget
    }
    
    if(length(fixed_pars) == 0) fixed_pars <- NULL
    
    gp_fit <- hetGP:::mleHomGP(X_fit, y_fit, covtype=map_kernel_name(kernel_name), known=fixed_pars, ...)
  },
  
  predict_package = function(X_new, output_idx, return_mean=TRUE, return_var=TRUE, 
                             return_cov=FALSE, return_cross_cov=FALSE, X_cross=NULL, 
                             include_nugget=TRUE, ...) {
    
    if(return_cov) X_prime <- X_new
    else X_prime <- NULL

    pred <- hetGP:::predict(gp_model[[output_idx]], X_new, xprime=X_prime)
    return_list <- list()
    if(return_mean) return_list$mean <- pred$mean
    if(return_cov) {
      return_list$cov <- pred$cov
      if(include_nugget) diag(return_list$cov) <- diag(return_list$cov) + pred$nugs
      return_list$var <- diag(return_list$cov)
    } else if(return_var) {
      return_list$var <- pred$sd2
      if(include_nugget) return_list$var <- return_list$var + pred$nugs
    }
    
    # Cross covariance.
    if(return_cross_cov) {
      pred_cross <- hetGP:::predict(gp_model[[output_idx]], X_new, xprime=X_cross)
      return_list$cross_cov <- pred_cross$cov
    }
    
    return(return_list)
  }, 
  
  # TODO: I should be careful to ensure X_train is aligned with the data that the homGP 
  # object stores. Currently these could become misaligned if X_train contains duplicates. 
  update_package = function(X_new, y_new, output_idx, update_hyperpar=FALSE, ...) {
    
    if(update_hyperpar) {
      return(hetGP:::update.homGP(gp_model[[output_idx]], Xnew=X_new, Znew=y_new, ...))
    }
    
    return(hetGP:::update.homGP(gp_model[[output_idx]], Xnew=X_new, Znew=y_new, maxit=0, ...))
  }
  
)


# -------------------------------------------------------------------------------
# gpKerGP: Class providing interface between calibration code and kergp package.
#
# Gaussian kernel parameterization:
#    v * exp{(x-y)^2 / ell^2}
# so the 1/2 factor is NOT included and lengthscale ell is on the same scale as 
# the inputs. Marginal variance v is on the squared scale. 
# Source: kergp function `kNormFun`. Somewhat confusingly the function 
# `k1FunGauss` which is used to define `kNormFun` (see the `kGauss ` function)
# does include the 1/2 factor. However, when this function is called from the 
# `kNormFun` function, the parameters (1.0,1.0) are passed, meaning that 
# the C code uses lengthscale/marginal variance parameters of 1.0. These 
# parameters are instead set directly in the R code. 
# -------------------------------------------------------------------------------

gpWrapperKerGP <- setRefClass(
  Class = "gpWrapperKerGP",
  contains = "gpWrapper"
)

gpWrapperKerGP$methods(
  
  initialize = function(X, Y, ...) {
    # X and Y should only be missing when the generator is being called due to 
    # a call to `$copy()`. 
    # TODO: think of better long-term solution. 
    if(missing(X) || missing(Y)) return(NULL)
    
    library(kergp)
    callSuper(X=X, Y=Y, lib="kergp", ...)
    
    # Defining the kernel (covariance function). kergp requires the kernel object 
    # to have names for each input variable that align with names in the X data 
    # when fitting and predicting. Currently, I am using a hack that adds a jitter
    # to the diagonal of the kernel matrix when the kernel matrix is square. 
    # This is not ideal - ideally `kergp` would make it easier to fix a jitter. 
    kernel_map <- function(kernel_name, X_fit, y_fit) {
      
      if(kernel_name == "Gaussian") {
        kernFun <- function(x1, x2, par) {
          K12 <- kergp:::kNormFun(x1, x2, par, k1FunGauss)
          
          # Hack: add jitter to diagonal if number of rows are equal. 
          if(nrow(x1) == nrow(x2)) K12 <- hetGP:::add_diag(K12, rep(default_nugget, nrow(x1)))
          return(K12)
        }
        
        ker <- covMan(kernel = kernFun,
                      hasGrad = TRUE,
                      acceptMatrix = TRUE,
                      d = X_dim,
                      par = c(rep(1, X_dim), 1),
                      parLower = rep(1e-8, X_dim + 1L),
                      parUpper = rep(Inf, X_dim + 1L),
                      parNames = c(X_names, "sigma2"),
                      label = "Gaussian kernel with jitter.")
      } else if(kernel_name == "Matern5_2") {
        .NotYetImplemented()
      } else if(kernel_name == "Matern3_2") {
        .NotYetImplemented()
      } else if(kernel_name == "Gaussian_plus_Quadratic") {
      
        kernFun <- function(x1, x2, par) {
          x1 <- as.matrix(x1)
          x2 <- as.matrix(x2)
          
          # Gaussian part. 
          K12_Gauss <- kergp:::kNormFun(x1, x2, par[1:(length(par)-1)], kergp::k1FunGauss)
          
          # Quadratic part. 
          affine_comb <- tcrossprod(x1, x2) + par[length(par)]
          K12_quad <- affine_comb^2
          attr(K12_quad, "gradient") <- list(quad_offset=2*affine_comb)
          
          # Add kernels. 
          K12 <- K12_Gauss + K12_quad
          attr(K12, "gradient") <- abind(attr(K12_Gauss, "gradient"), cst=attr(K12_quad, "gradient")$quad_offset, along=3)
          
          # Hack: add jitter to diagonal if number of rows are equal. 
          if(nrow(x1) == nrow(x2)) K12 <- hetGP:::add_diag(K12, rep(default_nugget, nrow(x1)))
          return(K12)
        }
        
        # Hyperparameters: lengthscales, marginal variance, additive cst in quadratic kernel.
        ker_par_names <- c(X_names, "sigma2", "quad_offset") 
        y_centered <- abs(y_fit-mean(y_fit))
        max_y <- max(abs(y_centered))
        qy <- quantile(y_centered, 0.2)
        lengthscale_bounds <- get_lengthscale_bounds(X_fit, p=c(0.1), include_one_half=FALSE)
        marg_var_max <- get_marginal_variance_bounds(max_y, p=0.99, return_variance=TRUE)
        marg_var_default <- get_marginal_variance_bounds(qy, p=0.90, return_variance=TRUE)
        par_lower <- setNames(c(lengthscale_bounds["min",], 1e-08, 0), ker_par_names)
        par_upper <- setNames(c(lengthscale_bounds["max",], marg_var_max, 0), ker_par_names)
        par_default <- setNames(c(lengthscale_bounds[3,], marg_var_default, 0), ker_par_names)
        
        ker <- covMan(kernel = kernFun,
                      acceptMatrix = TRUE, 
                      hasGrad = TRUE,
                      d = X_dim,
                      parNames = ker_par_names, 
                      parLower = par_lower,
                      parUpper = par_upper, 
                      label = "Gaussian plus quadratic kernel with jitter", 
                      par = par_default)
      } else {
        stop("gpWrapperKerGP does not support kernel ", kernel_name)
      }

      # Ensure the names attribute used by kergp agrees with the names attribute 
      # used by gpWrapper. 
      inputNames(ker) <- X_names
      
      return(ker)
    }
    
    # Defining the mean functions. 
    mean_func_map <- function(mean_func_name) {
      if(mean_func_name == "constant") {
        mean_func <- as.formula("y ~ 1")
      } else if(mean_func_name == "linear") {
        mean_func <- as.formula(paste0("y ~ ", paste(X_names, collapse=" + ")))
      } else if(mean_func_name == "quadratic") {
        mean_func <- as.formula(paste0("y ~ ", paste0("poly(", X_names, ", degree=2)", collapse=" + ")))
      } else {
        stop("gpWrapperKerGP does not support mean function: ", mean_func_name)
      }
      
      return(mean_func)
    }
    
    initFields(kernel_name_map=kernel_map, mean_func_name_map=mean_func_map) 
  }, 
  
  
  fit_package = function(X_fit, y_fit, kernel_name="Gaussian", mean_func_name="constant", 
                         estimate_nugget=TRUE, fixed_pars=list(), ...) {
    # `kergp` supports a variety of different optimization algorithms, which can be specified 
    # using the `...` arguments. 
    # For kergp, the nugget can be handled more on the prediction side than on the model fitting side. 
    # To fix a jitter, we can fit the model using `noise = FALSE`. Then for the returned GP object 
    # we can set `gp_fit$varNoise <- nugget`. Then when predicting want to use 
    # `predict(..., forceInterp=TRUE)` to include the nugget variance in the kernel matrix. 
    
    # kergp supports fixing known mean function coefficients `beta`. 
    if(isTRUE("beta" %in% names(fixed_pars))) {
      beta <- fixed_pars$beta
    } else {
      beta <- NULL
    }
    
    # TODO: temporarily do not support nugget variance estimation. 
    if(estimate_nugget) {
      stop("gpWrapperKerGP does not currently support estimation of the nugget variance")
    }

    # Fit GP. 
    gp_fit <- kergp::gp(map_mean_func_name(mean_func_name), data=data.frame(y=drop(y_fit), X_fit), 
                        inputs=X_names, cov=map_kernel_name(kernel_name, X_fit, y_fit), 
                        estim=TRUE, beta=beta, noise=FALSE, ...)
    
    # TODO: Commenting this out for now, as we are fixing a jitter within the kernel function. 
    # Set fixed nugget if not estimated.
    # if(!estimate_nugget) gp_fit$varNoise <- default_nugget
    
    return(gp_fit)
  },
  
  predict_package = function(X_new, output_idx, return_mean=TRUE, return_var=TRUE, 
                             return_cov=FALSE, return_cross_cov=FALSE, X_cross=NULL, 
                             include_nugget=TRUE, return_trend=TRUE, pred_type="SK", ...) {
    # For kergp, the predictive standard deviation "sd" includes the nugget variance 
    # if `force_interp` is TRUE; otherwise it is excluded. Specifically, `force_interp=FALSE`
    # returns the predictive dist f|f(X), while `force_interp=TRUE` returns the predictive dist 
    # y|y(X). Thus, `force_interp=TRUE` has two effects: (1) adding the nugget to the diagonal 
    # of the kernel matrix; and (2) adding the nugget to the predictive variance. 
    # If `return_trend` is TRUE, then the returned list will also contain element "trend" 
    # storing the predictions from the GP trend alone; i.e., the prior mean function 
    # evaluations.
    # `pred_type` can either be "SK" (simple kriging) or "UK" (universal kriging). The latter 
    # incorporates uncertainty due to estimation of the GP mean function coefficients. Note 
    # that if using "UK" then the predictive variance will not necessarily stabilize to the 
    # prior variance far away from the design points; in fact, the predictive variance can 
    # explode, causing issues when extrapolation is desired. 
    
    if(return_cross_cov) stop("`return_cross_cov` not yet implemented for gpKerGP class.")
    if(is.null(colnames(X_new))) colnames(X_new) <- X_names
    
    pred <- kergp:::predict.gp(gp_model[[output_idx]], newdata=X_new, 
                               forceInterp=FALSE, seCompute=return_var, 
                               covCompute=return_cov, type=pred_type, ...)
    
    return_list <- list()
    if(return_mean) return_list$mean <- pred$mean
    if(return_var) return_list$var <- (pred$sd)^2
    if(return_cov) return_list$cov <- pred$cov
    if(return_trend) return_list$trend <- pred$trend
    
    # Using the current nugget hack, the nugget variance will be added whenever the 
    # two matrices passed to the kernel have the same number of observations. 
    # TODO: need to check whether varVec includes this nugget or not. For now 
    # assuming it does; thus, the nugget will be added by default. 
    if(!include_nugget) {
      if(return_var) return_list$var <- return_list$var - default_nugget
      if(return_cov) return_list$cov <- hetGP:::fast_diag(return_list$cov, -rep(default_nugget, nrow(cov)))
    }

    return(return_list)
  }
  
)


