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
                 kernel_name_map="list", mean_func_name_map="list")
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
    colnames(X) <<- X_names
    colnames(Y) <<- Y_names
    
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
                     return_cross_cov=FALSE, X_cross=NULL, include_nugget=TRUE, ...) {
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
    
    # If covariance is requested, default to computing cov at inputs `X_new`.
    if(return_cross_cov && is.null(X_cross)) stop("`return_cross_cov` is TRUE but `X_cross` is NULL.") 
      
    # Predict for each independent GP. 
    pred_list <- vector(mode="list", length=Y_dim)
    for(i in seq_along(pred_list)) {
      pred_list[[i]] <- predict_package(X_new, i, return_mean, return_var, return_cov, 
                                        return_cross_cov, X_cross, include_nugget, ...)
    }

    names(pred_list) <- Y_names
    mean_pred <- do.call(cbind, lapply(pred_list, function(l) l$mean))
    var_pred <- do.call(cbind, lapply(pred_list, function(l) l$var))
    cov_pred <- abind(lapply(pred_list, function(l) l$cov), along=3)
    cross_cov_pred <- abind(lapply(pred_list, function(l) l$cross_cov), along=3)
    
    # Set negative variances to 0. 
    if(return_var || return_cov) {
      neg_var_idx <- (var_pred < 0)
      if(any(neg_var_idx)) {
        message("Thresholding negative variance predictions at 0.")
        var_pred[neg_var_idx] <- 0
        if(return_cov) diag(cov_pred)[neg_var_idx] <- 0
      }
    }
    
    # Return outputs to unnormalized scale. 
    if(normalize_output) {
      if(return_mean) mean_pred <- .self$normalize(mean_pred, inverse=TRUE)
      if(return_var) var_pred <- var_pred %*% diag(Y_std^2, nrow=Y_dim)
      if(return_cov) for(i in 1:Y_dim) cov_pred[,,i] <- Y_std[i]^2 * cov_pred[,,i]
      if(return_cross_cov) for(i in 1:Y_dim) cross_cov_pred[,,i] <- Y_std[i]^2 * cross_cov_pred[,,i]
    }
      
    return(list(mean=mean_pred, var=var_pred, cov=cov_pred, cross_cov=cross_cov_pred))
  }, 
  
  
  predict_parallel = function(X_new, return_mean=TRUE, return_var=TRUE, return_cov=FALSE,
                              return_cross_cov=FALSE, X_cross=NULL, include_nugget=TRUE) {
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
      if(is.null(pred_list$chol_cov)) pred_list$chol_cov <- abind(lapply(1:Y_dim, 
                                                                         function(i) t(chol(pred_list$cov[,,i]))), along=3)
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
  
  map_kernel_name = function(kernel_name) {
    return(kernel_name_map[[kernel_name]])
  },
  
  map_mean_func_name = function(mean_func_name) {
    return(mean_func_name_map[[mean_func_name]])
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
    # X and Y should only be missing when the generator is being called due to 
    # a call to `$copy()`. 
    # TODO: think of better long-term solution. 
    if(missing(X) || missing(Y)) return(NULL)
    
    
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


# -----------------------------------------------------------------------------
# gpKerGP: Class providing interface between calibration code and kergp package.  
# -----------------------------------------------------------------------------

gpKerGP <- setRefClass(
  Class = "gpKerGP",
  contains = "gpWrapper"
)

gpKerGP$methods(
  
  initialize = function(X, Y, ...) {
    # X and Y should only be missing when the generator is being called due to 
    # a call to `$copy()`. 
    # TODO: think of better long-term solution. 
    if(missing(X) || missing(Y)) return(NULL)
    
    library(kergp)
    
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




