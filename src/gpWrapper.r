#
# gpWrapper.r
# Class definitions for the `gpWrapper` class, and classes which inherit from 
# this class. All classes defined using reference classes. 
#
# Andrew Roberts
# 
# Depends: general_helper_functions.r, statistical_helper_functions.r,
#          gp_helper_functions.r
#

library(scoringRules)
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
#
# To standardize behavior across different GP packages that are "wrapped", the 
# default assumption here is that the predictive distribution of the model 
# is a GP, with standard predictive mean and covariance functions encoded by 
# the `predict` method. Note that this is NOT the default behavior of many R 
# GP packages out there. For example, packages like kergp and hetGP consider 
# the "ordinary kriging" approach of propagating uncertainty in the mean 
# function coefficients in the predictive variance. Other packages may 
# may marginalize the GP predictive distribution with respect to a marginal 
# variance parameter equipped with an inverse Gamma prior, resulting in a 
# Student-T - not Gaussian - predictive distribution. The gpWrapper class is 
# intended to be flexible enough to accomodate such features, but the default
# behavior should be the basic GP predictive distribution for consistency. 
#
# In particular, `gpWrapper` encapsulates a set of independent univariate GPs.
# This allows modeling a multi-output function. Each of the univariate GPs 
# is assumed to take the form
#
#    y(x) = f(x) + eps
#    f ~ GP(m, k)
#    eps ~ N(0, sig2)
#
# where `m(.)` and `k(.,.)` are the mean and covariance function (i.e., kernel),
# respectively. `sig2` is the "noise variance" for the iid Gaussian noise.
# In certain settings, `f(x)` is not latent; it is observed directly without 
# noise. In this case, the gpWrapper class attribute "noise" will be set to 
# FALSE, which causes `sig2` to be set to a small, near-zero value for 
# numerical stability. We refer to this small value as the "jitter". By default, 
# this value is provided by the "default_jitter", but the user can manually 
# override this default as well. If "noise" is TRUE, then "sig2" is by default
# viewed as another parameter to be estimated. In the gpWrapper class the 
# attribute `noise_var` defines the value of "sig2". "noise_var" is actually 
# a vector of length equal to the output dimension, storing one noise variance 
# per output variable.
#
# The `predict()` method is used to compute the GP predictive distribution 
# at a finite set of new input points. The argument `include_noise` controls 
# whether or not the predictions pertain to `y(x)` or the latent function 
# `f(x)`. This doesn't affect the means but adds `sig2` to predictive 
# variances.
# -----------------------------------------------------------------------------

gpWrapper <- setRefClass(
   Class = "gpWrapper", 
   fields = list(gp_model="ANY", lib="character", X="matrix", Y="matrix", 
                 X_dim="integer", Y_dim="integer", mean_func="ANY",
                 kernel="ANY", noise="logical", noise_var="numeric", 
                 default_jitter="numeric", fixed_pars="list", par_bounds="list",
                 scale_input="logical", normalize_output="logical", 
                 X_bounds="matrix", Y_mean="numeric", Y_std="numeric", 
                 X_names="character", Y_names="character", X_train="matrix", 
                 Y_train="matrix")
)

gpWrapper$methods(
  
  initialize = function(X, Y, scale_input=FALSE, normalize_output=FALSE, 
                        x_names=NULL, y_names=NULL, ...) {
    
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
               default_jitter=default_jitter)

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
  },
  
  set_gp_prior = function(kernel_name="Gaussian", mean_func_name="constant", 
                          include_noise=TRUE, fixed_hyperpars=NULL, 
                          kernel_bounds=NULL, mean_func_bounds=NULL,
                          noise_var_bounds=NULL, jitter=NULL, ...) {
    # This method specifies the GP prior distribution by setting the `kernel` 
    # and `mean_func` attributes, with default hyperparameters (this method does
    # not perform hyperparameter optimization). The kernel and mean function 
    # correspond to a single univariate GP. If there are multiple independent 
    # GPs, then the same kernel and mean function is used for each GP, but their 
    # hyperparameters will of course be optimized separately. Also, they 
    # are each assigned a separate noise variance in the `noise_var`
    # attribute. The specific type of the `kernel` and `mean_func` attributes 
    # will depend on the specific gpWrapper class. The methods 
    # `define_kernel()` and `define_mean_func()` are defined specially for 
    # each class that inherits from base gpWrapper. `fixed_hyperpars` is a 
    # list with elements "mean_func", "kernel", and "noise_var". The former two 
    # are themselves lists, with named elements that will depend on the specific 
    # gpWrapper class. "noise_var" is a numeric value or NULL, which can be 
    # used to specify a known noise variance. 
    # `par_bounds`, is a list with the same element names. The "kernel" and 
    # "mean_func" elements are each lists length `dim_Y` with names set to the 
    # corresponding output names `Y_names". Each element is a matrix with 
    # 4 columns with column names "par_name", "init", "lower", and "upper". 
    # "par_name" is the hyperparameter name and the remaining columns are used 
    # to specify initial values, lower bounds, and upper bounds. 
    # `par_bounds$noise_var` is a matrix with columns "y_name", "init", "lower", 
    # "upper" used to specify the analogous quantities for the noise variance of 
    # each output.
    
    # Set `fixed_pars` attribute.
    if(is.null(fixed_hyperpars)) {
      fixed_hyperpars <- list(kernel=list(), mean_func=list(), noise_var=NULL)
    } else {
      assert_that(all(c("kernel", "mean_func", 
                        "noise_var") %in% names(fixed_hyperpars)))
      if(length(noise_var)==1) {
        fixed_hyperpars$noise_var <- rep(fixed_hyperpars$noise_var, .self$Y_dim)
      }
    }
    
    assert_that(length(fixed_hyperpars$noise_var) %in% c(0, .self$Y_dim))
    fixed_pars <<- fixed_hyperpars
    
    # Whether or not noise is to be included, set a default jitter value.
    if(is.null(jitter)) {
      jitter <- sqrt(.Machine$double.eps)
    } else {
      jitter <- drop(jitter)
      assert_that(is.numeric(jitter) && (length(jitter)==1))
    }
    default_jitter <<- jitter
    
    noise <<- include_noise
    if(!include_noise) {
      # If not including noise term, then noise variance will be fixed at 
      # `default_jitter`.
      noise_var <<- rep(jitter, .self$Y_dim)
    } else if(include_noise && !is.null(fixed_pars$noise_var)) {
      # If noise is to be included, but variance parameter is fixed (not 
      # estimated).
      noise_var <<- .self$fixed_pars$noise_var
    } else {
      # If estimating noise variances.
      noise_var <<- rep(NA_real_, .self$Y_dim)
    }
    names(noise_var) <<- .self$Y_names
    
    # Set up mean and covariance functions (no hyperparameter optimization 
    # is performed here). 
    kernel <<- .self$define_kernel(kernel_name, ...)
    mean_func <<- .self$define_mean_func(mean_func_name, ...)
    
    # Set `par_bounds` attribute, which controls bounds/initialization values
    # for mean function and kernel hyperparameters, as well as the noise 
    # variance parameter.
    par_bounds <<- .self$define_par_bounds(kernel_bounds=kernel_bounds, 
                                           mean_func_bounds=mean_func_bounds,
                                           noise_var_bounds=noise_var_bounds,
                                           ...)
  },
  
  
  define_kernel = function(kernel_name, ...) {
    # This function is intended to be called from `set_gp_prior()`, and 
    # thus assumes that certain attributes (e.g., `fixed_pars`, `noise`)
    # have already been set. Must return a list with at minimum elements 
    # "name" (the kernel name), and "par_names" (a vector of names of the 
    # hyperparameters defining the kernel). "par_names" is used to set the 
    # official order of the kernel parameter names.

    stop("`define_kernel()` is implemented by each class that inherits from ",
         "the base gpWrapper class.")
  },
  
  define_mean_func = function(mean_func_name, ...) {
    # This function is intended to be called from `set_gp_prior()`, and 
    # thus assumes that certain attributes (e.g., `fixed_pars`, `noise`)
    # have already been set.
    
    stop("`define_mean_func()` is implemented by each class that inherits ",
         "from the base gpWrapper class.")
  },
  
  define_par_bounds = function(kernel_bounds=NULL, mean_func_bounds=NULL, 
                               noise_var_bounds=NULL, ...) {
    # Controls bounds/initialization values for kernel and mean function 
    # hyperparameters, as well as the noise variance parameter. See 
    # `get_noise_var_bounds()` and `get_mean_func_bounds()` for information
    # on the required formats for the arguments to this method.
    #
    # NOTE: get_mean_func_bounds() is not yet implemented, as neither 
    #       hetGP nor kergp support bounds or init values for the mean function
    #       parameters.
    #
    # The default method returns NULL, which does not set any bounds/init 
    # values. This method should typically be overloaded by classes inheriting 
    # from `gpWrapper`.
    
    noise_var_bounds <- .self$get_noise_var_bounds(noise_var_bounds, ...)
    kernel_bounds <- .self$get_kernel_bounds(kernel_bounds, ...)
    
    list(kernel=kernel_bounds, mean_func=NULL, noise_var=noise_var_bounds)
  },
  
  get_noise_var_bounds = function(noise_var_bounds=NULL, ...) {
    # Returns a matrix with rownames set to `Y_names` and colnames 
    # "init", "lower", "upper". Used to specify initialization values and 
    # lower/upper bounds for the noise variance of each 
    # output. Not set if `.self$noise` if FALSE (in which case the jitter 
    # will be used). This default method may be overloaded for certain classes
    # that inherit from gpWrapper.
    
    # In noiseless observation setting, no need to define bounds.
    if(!.self$noise) return(NULL)

    row_names <- .self$Y_names
    col_names <- c("init", "lower", "upper")
    
    # Set up empty matrix if missing.
    if(is.null(noise_var_bounds)) {
      noise_var_bounds <- matrix(nrow=length(row_names), ncol=length(col_names),
                                 dimnames=list(row_names, col_names))
    } else {
      assert_that(is.matrix(noise_var_bounds))
    }
    
    # Fill in missing columns.
    missing_cols <- setdiff(col_names, colnames(noise_var_bounds))
    if(length(missing_cols) > 0) {
      noise_var_bounds <- cbind(noise_var_bounds, 
                                matrix(NA, nrow=nrow(noise_var_bounds), 
                                      ncol=length(missing_cols),
                                      dimnames=list(row_names, missing_cols)))
      noise_var_bounds <- noise_var_bounds[,col_names]
    }
    
    # Fill in missing rows.
    missing_rows <- setdiff(row_names, rownames(noise_var_bounds))
    if(length(missing_rows) > 0) {
      new_mat <- matrix(NA, nrow=length(missing_rows), ncol=length(col_names),
                        dimnames=list(missing_rows, col_names))
      noise_var_bounds <- rbind(noise_var_bounds, new_mat)[row_names,]
    }
    
    # Fill in missing values with defaults derived from observed variability 
    # in design outputs.
    if(any(is.na(noise_var_bounds))) {
      defaults <- matrix(nrow=nrow(noise_var_bounds), 
                         ncol=ncol(noise_var_bounds), 
                         dimnames=list(row_names, col_names))
      y_train_sd <- apply(.self$Y_train, 2, sd)
      defaults[,"lower"] <- .self$default_jitter
      defaults[,"upper"] <- y_train_sd^2
      defaults[,"init"] <- (0.3 * y_train_sd)^2
      
      noise_var_bounds[is.na(noise_var_bounds)] <- defaults[is.na(noise_var_bounds)]
    }
  
    assert_that(all(as.vector(noise_var_bounds[,c("lower","upper","init")]) >= 0))
    return(noise_var_bounds[row_names, col_names])
  },
  
  get_kernel_bounds = function(kernel_bounds=NULL, ...) {
    # This method is intended to be called from `define_par_bounds()`. It 
    # accepts NULL or a list, with names set to a subset of `Y_names`. Each 
    # element of this list must be a matrix with rownames set to the kernel 
    # parameter names and colnames set to "init", "lower", and "upper". The 
    # argument can be used to specify only a subset subset of `Y_names`. Can 
    # also be used to specify bounds for only a subset of kernel 
    # hyperparameters, or only a subset of the columns for the bounds matrix. 
    # In either case, the bounds matrix need only contain the rows and/or 
    # columns necessary to set the desired bounds (with the row and column 
    # names always required). Missing rows and/or columns will be automatically 
    # added and filled with NAs. This method is intended to check that any 
    # specified information is in the right format, and fill in any information 
    # not provided with NAs. The list is then passed one element at a time 
    # (i.e., one output dimension at a time) to `get_default_kernel_bounds()`, 
    # which is a method that is defined on a class-by-class basis. The 
    # `get_kernel_bounds()` method may be overwritten for certain classes with 
    # special structure. 
    
    col_names <- c("init", "lower", "upper")
    row_names <- .self$kernel$par_names
    
    # Set up the `kernel_bounds` list. 
    if(is.null(kernel_bounds)) kernel_bounds <- list()
    assert_that(is.list(kernel_bounds))
    if(length(kernel_bounds) > 1) {
      assert_that(!is.null(names(kernel_bounds)))
      assert_that(all(names(kernel_bounds) %in% .self$Y_names))
    }
    
    # For each output, set up the matrix storing the bounds/init values, and 
    # call `get_default_kernel_bounds()` to populate any missing values in 
    # this matrix.
    for(j in 1:.self$Y_dim) {
      y_name <- .self$Y_names[j]
      mat <- kernel_bounds[[y_name]]
      
      # Set up empty matrix if missing.
      if(is.null(mat)) {
        mat <- matrix(nrow=length(row_names), ncol=length(col_names),
                      dimnames=list(row_names, col_names))
      } else {
        assert_that(is.matrix(mat))
      }
      
      # Fill in missing columns.
      missing_cols <- setdiff(col_names, colnames(mat))
      if(length(missing_cols) > 0) {
        mat <- cbind(mat, matrix(NA, nrow=nrow(mat), 
                                 ncol=length(missing_cols),
                                 dimnames=list(row_names, missing_cols)))
        mat <- mat[,col_names]
      }
      
      # Fill in missing rows.
      missing_pars <- setdiff(row_names, rownames(mat))
      if(length(missing_pars) > 0) {
        new_mat <- matrix(NA, nrow=length(missing_pars), ncol=length(col_names),
                          dimnames=list(missing_pars, col_names))
        mat <- rbind(mat, new_mat)[row_names,]
      }
      
      # Call `get_default_kernel_bounds()` to fill in missing bounds/init vals.
      kernel_bounds[[y_name]] <- .self$get_default_kernel_bounds(mat, j, ...)
    }
    
    return(kernel_bounds)
  },
  
  get_default_kernel_bounds = function(bounds_mat, output_idx, ...) {
    # Specified on a class by class basis. By default just returns 
    # `bounds_mat` as is.
    bounds_mat
  },
  
  gp_is_fit = function() {
    # TRUE if the GP model attribute `gp_model` is initialized. 
    class(.self$gp_model) != "uninitializedField"
  },
  
  prior_is_defined = function() {
    # TRUE if the attributes `mean_func`, `kernel`, `noise`, and `noise_var` 
    # are initialized. Note that the `noise_var` vector may contains NAs if the 
    # noise variances have not been fixed nor estimated yet, but the attribute
    # is still defined by `set_gp_prior()`.
    !any(c(class(.self$mean_func),
           class(.self$kernel),
           class(.self$noise),
           class(.self$noise_var)) %in% "uninitializedField")
  },
  
  get_summary_str = function(...) {
    summary_str <- "gpWrapper summary:\n"

    # GP prior summary.
    summary_str <- paste0(summary_str, .self$summarize_gp_prior())

    # Design information.
    summary_str <- paste0(summary_str, "\n-----> Design information:\n")
    summary_str <- paste0(summary_str, 
                          "X dimension: ", .self$X_dim, 
                          "\nY dimension: ", .self$Y_dim, 
                          "\nX names: ", paste0(.self$X_names, collapse=", "),
                          "\nY names: ", paste0(.self$Y_names, collapse=", "),
                          "\nNumber design points: ", nrow(.self$X),
                          "\nScale input: ", .self$scale_input, 
                          "\nNormalize output: ", .self$normalize_output, "\n")
    
    # Package-specific information. If the gpWrapper package-specific class 
    # has not implemented the method `get_summary_str_package()` then the 
    # default is for this method to just return NULL, which will not be printed
    # by `cat()`. 
    if(.self$gp_is_fit()) {
      for(i in seq_len(.self$Y_dim)) {
        
        # Package-specific summary.
        summary_str <- paste0(summary_str, "\n", 
                              .self$get_summary_str_package(i,...))
        
        # Noise variances are standardized across all packages.
        summary_str <- paste0(summary_str, "\n>>> Noise:\n")
        summary_str <- paste0(summary_str, "Std dev: ", 
                              sqrt(.self$noise_var[i]), "\n")
      }
    }
    
    return(summary_str)
  },
  
  summarize_gp_prior = function() {
    
    str <- "\n-----> GP Prior Summary:\n"
    
    # If prior is not defined, then there is nothing to summarize.
    if(!.self$prior_is_defined()) {
      str <- paste0(str, "GP prior defined: FALSE\n")
      return(str)
    }
    
    # Noise/jitter.
    str <- paste0(str, "GP prior defined: TRUE\n")
    str <- paste0(str, "Include noise: ", .self$noise, "\n")
    str <- paste0(str, "Default jitter: ", .self$default_jitter, "\n")
    
    if(class(.self$noise_var) == "uninitializedField") {
      str <- paste0(str, "Noise variance: not initialized")
    } else {
      str <- paste0(str, "Noise variance(s): ", 
                    paste(.self$noise_var, collapse=", "), "\n")
    }
    
    if(!.self$noise) {
      str <- paste0(str, "\t> Note: noise variance(s) fixed at `default_jitter` ",
                    "for numerical stability.\n")
    }
    
    # Mean function and kernel name.
    str <- paste0(str, "Mean function name: ", .self$mean_func$name, "\n")
    str <- paste0(str, "Kernel name: ", .self$kernel$name, "\n")
    
    # Fit/hyperparameters optimized. 
    gp_obj_exists <- .self$gp_is_fit()
    str <- paste0(str, "GP package object defined: ", gp_obj_exists, "\n")
    
    return(str)    
  },
  
  get_summary_str_package = function(...) {
    # Intended to be implemented by package-specific gpWrapper classes.
    NULL
  },
  
  summarize = function(...) {
    # Note that this function can be overwritten by package-specific gpWrapper
    # classes so that package-specific information is printed in addition to 
    # the generic gpWrapper info. 
    
    cat(.self$get_summary_str(...))
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
    if(any(.self$Y_std == 0)) {
      stop("`Y_std` is 0, so cannot normalize.")
    }
    
    if(inverse) {
      Ynew <- Ynew * matrix(Y_std, nrow=nrow(Ynew), ncol=Y_dim, byrow=TRUE) + 
              matrix(Y_mean, nrow=nrow(Ynew), ncol=Y_dim, byrow=TRUE)
    } else {
      Ynew <- (Ynew - matrix(Y_mean, nrow=nrow(Ynew), ncol=Y_dim, byrow=TRUE)) / 
              matrix(Y_std, nrow=nrow(Ynew), ncol=Y_dim, byrow=TRUE)
    }
    
    return(Ynew)
  },
  

  fit = function(...) {
    assert_that(.self$prior_is_defined(), 
                msg=paste("GP prior not defined. ",
                          "Call `set_gp_prior()` method prior to `fit()`.")) 
    
    fits_list <- vector(mode="list", length=.self$Y_dim)
    for(i in seq_along(fits_list)) {
      fits_list[[i]] <- .self$fit_package(X_fit=X_train, y_fit=Y_train[,i],
                                          output_idx=i, ...)
    }

    names(fits_list) <- .self$Y_names
    gp_model <<- fits_list
  }, 
  
  fit_package = function(X_fit, y_fit, output_idx, ...) {
    err_msg <- "A fit_package() method must be implemented for each class inheriting from 
                gpWrapper. gpWrapper does not implement its own fit_package() method. The
                fit() method should ultimately call fit_package()."
    stop(err_msg)
  },

  fit_parallel = function(...) {
    .NotYetImplemented()
  },
  
  predict_package = function(X_new, output_idx, return_mean=TRUE, 
                             return_var=TRUE, return_cov=FALSE, 
                             return_cross_cov=NULL, X_cross=NULL, 
                             include_noise=TRUE) {
    err_msg <- "A predict_package() method must be implemented for each class inheriting from 
                gpWrapper. gpWrapper does not implement its own predict_package() method. The
                predict() method should ultimately call predict_package()."
    stop(err_msg)
  }, 
  
  predict = function(X_new, return_mean=TRUE, return_var=TRUE, return_cov=FALSE, 
                     return_cross_cov=FALSE, X_cross=NULL, include_noise=TRUE, 
                     return_trend=TRUE, ...) {
    # Logic for all predict() functions:
    #   - `return_cov` refers to k(X_new, X_new) while `return_cross_cov` refers 
    #      to k(X_new, X_cross).
    #   - The former is always diagonal, and `include_noise==TRUE` will cause 
    #     noise variances to be added to the diagonal. The latter may or may 
    #     not be diagonal, and `include_noise` has no effect on the cross cov 
    #     matrix (noise variances will NOT be added in either case). 
    #   - If `return_cov==TRUE`, then the variances will always be returned as 
    #     well, by simply taking the diagonal of the covariance matrix. 
    #   - `return_cross_cov` has no effect on the other settings; can be thought 
    #      of as an optional add-on.
    #
    # This default predict method is overwritten in some special cases; e.g., 
    # `gpWrapperSum.` 
    #
    # Output dimensions (let m := nrow(X_new) and q := nrow(X_cross)):
    # "mean": matrix, (m, Y_dim)
    # "var": matrix, (m, Y_dim)
    # "cov": array, (m, m, Y_dim)
    # "cross_cov": array, (m, q, Y_dim).

    # If there is only one input, no covariance computation required.
    if(isTRUE(nrow(X_new)==1)) return_cov <- FALSE
    
    # Scale inputs, if required. 
    if(scale_input) {
      X_new <- .self$scale(X_new)
      if(return_cross_cov && !is.null(X_cross)) X_cross <- .self$scale(X_cross)
    }

    # If covariance is requested, default to computing scov at inputs `X_new`.
    if(return_cross_cov && is.null(X_cross)) {
      stop("`return_cross_cov` is TRUE but `X_cross` is NULL.") 
    }
      
    # Predict for each independent GP. 
    pred_list <- vector(mode="list", length=Y_dim)
    for(i in seq_along(pred_list)) {
      pred_list[[i]] <- .self$predict_package(X_new=X_new, output_idx=i, 
                                        return_mean=return_mean, 
                                        return_var=return_var, 
                                        return_cov=return_cov, 
                                        return_cross_cov=return_cross_cov, 
                                        X_cross=X_cross, 
                                        include_noise=include_noise, 
                                        return_trend=return_trend, ...)
    }

    names(pred_list) <- .self$Y_names
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
                              return_cross_cov=FALSE, X_cross=NULL, include_noise=TRUE,
                              return_trend=TRUE) {
    .NotYetImplemented()
  },
  
  
  sample = function(X_new, use_cov=FALSE, include_noise=TRUE, N_samp=1, 
                    pred_list=NULL, adjustment="none", bounds=c(0,Inf), ...) {
    # If `pred_list` is passed, it should have all the required components. 
    # Returns array of dimension (num input, num samp, Y_dim). The "adjustment"
    # argument can be used to truncate the Gaussian predictive distribution
    # to satisfy the bounds given by `bounds`. `adjustment = "truncated"` will 
    # sample from a truncated Gaussian, while `adjustment = "rectified"` will 
    # sample from a rectified Gaussian (i.e., a sample exceeding the upper bound
    # will be set to the upper bound, and likewise for the lower bound).

    n_input <- nrow(X_new)
    if(n_input < 2) use_cov <- FALSE
    assert_that(length(bounds)==2L)
    
    # Compute required predictive quantities if not already provided. 
    if(is.null(pred_list)) {
      pred_list <- .self$predict(X_new, return_mean=TRUE, return_var=!use_cov,  
                                 return_cov=use_cov, include_noise=include_noise)
    } else {
      if(use_cov) assert_that(!is.null(pred_list$cov))
    }

    # Compute lower Cholesky factors of the predictive covariance matrices. 
    # If not using predictive cov, Cholesky factors are diagonal with standard 
    # devs on diag. For truncated normal, `tmvnorm` doesn't accept Cholesky 
    # factor, so just pass covariance matrix. 
    if((adjustment=="truncated") && !use_cov) {
      pred_list$cov <- abind(lapply(1:Y_dim, function(i) diag(pred_list$var[,i],
                                                              nrow=n_input)), along=3)
    } else if((adjustment != "truncated") && use_cov) {
      if(is.null(pred_list$chol_cov)) pred_list$chol_cov <- abind(lapply(1:Y_dim, function(i) t(chol(pred_list$cov[,,i]))), along=3)
    } else if(adjustment != "truncated") {
      pred_list$chol_cov <- abind(lapply(1:Y_dim, function(i) diag(sqrt(pred_list$var[,i]), nrow=n_input)), along=3)
    }
    
    samp <- array(dim=c(nrow(X_new), N_samp, Y_dim))
    for(i in 1:Y_dim) {
      if(adjustment=="truncated") { # Zero left-truncated Gaussian. 
        samp[,,i] <- t(rtmvnorm(N_samp, mean=pred_list$mean[,i], 
                                sigma=pred_list$cov[,,i], 
                                lower=rep(bounds[1], nrow(X_new)),
                                upper=rep(bounds[2], nrow(X_new))))
      } else {
        samp[,,i] <- sample_Gaussian_chol(pred_list$mean[,i],
                                          as.matrix(pred_list$chol_cov[,,i]), 
                                          N_samp)
        if(adjustment=="rectified") {
          samp[,,i] <- pmax(samp[,,i], bounds[1])
          samp[,,i] <- pmin(samp[,,i], bounds[2])
        }
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

  calc_pred_func = function(func, type="pw", X_new=NULL, Y_new=NULL,  
                            pred_list=NULL, include_noise=TRUE,  
                            return_mean=TRUE, return_var=TRUE, 
                            return_cov=FALSE, ...) {
    # Returns the output of a function applied to `pred_list`. The function 
    # is applied separately for each output GP. The signature of the function 
    # should look like `func(pred_list, output_idx, ...)` or 
    # `func(pred_list, output_idx, y, ...)`, with output index specifying which 
    # of the output GPs the function should be applied to. The second 
    # signature is typically used for error measures, where `y` represents 
    # the true values. Alternatively, `func` may be a character in which 
    # case this method looks for a method called `calc_pw_<func>()`.
    # `type` can be "pw" (pointwise) or "agg" (aggregate). The former implies 
    # that `func` should operate on inputs point-by-point and thus return 
    # a vector of length `nrow(X_new)`. "agg" doesn't do any checks on the 
    # output size, but this typically means that `func` returns a single 
    # value. 

    # This is necessary to make these methods available in the current 
    # environment, since they are not explicitly included in the function body.
    usingMethods(calc_pw_entropy, calc_pw_crps, calc_pw_mae, calc_pw_mse, 
                 calc_pw_log_score, calc_agg_mah, calc_agg_log_score)
    
    # Select method to use.
    assert_that(is.function(func) || is.character(func))
    if(is.character(func)) {
      method_name <- paste("calc", type, func, sep="_")
      assert_that(method_name %in% names(.self), 
                  msg=paste0("Method not found: ", method_name))
      func <- .self[[method_name]]
    }
    
    # Compute required predictive quantities if not already provided. 
    if(is.null(pred_list)) {
      assert_that(return_mean || return_var || return_cov)
      pred_list <- .self$predict(X_new, return_mean=return_mean, 
                                 return_var=return_var, 
                                 include_noise=include_noise)
    }
    
    # Identify the number of test points. Only required for pointwise function 
    # evaluations to check that the output size is correct.
    n_test <- NULL
    if(type == "pw") {
      if(!is.null(pred_list$mean)) n_test <- nrow(pred_list$mean)
      else if(!is.null(pred_list$var)) n_test <- nrow(pred_list$var)
      else n_test <- dim(pred_list$cov)[1]
    }

    # Each independent GP is handled separately.
    l <- lapply(seq_len(.self$Y_dim), function(i) drop(func(pred_list, i, Y_new[,i], ...)))
    if(type == "pw") assert_that(all(sapply(l, length) == n_test))
    func_mat <- do.call(cbind, l)
    colnames(func_mat) <- .self$Y_names

    return(func_mat)
  },
  
  calc_pred_multi_func = function(func_list, type="pw", X_new=NULL, Y_new=NULL,  
                                  pred_list=NULL, return_list=FALSE, 
                                  include_noise=TRUE, return_mean=TRUE, 
                                  return_var=TRUE, return_cov=FALSE, ...) {
    # A wrapper around `calc_pred_func()` that allows `func` to be a list of 
    # functions (or method names). Note that `type` should not be a list; 
    # all of the functions should have the same return shape. If `return_list`
    # is TRUE, then the return type is a list of length `length(func_list)`. 
    # Otherwise, the list elements are combined into a single data.table with 
    # column names "<Y_names>", "func", where `<Y_names>` includes one column 
    # per output variable.
    
    # Compute required predictive quantities if not already provided. 
    if(is.null(pred_list)) {
      pred_list <- .self$predict(X_new, return_mean=return_mean, 
                                 return_var=return_var, 
                                 include_noise=include_noise)
    }
    
    # Return function evaluations for each function in list.
    l <- lapply(func_list, function(f) .self$calc_pred_func(f, type=type,  
                                                            X_new=X_new, Y_new=Y_new, 
                                                            pred_list=pred_list, ...))
  
    # Obtain names identifying each function.
    if(is.character(func_list)) func_names <- func_list
    else if(!is.null(names(func_list))) func_names <- names(func_list)
    else func_names <- paste0("f", seq_along(func_list))
    
    # Return as list, if requested.
    if(return_list) {
      names(l) <- func_names
      return(l)
    }

    # Otherwise stack into a single data.table.
    for(i in seq_along(l)) {
      l[[i]] <- as.data.table(l[[i]])
      l[[i]][, func := func_names[i]]
    }
    
    return(rbindlist(l, use.names=TRUE))
  },
  
  calc_pw_entropy = function(pred_list, output_idx, ...) {
    0.5 * log(2*pi*exp(1)*pred_list$var[,output_idx])
  },
  
  calc_pw_crps = function(pred_list, output_idx, y, ...) {
    assert_that(!is.null(y))
    scoringRules::crps_norm(y, mean=pred_list$mean[,output_idx], 
                            sd=sqrt(pred_list$var[,output_idx]))
  },
  
  calc_pw_mae = function(pred_list, output_idx, y, ...) {
    assert_that(!is.null(y))
    abs(pred_list$mean[,output_idx] - y)
  },
  
  calc_pw_mse = function(pred_list, output_idx, y, ...) {
    assert_that(!is.null(y))
    (pred_list$mean[,output_idx] - y)^2
  },
  
  calc_pw_log_score = function(pred_list, output_idx, y, ...) {
    assert_that(!is.null(y))
    scoringRules::logs_norm(y, mean=pred_list$mean[,output_idx],
                            sd=sqrt(pred_list$var[,output_idx]))
  },
  
  calc_agg_mah = function(pred_list, output_idx, y, ...) {
    assert_that(!is.null(y))
    stats::mahalanobis(y, center=pred_list$mean[,output_idx],
                       cov=pred_list$cov[,,output_idx])
  },
  
  calc_agg_log_score = function(pred_list, output_idx, y, ...) {
    # Note that this function is the multivariate version of 
    # `calc_pw_log_score()`, which does not consider the covariance.
    assert_that(!is.null(y))
    mvtnorm::dmvnorm(y, mean=pred_list$mean[,output_idx], 
                     sigma=pred_list$cov[,,output_idx], log=TRUE)
  },
  
  plot_pred_1d = function(X_new, include_noise=TRUE, include_interval=TRUE, 
                          interval_method="pm_std_dev", N_std_dev=1, CI_prob=0.9, 
                          pred_list=NULL, Y_new=NULL, plot_title=NULL, 
                          xlab=NULL, ylab=NULL, ...) {
    
    assert_that(X_dim==1, 
                msg=paste0("plot_pred_1d() requires 1d input space. X_dim = ", X_dim))
    
    # Compute required predictive quantities if not already provided. 
    if(is.null(pred_list)) {
      pred_list <- .self$predict(X_new, return_mean=TRUE, return_var=TRUE, 
                                 include_noise=include_noise)
    }
    
    # Default x-axis label.
    if(is.null(xlab)) xlab <- X_names

    plts <- vector(mode="list", length=Y_dim)
    for(i in 1:Y_dim) {
      
      # Default y-axis label. 
      if(is.null(ylab)) ylab <- Y_names[i]
      
      plts[[i]] <- plot_Gaussian_pred_1d(X_new=X_new[,1], pred_mean=pred_list$mean[,i], 
                                         pred_var=pred_list$var[,i], 
                                         include_interval=include_interval, 
                                         interval_method=interval_method, 
                                         N_std_dev=N_std_dev, CI_prob=CI_prob, 
                                         y_new=Y_new[,i], X_design=X[,1], 
                                         y_design=Y[,i], xlab=xlab, ylab=ylab, 
                                         plot_title=plot_title, ...)
    }
    
    return(plts)
    
  },
  
  plot_pred_2d = function(X_new, include_noise=TRUE, pred_list=NULL, Y_new=NULL) {
    # Returns a list of heatmaps for the predictive mean and standard deviation. 
    
    assert_that(X_dim==2, msg=paste0("plot_pred_2d() requires 1d input space. X_dim = ", X_dim))
    stop("gpWrapper$plot_pred_2d not yet implemented.")
  },
  
  
  plot_resid_hist = function(X_new, Y_new, include_noise=TRUE, CI_prob=0.9, 
                             normalize_resid=TRUE, pred_list=NULL, ...) {
    # If `norm_resid = TRUE` then the normalized reiduals (y-mean)/sd are 
    # plotted. Otherwise y-mean is plotted.
    # Compute required predictive quantities if not already provided. 
    # TODO: finish this method.
    
    if(!normalize_resid) .NotYetImplemented()
    
    if(is.null(pred_list)) {
      pred_list <- .self$predict(X_new, return_mean=TRUE, return_var=TRUE, 
                                 include_noise=include_noise)
    }
    
    CI_tail_prob <- 0.5 * (1-CI_prob)
    plts <- vector(mode="list", length=.self$Y_dim)
    
    for(i in 1:.self$Y_dim) {
      plot_title <- paste0("GP predictions: ", .self$Y_names[i])
      resids <- Y_new[,i] - pred_list$mean[,i]
      vars <- pred_list$var[,i]
      idx <- (vars > .self$default_jitter)
      resids <- resids[idx]
      vars <- vars[idx]
      resids <- resids / sqrt(vars)
      
      plts[[i]] <- ggplot(data.frame(x=resids)) + 
                    geom_histogram(aes(x=x), ...)
    }
    
    return(plts)
  },

    
  plot_pred = function(X_new, Y_new, include_CI=FALSE, include_noise=TRUE, 
                       CI_prob=0.9, pred_list=NULL, ...) {

    # Compute required predictive quantities if not already provided. 
    if(is.null(pred_list)) {
      pred_list <- .self$predict(X_new, return_mean=TRUE, return_var=include_CI, 
                                 include_noise=include_noise)
    }
    
    CI_tail_prob <- 0.5 * (1-CI_prob)
    plts <- vector(mode="list", length=.self$Y_dim)
    
    for(i in 1:.self$Y_dim) {
      plot_title <- paste0("GP predictions: ", .self$Y_names[i])
      y_mean <- pred_list$mean[,i]
      
      if(include_CI) {
        y_sd <- sqrt(pred_list$var[,i])
        CI_lower <- qnorm(CI_tail_prob, y_mean, y_sd)
        CI_upper <- qnorm(CI_tail_prob, y_mean, y_sd, lower.tail=FALSE)
      } else {
        CI_lower <- CI_upper <- NULL
      }
      
      plts[[i]] <- plot_true_pred_scatter(y_pred=y_mean, y_true=Y_new[,i], 
                                          include_CI=include_CI, 
                                          CI_lower=CI_lower,  
                                          CI_upper=CI_upper, 
                                          plot_title=plot_title, 
                                          xlab="Observed", ylab="Predicted", ...)
    }
    
    return(plts)
  },
  
  
  plot_1d_projection = function(x_names=NULL, n_points=100L, X_list=NULL, X_fixed=NULL, 
                                include_design=FALSE, include_interval=TRUE, 
                                interval_method="pm_std_dev", N_std_dev=1, CI_prob=0.9, 
                                include_noise=TRUE, ...) {
    # Plots GP predictions as a single input variable varies, with the remaining variables 
    # fixed at a single value. The default behavior of the method allows it to be called 
    # with no arguments; e.g., `plot_1d_projection()`. In this case, projections will be 
    # considered on all input dimensions. For each variable being varied, 
    # the remaining variables will by default be fixed at their midpoints with respect to 
    # the `X_bounds` class attribute. GP predictions will be computed with these variables fixed
    # and the variable "x1" varied at `n_points` (default 100) equally spaced points between 
    # its bound values in `X_bounds`. Optional arguments offer much more fine-tuned control; 
    # e.g. `X_list` controls the points at which the varied variables are plotted, and `X_fixed`
    # controls the values at which the fixed variables are fixed. One variable may be 
    # varied per plot, but multiple variables can be specified in `x_names`, in which multiple 
    # plots will be produced. Multiple combinations of fixed variables can also be specified 
    # in `X_fixed`, which will be plotted on the same plot. 
    #
    # Args:
    #    x_names: character, a subset of the `X_names` field determining variables to vary 
    #             (they will be varied one at a time, a different plot being produced for each).
    #             If NULL, set to all variables. 
    #    n_points: integer, the number of grid points at which predictions will be computed for 
    #              the varied variables. This is overwritten for variables with input points 
    #              explicitly provided in `X`.
    #    X_list: list, optionally provides the input points at which to compute predictions for 
    #            the varied variables. Names of the list elements should correspond to `x_names`. 
    #            Elements of the list should be numeric vectors or matrices with a single column. 
    #            Can contain elements corresponding to a strict subset of `x_names`, in which case
    #            the default input grid will be used for the remaining variables.
    #    X_fixed: matrix, determining the values at which the non-varied inputs will be fixed. 
    #             The columns of this matrix correspond to the input variables, and the 
    #             `colnames` attribute must be set to the corresponding variable names. Each  
    #             row of `X_fixed` will contribute one line to each produced plot; i.e., this 
    #             gives a way to investigate the effect of fixing different values for the 
    #             non-varied variables. Note that the code will select the columns corresponding
    #             to the fixed variables by name as needed, so `X_fixed` technically only requires
    #             columns for variables that will be used as fixed variables at some point. It is 
    #             typically easiest just to include all variables in this matrix. 
    #    include_design: logical, whether or not to plot the design points (projected onto 
    #                    the relevant dimension).
    #    include_interval, interval_method, N_std_dev, CI_prob: all used to plot a confidence
    #               interval in addition to mean predictions. Same behavior as in other plotting
    #               functions. 
    #
    # Returns:
    #    List of plots, containing `Y_dim` * `length(x_names)` plots in total; i.e., one plot 
    #    per varied variable per independent GP. The returned list has multiple levels, the first 
    #    being a list of length `Y_dim` with names set to the output variable names. Each of these
    #    elements contains a sublist of length `length(x_name)` corresponding to the projection 
    #    dimension (the variable being varied). Each plot will contain `nrow(X_fixed)` lines.
    #
    # TODO: would probably be nice to have helper function(s) that construct input grids that 
    # only vary one variable.
    # TODO: need to generalize code to allow `X_fixed` to have multiple rows; a helper function
    # will probably be best here as well.
    # TODO: allow ability to pass in true values.
    
    # Determine the variables that will be varied. 
    if(is.null(x_names)) {
      x_names <- .self$X_names
    } else {
      x_names <- unique(x_names)
      assert_that(all(is.element(x_names, .self$X_names)))
    }
    
    # Constructing input grids for varied parameters.
    n_vary <- length(x_names)
    if(is.null(X_list)) X_list <- list()
    for(x_name in x_names) {
      if(is.null(X_list[[x_name]])) {
        x_bounds <- .self$X_bounds[,x_name]
        X_list[[x_name]] <- seq(x_bounds[1], x_bounds[2], length.out=n_points)
      }
    }
    
    # If not provided, set values of fixed parameters to midpoints with respect to 
    # design bounds.
    if(is.null(X_fixed)) {
      X_fixed <- matrix(apply(.self$X_bounds, 2, mean), nrow=1)
      colnames(X_fixed) <- colnames(.self$X_bounds)
    } else {
      assert_that(nrow(X_fixed)==1, 
                  msg="Currently only supports the case where `X_fixed` has a single row.")
    }
    
    # Determine whether or not design points are included in plot.
    if(include_design) X_design <- .self$X
    else X_design <- NULL
    
    # Create plots.
    plts <- list()
    n_fixed <- .self$X_dim - 1L
    for(i in seq_len(n_vary)) {
      # Input matrix for prediction.
      x_name <- x_names[i]
      x_names_fixed <- setdiff(.self$X_names, x_name)
      x_grid <- matrix(X_list[[x_name]], ncol=1, dimnames=list(NULL,x_name))
      n_grid <- nrow(x_grid)
      X_pred <- cbind(x_grid, matrix(X_fixed[,x_names_fixed,drop=FALSE], ncol=n_fixed, 
                                     nrow=n_grid, byrow=TRUE, dimnames=list(NULL,x_names_fixed)))
      X_pred <- X_pred[, .self$X_names]
        
      # Compute GP predictions.
      pred_list <- .self$predict(X_pred, return_var=include_interval, include_noise=include_noise, ...)
      
      # Construct plots.
      for(j in seq_len(.self$Y_dim)) {
        y_name <- .self$Y_names[j]
        plts[[y_name]][[x_name]] <- plot_Gaussian_pred_1d(drop(x_grid), pred_mean=pred_list$mean[,j], 
                                                          pred_var=pred_list$var[,j],
                                                          include_interval=include_interval, 
                                                          interval_method=interval_method, N_std_dev=N_std_dev,
                                                          CI_prob=CI_prob, X_design=X_design[,x_name], 
                                                          y_design=.self$Y[,j], xlab=x_name, ylab=y_name)
      }
    }
    
    return(plts)
  },
  
  plot_loocv = function() {
    .NotYetImplemented()
  },
  
  sample_Gaussian_chol = function(mu, C_chol_lower, N_samp=1) {
    drop(mu) + C_chol_lower %*% matrix(rnorm(nrow(C_chol_lower)*N_samp), ncol=N_samp)
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
  
  calc_expected_exp_cond_var = function(X_eval, X_cond, include_noise=TRUE, log_scale=TRUE, ...) {
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
    # response at the evaluation locations. This leaves the predictive 
    # mean unchanged, but updates the predictive variance.
    gp_copy$update(X_cond, pred$mean, update_hyperpar=FALSE, ...)
    
    # Predict with the conditional ("cond") GP (i.e., the updated GP) at the  
    # evaluation locations. Convert to the exponentiated scale to obtain 
    # log-normal predictive quantities.
    pred_cond <- gp_copy$predict(X_eval, return_mean=TRUE, return_var=TRUE, ...)
    lnp_cond_log_var <- matrix(nrow=N_eval, ncol=.self$Y_dim)
    for(j in seq_len(.self$Y_dim)) {
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


# ------------------------------------------------------------------------------
# gpWrapperSum
# 
# This class defines a multi-output GP constructed via a weighted sum of 
# non-random basis vectors, where the weights are modeled as independent GPs.
# i.e., it defines a model of the form G: R^d -> R^p by 
#
# G(u) = g + sum_{j=1}^{r} w_j(u)b_j + eps =: g + Bw(u)
# w_j ~ GP(mu_j, k_j)
# eps ~ N(0, sig2 * I)
#
# where all of the `w_j` and `eps` are pairwise independent, `g` and 
# `b_1, ..., b_r` are constant vectors in R^p. `eps` is a zero-mean Gaussian 
# noise term with constant (not dependent on `u`) variance `sig2`. `B` is 
# the matrix with columns set to the `b_j`, and `w(u)` is the column vector 
# with entries `w_j(u)`. `w(u)` is a multi-output GP consisting of independent 
# GPs with `w ~ GP(mu, k)`, with `mu` and `k` constructed from the `mu_j` and 
# `k_j`, respectively. Note that this model implies that `G(u)` is a 
# multi-output GP G ~ GP(m, c), where
# 
# m(u) = E[G(u)] = g + B*mu(u) = g + sum_{j=1}^{r} mu_j(u) b_j
# c(u,v) = Cov[G(u),G(v)] = B*k(u,v)*B^T = sum_{j=1}^{r} k_j(u,v) b_j b_j^T
#
# This class is implemented by inheriting from `gpWrapper`, but it does not 
# implement a `fit()` method. Instead, the GP `w(u)` is assumed to have already 
# been fit via some `gpWrapper` class. The fit GP is then passed to 
# `gpWrapperSum`. This class essentially acts as a wrapper around a set of 
# independent GPs, and its `predict()` method leverages the predictions from 
# the underlying GPs to construct `m(u)` and `c(u,v)`. The class fields
# `gp_model`, `B`, shift`, and `sig2` corresponds to `w`, `B`, `g`, and `sig2` 
# respectively, in the above notation. 
#
# Note that the primary reason this class inherits from `gpWrapper` is to allow 
# it to be compatible with downstream algorithms that expect a `gpWrapper` 
# object.
#
# See this blog post for more details on 
# this basis function GP model:
# https://arob5.github.io/blog/2024/06/25/basis-func-emulator/
# ------------------------------------------------------------------------------

gpWrapperSum <- setRefClass(
  Class = "gpWrapperSum",
  contains="gpWrapper",
  fields=list(B="matrix", shift="numeric", sig2="numeric")
)

gpWrapperSum$methods(
  
  initialize = function(gp, B, shift=NULL, sig2=NULL) {
    # `gp` should only be missing when the generator is being called due to 
    # a call to `$copy()`. 
    # TODO: think of better long-term solution. 
    if(missing(gp) || missing(B)) return(NULL)
    
    # Argument checking.
    assert_that(inherits(gp, "gpWrapper"),
                msg="`gp` must inherit from `gpWrapper` class.")
    assert_that(is.matrix(B))
    
    if(!is.null(shift)) {
      assert_that(nrow(B) == length(drop(shift)))
    }
    
    if(!is.null(sig2)) assert_that(sig2 > 0)
      
    # Initialize fields not inherited from `gpWrapper`, or not set by 
    # `gpWrapper` constructor.
    initFields(gp_model=gp, B=B, shift=drop(shift), sig2=sig2)
      
    # Output dimension corresponds to the dimension of `B`, not to be confused 
    # with `gp$Y_dim`, which is the dimension of `w(u)`. Here we define the 
    # design set (X,Y) by transforming the design of the underlying GP object.
    # These matrices don't really need to be computed, so could improve 
    # efficiency in the future by not storing them, and only computing them 
    # if needed.
    out_dim <- nrow(B)
    B_names <- paste0("b", 1:out_dim)
    Y_design <- tcrossprod(.self$gp_model$Y, .self$B)
    if(!is.null(.self$shift)) Y_design <- add_vec_to_mat_rows(.self$shift, Y_design)
    
    # Call `gpWrapper` constructor so that `Y_names` and `Y_dim` will be 
    # correctly set to the output dimension of the basis vectors.
    callSuper(X=.self$gp_model$X, Y=Y_design, normalize_output=FALSE, 
              scale_input=FALSE, x_names=.self$gp_model$X_names, 
              y_names=B_names)
  }, 
  
  predict = function(X_new, return_mean=TRUE, return_var=TRUE, return_cov=FALSE, 
                     return_cross_cov=FALSE, X_cross=NULL, include_noise=TRUE, 
                     return_trend=TRUE, pred_list=NULL, include_output_cov=FALSE, 
                     ...) {
    # `pred_list` is a GP prediction list corresponding to `gp_model`, the 
    # set of independent GPs. Given that the gpWrapperSum GP is a multi-output 
    # GP, there are two notions of covariance here: covariance across input 
    # values, and covariance across output values. By default, this method 
    # excludes the latter. This is done for two reasons: (1) to align with the 
    # behavior of `predict()` for other gpWrapper classes; and (2) to avoid 
    # computing with very large matrices when not required. The different 
    # return options are summarized below. Let `M := nrow(X_new)`, 
    # `Q := nrow(X_cross)`, and `P := dim_Y`.
    #
    # Returns:
    # list, with different potential elements:
    # "mean": included if `return_mean = TRUE`. matrix of shape (M,P), where 
    #         the (m,p) element is E[G(x_m)_p].
    # "trend": included if `return_trend = TRUE`. Same shape as "mean". The 
    #          GP prior mean predictions.
    # "var": included if `return_var = TRUE`. matrix of shape (M,P), where the 
    #        (m,p) element is Var[G(x_m)_p].
    #
    # The elements "cov" and "cross_cov" will depend on the value of 
    # `include_output_cov`. At present, no output cross-covariance computation
    # is supported.
    #
    # If `include_output_cov = FALSE` (the default): 
    # "cov": included if `return_cov = TRUE`. matrix of shape (M,M,P)
    # "cross_cov": included if `return_cross_cov = TRUE`. matrix of shape 
    #              (M,Q,P).
    #
    # If `include_output_cov = TRUE`:
    # "cov": included if `return_cov = TRUE`. matrix of shape (M,P,P), where the
    #        (m,,) element contains the covariance matrix Cov[G(x_m)].
    # "cross_cov": not yet supported.
    
    n_input <- nrow(X_new)
    if((n_input==1) && return_cov) return_var <- TRUE
    
    return_list <- list()
    
    # Predict using independent GPs, if required.
    if(is.null(pred_list)) {
      pred_list <- .self$gp_model$predict(X_new, return_mean=return_mean, 
                                          return_var=return_var, 
                                          return_cov=return_cov,
                                          return_cross_cov=return_cross_cov,
                                          X_cross=X_cross, 
                                          include_noise=include_noise,
                                          return_trend=return_trend, ...)
    }
    
    # Compute mean predictions.
    if(return_mean) {
      mu <- tcrossprod(pred_list$mean, .self$B)
      if(!is.null(.self$shift)) mu <- add_vec_to_mat_rows(.self$shift, mu)
      return_list$mean <- mu
    }
    
    # Compute trend (prior mean) predictions.
    if(return_trend) {
      trend <- tcrossprod(pred_list$trend, .self$B)
      if(!is.null(.self$shift)) trend <- add_vec_to_mat_rows(.self$shift, trend)
      return_list$trend <- trend
    }
    
    # Compute variance predictions.
    if(return_var) {
      vars <- tcrossprod(pred_list$var, (.self$B)^2)
      if(!is.null(.self$sig2)) vars <- vars + sig2
      return_list$var <- vars
    }
    
    # Compute covariances and cross-covariances. Predictive covariances will 
    # either be with respect to the output space or input space.
    if(return_cov || return_cross_cov) {
      if(include_output_cov) {
        return_list_covs <- .self$predict_output_cov(X_new, return_cov=FALSE, 
                                                     return_cross_cov=FALSE, 
                                                     X_cross=NULL, 
                                                     include_noise=TRUE,
                                                     pred_list=pred_list)
      } else {
        return_list_covs <- .self$predict_input_cov(X_new, return_cov=FALSE, 
                                                    return_cross_cov=FALSE, 
                                                    X_cross=NULL, 
                                                    include_noise=TRUE,
                                                    pred_list=pred_list)
      }
      
      return_list <- c(return_list, return_list_covs)
    }
    
    return(return_list)
  }, 
  
  predict_input_cov = function(X_new, return_cov=FALSE, return_cross_cov=FALSE, 
                               X_cross=NULL, include_noise=TRUE, 
                               pred_list=NULL, ...) {
    # Compute predictive covariances and cross-covariances, where 
    # these quantities are interpreted as covariances with respect to the 
    # input ("x") variables. No output covariances are computed here. See 
    # `predict()` for specifics on the computed quantities.

    n_input <- nrow(X_new)
    return_var <- ifelse((n_input==1) && return_cov, TRUE, FALSE)
    return_list <- list()
    
    # Predict using independent GPs, if required.
    if(is.null(pred_list)) {
      pred_list <- .self$gp_model$predict(X_new, return_mean=FALSE, 
                                          return_var=return_var, 
                                          return_cov=return_cov,
                                          return_cross_cov=return_cross_cov,
                                          X_cross=X_cross, 
                                          include_noise=include_noise, ...)
    }
    
    # Compute covariance predictions.
    if(return_cov) {
      covs <- array(0, dim=c(n_input, n_input, .self$dim_Y))
      
      if(n_input==1) {
        covs[1,1,] <- return_list$var[1,]
      } else {
        for(p in seq_len(.self$dim_Y)) {
          for(j in 1:nrow(.self$B)) {
            cov[,,p] <- cov[,,p] + (.self$B[p,j])^2 * pred_list$cov[,,j,drop=FALSE]
          }
        }
      }
      
      return_list$cov <- covs
    }
    
    # Compute cross-covariance predictions.
    if(return_cross_cov) {
      cross_covs <- array(0, dim=c(n_input, nrow(X_cross), .self$dim_Y))
      for(p in seq_len(.self$dim_Y)) {
        for(j in 1:nrow(.self$B)) {
          cross_covs[,,p] <- cross_covs[,,p] + (.self$B[p,j])^2 * pred_list$cross_cov[,,j,drop=FALSE]
        }
      }
      
      return_list$cross_cov <- cross_covs
    }
    
    return(return_list)
  },
  
  
  predict_output_cov = function(X_new, return_cov=FALSE, return_cross_cov=FALSE, 
                                X_cross=NULL, include_noise=TRUE, 
                                pred_list=NULL, ...) {
    # Compute predictive covariances and cross-covariances, where 
    # these quantities are interpreted as covariances with respect to the 
    # output variables. No input covariances are computed here. See 
    # `predict()` for specifics on the computed quantities. At present, this 
    # method does not support computation of cross-covariances. For covariances,
    # returns array of shape (M,P,P), where M=nrow(X_new) and P=Y_dim.
    # The (i,,) element of this array contains the PxP covariance matrix 
    # Cov[G(x_i), G(x_i)], where x_i is the ith row of X_new.
    
    return_list <- list()
    
    # Predict using independent GPs, if required.
    if(is.null(pred_list)) {
      pred_list <- .self$gp_model$predict(X_new, return_mean=FALSE, 
                                          return_var=return_cov, 
                                          return_cov=FALSE,
                                          return_cross_cov=return_cross_cov,
                                          X_cross=X_cross, 
                                          include_noise=include_noise, ...)
    }
    
    # Loop input-by-input to compute the (cross) covariance matrices.
    n_input <- nrow(X_new)
    
    if(return_cov) {
      cov_arr <- array(NA, dim=c(n_input, .self$Y_dim, .self$Y_dim))
      for(i in seq_len(n_input)) {
        # Compute covariance predictions.
        if(return_cov) {
          cov_arr[i,,] <- tcrossprod(mult_vec_with_mat_rows(pred_list$var[i,], .self$B), .self$B)
          diag(cov_arr[i,,]) <- diag(cov_arr[i,,]) + .self$sig2
        }
      }
      return_list$cov <- cov_arr
    }

    # Compute cross-covariance predictions.
    if(return_cross_cov) {
      message("gpWrapperSum does not currently support cross-covariances ",
              "across output variables.")
    }
    
    return(return_list)
  },
  
  
  sample = function(X_new, use_cov=FALSE, include_noise=TRUE, N_samp=1, 
                    pred_list=NULL, ...) {
    # Returns samples from G(x) at different input points x. This simply 
    # requires sampling the underlying independent GPs w(x) and using these 
    # weight samples to compute linear combinations of the basis vectors.
    # Note that `use_cov` here refers to covariance across different inputs `x`;
    # output covariances are always included.
    #
    # Returns array of dimension (nrow(X_new), N_samp, Y_dim). 
    
    # Sample weights (independent GPs).
    w_samp <- .self$gp_model$sample(X_new, use_cov=use_cov, 
                                    include_noise=include_noise,
                                    N_samp=N_samp, pred_list=pred_list, ...)
    n_input <- dim(w_samp)[1]

    # Compute linear combination of weight samples to produce samples from G.
    G_samp <- array(NA, dim=c(n_input, N_samp, .self$Y_dim))
    for(i in 1:n_input) {
      W <- w_samp[i,,]
      G <- tcrossprod(W, .self$B) # Works for N_samp=1 case when W is numeric().
      
      if(!is.null(.self$shift)) G <- add_vec_to_mat_rows(.self$shift, G)
      if(!is.null(.self$sig2)) {
        Z <- matrix(rnorm(prod(dim(G)), sd=sqrt(.self$sig2)), 
                    nrow=nrow(G), ncol=ncol(G))
        G <- G + Z
      }
      
      G_samp[i,,] <- G
    }
    
    return(G_samp)
  }, 
  
  get_summary_str = function(...) {
    
    # gpWrapperSum attributes.
    summary_str <- "----- gpWrapperSum -----\n"
    summary_str <- paste0(summary_str, "Output dimension: ", .self$Y_dim, "\n")
    summary_str <- paste0(summary_str, "Number independent GPs: ", ncol(.self$B), "\n")
    summary_str <- paste0(summary_str, "Include shift: ", !is.null(.self$shift), "\n")
    summary_str <- paste0(summary_str, "Include noise: ", !is.null(.self$sig2), "\n")
    
    # Append attributes for underlying GPs.
    summary_str <- paste0(summary_str, 
                          "\n----- Underlying independent GPs -----\n\n")
    summary_str <- paste0(summary_str, .self$gp_model$get_summary_str(...))
  }, 
  
  scale = function(Xnew, inverse=FALSE) {
    # Call the underlying gpWrapper scale method.
    
    .self$emulator_model$scale(Xnew, inverse=inverse)
  }, 
  
  normalize = function(Ynew, inverse=FALSE) {
    # Call the underlying gpWrapper normalize method.
    
    .self$emulator_model$normalize(Ynew, inverse=inverse)
  }, 
  
  fit = function(...) {
    stop("gpWrapperSum does not implement its own `fit()` method.")
  },
  
  update = function(X_new, Y_new, update_hyperpar=FALSE, ...) {
    # Call the underlying gpWrapper normalize method.
    
    .self$emulator_model$update(X_new, Y_new, update_hyperpar=update_hyperpar)
  }, 
  
  calc_pred_func = function(...) {
    .NotYetImplemented()
  }, 
  
  plot_pred = function(X_new, Y_new=NULL, include_CI=TRUE, include_noise=TRUE, 
                       CI_method="CI", CI_prob=0.9, N_std_dev=1, 
                       pred_list=NULL, g_true=NULL, g_func=NULL, ...) {
    # This method is overloaded to plot the predictive distribution of G(x) as 
    # a function of the output index; i.e., plot pairs of the form (p, G(x)_p),
    # where `x` is a fixed input. In particular, plots the mean 
    # (p, E[G(x)_p]). Optionally also plots the true trajectory (p, G*(x)_p)
    # for comparison, as well as a credible interval around (p, E[G(x)_p]).
    # Only allows plotting for a single input `x`, so requires
    # `nrow(X_new)==1`. "CI_method", "CI_prob", and "N_std_dev" are passed 
    # to `plot_Gaussian_plot_1d()` to control the type of confidence interval
    # that is plotted. `g_true` can optionally be passed as a numeric vector or 
    # one-row/column matrix, in which case the trajectory `g_true`
    # will also be plotted. Alternatively, a function `g_func` can be provided 
    # and the "true" trajectory will be computed by calling `g_func(X_new)`.
    
    assert_that(nrow(X_new)==1, msg="plot_pred() requires `nrow(X_new)=1`")
    
    pred_list_g <- .self$predict(X_new, return_var=include_CI, 
                                 include_noise=include_noise, 
                                 pred_list=pred_list, ...)
    
    # Determine if "true" trajectory will be plotted.
    if(!is.null(g_true) || !is.null(g_func)) {
      if(is.null(g_true)) g_true <- g_func(X_new)
      assert_that(length(drop(g_true)) == .self$Y_dim)
    }
    
    plot_Gaussian_pred_1d(seq_len(.self$Y_dim), drop(pred_list_g$mean), 
                          pred_var=drop(pred_list_g$var),
                          y_new=g_true, include_interval=include_CI, 
                          interval_method=CI_method,
                          N_std_dev=N_std_dev, CI_prob=CI_prob,
                          xlab="output index", ylab="G", ...)
  },
  
  plot_pred_lin_proj = function() {
    .NotYetImplemented()
  },
  
  plot_pred_output = function(output_idx) {
    .NotYetImplemented()
  },
  
  plot_pred_mean = function() {
    .NotYetImplemented()
  },
  
  plot_samp = function(X_new, use_cov=FALSE, include_noise=TRUE, N_samp=1,
                       pred_list=NULL, g_true=NULL, g_func=NULL, ...) {
    # A helper function that wraps around the `sample()` method to plot sample 
    # trajectories of G(x) at different x values given in `X_new`. This method 
    # operates in one of two different modes. `N_samp` is the number of samp
    # per input. It is only allowed to be greater than 1 if `nrow(X_new)=1`.
    # In this case, the plot contains samples from G(x) at a single fixed input 
    # value `x`. If `nrow(X_new) > 1` then `N_samp` must be 1, in which case 
    # the function plots one sample each from `G(x_1), ..., G(x_M)`. `g_true`
    # can optionally be passed as a numeric vector or one-row/column matrix 
    # in the case that `nrow(X_new) = 1`, in which case the trajectory `g_true`
    # will also be plotted. Alternatively, a function `g_func` can be provided 
    # and the "true" trajectory will be computed by calling `g_func(X_new)`.
    
    assert_that(!((nrow(X_new)>1) && (N_samp>1)), 
                msg="`plot_samp()` does not allow `nrow(Xnew) > 1` and `N_samp > 1`.")
    
    G_samp <- .self$sample(X_new, use_cov=use_cov, N_samp=N_samp, 
                           pred_list=pred_list)
    
    # Prepare input for passing to `plot_curves_1d_helper()`.
    if(N_samp > 1) {
      G_samp <- G_samp[1,,]
      title <- "G(x) samples (single input)"
      
      if(!is.null(g_true) || !is.null(g_func)) {
        if(is.null(g_true)) g_true <- g_func(X_new)
        assert_that(length(drop(g_true)) == .self$Y_dim)
      }
    } else {
      G_samp <- G_samp[,1,]
      title <- "G(x) samples (one sample per input)"
      g_true <- NULL
    }
    
    if(!is.matrix(G_samp)) G_samp <- matrix(G_samp, nrow=1)
    G_samp <- t(G_samp)
    
    plot_curves_1d_helper(seq_len(.self$Y_dim), G_samp, df_by=NULL, 
                          plot_title=title, xlab="output index", 
                          ylab="G", y_new=g_true, legend=FALSE)
  }
  
)



# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# gpWrapperHet: Class providing interface between calibration code and hetGP 
#               package.
#
# hetGP focuses on GPs with standard stationary kernels (Gaussian and Metérn),
# and only allows a constant mean function. It offers methods for fast GP 
# updates (i.e., fast conditioning on new points). hetGP also offers 
# functionality for GP modeling with a heteroskedastic noise term, though 
# this functionality is not currently supported by gpWrapperHet.
# For the purposes of achieving a closed-form estimate of the marginal variance, 
# hetGP assumes a GP model of the form:
#
# y(x) = f(x) + eps
# f ~ GP(beta0, nu_hat * k(.,.))
# eps ~ N(0, nu_hat * g)
#
# Thus, the noise variance in this case is actually `noise_var := nu_hat * g`.
# `nu_hat` is the marginal variance, and `beta0` is a mean constant. In the 
# noiseless setting, this means that we cannot simply fix `noise_var` to 
# `jitter_default`, as the noise variance is also a function of the unknown 
# parameter `nu_hat`. The approach taken here is to instead set 
# `g := jitter_default`, which is reasonable as long as `nu_hat` is not to 
# big. Then after hyperparameter optimization we set the attribute 
# `noise_var := nu_hat + g = nu_hat + jitter_default`, using the 
# optimized marginal variance `nu_hat`.
#
# A second difference from the standard gpWrapper assumptions is that, if 
# `beta0` is estimated, then hetGP sets a "trendtype" parameter to "OK"
# (ordinary kriging), which inflates the GP predictive variance at prediction 
# time to account for the uncertainty in the `beta0` estimation. The default
# behavior of gpWrapper is not to include such variance inflation. Thus, 
# the `fit()` method of this class has to use a bit of a hack in setting 
# "trendtype := SK" (simple kriging - no variance inflation) post-hoc after 
# hyperparameter optimization. It would be better to change this pre-MLE but it 
# doesn't appear that hetGP provides an easy way to do this. This is not ideal 
# since it means that hetGP uses "OK" during hyperparameter optimization, and 
# then "SK" for prediction. If `beta0` is not estimated, then this is not an 
# issue. 
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
    callSuper(X=X, Y=Y, lib="hetGP", ...)
  }, 
  
  define_kernel = function(kernel_name, ...) {
    # For gpWrapperHet, we define the `kernel` attribute to be a list with 
    # element "name" (the kernel name, using the hetGP conventions).

    l <- list()
    
    # Map to hetGP kernel name.
    if(kernel_name == "Gaussian") l$name <- "Gaussian"
    else if(kernel_name == "Materm5_2") l$name <- "Matern5_2"
    else if(kernel_name == "Matern3_2") l$name <- "Matern3_2"
    else stop("gpWrapperHet does not support kernel: ", kernel_name)

    return(l)
  },
  
  define_mean_func = function(mean_func_name, ...) {
    # gpWrapperHet defines the attribute `mean_func` to be a list with 
    # element "name" set to the name of the mean function. hetGP only allows 
    # a constant mean, so an error is thrown for any other mean function.
    
    l <- list()
    
    if(mean_func_name == "constant") l$name <- mean_func_name 
    else stop("gpWrapperHet does not support mean function: ", mean_func_name)
    
    return(l)
  },
  
  define_par_bounds = function(kernel_bounds=NULL, mean_func_bounds=NULL, 
                               noise_var_bounds=NULL, ...) {
    # hetGP has a method for choosing default bounds for hyperparameters; 
    # if `hyperpar_bounds` is not explicitly provided, we fall back on this 
    # hetGP method.
    
    if(is.null(hyperpar_bounds)) {
      message("No bounds/initialization values explicitly defined for ",
              "kernel and mean function hyperparameters. Falling back on ",
              "hetGP's method for defining bounds/init values.")
      return(NULL)
    } else {
      message("gpWrapperHet does not yet support explicit setting of ",
              " hyperparameter bounds. hetGP's default bounds will be used.")
      return(NULL)
    }
  },
  
  fit_package = function(X_fit, y_fit, output_idx, ...) {
    # We opt not to define bounds on the kernel hyperparameters here, as hetGP's 
    # `auto_bounds()` function tends to work pretty well. If desired, explicit 
    # bounds can be passed by specifying "lower" and "upper" arguments in `...`.
    
    # Convert `fixed_pars` to hetGP format.
    hetgp_fixed_pars <- c(.self$fixed_pars$kernel, .self$fixed_pars$mean_func)
    
    # Deal with noiseless case. See class description for why this slightly 
    # differs from the typical approach.
    if(!.self$noise) {
      hetgp_fixed_pars$g <- .self$default_jitter
    }
    
    if(length(hetgp_fixed_pars) == 0) hetgp_fixed_pars <- NULL
    
    gp_fit <- hetGP:::mleHomGP(X_fit, y_fit, covtype=.self$kernel$name, 
                               known=hetgp_fixed_pars, ...)
    
    # Change the trendtype from "ordinary kriging" to "simple kriging" to align 
    # with the default gpWrapper behavior of not inflating the predictive 
    # variance due to mean function estimation.
    gp_fit$trendtype <- "SK"
    
    # Update the noise variance attribute based on the hyperparameter 
    # optimization. See class comments for notes on the parameterization used 
    # here.
    .self$noise_var[output_idx] <<- gp_fit$nu_hat * gp_fit$g
    
    return(gp_fit)
  },
  
  predict_package = function(X_new, output_idx, return_mean=TRUE, 
                             return_var=TRUE, return_cov=FALSE, 
                             return_cross_cov=FALSE, X_cross=NULL, 
                             include_noise=TRUE, ...) {
    
    if(return_cov) X_prime <- X_new
    else X_prime <- NULL

    pred <- hetGP:::predict(gp_model[[output_idx]], X_new, xprime=X_prime)
    return_list <- list()
    if(return_mean) return_list$mean <- pred$mean
    if(return_cov) {
      return_list$cov <- pred$cov
      if(include_noise) diag(return_list$cov) <- diag(return_list$cov) + pred$nugs
      return_list$var <- diag(return_list$cov)
    } else if(return_var) {
      return_list$var <- pred$sd2
      if(include_noise) return_list$var <- return_list$var + pred$nugs
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
  }, 
  
  get_summary_str_package = function(idx, ...) {
    # This method returns a string summarizing the GP stored in
    # `.self$gp_model[[idx]]`. The `gpWrapper` method `get_summary_str()`
    # loops over the GP list to assemble these strings for all GPs.
    if(!.self$gp_is_fit()) return(NULL)
    gp_het <- .self$gp_model[[idx]]

    summary_str <- paste0("\n-----> gpWrapperHet GP ", idx, "; output = ", .self$Y_names[idx], ":\n")

    # Mean function.
    summary_str <- paste0(summary_str, "\n>>> Mean function:\n")
    summary_str <- paste0(summary_str, "Mean constant: ", gp_het$beta0, "\n")
    summary_str <- paste0(summary_str, "Trend type: ", gp_het$trendtype, "\n")

    # Kernel.
    summary_str <- paste0(summary_str, "\n>>> Kernel:\n")
    summary_str <- paste0(summary_str, "Lengthscales:\n")
    ls <- sqrt(gp_het$theta)
    for(i in seq_along(ls)) {
      summary_str <- paste0(summary_str, "\t", names(ls)[i], ": ", ls[i], "\n")
    }
    summary_str <- paste0(summary_str, "Marginal std dev: ", sqrt(gp_het$nu_hat), "\n")

    return(summary_str)
  }, 
  
  plot_corr_summary = function(idx, types="default") {
    # Options for `types` are "default" (using the lengthscale value 
    # is actually used in the correlation function - either the default 
    # or the optimized value), "lower" (use lower bound lengthscale), 
    # "upper" (use upper bound lengthscale).
    
    .NotYetImplemented()
  },
  
  get_lengthscale_summary = function(idx) {
    # Lengthscale values are returned on the normalized scale [0,1], such that 
    # 0 corresponds to `X_bounds[0,idx]` and 1 corresponds to 
    # `X_bounds[1,idx]` (the minimum and maximum design point value in input 
    # dimension `idx`). Thus, values outside of [0,1] indicate 
    .NotYetImplemented()
  }
  
)


# -------------------------------------------------------------------------------
# gpWrapperKerGP: Class providing interface between calibration code and kergp 
#                 package.
#
# The kergp provides access to a large number of kernels, including the ability 
# to write user-defined kernels and combine kernels via the typical algebraic
# operations. It also allows for linear models to be specified for the GP prior 
# mean using R formulas, analogously to the way it is done in `lm()`. The 
# flexibility in mean function/kernel specification is more that I am aware of 
# any other R package offering, but there are some quirks to be aware of, which 
# are discussed below.
#
# Hyperparameter optimization:
# The optimization approach allows interval bound constraints for the kernel 
# hyperparameters and noise variance, but not for the coefficients of the 
# linear model defining the prior mean. The reason for this is that the 
# optimization is with respect to the concentrated log-likelihood, where the 
# mean coefficients have been "concentrated out". Conditional on the kernel 
# hyperparameters, the mean coefficients can be optimized in closed-form, and 
# their esimtate simply takes the form of a generalized least squares model. 
# This is nice, unless you want to specify constraints on the coefficients.
# Also, it does not appear that kergp provides the ability to fix a known 
# noise variance (e.g., a small jitter) during the optimization. To get around 
# this, this class sets lower/upper bounds on the `varNoise` parameter to the 
# desired fixed value, which seems to work fine. The optimized `varNoise` 
# parameter is always added to the diagonal of the kernel matrix.
# Unfortunately, kergp currently provides no way to fix kernel hyperparameters 
# during the optimization. The documentation lists an argument "parFixed" to 
# the `mle` method, but notes that it is not yet implemented. The mean function
# coefficients can be fixed.
# 
# Prediction:
# The behavior fo the predict method is a bit odd to me. The method has an 
# argument "forceInterp" that does two things:
# (1.) adds `varNoise` to predictive variances.
# (2.) adds a "white-noise" kernel to the cross covariance k(X,Xnew).
#
# The white noise kernel takes the value `varNoise` when its inputs are the 
# same (within a small threshold) and otherwise returns zero. I don't understand
# why we would want to do this, or why the `forceInterp` argument does both 
# of these things together (the documentation says this argument will likely 
# be removed in the future). To avoid dealing with this, this class always 
# sets `forceInterp = FALSE` and optionally adds the noise variance post-hoc 
# to control whether predictive variances are computed for f(x) or y(x).
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
  }, 
  
  define_kernel = function(kernel_name, ...) {
    # The `kernel` attribute is defined to be a list with elements "name" and 
    # "object", the latter being the kergp kernel object of class "covMan".
    # Some kernel names have been standardized, so that a string can be passed 
    # for "kernel_name" and the kernel will be automatically created.
    # Alternatively, `kernel_name` can be a "covMan" kernel object, 
    # in which case `kernel$name` will be set to "user-specified kernel".
    # Note that kergp requires the kernel object to have names for each input 
    # variable that align with names in the X data when fitting and predicting.
    # Hyperparameters for a `covMan` object "obj" can be extracted via 
    # `coef(obj)`. Hyperparameter names can be accessed via 
    # `attr(obj, "kernParNames")`See `attributes(obj)` for other class 
    # attributes (parameter bounds, etc.). kergp only enforces the most basic 
    # bounds for hyperparameters by default (e.g., non-negativity for 
    # lengthscale parameters), so setting bounds is often quite important.
    
    if(kernel_name == "Gaussian") {
      # Hyperparameters are the lengthscales (on per input dimension), and 
      # the marginal variance.
      ker <- kergp::kGauss(d=.self$X_dim)
      attr(ker, "kernParNames") <- c(paste0("ell_", .self$X_names), "marg_var")
    } else if(kernel_name == "Quadratic") {
      # Hyperparameters are "quad_center" and "quad_offset", which are the 
      # constants `c` and `a` in the expression (<x-a,z-a> + c)^2, 
      # respectively. `c` is a scalar and `a` is a vector of length `X_dim`.
      
      #
      # TODO: below is incomplete; need to finish.
      #
      .NotYetImplemented()
      
      # No pre-defined kergp kernel of this type, so define our own.
      kernFun <- function(x1, x2, par) {
        x1 <- as.matrix(x1)
        x2 <- as.matrix(x2)
        a <- par[1:(length(par)-1)]

        affine_comb <- tcrossprod(add_vec_to_mat_rows(-a, x1), 
                       add_vec_to_mat_rows(-a, x2)) + par[length(par)]
        K12_quad <- affine_comb^2
        attr(K12_quad, "gradient") <- list(quad_center=-1, 
                                           quad_offset=2*affine_comb)
 
      }
      
    } else if(kernel_name == "Gaussian_plus_Quadratic") {
      # Hyperparameters are the lengthscales (one per input dimension), the 
      # the marginal variance, and "quad_offset", which is the constant c in 
      # the expression (<x,z> + c)^2.
      
      # No pre-defined kergp kernel of this type, so define our own.
      kernFun <- function(x1, x2, par) {
        x1 <- as.matrix(x1)
        x2 <- as.matrix(x2)
        
        # Gaussian part. 
        K12_Gauss <- kergp:::kNormFun(x1, x2, par[1:(length(par)-1)], 
                                      kergp::k1FunGauss)
        
        # Quadratic part. 
        affine_comb <- tcrossprod(x1, x2) + par[length(par)]
        K12_quad <- affine_comb^2
        attr(K12_quad, "gradient") <- list(quad_offset=2*affine_comb)
        
        # Add kernels. 
        K12 <- K12_Gauss + K12_quad
        attr(K12, "gradient") <- abind(attr(K12_Gauss, "gradient"), 
                                       cst=attr(K12_quad, "gradient")$quad_offset, 
                                       along=3)
        return(K12)
      }
      
      ker <- covMan(kernel = kernFun,
                    acceptMatrix = TRUE, 
                    hasGrad = TRUE,
                    d = .self$X_dim,
                    parNames = c(paste0("ell_", .self$X_names), 
                                 "marg_var", "quad_offset"), 
                    label = "Gaussian plus quadratic kernel")
    } else if(inherits(kernel_name, "covMan")) {
      ker <- kernel_name
      kernel_name <- "user-specified kernel"
    } else {
      stop("gpWrapperKerGP does not support kernel ", kernel_name)
    }
    
    # Ensure the names attribute used by kergp agrees with the names attribute 
    # used by gpWrapper. 
    inputNames(ker) <- X_names
    
    # Define kernel hyperparameters.
    par_names <- attr(ker, "kernParNames")
    
    return(list(name=kernel_name, object=ker, par_names=par_names))
  },
  
  
  define_mean_func = function(mean_func_name, ...) {
    # kergp allows specification of a mean function by an R formula (analogous 
    # to the way it is done in `lm()` to fit a linear model). The `mean_func`
    # attribute for gpWrapperKer is defined to be a list with elements "name",
    # "formula", and "par_names". Some mean function names have been 
    # standardized, so that a string can be passed for "mean_func_name" and the 
    # formula will be automatically created. Alternatively, `mean_func_name` can
    # be a formula, in which case the "name" will be set to 
    # "user-defined formula".
    
    l <- list()
    
    if(mean_func_name == "constant") {
      l$name <- "constant"
      l$formula <- as.formula("y ~ 1")
    } else if(mean_func_name == "linear") {
      l$name <- "linear"
      l$formula <- as.formula(paste0("y ~ ", paste(X_names, collapse=" + ")))
    } else if(mean_func_name == "quadratic") {
      l$name <- "quadratic (no cross terms)"
      l$formula <- as.formula(paste0("y ~ ", 
                                     paste0("poly(", .self$X_names, ", degree=2)", 
                                            collapse=" + ")))
    } else if(inherits(mean_func_name,"formula")) {
      l$name <- "user-specified formula"
      l$formula <- mean_func_name
    } else {
      stop("gpWrapperKerGP does not support mean function: ", mean_func_name)
    }
    
    return(l)
  },
  
  get_default_kernel_bounds = function(bounds_mat, output_idx, ...) {
    # Set default bounds and initial values for kernel hyperparameters that have
    # not been provided these values explicitly. A procedure for defining 
    # default bounds is only defined for certain kernels. For example, 
    # no default method is defined for user-defined kernels, as the 
    # interpretation of the kernel hyperparameters is not known a priori.
    # Note that kergp uses the following parameterization for the Gaussian
    # kernel: v * exp{(x-y)^2 / ell^2}.
    #
    # Logic for defining defaults:
    #
    # For marginal variance parameters:
    #   Fit a linear model using the formula `mean_func$formula`, then use the 
    #   residuals from this model to inform the default bounds. If the 
    #   mean function is "constant" then no linear model is fit and the 
    #   `Y_train` data is used in place of the residuals. The residuals are 
    #   then passed to the helper function `get_marginal_variance_bounds()`
    #   in `gp_helper_functions.r`.
    #
    # For lengthscale parameters:
    #   Utilizes the set of pairwise distances calculated between all design 
    #   inputs ("x" values) to set bounds that avoid lengthscales that are 
    #   well above or below the minimum/maximum observed pairwise input 
    #   distance. Details are described in the helper function 
    #   `get_lengthscale_bounds()` in `gp_helper_functions.r`. Note that the 
    #   current convention is to name lengthscale parameters as
    #   "ell_<x_name>", where `x_name` is the respective name in `X_names`.
    
    kernel_name <- .self$kernel$name
    mean_func_name <- .self$mean_func$name
    col_names <- c("init", "lower", "upper")

    # Lengthscale bounds: currently just applies to "Gaussian" kernel, but this
    # could also apply to Matérn. For now, using default arguments in 
    # `get_lengthscale_bounds()`, and setting init lengthscale to the 15th 
    # percentile of the pairwise distance distribution. Note that doing this 
    # for each independent GP is unneccesary, but leaving it for now.
    if(kernel_name %in% c("Gaussian", "Gaussian_plus_Quadratic")) {
      
      # Lengthscale bounds and defaults.
      ls <- get_lengthscale_bounds(.self$X_train, dim_by_dim=FALSE,
                                   include_one_half=FALSE, 
                                   convert_to_square=FALSE)
      ell_names <- paste0("ell_", .self$X_names)
      ell_bounds <- t(ls$ell_bounds)[.self$X_names,,drop=FALSE]
      ell_init <- 0.2*ell_bounds[,"lower"] + 0.8*ell_bounds[,"upper"]
      ell_bounds <- cbind(init=ell_init, ell_bounds)
  
      bounds_mat[ell_names, col_names] <- ifelse(is.na(bounds_mat[ell_names, col_names]), 
                                                 ell_bounds,
                                                 bounds_mat[ell_names, col_names])
      
      # Marginal variance bounds and defaults.
      if(any(is.na(bounds_mat["marg_var",]))) {
        beta <- .self$fixed_pars$mean_func$beta
        df <- data.frame(y=.self$Y_train[,output_idx], .self$X_train)
        
        if(is.null(beta)) {
          lm_fit <- lm(.self$mean_func$formula, data=df)
          lm_resid <- resid(lm_fit)
          p <- 0.7
        } else {
          mf <- model.frame(formula, data=df)
          basis_mat <- model.matrix(.self$mean_func$formula, data=mf)
          lm_resid <- .self$Y_train[,output_idx] - drop(basis_mat %*% beta)
          p <- 0.99
        }
        
        if(is.na(bounds_mat["marg_var","init"])) {
          bounds_mat["marg_var","init"] <- max(var(lm_resid), 
                                               sqrt(.Machine$double.eps))
        }
        
        if(is.na(bounds_mat["marg_var","lower"])) {
          bounds_mat["marg_var","lower"] <- sqrt(.Machine$double.eps)
        }
        
        if(is.na(bounds_mat["marg_var","upper"])) {
          v <- get_marginal_variance_bounds(max(abs(lm_resid)), p=p, 
                                            return_variance=TRUE)
          bounds_mat["marg_var","upper"] <- v
        }
      }
    }
    
    # Defining bound for the "c" parameter in the quadratic portion of the kernel.
    if(kernel_name=="Gaussian_plus_Quadratic") {
      if(is.na(bounds_mat["quad_offset","lower"])) {
        bounds_mat["quad_offset","lower"] <- sqrt(.Machine$double.eps)
      }
    }
    
    # Print warning if any hyperparameters are missing bounds.
    if(any(is.na(bounds_mat[,c("lower","upper")]))) {
       message("Some kernel hyperparameters do not have bounds set, output = ",
               .self$Y_names[output_idx], 
               " This implies there is no default procedure set up to define",
               " bounds for kernel ", kernel_name, ". kergp does not",
               " provide defaults so this may cause issues.")
    }
    
    # Include bounds in kergp kernel object.
    ker <- .self$kernel$object
    assert_that(all(attr(ker, "kernParNames") == rownames(bounds_mat)))
    
    attr(ker, "par") <- ifelse(is.na(bounds_mat[,"init"]), 
                               attr(ker, "par"), bounds_mat[,"init"])
    attr(ker, "parLower") <- ifelse(is.na(bounds_mat[,"lower"]), 
                                    attr(ker, "parLower"), bounds_mat[,"lower"])
    attr(ker, "parUpper") <- ifelse(is.na(bounds_mat[,"upper"]), 
                                    attr(ker, "parUpper"), bounds_mat[,"upper"])
    kernel$object <<- ker 
    
    
    return(bounds_mat)  
  },
  

  fit_package = function(X_fit, y_fit, output_idx, ...) {
    # `kergp` supports a variety of different optimization algorithms, which can 
    # be specified using the `...` arguments. 
    # For kergp, the nugget can be handled more on the prediction side than on 
    # the model fitting side. To fix a jitter, we can fit the model using 
    # `noise = FALSE`. Then for the returned GP object we can set 
    # `gp_fit$varNoise <- nugget`. Then when predicting want to use 
    # `predict(..., forceInterp=TRUE)` to include the nugget variance in the 
    # kernel matrix. By passing and integer argument `multistart = k` via `...`
    # the optimization will be run from `k` different initializations.
    # Bounds on the kernel hyperparameters are included as part of the 
    # `kernel$object` object, which are set when setting the prior.

    y_name <- .self$Y_names[output_idx]
    noise_var_bounds <- .self$par_bounds$noise_var[y_name,]
    
    # kergp supports fixing known mean function coefficients `beta`. Will be 
    # NULL if no fixed beta was provided.
    beta <- .self$fixed_pars$mean_func[[y_name]]$beta

    # When the noise variance is fixed, we use a hack to fix in during kergp's 
    # optimization: define the upper and lower bounds to be the desired value.
    fixed_noise_var <- .self$fixed_pars$noise_var[y_name]
    if(!.self$noise) fixed_noise_var <- .self$default_jitter
    if(!is.null(fixed_noise_var)) {
      noise_var_bounds["init"] <- fixed_noise_var
      noise_var_bounds["lower"] <- fixed_noise_var
      noise_var_bounds["upper"] <- fixed_noise_var
    }

    # Fit GP. 
    gp_fit <- kergp::gp(.self$mean_func$formula, 
                        data=data.frame(y=drop(y_fit), X_fit), 
                        inputs=.self$X_names, 
                        cov=.self$kernel$object, 
                        estim=TRUE, beta=beta, noise=TRUE, 
                        varNoiseIni=noise_var_bounds["init"],
                        varNoiseLower=noise_var_bounds["lower"],
                        varNoiseUpper=noise_var_bounds["upper"],
                        ...)
    
    # Store noise variance.
    noise_var[output_idx] <<- gp_fit$varNoise
    
    return(gp_fit)
  },
  
  predict_package = function(X_new, output_idx, return_mean=TRUE, 
                             return_var=TRUE, return_cov=FALSE, 
                             return_cross_cov=FALSE, X_cross=NULL, 
                             include_noise=TRUE, return_trend=TRUE, 
                             pred_type="SK", ...) {
    # See class notes above for information on `force_interp`; we set this 
    # FALSE here. To align with gpWrapper conventions, the default behavior 
    # is to call kergp's predict function with `type = "SK"` (simple kriging)
    # and `biasCorrect = FALSE` (the latter is already kergp's default). 
    # The former can be overwritten by passing `pred_type="UK"` (universal 
    # kriging). The latter can be overwritten by including `biasCorrect=TRUE`
    # in `...`. The kergp `predict` method returns a list with elements 
    # "sdSK" and "sd" regardless of these specifications. If using "UK" then 
    # "sd" is the inflated variance. If "SK", then both "sd" and "sdSK" are 
    # the same. Since we are setting `force_interp = FALSE` then kergp will 
    # return predictive variances for f(x), not y(x). Therefore, if 
    # `include_noise=TRUE`, we add the noise variance post-hoc.

    if(return_cross_cov) {
      stop("`return_cross_cov` not yet implemented for gpKerGP class.")
    }
    if(is.null(colnames(X_new))) colnames(X_new) <- .self$X_names
    
    pred <- kergp:::predict.gp(.self$gp_model[[output_idx]], newdata=X_new, 
                               forceInterp=FALSE, seCompute=return_var, 
                               covCompute=return_cov, type=pred_type, ...)
    
    return_list <- list()
    if(return_mean) return_list$mean <- pred$mean
    if(return_var) return_list$var <- (pred$sd)^2
    if(return_cov) return_list$cov <- pred$cov
    if(return_trend) return_list$trend <- pred$trend
    
    # If requested, convert predictive variances of f(x) to predictive variances
    # of y(x).
    if(include_noise) {
      if(return_var) return_list$var <- return_list$var + .self$noise_var[output_idx]
      if(return_cov) diag(return_list$cov) <- diag(return_list$cov) + .self$noise_var[output_idx] 
    }

    return(return_list)
  }, 
  
  update_package = function(X_new, y_new, output_idx, update_hyperpar=FALSE, ...) {
    # NOTE: this update is O(N^3), where N is the total number of design points 
    # (including the old design points). This is due to the fact that kergp 
    # stores the Cholesky factor and some other quantities related to the QR 
    # decomposition to form predictions. There is no obvious way to update the 
    # existing quantities other than just re-computing them at the union of the 
    # old and new design points. 
    if(update_hyperpar) {
      stop("kergp update with hyperparameter update not yet implemented.")
    }
    
    gp_obj <- .self$gp_model[[output_idx]]
    
    # Update design. Note that this assumes that this function is called via the 
    # `update()` method, meaning that `.self$X_train` and `.self$Y_train` have 
    # already been updated. 
    gp_obj$X <- .self$X_train
    gp_obj$y <- .self$Y_train[,output_idx]

    # Update basis functions for the GP mean.
    df <- data.frame(y=gp_obj$y, gp_obj$X)
    mean_func_formula <- .self$mean_func$formula
    mf <- model.frame(mean_func_formula, data=df)
    gp_obj$F <- model.matrix(mean_func_formula, data=mf)

    # Call kergp's generalized least squares method, which computes the betaHat 
    # estimate conditional on fixed kernel. However, here we pass the current 
    # value of betaHat so that it will not be re-estimated. The point of calling 
    # this is thus to simply compute intermediate quantities (e.g., Cholesky 
    # factor of covariance). 
    gls_list <- kergp::gls(object=gp_obj$covariance, y=gp_obj$y, X=gp_obj$X, 
                           F=gp_obj$F, beta=gp_obj$betaHat, 
                           varNoise=.self$noise_var)
    
    # Update quantities computed by the GLS function. 
    gls_quantities <- c("L","eStar","sseStar","FStar","RStar")
    gp_obj[gls_quantities] <- gls_list[gls_quantities]

    return(gp_obj)
  }, 
  
  get_default_hyperpar_bounds = function(mean_func_name, kernel_name, X_fit, y_fit, ...) {
    hyperpar_list <- list()

    # Statistics summarizing observed response - use to empirically determine 
    # reasonable default hyperparamter bounds. 
    y_centered <- abs(y_fit-mean(y_fit))
    max_y <- max(abs(y_centered))

    if((mean_func_name=="constant") && (kernel_name=="Gaussian")) {

      # Gaussian kernel hyperparameters.
      ker_par_names <- c(.self$X_names, "sigma2") 
      ls <- get_lengthscale_bounds(X_fit, p_extra=0.15, dim_by_dim=FALSE,
                                   include_one_half=FALSE, convert_to_square=FALSE)
      marg_var_max <- get_marginal_variance_bounds(max_y, p=0.99, return_variance=TRUE)
      qy <- quantile(y_centered, 0.2)
      marg_var_default <- get_marginal_variance_bounds(qy, p=0.90, return_variance=TRUE)
      hyperpar_list$lower <- setNames(c(ls["lower",], 1e-08), ker_par_names)
      hyperpar_list$upper <- setNames(c(ls["upper",], marg_var_max), ker_par_names)
      hyperpar_list$default <- setNames(c(ls[3,], marg_var_default), ker_par_names)
    } else if((mean_func_name=="quadratic") && (kernel_name=="Gaussian")) {
      ker_par_names <- c(X_names, "sigma2")
      ls <- get_lengthscale_bounds(X_fit, p_extra=0.15, dim_by_dim=FALSE,
                                   include_one_half=FALSE, convert_to_square=FALSE)
      rng_y <- diff(range(y_fit))
      marg_var_max <- get_marginal_variance_bounds(rng_y/4, p=0.80, return_variance=TRUE)
      marg_var_default <- get_marginal_variance_bounds(rng_y/8, p=0.80, return_variance=TRUE)
      
      hyperpar_list$lower <- setNames(c(ls["lower",], 1e-08), ker_par_names)
      hyperpar_list$upper <- setNames(c(ls["upper",], marg_var_max), ker_par_names)
      hyperpar_list$default <- setNames(c(ls[3,], marg_var_default), ker_par_names)
      
    } else if((mean_func_name=="constant") && (kernel_name=="Gaussian_plus_Quadratic")) {
      # Hyperparameters: lengthscales, marginal variance, additive cst in quadratic kernel.
      ker_par_names <- c(X_names, "sigma2", "quad_offset")
      ls <- get_lengthscale_bounds(X_fit, p_extra=0.15, dim_by_dim=FALSE,
                                   include_one_half=FALSE, convert_to_square=FALSE)
      rng_y <- diff(range(y_fit))
      marg_var_max <- get_marginal_variance_bounds(rng_y/4, p=0.80, return_variance=TRUE)
      marg_var_default <- get_marginal_variance_bounds(rng_y/8, p=0.80, return_variance=TRUE)
      
      hyperpar_list$lower <- setNames(c(ls["lower",], 1e-08, 0), ker_par_names)
      hyperpar_list$upper <- setNames(c(ls["upper",], marg_var_max, sqrt(marg_var_max)), ker_par_names)
      hyperpar_list$default <- setNames(c(ls[3,], marg_var_default, sqrt(marg_var_default)), ker_par_names)
    } else {
      message("No default kergp hyperparameters defined for mean function <", 
              mean_func_name, "> and kernel <", kernel_name, ">")
      hyperpar_list <- NULL
    }
    
    return(hyperpar_list)
  }, 
  
  get_summary_str_package = function(idx, ...) {
    # This method returns a string summarizing the GP stored in
    # `.self$gp_model[[idx]]`. The `gpWrapper` method `get_summary_str()`
    # loops over the GP list to assemble these strings for all GPs.
    if(!.self$gp_is_fit()) return(NULL)
    gp_kergp <- .self$gp_model[[idx]]
    
    summary_str <- paste0("\n-----> gpWrapperKerGP GP ", idx, "; output = ", 
                          .self$Y_names[idx], ":\n")
    
    # Mean function.
    mean_coefs <- gp_kergp$betaHat
    summary_str <- paste0(summary_str, "\n>>> Mean function:\n")
    summary_str <- paste0(summary_str, "Estimated trend: ", 
                          !gp_kergp$trendKnown, "\n")
    summary_str <- paste0(summary_str, "Mean function coefs:\n")
    for(i in seq_along(mean_coefs)) {
      summary_str <- paste0(summary_str, "\t", names(mean_coefs)[i], ": ", 
                            mean_coefs[i], "\n")
    }
    
    # Kernel.
    # TODO: need to standardize parameter names to be able to extract the 
    # correct parameters from the parameter vector. And need to account for  
    # different kernel types like quadratic. 
    cov_obj <- gp_kergp$covariance
    summary_str <- paste0(summary_str, "\n>>> Kernel:\n")
    summary_str <- paste0(summary_str, "Kernel description: ", attr(cov_obj,"label"), "\n")
    summary_str <- paste0(summary_str, "Lengthscales:\n")
    ls <- attr(cov_obj, "par")[.self$X_names]
    for(i in seq_along(ls)) {
      summary_str <- paste0(summary_str, "\t", names(ls)[i], ": ", ls[i], "\n")
    }
    summary_str <- paste0(summary_str, "Marginal std dev: ", sqrt(attr(cov_obj, "par")["sigma2"]), "\n")
    
    # Nugget/Noise variance.
    summary_str <- paste0(summary_str, "\n>>> Nugget/Noise:\n")
    summary_str <- paste0(summary_str, "Std dev: ", sqrt(.self$default_jitter), "\n")
    
    return(summary_str)
  }
  
)


