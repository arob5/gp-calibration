#
# gp.r
# Class definitions for the gp, gpInd, and gpLik classes. All classes defined 
# using reference classes. 
#
# Andrew Roberts
# 

library(assertthat)

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
                 X_train="matrix", Y_train="matrix", nugget="numeric")
)

gpWrapper$methods(
  
  initialize = function(X, Y, scale_input=FALSE, normalize_output=FALSE, 
                        x_names=NULL, y_names=NULL, 
                        nugget=sqrt(.Machine$double.eps), ...) {
    
    initFields(X=X, Y=Y, X_dim=ncol(X), Y_dim=ncol(Y), scale_input=scale_input, 
               normalize_output=normalize_output, X_bounds=apply(X, 2, range))
    
    if(normalize_output) {
      initFields(Y_mean=colMeans(Y), Y_std=apply(Y, 2, sd))
      initFields(Y_train=.self$normalize(Y))
    } else {
      initFields(Y_train=Y)
    }
    
    if(scale_input) {
      initFields(X_train=.self$scale(X))
    } else {
      initFields(X_train=X)
    }
    
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
  
  fit_package = function(output_idx, kernel="Gaussian", mean_func="constant", fixed_pars=list(), ...) {
    err_msg <- "A fit_package() method must be implemented for each class inheriting from 
                gpWrapper. gpWrapper does not implement its own fit_package() method. The
                fit() method should ultimately call fit_package()."
    stop(err_msg)
  },
  
  fit = function(kernel="Gaussian", mean_func="constant", fixed_pars=list(), ...) {
    fits_list <- vector(mode="list", length=Y_dim)
    for(i in seq_along(fits_list)) fits_list[[j]] <- fit_package(i, kernel, mean_func, fixed_pars, ...)
    gp_model <<- fits_list
  }, 
  
  raw_fit_parallel = function(kernel="Gaussian", mean_func="constant", fixed_pars=list(), ...) {
    .NotYetImplemented()
  },
  
  predict_package = function(X_new, mean=TRUE, var=TRUE, cov=FALSE, X_cov=X_new, include_nugget=TRUE) {
    err_msg <- "A predict_package() method must be implemented for each class inheriting from 
                gpWrapper. gpWrapper does not implement its own predict_package() method. The
                predict() method should ultimately call predict_package()."
    stop(err_msg)
  }, 
  
  predict = function(X_new, mean=TRUE, var=TRUE, cov=FALSE, X_cov=X_new, include_nugget=TRUE) {
    pred_list <- vector(mode="list", length=Y_dim)
    for(i in seq_along(pred_list)) pred_list[[j]] <- predict_package(X_new, mean, var, cov, X_cov, include_nugget)
    names(pred_list) <- Y_names
      
    return(pred_list)
  }, 
  
  predict_parallel = function(X_new, mean=TRUE, var=TRUE, cov=FALSE, X_cov=X_new, include_nugget=TRUE) {
    .NotYetImplemented()
  }

)

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
  contains = "gpWrapper",
  fields = list(lib="character", valid_kernels="character", valid_means="character")
)

gpWrapperHet$methods(
  
  initialize = function(gp_model=NULL, ...) {
    initFields(valid_kernels=c("Gaussian", "Matern5_2", "Matern3_2"), 
               valid_mean_funcs=c("constant"))
    callSuper(lib="hetGP", gp_model=gp_model, ...)
  }, 
  
  map_kernel_name = function(kernel) {
    if(kernel == "Gaussian") return(invisible("Gaussian"))
    else if(kernel == "Matern5_2") return(invisible("Matern5_2"))
    else if(kernel == "Matern3_2") return(invisible("Matern3_2"))
    else stop(paste0("hetGP only supports kernels ",  paste0(valid_kernels, collapse=", ")))
  },
  
  map_mean_func_name = function(mean_func) {
    if(mean_func == "constant") return(invisible("beta0"))
    else stop(paste0("hetGP only supports mean functions ",
                     paste0(valid_mean_funcs, collapse=", ")))
  },
  
  map_fixed_pars = function(fixed_pars) {
    # Valid element names for `fixed_pars`: "mean_func", "kernel", "nugget"
  },
  
  fit_package = function(output_idx, kernel="Gaussian", mean="constant", fixed_pars=list(), ...) {
    gp_fits <- mleHomGP(X_train, Y_train[,output_idx], covtype=map_kernel_name(kernel), known=fixed_pars, ...)
  },
  
  predict_package = function(X_new, mean=TRUE, var=TRUE, cov=FALSE, X_cov=X_new, include_nugget=TRUE) {
    if(cov) {
      X_prime <- X_cov
    } else {
      X_prime <- NULL
    }
    
    pred <- predict(gp_model, X_new, xprime=X_prime)
    return_list <- list()
    if(mean) return_list$mean <- pred$mean
    if(var) {
      return_list$var <- pred$sd2
      if(include_nugget) return_list$var <- return_list$var + pred$nugs
    }
    
    if(cov) {
      return_list$cov <- pred$cov
      # TODO: need to check that this is the right. 
      if(include_nugget) diag(return_list$cov) <- diag(return_list$cov) + pred$nugs 
    }
    
    return(return_list)
  }
  
)




#
# TEST
#

# Testing base class. 
X <- matrix(seq(10,20,length.out=5), ncol=1)
Y <- X^2 + 0.2*matrix(rnorm(nrow(X)), ncol=1)

gp <- gpWrapper(X,Y, normalize_output=TRUE, scale_input=TRUE)
gp$field("lib")
gp$field("normalize_output")
gp$field("scale_input")
gp$field("X_names")
gp$field("Y_names")

gp$field("Y")
gp$field("Y_mean")
gp$field("Y_std")
gp$field("Y_train")
gp$field("X")
gp$field("X_bounds")
gp$field("X_train")

gp$fit()


# Testing hetGP wrapper.
gpHet <- gpWrapperHet()

gpHet$field("lib")
gpHet$field("model")




