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
                 input_bounds="matrix", Y_mean="numeric", Y_std="numeric",
                 X_names="character", Y_names="character",
                 X_std="matrix", Y_norm="matrix")
)

gpWrapper$methods(
  
  initialize = function(X, Y, scale_input=FALSE, normalize_output=FALSE, ...) {

    initFields(X_dim=ncol(X), Y_dim=ncol(Y), scale_input=scale_input, 
               normalize_output=normalize_output)
    
    print(.self$X_dim)
    
    if(normalize_output) {
      initFields(Y_mean=colMeans(Y), Y_std=apply(Y, 2, sd))
      initFields(Y_norm=.self$normalize(Y))
    }
    
    initFields(X=X, Y=Y)
    initFields(...)
  },
  
  scale = function(Xnew, inverse=FALSE) {
    if(inverse) {
      Xnew <- Xnew %*% diag(input_bounds[2,] - input_bounds[1,], X_dim) + 
              matrix(input_bounds[1,], nrow=nrow(Xnew), ncol=X_dim, byrow=TRUE)
    } else {
      Xnew <- (Xnew - matrix(input_bounds[1,], nrow=nrow(Xnew), ncol=X_dim, byrow=TRUE)) %*% 
              diag(1/(input_bounds[2,] - input_bounds[1,]), X_dim)
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
    callSuper(lib="hetGP", gp_model=gp_model, ...)
  }, 
  
  predict = function(X_new, mean=TRUE, var=TRUE, cov=FALSE, X_cov=X_new, include_nug=TRUE) {
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
      if(include_nug) return_list$var <- return_list$var + pred$nugs
    }
    
    if(cov) return_list$cov <- pred$cov
    
    return(return_list)
  }
  
)




#
# TEST
#

# Testing base class. 
X <- matrix(seq(0,1,length.out=5), ncol=1)
Y <- X^2 + 0.2*matrix(rnorm(nrow(X)), ncol=1)

gp <- gpWrapper(X,Y, normalize_output=TRUE)
gp$field("lib")
gp$field("X")
gp$field("Y")
gp$Y_mean
gp$Y_std
gp$Y_norm


# Testing hetGP wrapper.
gpHet <- gpWrapperHet()

gpHet$field("lib")
gpHet$field("model")




