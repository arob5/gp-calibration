#
# gp.r
# Class definitions for the gp, gpInd, and gpLik classes. All classes defined 
# using reference classes. 
#
# Andrew Roberts
# 

# -----------------------------------------------------------------------------
# gp: the base Gaussian Process (GP) class. 
# -----------------------------------------------------------------------------

gpWrapper <- setRefClass(
   Class = "gpWrapper", 
   fields = list(gp_model="ANY", lib="character", X="matrix", y="numeric")
)

# -----------------------------------------------------------------------------
# gpHet: Class providing interface between calibration code and hetGP package.  
# -----------------------------------------------------------------------------

gpWrapperHet <- setRefClass(
  Class = "gpWrapperHet",
  contains = "gpWrapper",
  fields = list(lib="character")
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

gp <- gpWrapperHet()

gp$field("lib")
gp$field("model")




