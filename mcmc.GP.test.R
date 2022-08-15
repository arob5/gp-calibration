# mcmc.GP.test.R
# Contains function mcmc.GP.test(), which is a simplified version of the PEcAn
# function mcmc.GP() to be used for testing purposes and comparison to 
# alternative algorithms. 
#
# Andrew Roberts

library(mvtnorm)

mcmc.GP.test <- function(gp.obj) {
  
  
}


sample.SS <- function(gp.obj, x.pred, x.curr = NULL) {
  gp.pred <- predict_gp(as.matrix(x.pred), gp.obj)
  
  repeat {
    SS <- rnorm(1, gp.pred$mean, sqrt(gp.pred$var))
    if(SS >= 0) break
  }
  
  return(SS)
  
}

