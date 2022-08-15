# mcmc.GP.test.R
# Contains function mcmc.GP.test(), which is a simplified version of the PEcAn
# function mcmc.GP() to be used for testing purposes and comparison to 
# alternative algorithms. 
#
# Andrew Roberts

library(TruncatedNormal)
library(mvtnorm)

# TODO: generalize so that each dimension can be updated independently 
mcmc.GP.test <- function(gp.obj, n.itr) {
  for(itr in seq(1, n.itr)) {
    # TODO: adapt proposal variance
    
    # Propose new calibration parameters
    # TODO: define rng
    u.new <- TruncatedNormal::rtmvnorm(1, mu = c(u.curr), sigma = cov.proposal, lb = rng[,1], ub = rng[,2])
    
    
  }
  
}


sample.SS <- function(gp.obj, X.pred) {

  gp.pred <- predict_gp(X.pred, gp.obj, pred.cov = TRUE)

  repeat {
    SS <- rmvnorm(1, mean = gp.pred$mean, sigma = gp.pred$cov)
    if(all(SS >= 0)) break
  }
  
  return(SS)
  
}

