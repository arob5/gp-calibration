# mcmc.GP.test.R
# Contains function mcmc.GP.test(), which is a simplified version of the PEcAn
# function mcmc.GP() to be used for testing purposes and comparison to 
# alternative algorithms. 
#
# Andrew Roberts


mcmc.GP.test <- function(gp.obj) {
  
  L <- t(chol(K(X_obs, X_obs, rho, alpha) + diag(rep(eps, nrow(X_obs)))))
  
  # Add nugget
  if(add.nugget) {
    if(sigma == 0) {
      eps <- sqrt(.Machine$double.eps)
    } else {
      eps <- sigma^2
    }
    
    K <- K + diag(rep(eps, nrow(X1)))
  }
  
  
}

