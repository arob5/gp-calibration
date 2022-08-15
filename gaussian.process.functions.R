# gaussian.process.functions.R
# Contains R implementation of the functions in gaussian_process_functions.stan.
# These functions are used for testing the Stan functions as well as for running 
# tests in R without having to compile Stan code.
#
# Andrew Roberts
# Dependencies: stan.helper.functions.R

cov_exp_quad <- function(X1, X2 = NULL, rho, alpha) {
  
  if(is.null(X2)) X2 <- X1

  # Create correlation matrix
  C <- matrix(NA, nrow = nrow(X1), ncol = nrow(X2))
  weights <- 1 / rho^2
  
  for(i in 1:nrow(X1)) {
    for(j in 1:nrow(X2)) {
      C[i, j] <- sum(weights * (X1[i, ] - X2[j, ])^2)
    }
  }
  
  return(alpha^2 * exp(-0.5 * C))
  
}

create.gp.obj <- function(gp, library, X.obs, y.obs) {
  
  # Map GP parameters to common parameterization
  gp.obj <- create.gp.params.list(gp, library)
  
  # Define nugget
  if(gp.obj$gp_sigma == 0) {
    gp.obj$nugget <- sqrt(.Machine$double.eps)
  } else {
    gp.obj$nugget <- gp.obj$gp_sigma^2
  }
  
  # Design matrix/knots
  gp.obj$X.obs <- X.obs
  gp.obj$N <- nrow(X.obs)
  
  # Cholesky Factor
  gp.obj$L <- t(chol(cov_exp_quad(gp.obj$X.obs, NULL, gp.obj$gp_rho, gp.obj$gp_alpha) + diag(rep(gp.obj$nugget, gp.obj$N))))
  
  # Pre-computed term used in predictive mean calculation
  gp.obj$K.inv.y <- solve(t(gp.obj$L), solve(gp.obj$L, y.obs - gp.obj$gp_mean))
  
  return(gp.obj)
  
}


predict_gp <- function(X.pred, gp.obj) {
  k.Xx <- cov_exp_quad(gp.obj$X.obs, X.pred, gp.obj$gp_rho, gp.obj$gp_alpha)
  pred.list <- list()
  pred.list[["mean"]] <- predict_mean(k.Xx, gp.obj$K.inv.y, gp.obj$gp_mean)
  pred.list[["var"]] <- predict_var(X.pred, k.Xx, gp.obj$L, gp.obj$gp_rho, gp.obj$gp_alpha, gp.obj$nugget)
  
  return(pred.list)
}


predict_mean <- function(k.Xx, K.inv.y, mu) {

  return(mu + t(k.Xx) %*% K.inv.y)
 
}


predict_var<- function(X.pred, k.Xx, L, rho, alpha, nugget) {
  
  pred.vars <- vector("numeric", nrow(X.pred))
  for(i in seq_along(pred.vars)) {
    k.x <- cov_exp_quad(X.pred[i,,drop=FALSE], NULL, rho, alpha) + nugget
    v <- solve(L, k.Xx[,i,drop=FALSE])
    pred.vars[i] <- k.x - sum(v^2)
  }
  
  return(pred.vars)
}










