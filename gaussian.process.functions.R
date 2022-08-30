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


create.gp.obj <- function(gp, library, X, y, nugget.override = NULL) {
  
  # Map GP parameters to common parameterization
  gp.obj <- create.gp.params.list(gp, library)
  
  # Define nugget
  if(is.null(nugget.override)) {
    if(gp.obj$gp_sigma == 0) {
      gp.obj$nugget <- sqrt(.Machine$double.eps)
    } else {
      gp.obj$nugget <- gp.obj$gp_sigma^2
    }
  } else {
    gp.obj$nugget <- nugget.override
  }
  
  # Ordinary kriging implies that a plug-in estimate is used for the constant GP mean.
  # In this case the predictive variance equations is modified to account for this 
  # additional uncertainty. 
  if(library == "hetGP") {
    gp.obj$ordinary.kriging <- isTRUE(gp$trendtype == "OK")
  } else {
    gp.obj$ordinary.kriging <- FALSE
  }
  
  # Design matrix/knots
  gp.obj$X <- X
  gp.obj$N <- nrow(X)
  
  # Cholesky Factor
  gp.obj$L <- t(chol(cov_exp_quad(gp.obj$X, NULL, gp.obj$gp_rho, gp.obj$gp_alpha) + diag(rep(gp.obj$nugget, gp.obj$N))))
  
  # Pre-computed term used in predictive mean calculation
  y <- as.matrix(y, ncol = 1)
  gp.obj$K.inv.y <- solve(t(gp.obj$L), solve(gp.obj$L, y - gp.obj$gp_mean))
  
  return(gp.obj)
  
}


predict_gp <- function(X.pred, gp.obj, pred.cov = FALSE) {
  k.Xx <- cov_exp_quad(gp.obj$X, X.pred, gp.obj$gp_rho, gp.obj$gp_alpha)
  pred.list <- list()
  pred.list[["mean"]] <- predict_mean(k.Xx, gp.obj$K.inv.y, gp.obj$gp_mean)
  
  if(pred.cov) {
    pred.list[["cov"]] <- predict_cov_mat(X.pred, k.Xx, gp.obj$L, gp.obj$gp_rho, gp.obj$gp_alpha, gp.obj$nugget, gp.obj$ordinary.kriging)
  } else {
    pred.list[["var"]] <- predict_var(X.pred, k.Xx, gp.obj$L, gp.obj$gp_alpha, gp.obj$nugget, gp.obj$ordinary.kriging)
  }
  
  return(pred.list)
}


predict_mean <- function(k.Xx, K.inv.y, mu) {

  return(mu + t(k.Xx) %*% K.inv.y)
 
}


predict_var<- function(X.pred, k.Xx, L, alpha, nugget, ordinary.kriging) {
  
  if(ordinary.kriging) {
    L.inv.1 <- solve(L, rep(1, nrow(L)))
    K.inv.sum <- sum(L.inv.1^2)
  }
  
  pred.vars <- vector("numeric", nrow(X.pred))
  for(i in seq_along(pred.vars)) {
    v <- solve(L, k.Xx[,i,drop=FALSE])
    pred.vars[i] <- alpha^2 + nugget - sum(v^2)
    
    if(ordinary.kriging) {
       pred.vars[i] <- pred.vars[i] + (1 - sum(v * L.inv.1))^2 / K.inv.sum
    }
    
  }
  
  return(pred.vars)
}

predict_cov_mat <- function(X.pred, K.Xx, L, rho, alpha, nugget, ordinary.kriging) {
  if(ordinary.kriging) {
    L.inv.1 <- solve(L, rep(1, nrow(L)))
    K.inv.sum <- sum(L.inv.1^2)
  }
  
  N.pred <- nrow(X.pred)
  V <- solve(L, K.Xx)
  K.x <- cov_exp_quad(X.pred, NULL, rho, alpha) + diag(rep(nugget, N.pred), nrow = N.pred)
  
  if(ordinary.kriging) {
    return(K.x - crossprod(V) + crossprod(1 - crossprod(V, h)) / K.inv.sum)
  }
  
  return(K.x - crossprod(V))
}









