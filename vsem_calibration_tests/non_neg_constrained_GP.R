# non_neg_contrained_GP.r
#
# An implementation of the constrained MLE optimization in the paper 
# "Nonnegativity-enforced Gaussian Process Regression" (Pensoneault et al, 2020).
# This relies on functions from the package "hetGP" for computing the log marginal 
# likelihood and its gradient. 
#
# Andrew Roberts

library(hetGP)


# ----------------------------------------------------------------------------------
# Unconstrained optimization:
#    This is for testing purposes, and making sure I understand the behavior of the 
#    hetGP functions. I try to reproduce the hetGP unconstrained MLE, but simplify 
#    things by removing a lot of the options/customization. 
# ----------------------------------------------------------------------------------

# TODO: should map params to correct parameterization in this function, so it doesn't have to be done both in the 
# marg llik and gradient functions. 

mleHomGP_test <- function(X, Z, lower = NULL, upper = NULL, known = NULL,
                          noiseControl = list(g_bounds = c(sqrt(.Machine$double.eps), 1e2)),
                          init = NULL,
                          covtype = c("Gaussian", "Matern5_2", "Matern3_2"),
                          maxit = 100, eps = sqrt(.Machine$double.eps), settings = list(return.Ki = TRUE, factr = 1e7)){
  # Here I have taken the source code of `mleHomGP()` and modified it, removing many of the bells and whistles. 
  # `mleHomGP` uses a numerical optimizer to optimize the lengthscales (theta) and nugget (g). The marginal variance 
  # (nu) and (if performing ordinary kriging) the constant mean (beta0) MLE are available in closed-form. 
  # In the constrained case we cannot use these closed-form plug-in estimates since they may violate the constraints. 
  # Therefore, the biggest change I make to this function is to integrate `nu` into the numerical optimization as well; 
  # we can then compare to `mleHomGP` to see how the result compares to the plug-in estimate. For now, I am not altering 
  # how bet0 is treated and only simple kriging should be used with fixed beta0. 

  # Unlike in the original code, I assume there are no replicates in `X`. We still need to set 
  # `mult` in order to make this work in the hetGP functions. 
  if(is.null(dim(X))) X <- matrix(X, ncol = 1)
  if(nrow(X) != length(Z)) stop("Dimension mismatch between Z and X")
  X0 <- X
  Z0 <- Z
  mult <- rep(1, length(Z))

  # Bounds on lengthscale parameters `theta`. 
  covtype <- match.arg(covtype)
  
  if(is.null(lower) || is.null(upper)){
    auto_thetas <- auto_bounds(X = X0, covtype = covtype)
    if(is.null(lower)) lower <- auto_thetas$lower
    if(is.null(upper)) upper <- auto_thetas$upper
    if(is.null(known[["theta"]]) && is.null(init$theta)) init$theta <- sqrt(upper * lower)
  }
  if(length(lower) != length(upper)) stop("upper and lower should have the same size")
  
  # Save time to train model
  tic <- proc.time()[3]
  
  if(is.null(settings$return.Ki)) settings$return.Ki <- TRUE
  if(is.null(noiseControl$g_bounds)) noiseControl$g_bounds <- c(sqrt(.Machine$double.eps), 1e2)
  
  g_min <- noiseControl$g_bounds[1]
  g_max <- noiseControl$g_bounds[2]
  
  beta0 <- known$beta0
  
  N <- length(Z)
  n <- nrow(X0)
  
  if(is.null(n))
    stop("X0 should be a matrix. \n")
  
  if(is.null(known[["theta"]]) && is.null(init$theta)) init$theta <- 0.9 * lower + 0.1 * upper # useful for mleHetGP

  # OK = Ordinary Kriging, estimates constant mean function which impacts the predictive variance calculation. 
  # SK = Simple Kriging, fixed constant mean function. 
  trendtype <- 'OK'
  if(!is.null(beta0)) trendtype <- 'SK'
    
  ## General definition of fn and gr
  fn <- function(par, X0, Z0, Z, mult, beta0, theta, g, env){
    idx <- 1 # to store the first non used element of par
    
    if(is.null(theta)){
      theta <- par[1:length(init$theta)]
      idx <- idx + length(init$theta)
    }
    
    if(is.null(g)){
      g <- par[idx]
    }
    
    loglik <- logLikHom(X0 = X0, Z0 = Z0, Z = Z, mult = mult, theta = theta, g = g, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    
    if(!is.null(env) && !is.na(loglik)){
      if(is.null(env$max_loglik) || loglik > env$max_loglik){
        env$max_loglik <- loglik
        env$arg_max <- par
      }
    } 
    
    return(loglik)
  }
  
  gr <- function(par, X0, Z0, Z, mult, beta0, theta, g, env){
    idx <- 1
    components <- NULL
    
    if(is.null(theta)){
      theta <- par[1:length(init$theta)]
      idx <- idx + length(init$theta)
      components <- "theta"
    }
    
    if(is.null(g)){
      g <- par[idx]
      components <- c(components, "g")
    }
    return(dlogLikHom(X0 = X0, Z0 = Z0, Z = Z, mult = mult, theta = theta, g = g, beta0 = beta0, covtype = covtype, eps = eps,
                      components = components, env = env))
  }
  
  ## Both lengthscales `theta` and nugget `g` known; no optimization needed. 
  envtmp <- environment()
  if(!is.null(known$g) && !is.null(known[["theta"]])){
    theta_out <- known[["theta"]]
    g_out <- known$g
    out <- list(value = logLikHom(X0 = X0, Z0 = Z0, Z = Z, mult = mult, theta = theta_out, g = g_out, beta0 = beta0, covtype = covtype, eps = eps),
                message = "All hyperparameters given", counts = 0, time = proc.time()[3] - tic)
  } else { # Run optimization. 
    parinit <- lowerOpt <- upperOpt <- NULL
    if(is.null(known[["theta"]])){
      parinit <- init$theta
      lowerOpt <- c(lower)
      upperOpt <- c(upper)
    }
    if(is.null(known$g)){
      parinit <- c(parinit, init$g)
      lowerOpt <- c(lowerOpt, g_min)
      upperOpt <- c(upperOpt, g_max)
    }
    
    out <- try(optim(par = parinit, fn = fn, gr = gr, method = "L-BFGS-B", lower = lowerOpt, upper = upperOpt, theta = known[["theta"]], g = known$g,
                     X0 = X0, Z0 = Z0, Z = Z, mult = mult, beta0 = beta0,
                     control = list(fnscale = -1, maxit = maxit, factr = settings$factr, pgtol = settings$pgtol), env = envtmp))
    ## Catch errors when at least one likelihood evaluation worked
    if(is(out, "try-error"))
      out <- list(par = envtmp$arg_max, value = envtmp$max_loglik, counts = NA,
                  message = "Optimization stopped due to NAs, use best value so far")
    
    if(is.null(known$g)) g_out <- out$par[length(out$par)] else g_out <- known$g
    if(is.null(known[["theta"]])) theta_out <- out$par[1:length(init$theta)] else theta_out <- known[["theta"]]
    
  }
  
  Ki <- chol2inv(chol(add_diag(cov_gen(X1 = X0, theta = theta_out, type = covtype), eps + g_out/ mult)))
  
  if(is.null(beta0)) beta0 <- drop(colSums(Ki) %*% Z0 / sum(Ki))
    
  psi_0 <- drop(crossprod(Z0 - beta0, Ki) %*% (Z0 - beta0))
  
  nu <- 1/N * ((crossprod(Z - beta0) - crossprod((Z0 - beta0) * mult, Z0 - beta0))/g_out + psi_0)
  
  # # Get hessian of our cost function/
  # if (is.null(known[["theta"]])) {
  #   # Jacobian is more precise numerically but doesn't seem to work for some reason
  #   #hess <- jacobian(func = gr, x = out$par, theta = known[["theta"]], g = known$g,
  #   #             X0 = X0, Z0 = Z0, Z = Z, mult = mult, beta0 = beta0, env = envtmp)
  #   fwrap <- function(par, ...) fn(c(par, out$par[length(out$par)]), ...)
  #   hess <- hessian(func = fwrap, x = out$par[1:(length(out$par)-1)], theta = known[["theta"]], g = known$g,
  #                   X0 = X0, Z0 = Z0, Z = Z, mult = mult, beta0 = beta0, env = NULL)
  # } else {
  #   hess <- NULL
  # }
  
  res <- list(theta = theta_out, g = g_out, nu_hat = as.numeric(nu), ll = out$value, nit_opt = out$counts,
              beta0 = beta0, trendtype = trendtype, covtype = covtype, msg = out$message, eps = eps,
              X0 = X0, Z0 = Z0, Z = Z, mult = mult, call = match.call(),
              used_args = list(lower = lower, upper = upper, known = known, noiseControl = noiseControl),
              time = proc.time()[3] - tic) # hess = hess)
  
  if(settings$return.Ki) res <- c(res, list(Ki = Ki))
  
  class(res) <- "homGP"
  return(res)
}


map_llik_parameterization <- function(param_val, param = NA_character_) {
  
  if(param_val == "sqrt") {
    param_val <- param_val^2 
  } else if(param_val == "log") {
    param_val <- exp(2*param_val)
  } else if(param_val == "default") {
    return(param_val)
  } else {
    stop("<param> must be either 'default', 'sqrt', or 'log'.")
  }

  return(param_val)
    
}


marg_llik_GP <- function(X, y, l2, eta2, tau2, beta0 = NULL, covtype = "Gaussian", eps = sqrt(.Machine$double.eps), env = NULL, 
                         param = c(l = "default", eta = "default", tau = "default")) {
  # The default parameterization here (as indicated by the argument names) is that 
  # l2 = [l_1^2, ..., l_D^2] are the lengthscales, eta2 = eta^2 is the nugget variance, and tau2 = tau^2 is the 
  # marginal variance. The `param` argument can be set to "sqrt" to instead use the parameterization 
  # `phi = [l_1, ..., l_D, tau, eta]` and to "log" to use `phi = [log(l_1), ..., log(l_D), log(tau), log(eta)]`. In these 
  # cases the arguments `l2`, `eta2`, `tau2` (despite their names) are interpreted as using the specified parameterization. 
  # They are thus adjusted to make sure the calculation is correct. Note that different parameterizations can be used 
  # for different parameters. 
  # `param` options: "sqrt", "identity", "log"
  
  # Specify correct parameterization of likelihood. 
  l2 <- map_llik_parameterization(l2, param = param$l)
  eta2 <- map_llik_parameterization(eta2, param = param$eta) 
  tau2 <- map_llik_parameterization(tau2, param = param$tau) 
  
  N <- length(y)
  
  # Temporarily store Cholesky transform of K in Ki
  C0 <- hetGP::cov_gen(X1, theta = l2, type = covtype)
  if(!is.null(env)) env$C0 <- C0
  
  C_inv <- chol(C0 + diag(eps + eta2, nrow = nrow(C0))) # Not actually the inverse, but no need to allocate more memory for chol(C). 
  ldetC <- -2 * sum(log(diag(C_inv)))                   # log determinant from Cholesky factor.
  C_inv <- chol2inv(C_inv)                              # Now this is actually the inverse of C. 
  
  if(!is.null(env)) env$C_inv <- C_inv
  if(is.null(beta0)) beta0 <- drop(colSums(C_inv) %*% y / sum(C_inv))
  
  # Compute log likelihood. 
  quad_form_term <- drop(crossprod(y - beta0, C_inv) %*% (y - beta0))
  
  loglik <- -0.5*N*log(2*pi) - 0.5*N*log(tau2) - 0.5*ldetC - 0.5 * (1/tau2) * quad_form_term
  
}


grad_marg_llik_GP <- function(X, y, l2, eta2, tau2, beta0 = NULL, covtype = "Gaussian",
                              eps = sqrt(.Machine$double.eps), components = c("l2", "eta2", "tau2"), env = NULL, 
                              param = c(l = "default", eta = "default", tau = "default")) {
  
  # Specify correct parameterization of likelihood. 
  l2 <- map_llik_parameterization(l2, param = param$l)
  eta2 <- map_llik_parameterization(eta2, param = param$eta) 
  tau2 <- map_llik_parameterization(tau2, param = param$tau) 
  
  N <- length(y)
  
  # Get matrices C0 and inverse of C if they have already been computed, else 
  # compute them from scratch. 
  if(!is.null(env)) {
    C0 <- env$C0
    C_inv <- env$C_inv
  } else {
    C0 <- cov_gen(X1 = X, theta = l2, type = covtype)
    C_inv <- chol2inv(chol(C0 + diag(eps + eta2, nrow = nrow(C0))))
  }
  
  if(is.null(beta0)) beta0 <- drop(colSums(C_inv) %*% y / sum(C_inv))
  
  y <- y - beta0
  
  C_inv_y <- C_inv %*% y # to avoid recomputing  
  
  grad_l2 <- grad_eta2 <- grad_tau2 <- NULL
  
  # First component, derivative with respect to lengthscale `l2`.  
  if("l2" %in% components){
    grad_l2 <- rep(NA, length(l2))
    
    if(length(theta)==1) {
      dC_dl2 <- partial_cov_gen(X1 = X, theta = l2, type = covtype, arg = "theta_k") * C0
      grad_l2 <- -0.5 * (1/tau2) * crossprod(C_inv_y, dC_dl2) %*% C_inv_y - 0.5 * trace_sym(C_inv, dC_dl2) # `trace_sym` extracts diagonal of matrix product. 
    } else {
      for(d in 1:length(l2)) { # One dimension at a time. 
        dC_dl2 <- partial_cov_gen(X1 = X[,d, drop = FALSE], theta = l2[d], type = covtype, arg = "theta_k") * C0
        grad_l2[d] <- -0.5 * (1/tau2) * crossprod(C_inv_y, dC_dl2) %*% C_inv_y - 0.5 * trace_sym(C_inv, dC_dl2)
      }
    } 
  }
  
  # Second component derivative with respect to nugget `eta2`. 
  if("eta2" %in% components) {
    grad_eta2 <- 0.5 * (1/tau2) * sum(C_inv_y^2) - 0.5 * sum(diag(C_inv))
  }
  
  # Third component derivative with respect to marginal variance `tau2`. 
  if("tau2" %in% components) {
    grad_tau2 <- 0.5 * (1/tau2^2) * drop(crossprod(y, C_inv_y)) - 0.5 * N / tau2
  }
  
  
  return(c(grad_l2, grad_eta2, grad_tau2))
  
}

















## General definition of fn and gr
marg_llik_GP <- function(par, X0, Z0, Z, mult, beta0, theta, g, env){
  idx <- 1 # to store the first non used element of par
  
  if(is.null(theta)){
    theta <- par[1:length(init$theta)]
    idx <- idx + length(init$theta)
  }
  
  if(is.null(g)){
    g <- par[idx]
  }
  
  loglik <- logLikHom(X0 = X0, Z0 = Z0, Z = Z, mult = mult, theta = theta, g = g, beta0 = beta0, covtype = covtype, eps = eps, env = env)
  
  if(!is.null(env) && !is.na(loglik)){
    if(is.null(env$max_loglik) || loglik > env$max_loglik){
      env$max_loglik <- loglik
      env$arg_max <- par
    }
  } 
  
  return(loglik)
}


marg_llik_GP_gradient <- function(par, X0, Z0, Z, mult, beta0, theta, g, env){
  idx <- 1
  components <- NULL
  
  if(is.null(theta)){
    theta <- par[1:length(init$theta)]
    idx <- idx + length(init$theta)
    components <- "theta"
  }
  
  if(is.null(g)){
    g <- par[idx]
    components <- c(components, "g")
  }
  return(dlogLikHom(X0 = X0, Z0 = Z0, Z = Z, mult = mult, theta = theta, g = g, beta0 = beta0, covtype = covtype, eps = eps,
                    components = components, env = env))
}







