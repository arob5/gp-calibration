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
# marg llik and gradient functions. I haven't yet dealt with the paramterization in this function. Also need to 
# deal with bounds on tau2. And to change the bounds based on the parameterization. 
# - This function returns a class of type `homGP`. I need to make sure that I map parameters to correct parameterization 
#   before creating this object. 

mleHomGP_test <- function(X, y, lower_l2 = NULL, upper_l2 = NULL, known = NULL,
                          noiseControl = list(eta2_bounds = c(sqrt(.Machine$double.eps), 1e2)),
                          init = NULL,
                          covtype = c("Gaussian", "Matern5_2", "Matern3_2"),
                          maxit = 100, eps = sqrt(.Machine$double.eps), settings = list(return_C_inv = TRUE, factr = 1e7), 
                          param = c(l = "default", eta = "default", tau = "default")){
  # Here I have taken the source code of `mleHomGP()` and modified it, removing many of the bells and whistles. 
  # `mleHomGP` uses a numerical optimizer to optimize the lengthscales (theta) and nugget (g). The marginal variance 
  # (nu) and (if performing ordinary kriging) the constant mean (beta0) MLE are available in closed-form. 
  # In the constrained case we cannot use these closed-form plug-in estimates since they may violate the constraints. 
  # Therefore, the biggest change I make to this function is to integrate `nu` into the numerical optimization as well; 
  # we can then compare to `mleHomGP` to see how the result compares to the plug-in estimate. For now, I am not altering 
  # how bet0 is treated and only simple kriging should be used with fixed beta0. 

  # Unlike in the original code, I assume there are no replicates in `X`. 
  if(is.null(dim(X))) X <- matrix(X, ncol = 1)
  if(nrow(X) != length(y)) stop("Dimension mismatch between y and X")

  # Settings for lengthscale parameters `l2`: bounds and initial values.  
  covtype <- match.arg(covtype)
  
  if(is.null(lower_l2) || is.null(upper_l2)) {
    auto_lengthscales <- hetGP:::auto_bounds(X = X, covtype = covtype)
    if(is.null(lower_l2)) lower_l2 <- auto_lengthscales$lower
    if(is.null(upper_l2)) upper_l2 <- auto_lengthscales$upper
    if(is.null(known[["l2"]]) && is.null(init$l2)) init$l2 <- sqrt(upper_l2 * lower_l2)
  } else {
    if(is.null(known[["l2"]]) && is.null(init$l2)) init$l2 <- 0.9 * lower_l2 + 0.1 * upper_l2 # useful for mleHetGP
  }
  if(length(lower_l2) != length(upper_l2)) stop("`upper_l2` and `lower_l2` should have the same size")
  
  # TODO: TEMP
  
  
  # Save time to train model
  tic <- proc.time()[3]
  
  # Settings for nugget variance `eta2`.
  if(is.null(noiseControl$eta2_bounds)) noiseControl$eta2_bounds <- c(sqrt(.Machine$double.eps), 1e2)
  eta2_min <- noiseControl$eta2_bounds[1]
  eta2_max <- noiseControl$eta2_bounds[2]
  
  # Constant GP mean.  
  beta0 <- known$beta0
  
  if(is.null(settings$return_C_inv)) settings$return_C_inv <- TRUE
  N <- nrow(X)
  if(is.null(N)) stop("X should be a matrix. \n")
  
  # OK = Ordinary Kriging, estimates constant mean function which impacts the predictive variance calculation. 
  # SK = Simple Kriging, fixed constant mean function. 
  trendtype <- "OK"
  if(!is.null(beta0)) trendtype <- "SK"
    
  ## General definition of fn and gr
  fn <- function(par, X, y, beta0, l2, eta2, tau2, env, parameterization) {
    idx <- 1 # to store the first non used element of par
    
    if(is.null(l2)) {
      l2 <- par[1:length(init$l2)]
      idx <- idx + length(init$l2)
    }
    if(is.null(eta2)){
      eta2 <- par[idx]
      idx <- idx + 1
    }
    if(is.null(tau2)) {
      tau2 <- par[idx]
    }
    
    # GP log marginal likelihood evaluation. 
    llik_marg <- marg_llik_GP(X=X, y=y, l2=l2, eta2=eta2, tau2=tau2, beta0=beta0, covtype=covtype, eps=eps, env=env, parameterization=parameterization)
    
    if(!is.null(env) && !is.na(llik_marg)) {
      if(is.null(env$max_llik) || llik_marg > env$max_llik){
        env$max_llik <- llik_marg
        env$arg_max <- par
      }
    } 
    
    return(llik_marg)
  }
  
  gr <- function(par, X, y, beta0, l2, eta2, tau2, env, parameterization) {
    idx <- 1
    components <- NULL
    
    if(is.null(l2)){
      l2 <- par[1:length(init$l2)]
      idx <- idx + length(init$l2)
      components <- "l2"
    }
    
    if(is.null(eta2)){
      eta2 <- par[idx]
      idx <- idx + 1
      components <- c(components, "eta2")
    }
    
    if(is.null(tau2)) {
      tau2 <- par[idx]
      components <- c(components, "tau2")
    }
    
    return(grad_marg_llik_GP(X=X, y=y, l2=l2, eta2=eta2, tau2=tau2, beta0=beta0, covtype=covtype,
                             eps = eps, components=components, env = env, parameterization = parameterization))
  }
  
  envtmp <- environment()
  if(!is.null(known$l2) && !is.null(known$eta2) && !is.null(known$tau2)) { # All parameters `l2`, `tau2`, and `eta2` known; no optimization needed. 
    l2_out <- known$l2
    eta2_out <- known$eta2
    tau2_out <- known$tau2
    
    out <- list(value = logLikHom(X = X, y = y, l2 = l2_out, eta2 = eta2_out, tau2 = tau2_out, beta0 = beta0, covtype = covtype, eps = eps),
                message = "All hyperparameters given", counts = 0, time = proc.time()[3] - tic)
    
  } else { # Run optimization. 
    parinit <- lowerOpt <- upperOpt <- NULL
    if(is.null(known$l2)) {
      parinit <- init$l2
      lowerOpt <- c(lower_l2)
      upperOpt <- c(upper_l2)
    }
    if(is.null(known$eta2)) {
      parinit <- c(parinit, init$eta2)
      lowerOpt <- c(lowerOpt, eta2_min)
      upperOpt <- c(upperOpt, eta2_max)
    }
    if(is.null(known$tau2)) {
      parinit <- c(parinit, init$tau2)
      # TODO: Add param init info for tau2. Look at my final deep GP paper for idea on prior for tau2. 
    }
    
    
    # --------------------- TODO: TEMP
    out_temp <- try(optim(par = as.numeric(init), fn = fn, gr = gr, method = "L-BFGS-B", l2 = known$l2, 
                          eta2 = known$eta2, tau2 = known$tau2, X = X, y = y, beta0 = beta0, parameterization = parameterization,
                          control = list(fnscale = -1, maxit = maxit, factr = settings$factr, pgtol = settings$pgtol), env = envtmp))
    # ---------------------
    
  
    out <- try(optim(par = parinit, fn = fn, gr = gr, method = "L-BFGS-B", lower = lowerOpt, upper = upperOpt, l2 = known$l2, 
                     eta2 = known$eta2, tau2 = known$tau2, X = X, y = y, beta0 = beta0,
                     control = list(fnscale = -1, maxit = maxit, factr = settings$factr, pgtol = settings$pgtol), env = envtmp))
    
    # Catch errors when at least one likelihood evaluation worked
    if(is(out, "try-error")) {
      out <- list(par = envtmp$arg_max, value = envtmp$max_loglik, counts = NA,
                  message = "Optimization stopped due to NAs, use best value so far")
    }
    
    if(is.null(known$l2)) l2out <- out$par[1:length(init$l2)] else l2_out <- known$l2
    if(is.null(known$eta2)) eta2_out <- out$par[length(out$par)-1] else eta2_out <- known$eta2
    if(is.null(known$tau2)) tau2_out <- out$par[length(out$par)] else tau2_out <- known$tau2
    
  }
  
  # Compute estimate of constant mean. 
  if(is.null(beta0)) {
    C_inv <- chol2inv(chol(add_diag(cov_gen(X1 = X, theta = l2_out, type = covtype), eps + eta2_out)))
    beta0 <- drop(colSums(C_inv) %*% y / sum(C_inv))
  } else {
    C_inv <- NULL
  }
  
  out_list <- list(l2 = l2_out, eta2 = eta2_out, tau2 = tau2_out, llik_marg = out$value, nit_opt = out$counts,
                   beta0 = beta0, trendtype = trendtype, covtype = covtype, msg = out$message, eps = eps,
                   X = X, y = y, call = match.call(),
                   used_args = list(lower_l2 = lower_l2, upper_l2 = upper_l2, known = known, noiseControl = noiseControl),
                   time = proc.time()[3] - tic) 
  
  if(settings$return_C_inv) {
    if(is.null(C_inv)) C_inv <- chol2inv(chol(add_diag(cov_gen(X1 = X, theta = l2_out, type = covtype), eps + eta2_out)))
    out_list <- c(out_list, list(C_inv = C_inv))
  }
    
  # Create `homGP` object, ensuring that the parameters are mapped to the parameterization used by `hetGP` package. 
  homGP_obj <- convert_to_homGP_obj(out_list, parameterization)
  
  return(list(out_list = out_list, homGP_obj = homGP_obj))
  
}


convert_to_homGP_obj <- function(mle_out_list, parameterization = c(l2 = "default", eta2 = "default", tau2 = "default")) {
  
  # Convert to paramterization used by `hetGP`. 
  l2 <- map_llik_parameterization(mle_out_list$l2, param = parameterization$l2)
  eta2 <- map_llik_parameterization(mle_out_list$eta2, param = parameterization$eta2) 
  tau2 <- map_llik_parameterization(mle_out_list$tau2, param = parameterization$tau2)
  
  # Assemble list to be converted to `homGP` object. 
  homGP_list <- list(theta = l2, g = eta2, nu_hat = tau2, ll = mle_out_list$llik_marg, nit_opt = mle_out_list$nit_opt,
                     beta0 = mle_out_list$beta0, trendtype = mle_out_list$trendtype, covtype = mle_out_list$covtype, 
                     msg = mle_out_list$message, eps = mle_out_list$eps,
                     X0 = mle_out_list$X, Z0 = mle_out_list$y, Z = mle_out_list$y, mult = rep(1, length(mle_out_list$y)), 
                     call = mle_out_list$call, used_args = mle_out_list$used_args,
                     time = mle_out_list$time, Ki = mle_out_list$C_inv)
  class(homGP_list) <- "homGP"
  
  return(homGP_list)
  
}


map_llik_parameterization <- function(param_val, parameterization = NA_character_) {
  # This function will always return the "squared parameterization" (i.e. l^2, tau^2, eta^2). 
  # It converts the alternative parameterizations [l, tau, eta] and [log(l), log(tau), log(eta)]
  # to this squared parameterization. 
  
  if(parameterization == "sqrt") {
    param_val <- param_val^2 
  } else if(parameterization == "log") {
    param_val <- exp(2*param_val)
  } else if(parameterization == "default") {
    return(param_val)
  } else {
    stop("<parameterization> must be either 'default', 'sqrt', or 'log'.")
  }

  return(param_val)
    
}


marg_llik_GP <- function(X, y, l2, eta2, tau2, beta0 = NULL, covtype = "Gaussian", eps = sqrt(.Machine$double.eps), env = NULL, 
                         parameterization = c(l2 = "default", eta2 = "default", tau2 = "default")) {
  # The default parameterization here (as indicated by the argument names) is that 
  # l2 = [l_1^2, ..., l_D^2] are the lengthscales, eta2 = eta^2 is the nugget variance, and tau2 = tau^2 is the 
  # marginal variance. The `param` argument can be set to "sqrt" to instead use the parameterization 
  # `phi = [l_1, ..., l_D, tau, eta]` and to "log" to use `phi = [log(l_1), ..., log(l_D), log(tau), log(eta)]`. In these 
  # cases the arguments `l2`, `eta2`, `tau2` (despite their names) are interpreted as using the specified parameterization. 
  # They are thus adjusted to make sure the calculation is correct. Note that different parameterizations can be used 
  # for different parameters. 
  # `param` options: "sqrt", "identity", "log"
  
  # Specify correct parameterization of likelihood. 
  l2 <- map_llik_parameterization(l2, parameterization = parameterization$l2)
  eta2 <- map_llik_parameterization(eta2, parameterization = parameterization$eta2) 
  tau2 <- map_llik_parameterization(tau2, parameterization = parameterization$tau2) 
  
  N <- length(y)
  
  # Temporarily store Cholesky transform of K in Ki
  C0 <- hetGP::cov_gen(X, theta = l2, type = covtype)
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
                              parameterization = c(l = "default", eta = "default", tau = "default")) {
  
  # Specify correct parameterization of likelihood. 
  l2 <- map_llik_parameterization(l2, parameterization = parameterization$l)
  eta2 <- map_llik_parameterization(eta2, parameterization = parameterization$eta) 
  tau2 <- map_llik_parameterization(tau2, parameterization = parameterization$tau) 
  
  N <- length(y)
  
  # Get matrices C0 and inverse of C if they have already been computed, else 
  # compute them from scratch. In MLE these are typically stored in `marg_llik_GP()`
  # so they need not be recomputed here. 
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
    
    if(length(l2)==1) {
      dC_dl2 <- hetGP:::partial_cov_gen(X1 = X, theta = l2, type = covtype, arg = "theta_k") * C0
      grad_l2 <- -0.5 * (1/tau2) * crossprod(C_inv_y, dC_dl2) %*% C_inv_y - 0.5 * hetGP:::trace_sym(C_inv, dC_dl2) # `trace_sym` extracts diagonal of matrix product. 
    } else {
      for(d in 1:length(l2)) { # One dimension at a time. 
        dC_dl2 <- hetGP:::partial_cov_gen(X1 = X[,d, drop = FALSE], theta = l2[d], type = covtype, arg = "theta_k") * C0
        grad_l2[d] <- -0.5 * (1/tau2) * crossprod(C_inv_y, dC_dl2) %*% C_inv_y - 0.5 * hetGP:::trace_sym(C_inv, dC_dl2)
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




















