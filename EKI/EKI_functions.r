#
# EKI_functions.r
# Functions related to Ensemble Kalman Inversion (EKI) and Ensemble Kalman 
# Sampling (EKS).
#
# Andrew Roberts
#

# -----------------------------------------------------------------------------
# Data Generation and Helper Functions 
# -----------------------------------------------------------------------------

gen_linear_Gaussian_data <- function(G, mu0, Sig0, Sig_eps, theta_true=NULL) {
  
  # Parameter and data dimension. 
  d <- ncol(G)
  n <- nrow(G)
  
  # Cholesky factors and precision matrices. 
  L0 <- t(chol(Sig0))
  L_eps <- t(chol(Sig_eps))
  Sig0_inv <- chol2inv(t(L0))
  Sig_eps_inv <- chol2inv(t(L_eps))
  
  # Generate data. 
  if(is.null(theta_true)) theta_true <- rep(0, d)
  y_true <- G %*% theta_true
  y <- y_true ++ L_eps %*% matrix(rnorm(n), ncol=1)
  
  # True posterior moments. 
  Sig <- chol2inv(chol(crossprod(G, Sig_eps_inv %*% G) + Sig0_inv))
  mu <- Sig %*% (crossprod(G, Sig_eps_inv %*% y) + Sig0_inv %*% mu0)
  
  return(list(theta_true=theta_true, y_true=y_true, y=y, mu=mu, Sig=Sig))
  
}


calc_KL_div_Gaussian <- function(m1, m2, C1, C2, L1=NULL, L2=NULL) {
  # Computes the KL divergence between two Gaussians. `m1` and `C1` are taken to 
  # be the mean vector and covariance matrix of the first entry in the KL 
  # divergence (the one that the expectation is taken with respect to). 
  
  if(is.null(L1)) L1 <- t(chol(C1))
  if(is.null(L2)) L2 <- t(chol(C2))
  C2_inv <- chol2inv(t(L2))
  d <- nrow(C1)
  
  term1 <- sum(log(diag(L2)) - log(diag(L1)))
  term2 <- 0.5 * sum(forwardsolve(L2, m1-m2)^2)
  term3 <- 0.5 * sum(C2_inv*C1)
  
  term1 + term2 + term3 - 0.5*d
  
}


compute_running_err <- function(samp, mu_true, cov_true) {
  
  mean_err <- vector(mode="numeric", length=nrow(samp))
  cov_err <- vector(mode="numeric", length=nrow(samp))
  KL_div <- vector(mode="numeric", length=nrow(samp))
  
  mu_true <- drop(mu_true)
  L_true <- t(chol(cov_true))
  
  mean_curr <- colMeans(samp[1:2,,drop=FALSE])
  cov_curr <- cov(samp[1:2,,drop=FALSE])
  
  for(i in seq(3,nrow(samp))) {
    
    # Update mean and covariance. 
    cov_curr <- tcrossprod(samp[i,]-mean_curr)/i + (i-2)/(i-1)*cov_curr
    mean_curr <- samp[i,]/i + ((i-1)/i)*mean_curr

    # Compute error measures. 
    mean_err[i] <- sqrt(sum((mean_curr - mu_true)^2))
    cov_err[i] <- sqrt(sum((cov_curr - cov_true)^2))
    KL_div[i] <- calc_KL_div_Gaussian(mu_true, mean_curr, cov_true, cov_curr, L1=L_true)
    
  }
  
  return(list(mean=mean_err, cov=cov_err, KL=KL_div))
  
}


# -----------------------------------------------------------------------------
# EKI and Related Algorithms 
# -----------------------------------------------------------------------------

run_EKI_finite_time <- function(computer_model_data, theta_prior_params, Sig_eps, total_steps, 
                                N_ensemble, L_eps=NULL, perturb_data=TRUE) {
  # Runs a version of Ensemble Kalman Sampling (EKI) that transforms prior 
  # samples to approximate posterior samples in a fixed finite number of time steps.
  # This is accomplished by approximating the likelihood tempering/Sequential Monte 
  # Carlo approach to Bayesian inverse problems with ensemble Kalman updates. The 
  # initial ensemble is drawn from the prior distribution. Note that this procedure 
  # is different from alternative EKI algorithms which seek to optimize or sample 
  # in infinite-time, relying on ergodicity. 
  
  ensemble_list <- vector(mode="list", length=total_steps+1)
  if(is.null(L_eps)) L_eps <- t(chol(Sig_eps))
  
  # Initial ensemble. 
  ensemble_list[[1]] <- do.call(rbind, lapply(seq_len(N_ensemble), sample_prior_theta(theta_prior_params)))

  for(itr in 2:(total_steps+1)) {
    ensemble_list[[itr]] <- run_EKI_one_step(computer_model_data, total_steps, ensemble_list[[itr-1]], 
                                             N_ensemble, Sig_eps, L_eps, perturb_data)
  }
  
  return(ensemble_list)
  
}


run_EKI_one_step <- function(computer_model_data, total_steps, ensemble_mat, 
                             N_ensemble, Sig_eps, L_eps=NULL, perturb_data=TRUE) {
  # For now this only works for a single objective; i.e. run_computer_model() should return 
  # 
  # Args:
  #    ensemble_mat: JxD (each row is an particle/ensemble member)
  
  if(length(computer_model_data$output_vars) > 1) {
    stop("EKI function currently only works with single output variable.")
  }
  
  # Prediction Step.
  G_mat <- do.call(cbind, run_computer_model(theta_vals=ensemble_mat, computer_model_data))
  C_G <- cov(G_mat)
  C_uG <- crossprod(t(ensemble_mat), G)
  
  # Analysis Step. 
  if(is.null(L_eps)) L_eps <- t(chol(Sig_eps))
  N <- nrow(L_eps)
  if(perturb_data) {
    psi <- sqrt(total_steps) * L_eps %*% matrix(rnorm(N*N_ensemble), nrow=N, ncol=N_ensemble)
  } else {
    psi <- 0
  }
  
  Y <- drop(computer_model_data$data_obs) + psi
  L <- t(chol(C_G + total_steps*Sig_eps))
  D_mat <- backsolve(t(L), forwardsolve(L, Y-G_mat))
  ensemble_mat_updated <- ensemble_mat + t(C_uG %*% D_mat)
  
  return(ensemble_mat_updated)
  
}


run_multiscale_inversion <- function(forward_model, y, Sig_eps, m0, Sig0, N_itr, time_step,
                                     sigma, delta, N_ensemble, theta_init=NULL) {
  # Implements the multiscale Bayesian inversion algorithm from the paper "Derivative-Free
  # Bayesian Inversion Using Multiscale Dynamics" (Pavliotis, Stuart, and Vaes). This 
  # algorithm discretizes a system of slow-fast SDEs to approximately sample from a 
  # Bayesian posterior. The system involves an ensemble of "fast explorers", which 
  # guide a single "distinguished" particle. The discretized trajectory of the distinguished
  # particle (called u here) are the approximate samples. This algorithm assumes a 
  # Gaussian prior and likelihood. It is intended to approximately sample non-Gaussian 
  # unimodal posteriors of non-linear inverse problems. In the refinement limit as 
  # the algorithm parameters sigma and delta go to 0, the algorithm should sample 
  # from the true posterior in infinite time. 
  #
  # Args:
  #    computer_model_data: the list containing the forward model and observed data. 
  #    m0: matrix of shape (d,1), where d is the dimension of the parameter in the 
  #        inverse problem. This is the prior mean. 
  #    Sig0: matrix of shape (d,d), the prior covariance.
  #    N_itr: integer, the number of time steps to run the numerical time-stepping scheme. 
  #    time_step: numeric, the constant size of the time steps. 
  #    sigma: the noise parameter in the algorithm. 
  #    delta: numeric, the refinement parameter.
  #    N_ensemble: the number of "fast explorers". 
  #   
  # Returns:
  #    matrix of shape (N_itr, d). Each row is an approximate sample from the posterior. 
  #    Early samples should typically be dropped as burn-in. 

  L0 <- t(chol(Sig0))
  d <- nrow(Sig0)
  
  fast_ensemble <- matrix(rnorm(d*N_ensemble), nrow=d, ncol=N_ensemble)
  u_samp <- matrix(nrow=N_itr, ncol=d)
  
  # Set initial condition. 
  if(is.null(theta_init)) {
    theta_init <- m0 + L0 %*% matrix(rnorm(d), ncol=1)
  }
  u_samp[1,] <- drop(theta_init)
  
  for(itr in seq(2, N_itr)) {
    # Forward model evaluations. 
    g <- forward_model(matrix(u_samp[itr-1,], ncol=1))
    G <- forward_model(u_samp[itr-1,] + sigma*fast_ensemble)
    
    # g <- run_computer_model(theta_vals=u_samp[itr-1,], computer_model_data)
    # G <- do.call(cbind, run_computer_model(theta_vals=t(u_samp[itr-1,] + sigma*fast_ensemble), computer_model_data))
    
    # u (distinguished particle) update.
    C_hat <- cov(t(fast_ensemble))
    u_samp[itr,] <- u_samp[itr-1,] - 
                    time_step*drop(fast_ensemble %*% crossprod(G-drop(g), g-y) / (N_ensemble*sigma)) - 
                    time_step*drop(C_hat %*% backsolve(t(L0), forwardsolve(L0, u_samp[itr-1,]-drop(m0)))) + 
                    sqrt(2*time_step) * drop(crossprod(chol(C_hat), matrix(rnorm(d))))
    
    # Fast particles update. 
    fast_ensemble <- exp(-time_step/delta^2)*fast_ensemble + 
                     sqrt(1-exp(-2*time_step/delta^2))*matrix(rnorm(d*N_ensemble), nrow=d, ncol=N_ensemble)
    
  }
  
  return(u_samp)
  
}


# -----------------------------------------------------------------------------
# Baseline MCMC Algorithm for Comparison
# -----------------------------------------------------------------------------

rwmh_Gaussian <- function(computer_model_data, sig2_eps, m0, Sig0, theta_init=NULL, N_mcmc=50000, 
                          adapt_frequency=1000, accept_rate_target=0.24, proposal_scale_decay=0.7,
                          proposal_scale_multiplier=1, Cov_prop_init_diag=NULL, adapt_init_threshold=3) {
  # MCMC random walk Metropolis-Hatings implementation for Bayesian inversion  
  # of a potentially non-linear inverse with Gaussian likelihood and Gaussian prior. 
  # Currently assumes a multiplicative Gaussian likelihood (single output), but should
  # generalize this to Gaussian likelihood with general covariance structure. Also 
  # currently assumes known noise variance; should also generalize this by adding 
  # a Gibbs step. 
  #
  # Args:
  #    computer_model_data: list, standard computer model data list.
  #    m0: matrix of dimension (d,1) where d is the number of parameters. The prior mean. 
  #    Sig0: matrix of dimension (d,d), the prior covariance. 
  #
  # Returns:
  #    list, with named elements "theta" and "Cov_prop". The former is the matrix of samples;
  #    the latter is the proposal covariance at the final iteration. 
  
  if(length(computer_model_data$output_vars) > 1) {
    stop("rwmh_Gaussian() function currently only works with single output variable.")
  }
  
  # Dimension of parameter space.
  d <- length(computer_model_data$pars_cal_names)
  
  # Objects to store samples.
  theta_samp <- matrix(nrow=N_mcmc, ncol=d)
  colnames(theta_samp) <- computer_model_data$pars_cal_names

  # Set initial condition. 
  L0 <- t(chol(Sig0))
  if(is.null(theta_init)) {
    theta_init <- m0 + L0 %*% matrix(rnorm(d), ncol=1)
  }
  theta_samp[1,] <- theta_init
  
  # Proposal covariance setup. 
  if(is.null(Cov_prop_init_diag)) Cov_prop_init_diag <- (2.4)^2 / d
  Cov_prop <- diag(Cov_prop_init_diag, nrow=d)
  L_prop <- t(chol(Cov_prop))
  log_scale_prop <- 0
  effective_log_scale_prop <- log_scale_prop + 0.5*log(Cov_prop_init_diag)
  accept_count <- 0
  samp_mean <- theta_init
  
  # Variables to store partial computations. 
  lprior_curr <- -0.5 * sum(forwardsolve(L0, theta_init-drop(m0))^2)
  llik_curr <- llik_product_Gaussian(computer_model_data, vars_obs=sig2_eps, theta_vals=theta_init, 
                                     normalize=FALSE, na.rm=TRUE)
  lpost_curr <- lprior_curr + llik_curr
  
  for(itr in seq(2, N_mcmc)) {
    # theta proposal.
    theta_prop <- theta_samp[itr-1,] + (exp(log_scale_prop) * L_prop %*% matrix(rnorm(d), ncol = 1))[,1]
    
    # Calculate log_prior and log-likelihood for proposed theta.
    llik_prop <- llik_product_Gaussian(computer_model_data, vars_obs=sig2_eps, theta_vals=theta_prop, 
                                       normalize=FALSE, na.rm=TRUE)
    lprior_prop <- -0.5 * sum(forwardsolve(L0, theta_prop-drop(m0))^2)
    lpost_prop <- lprior_prop + llik_prop
    
    # Accept-Reject step. 
    alpha <- min(1.0, exp(lpost_prop - lpost_curr))
    if(runif(1) <= alpha) {
      theta_samp[itr,] <- theta_prop
      lpost_curr <- lpost_prop 
      accept_count <- accept_count + 1 
    } else {
      theta_samp[itr,] <- theta_samp[itr-1,]
    }
    
    # Adapt proposal covariance matrix and scaling term.
    adapt_list <- adapt_cov_prop(Cov_prop, L_prop, log_scale_prop, theta_samp, itr, accept_count, alpha, samp_mean, 
                                 effective_log_scale_prop, adapt_frequency, accept_rate_target, proposal_scale_decay, 
                                 proposal_scale_multiplier, adapt_init_threshold)
    Cov_prop <- adapt_list$C
    L_prop <- adapt_list$L
    log_scale_prop <- adapt_list$log_scale
    effective_log_scale_prop <- adapt_list$effective_log_scale
    samp_mean <- adapt_list$samp_mean
    accept_count <- adapt_list$accept_count
  }
  
  return(list(theta=theta_samp, Cov_prop=Cov_prop))
  
}




