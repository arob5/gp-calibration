#
# EKI_functions.r
# Functions related to Ensemble Kalman Inversion (EKI) and Ensemble Kalman 
# Sampling (EKS).
#
# Andrew Roberts
#

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


run_multiscale_inversion <- function(computer_model_data, m0, Sig0, N_itr, time_step,
                                     sigma, delta, N_ensemble) {
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
  
  if(length(computer_model_data$output_vars) > 1) {
    stop("Multiscale inversion function currently only works with single output variable.")
  }
  
  L0 <- t(chol(Sig0))
  d <- length(computer_model_data$pars_cal_names)
  
  fast_ensemble <- matrix(rnorm(d*N_ensemble), nrow=d, ncol=N_ensemble)

  y <- computer_model_data$data_obs
  u_samp <- matrix(nrow=N_itr, ncol=d)
  u_samp[1,] <- drop(m0 + L0 %*% matrix(rnorm(d), ncol=1))
  
  for(itr in seq(2, N_itr)) {
    # Forward model evaluations. 
    g <- run_computer_model(theta_vals=u_samp[itr-1,], computer_model_data)
    G <- do.call(cbind, run_computer_model(theta_vals=t(u_samp[itr-1,] + sigma*fast_ensemble), computer_model_data))
    
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













