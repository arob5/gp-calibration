# ens_Kalman_inversion.r
# Implements methods for the approximate solution of Bayesian inverse problems
# using methodology based on the ensemble Kalman filter (EnKF). 
#
# Authors: Meng Lai and Andrew Roberts
#
# Depends: 
#    general_helper_functions.r

compute_enkf_update <- function(U, y, Y, C_uy, C_y=NULL, L_y=NULL) {
  # Applies the standard EnKF update (analysis step) to an ensemble of particles 
  # `U`. No sample means or covariances are estimated in this function; this 
  # function simply evaluates the update equation given the requisite covariances
  # as arguments. 
  #
  # Args:
  #    U: matrix of shape (J,D), where J = number of ensemble members and 
  #       D = dimension of each ensemble member.
  #    y: numeric vector of length `P`, where `P` is the observation dimension. 
  #       This vector is the observation being conditioned on. 
  #    Y: matrix of shape (J,P), the ensemble of "simulated observations" that 
  #       are subtracted from the true observation `y` in the update equation.
  #    C_uy: matrix of shape (D,P), the cross-covariance matrix used in the 
  #          update equation. 
  #    C_y: matrix of shape (P,P), the "y" part of the covariance. Optional if 
  #         its Cholesky factor `L_y` is provided. 
  #    L_y: matrix of shape (P,P), the lower Cholesky factor of `C_y`. If NULL, 
  #         will be computed using `C_y`. 
  #
  # Returns:
  #    matrix of shape (J,D), the updated ensemble. 
  
  assert_that(!(is.null(C_y) && is.null(L_y)))
  
  # Compute lower Cholesky factor.
  if(is.null(L_y)) L_y <- t(chol(C_y))
  
  # Update ensemble.
  Y_diff <- add_vec_to_mat_rows(y,-Y)
  U + t(C_uy %*% backsolve(t(L_y), forwardsolve(L_y, t(Y_diff))))
}


run_eki_step <- function(U, y, G, Sig) {
  # Computes one iteration of Ensemble Kalman inversion, which involves 
  # computing the sample mean and covariance estimates using the current 
  # ensembles, and then calling `compute_enkf_update()`. This function assumes
  # that the forward model has already run at the ensemble members `U`, with 
  # the corresponding outputs provided by `G`. At present, this function 
  # assumes that the inverse problem has a Gaussian likelihood with covariance 
  # `Sig`, but this can potentially be generalized in the future. To be explicit,
  # the current assumption is an additive Gaussian noise model, with independent
  # noise eps ~ N(0,Sig). 
  # 
  # Args:
  #    U: matrix of shape (J,D), where J = number of ensemble members and 
  #       D = dimension of each ensemble member.
  #    y: numeric vector of length `P`, where `P` is the observation dimension. 
  #       This vector is the observation being conditioned on.
  #    G: matrix of shape (J,P), storing the results of the forward model runs 
  #       evaluated at the ensemble members `U`. 
  #    Sig: matrix of shape (P,P), the covariance matrix for the Gaussian 
  #         likelihood. 
  #
  # Returns:
  #  list, with elements:
  #     `ens`: the updated "u" ensemble, stored in a (J,D) matrix.
  #     `m_u`: the estimated mean for the "u" part of the joint Gaussian 
  #            approximation.
  #     `m_y`: the estimated mean for the "y" part. 
  #     `C_u`: the estimated covariance for the "u" part. 
  #     `C_y`: the estimated covariance matrix for the "y" part.
  #     `L_y`: the lower Cholesky factor of `C_y`. 
  #     `C_uy`: the estimated cross-covariance matrix.
  #
  # The means and covariances define the joint Gaussian approximation implicit 
  # in the EnKF update.
  
  # Estimate means.
  m_u <- colMeans(U)
  m_y <- colMeans(G)
  
  # Estimate covariances.
  C_u <- cov(U)
  C_y <- cov(G) + Sig
  C_uy <- cov(U,G)
  L_y <- t(chol(C_y))
  
  # Generate "simulated observations".
  P <- ncol(Sig)
  J <- nrow(U)
  eps <- matrix(rnorm(P*J),nrow=J,ncol=P) %*% chol(Sig)
  Y <- G + eps
  
  # Compute EnKF update.
  U_updated <- compute_enkf_update(U, y, Y, C_uy, L_y=L_y)
  
  # Return list.
  list(ens=U_updated, m_u=m_u, m_y=m_y, C_u=C_u, C_y=C_y,
       L_y=L_y, C_uy=C_uy)
}

calc_lin_Gauss_cond_moments <- function(G_fwd, y, m0, C0=NULL, Sig=NULL, 
                                        L0=NULL, L_Sig=NULL) {

  # Relative size of parameter and observation dimensions determine which 
  # formulae will be utilized to compute the posterior moments.
  p <- ncol(G_fwd)
  d <- nrow(G_fwd)
  big_param_dim <- (d > p)
  
  if(big_param_dim) return(calc_lin_Gauss_cond_moments_par(G_fwd, y, m0, C0=NULL, Sig=NULL, 
                                                           L0=NULL, L_Sig=NUL))
  calc_lin_Gauss_cond_moments_obs(G_fwd, y, m0, C0=NULL, Sig=NULL, L0=NULL, L_Sig=NULL)
}

calc_lin_Gauss_cond_moments_par <- function(G_fwd, y, m0, C0=NULL, Sig=NULL, 
                                            L0=NULL, L_Sig=NULL) {
  .NotYetImplemented()
}

calc_lin_Gauss_cond_moments_obs <- function(G_fwd, y, m0, C0=NULL, Sig=NULL, 
                                            L0=NULL, L_Sig=NULL) {
  
  assert_that(!(is.null(C0) && is.null(L0)))
  assert_that(!(is.null(Sig) && is.null(L_Sig)))
  
  # Calculations require the prior and noise covariances. Could consider 
  # alternative approaches here; e.g., using SVD. 
  if(is.null(Sig)) Sig <- tcrossprod(L_sig)
  if(is.null(C0)) C0 <- tcrossprod(L0)
  
  # Linear solves.
  L <- t(chol(tcrossprod(G%*%C0, G) + Sig))
  L_inv_G <- forwardsolve(L, G)
  L_inv_y <- forwardsolve(L, y)
  
  # Conditional mean and covariance.
  m_cond <- m0 + tcrossprod(C0,G) %*% backsolve(t(L), L_inv_y - L_inv_G%*%m0)
  C_cond <- C0 - crossprod(L_inv_G %*% C0)
  
  return(list(mean=drop(m_cond), cov=C_cond))
}

