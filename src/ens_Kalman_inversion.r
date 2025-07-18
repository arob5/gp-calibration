# ens_Kalman_inversion.r
# Implements methods for the approximate solution of Bayesian inverse problems
# using methodology based on the ensemble Kalman filter (EnKF). 
#
# Authors: Meng Lai and Andrew Roberts
#
# Depends: 
#    general_helper_functions.r, seq_design.r

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


run_eki <- function(y, fwd, Sig, n_itr=1L, par_prior=NULL, U0=NULL, G0=NULL,
                    transform_pars=TRUE, par_map=NULL, design_method="LHS", 
                    n_ens=NULL) {
  # Runs Ensemble Kalman inversion (EKI) for `n_itr` iterations. This extends 
  # the one-step update implemented by `run_eki_step()` via a likelihood 
  # tempering approach. At present, this function assumes that the inverse 
  # problem has a Gaussian likelihood with covariance `Sig`, but this can 
  # potentially be generalized in the future. To be explicit, the current 
  # assumption is an additive Gaussian noise model, with independent noise 
  # eps ~ N(0,Sig). Each iteration of EKI requires evaluation of the forward 
  # model `fwd`. The final ensemble returned by this function provides a 
  # Monte Carlo approximation of the posterior p(u|y). For linear Gaussian 
  # inverse problems, these samples will be exactly distributed according to 
  # p(u|y), otherwise this represents a pure approximation. Since this EnKF
  # based method relies on Gaussian approximations, it is typically advisable 
  # to transform the parameter ensemble members to an unconstrained space 
  # prior to performing the EnKF updates. The arguments `transform_pars` and 
  # `par_map` allow for such transformations. See details below.
  # 
  # Args:
  #    y: numeric vector of length `P`, where `P` is the observation dimension. 
  #       This vector is the observation being conditioned on.
  #    fwd: function, representing the forward model. Must be vectorized so 
  #         that it can accept a matrix with each row representing a different 
  #         parameter vector inputs. It should return a matrix with rows 
  #         corresponding to the different inputs, and number of columns equal 
  #         to `P`, the dimension of the observation space.
  #    Sig: matrix of shape (P,P), the covariance matrix for the Gaussian 
  #         likelihood.
  #    n_itr: integer, the number of iterations the algorithm will be run.
  #    par_prior: data.frame storing the prior distributions for the parameters.
  #               This is required if `U0` is not provided, or if 
  #               `transform_pars` is TRUE but `par_map` is FALSE.
  #    U0: matrix, initial parameter ensemble of dimension (J,D), where J is 
  #        the number  of ensemble members and D the parameter space dimension. 
  #        If not  provided will be sampled from prior.
  #    G0: matrix, representing the output of `fwd(U0)`, the model outputs based
  #        on the initial ensemble. Optional.
  #    transform_pars: if TRUE, will use `par_map()` (or the default transport
  #                    map) to transform the ensemble members to an unconstrained 
  #                    space prior to the application of the EnKF analysis step. 
  #                    The inverse transformation is then applied before executing
  #                    the next round of forward model evaluations.
  #    par_map: function, representing the parameter transformation map (i.e., 
  #             transport map) and its inverse. See comments on the return value
  #             of `get_default_par_map()` for the requirements on this function.
  #    design_method: character, the sampling method used to generate the 
  #                   initial ensemble; only used if `U0` is NULL. See 
  #                   `get_batch_design()` for the different sampling options.
  #    n_ens: integer, the number of ensemble members J. Only required if 
  #           `U0` is NULL. 
  #
  # Returns:
  # list, with elements:
  #    `U`: matrix, the final JxD parameter ensemble in the untransformed 
  #         (original) space.
  #    `par_map`: function, the transport map used in the algorithm. NULL if 
  #               `transform_pars` is FALSE.
  #    `eki_list`: list of length `n_itr`. The kth element is the list returned 
  #                by the call to `run_eki_step()` at the kth iteration of the 
  #                algorithm.

  # If initial ensemble is not provided, sample from prior.
  if(is.null(U0)) {
    assert_that(!is.null(n_ens) && !is.null(par_prior))
    U <- get_batch_design(design_method, N_batch=n_ens, prior_params=par_prior)
  } else {
    assert_that(is.matrix(U0))
    n_ens <- nrow(U0)
    U <- U0
  }
  
  # Construct default transport map, if not explicitly provided.
  if(transform_pars && is.null(par_map)) {
    par_map <- get_default_par_map(par_prior)
  }
  
  # List storing intermediate outputs.
  eki_list <- vector(mode="list", length=n_itr)
  
  for(k in 1:n_itr) {
    
    # Run forward model. On first iteration, don't run if initial forward 
    # model evaluations have been explicitly passed in args.
    if((k==1L) && !is.null(G0)) G <- G0
    else G <- fwd(U)
    
    # Transform ensemble members to unconstrained space.
    if(transform_pars) U <- par_map(U)
    
    # EnKF update with tempered likelihood.
    Sig_scaled <- n_itr * Sig
    eki_step_list <- run_eki_step(U=U, y=y, G=G, Sig=Sig_scaled)
    U <- eki_step_list$U
    eki_list[[k]] <- eki_step_list
    
    # Map parameters back to original space.
    if(transform_pars) U <- par_map(U, inverse=TRUE)
  }
  
  return(list(U=U, par_map=par_map, eki_list=eki_list))
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
  # noise eps ~ N(0,Sig). Note that no parameter transformations are performed
  # here; see `run_eki()` for such transformations.
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
  #     `U`: the updated "u" ensemble, stored in a (J,D) matrix.
  #     `G`: the argument `G`.
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
  list(U=U_updated, G=G, m_u=m_u, m_y=m_y, C_u=C_u, C_y=C_y,
       L_y=L_y, C_uy=C_uy)
}


# ------------------------------------------------------------------------------
# Parameter Transformations
# ------------------------------------------------------------------------------

get_default_par_map <- function(par_prior) {
  # Given a prior distribution object `par_prior`, returns a function 
  # representing a map from the original prior space to a new (typically 
  # unconstrained) space. The function also contains an argument allowing 
  # it to also represent the inverse of this map. This is a convenience 
  # function returning a default map, where each parameter is treated 
  # independently and the default univariate maps for each parameter are 
  # determined by `get_default_dist_map()`.
  #
  # Args:
  #    par_prior: data.frame, the prior distribution object.
  # 
  # Returns:
  # A function with arguments "inputs" and "inverse". The latter is an argument 
  # controlling whether the function evaluates the transport map (from u to phi)
  # or its inverse (from phi to u). The first argument "inputs" is a matrix 
  # with column names set to the parameter names, and each row representing a
  # different parameter value. The function returns a matrix of the same 
  # shape, representing the application of the transport map or its inverse 
  # to each row of "inputs".
  #
  # Note: 
  # Be careful with scoping issues for the closure in this function. Everything
  # seems to be working now, but before had an issue that every element of 
  # `par_map_list` was being assigned the prior parameters for the last 
  # parameter.
  
  # Construct univariate map for each parameter.
  n_par <- nrow(par_prior)
  par_names <- par_prior$par_name
  par_map_list <- lapply(1:n_par, function(j) get_default_dist_map(par_prior$dist[j],
                                                                   par_prior$param1[j],
                                                                   par_prior$param2[j],
                                                                   par_prior$bound_lower[j], 
                                                                   par_prior$bound_upper[j]))
  names(par_map_list) <- par_names

  # Construct multivariate transport map by combining the univariate maps for 
  # each parameter.
  par_map <- function(inputs, inverse=FALSE) {
    outputs <- lapply(colnames(inputs), 
                      function(par_name) par_map_list[[par_name]](inputs[,par_name], inverse))
    outputs <- do.call(cbind, outputs)
    colnames(outputs) <- colnames(inputs)
    return(outputs)
  }
  
  return(par_map)
}


get_default_dist_map <- function(dist_name, param1=NULL, param2=NULL, 
                                 bound_lower=NULL, bound_upper=NULL) {
  # A convenience function that returns a default transport map and its 
  # inverse for a given distribution. The defaults here are chosen to 
  # transform random variables to N(0,1) Gaussians, using the inverse
  # transform method. This is a reasonable default for methods relying 
  # on the EnKF (which relies on Gaussian approximations), but note that 
  # this only guarantees that each transformed parameter will be marginally 
  # Gaussian. Other multivariate transformations may be preferable in 
  # certain applications. The variable names below use "u" to refer to the 
  # original parameter, and "phi" to refer to the transformed/unconstrained
  # parameter.
  #
  # Args:
  #    dist_name: character, character string defining the distribution, 
  #               aligning with the naming conventions used in the typical
  #               `par_prior` object.
  #    Remaining arguments correspond to the other columns in the `par_prior`
  #    object, and specify the parameters defining the specific distribution 
  #    within the parameterized distribution family specified by `dist_name`.
  #
  # Returns:
  # function, with arguments "inputs" and "inverse". The latter is an argument 
  # controlling whether the function evaluates the transport map (from u to phi)
  # or its inverse (from phi to u). The former is a numeric vector of inputs 
  # (either u or phi) that will be fed through the univariate map.

  map_list <- list()
  
  if(dist_name == "Gaussian") {
    # Identity maps.
    map_u_to_phi <- function(u) u
    map_phi_to_u <- function(phi) phi
  } else if(dist_name == "Uniform") {
    map_u_to_phi <- function(u) qnorm(punif(u, param1, param2))
    map_phi_to_u <- function(phi) qunif(pnorm(phi), param1, param2)
  } else if(dist_name == "Truncated_Gaussian") {
    map_u_to_phi <- function(u) qnorm(truncnorm::ptruncnorm(u, a=bound_lower, 
                                                            b=bound_upper, 
                                                            mean=param1, 
                                                            sd=param2))
    map_phi_to_u <- function(phi) truncnorm::qtruncnorm(pnorm(phi), a=bound_lower, 
                                                        b=bound_upper, 
                                                        mean=param1, 
                                                        sd=param2)
  } else {
    stop("Prior distribution ", dist_name, " not supported.")
  }
  
  # Define the forward and inverse map.
  par_map <- function(inputs, inverse=FALSE) {
    if(inverse) return(map_phi_to_u(inputs))
    map_u_to_phi(inputs)
  }
  
  return(par_map)
}


logit_map <- function(theta, a, b) {
  if (length(theta) != length(a) || length(theta) != length(b)) {
    stop("vector length differ")
  }
  if (any(theta < a | theta > b)) {
    stop("param outbound")
  }
  log((theta - a) / (b - theta))
}


inv_logit <- function(phi, a, b) {
  if (length(phi) != length(a) || length(phi) != length(b)) {
    stop("vector length differ")
  }
  a + (b - a) / (1 + exp(-phi))
}

### each row of dataframe

logit_map_mat <- function(df, a, b) {
  mapped <- apply(df, 1, function(row) {logit_map(row, a, b)}) 
  t(mapped)
}

inv_logit_mat <- function(df, a, b) {
  mapped <- apply(df, 1, function(row) {inv_logit(row, a, b)}) 
  t(mapped)
}


# ------------------------------------------------------------------------------
# Linear Gaussian model helper functions 
# ------------------------------------------------------------------------------

calc_lin_Gauss_cond_moments <- function(G_fwd, y, m0, C0=NULL, Sig=NULL, 
                                        L0=NULL, L_Sig=NULL) {
  # Relative size of parameter and observation dimensions determine which 
  # formulae will be utilized to compute the posterior moments.
  p <- ncol(G_fwd)
  d <- nrow(G_fwd)
  big_param_dim <- (d > p)
  
  if(big_param_dim) {
    calc_lin_Gauss_cond_moments_par(G_fwd, y, m0, C0, Sig, L0, L_Sig)
  } else {
    calc_lin_Gauss_cond_moments_obs(G_fwd, y, m0, C0, Sig, L0, L_Sig)
  }
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


################
# Diagnostics

run_eki_diagnostics <- function(eki_output, U, mcmc_samp, num_margin_pairs=2, prec_trans_fun=NULL) {
  # Args:
  #    U: matrix of shape (J,D), where J = number of ensemble members and 
  #       D = dimension of each ensemble member.
  #    eki_output:       list object that is the output of eki update. 
  #    mcmc_samp:        baselins sample exclusing burn-in period.
  #    num_margin_pairs: number of random plots for 2-d marginals.
  #    prec_trans_fun:   function applied to transform precision matrix. 
  #
  # Returns:
  #    a list containing diagnostic plots
  
  U_new <- eki_output$U
  eki_list <- eki_output$eki_list
  # transform U_new to gaussian space
  U_new_gauss <- eki_output$par_map(U_new)
  
  samp_dt01 <- append_samples_mat(mcmc_samp, U_new, param_type="par",  
                                  test_label="eki", chain_idx=1L)
  samp_dt02 <- append_samples_mat(samp_dt01, U, param_type="par", test_label="prior")
  
  #######################
  ##     Posterior     ##
  #######################
  # 1d marginals of the posterior ensemble
  kde_plts <- get_1d_kde_plots(samp_dt02, test_label_baseline="mcmc")
  # 2d marginals of the posterior ensemble
  samp_f2 <- dcast(mcmc_samp, itr ~ param_name, value.var = "sample")
  samp_mcmc2 <- samp_f2[,inv_prob$par_names]
  
  post_marg_2d_unbd <- plot_density_pairs(U_new_gauss, num_samp = num_margin_pairs)
  ### constrained space U, with baseline
  post_marg_2d_bd <- plot_density_pairs(U_new, num_samp = num_margin_pairs, baseline_df = samp_mcmc2)
  
  #######################
  ##       Prior       ##
  #######################
  G <- fwd(U)
  U_gauss <- eki_output$par_map(U)
  samp_dtp1 <- append_samples_mat(samp_filt, U_gauss, param_type="par",  
                                  test_label="prior_gauss", chain_idx=1L)
  samp_dtp2 <- append_samples_mat(samp_dtp1, U, param_type="par", test_label="prior")
  
  prior_kde_plts <- get_1d_kde_plots(samp_dtp2, test_label_baseline="mcmc")
  ##2d
  pri_marg_2d_unbd <- plot_density_pairs(U_gauss, num_samp = num_margin_pairs)
  pri_marg_2d_bd <- plot_density_pairs(U, num_samp = num_margin_pairs)
  
  # Heatmap of covariance matrix cov[u] (in transformed space)
  corr_u <- cor(U_gauss)
  hm_corr_U <- plot_corr_heatmap(corr_u, triangle = "upper", plt_title = "correlation")
  
  # prior predictive
  #1. Heatmap of covariance matrix cov[y]
  cov_y <- cov(G) + Sig
  corr_y <- get_corr(cov_y)
  hm_corr_y <- plot_corr_heatmap(corr_y, triangle = "upper", plt_title = "correlation")
  
  #2. Heatmap of precision matrix cov[y]^{-1}
  prec <- chol2inv(chol(cov_y))
  part_cor_y <- get_part_corr(prec)
  hm_part_cor_y <- plot_corr_heatmap(part_cor_y, triangle = "upper", plt_title = "partial correlation")
  # plot_corr_heatmap(prec, triangle = "upper", plt_title = "Precision Matrix")
  
  #3. If sparse estimate of precision matrix is used (e.g., banded precision), then also plot heatmap of this.
  hm_part_cor_y_spar <- NULL
  if(!is.null(prec_trans_fun)) {
    alt_prec <- prec_trans_fun(prec)
    part_cor_y_spar <- get_part_corr(alt_prec)
    hm_part_cor_y_spar <- plot_corr_heatmap(part_cor_y_spar, triangle = "upper", plt_title = "partial correlation sparse")
  }
  
  # 1. Heatmap of cross-covariance matrix cov[u,y] (where u samples are in transformed space)
  cov_uy <- cov(U,G)
  cor_uy <- cor(U,G)
  C_uy <- cor(U,G)  
  hm_corr_uy <- plot_cross_cov(C_uy, corr = TRUE)
  
  list(post_ens_kde = kde_plts, post_ens_kde2d_unbd = post_marg_2d_unbd, post_ens_kde2d_bd = post_marg_2d_bd,
       pri_ens_kde = prior_kde_plts, pri_ens_kde2d_unbd = pri_marg_2d_unbd, 
       pri_ens_kde2d_bd = pri_marg_2d_bd, pri_corhm_U = hm_corr_U,
       pri_corhm_y = hm_corr_y, pri_partcorhm_y = hm_part_cor_y, 
       pri_partcorhm_y_spar = hm_part_cor_y_spar, pri_corhm_uy = hm_corr_uy)
}

plot_cross_cov <- function(mat, corr=TRUE) {
  # helper funciton to plot cross covariance
  C_melted <- melt(mat)
  pltname <- "Cross-Correlation"
  if (corr==FALSE) pltname <- "Cross-Covariance"
  
  ggplot(C_melted, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, space = "Lab",
                         name = "rho") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8)) +
    labs(x = "", y = "", title = pltname)
}

get_part_corr <- function(prec_mat) {
  # input: precision matrix
  # output: partial correlation matrix
  if (!is.matrix(prec_mat) || nrow(prec_mat) != ncol(prec_mat)) {
    stop("need square precision matrix")
  }
  # norm and negate
  diag_vec <- diag(prec_mat)
  norm_fac <- sqrt(outer(diag_vec, diag_vec, "*"))
  norm_prec <- prec_mat/norm_fac
  corr_mat  <- -norm_prec
  diag(corr_mat) <- 1  
  return(corr_mat)
}

get_corr <- function(cov_mat) {
  # input: covariance matrix
  # output: correlation matrix
  if (!is.matrix(cov_mat) || nrow(cov_mat) != ncol(cov_mat)) {
    stop("need square cov matrix")
  }
  diag_vec <- diag(cov_mat)
  norm_fac <- sqrt(outer(diag_vec, diag_vec, "*"))
  corr <- cov_mat/norm_fac
  return(corr)
}

plot_corr_heatmap <- function(cormat, triangle = "upper", output_file = NULL,
                              plt_title = NULL, mid=0) {
  # helper function to plot covaraince matrix
  if (triangle == "upper") {
    cormat[lower.tri(cormat,diag = TRUE)] <- NA
    matrix_tri <- cormat
  } else if (triangle == "lower") {
    cormat[upper.tri(cormat,diag = TRUE)] <- NA
    matrix_tri <- cormat
  } else {
    stop("Invalid triangle argument. Use 'upper' or 'lower'.")
  }
  
  melted_cormat <- melt(matrix_tri, na.rm = TRUE)
  
  heatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, space = "Lab", 
                         name = "rho") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, 
                                 size = 12, hjust = 1),
      axis.text.y = element_text(size = 12)
    ) +
    coord_fixed() + ggtitle(plt_title) +
    xlab("") +
    ylab("")
  
  if (!is.null(output_file)) {
    ggsave(output_file, plot = heatmap, width = 8, height = 6)
    message("Heatmap saved to: ", output_file)
  } else {
    print(heatmap)
  }
}
