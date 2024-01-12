#
# gp_mcmc_functions.r
# MCMC algorithms utilizing a Gaussian process (GP) emulator to accelerate inference. 
#
# Andrew Roberts
#

# -----------------------------------------------------------------------------
# Code Standards
# -----------------------------------------------------------------------------

# All GP-MCMC functions must satisfy the following requirements: 
#    - Function name has format "mcmc_calibrate_<algorithm-name>".
#    - Required arguments: "computer_model_data", "theta_prior_params", "emulator_info_list", "..."; 
#                          the ellipses is required as it allows for different algorithms to have 
#                          additional unique arguments not required by other algorithms. 
#    - Optional arguments taken by all GP-MCMC functions: "theta_init", "sig2_eps_init", 
#      "learn_sig_eps", "sig_eps_prior_params", "N_itr", "adapt_cov", "adapt_scale", 
#      "adapt_frequency", "accept_rate_target", "proposal_scale_decay", "proposal_scale_multiplier", 
#      "cov_prop_init", "adapt_init_threshold". 
#    - GP-MCMC functions can also have additional arguments beyond the required ones above. 
#    - Return type: list, with potential (not required) named elements "theta", "sig_eps", "Cov_prop", 
#                   "SSR", "cov_prop_scale". 

# Description of required arguments: 
#    - computer_model_data: list, see file "inverse_problem_functions.r" for description. 
#    - theta_prior_params: data.frame, see file "inverse_problem_functions.r" for description. Note 
#                          that for the GP-MCMC algorithms, the prior is typically truncated to prevent
#                          posterior exploration outside of the bounds determined by the design points. 
#                          Thus, "theta_prior_params" will often be a truncated version of the original 
#                          prior.
#    - emulator_info_list: list, see file "gp_emulator_functions.r" for description. 
#    - theta_init: numeric(D), the initial calibration parameter value for MCMC. 
#    - sig2_eps_init: numeric(P), the initial value of the vector of observation variances for MCMC. 
#                     If "learn_sig_eps" is FALSE, this is a required argument and its value is used 
#                     as the fixed observation variance values. 
#    - learn_sig_eps: logical, if TRUE "sig2_eps" is treated as an unknown parameter and MCMC samples
#                     over the joint state space (theta, sig2_eps), with independent inverse Gamma 
#                     priors placed on the variances. 
#   - sig_eps_prior_params: list, see file "inverse_problem_functions.r" for description. 
#   - N_itr: integer, the number of MCMC iterations. 
#   - cov_prop_init: matrix of dimension (D,D), the initial proposal covariance matrix, which will also 
#                    be scaled by a global scaling constant; see "adapt_cov_prop()". If not adaptation 
#                    is performed, then this is the fixed proposal covariance. 
#   - All arguments with prefix "adapt_" concern the adaptation of the covariance proposal. See function 
#     "adapt_cov_prop()". 

# Description of elements in returned list: 
#    - theta: matrix of dimension (N_itr,D) containing the "theta" MCMC samples (one per row). 
#    - sig_eps: matrix of dimension (N_itr,P) containing the "sig2_eps" MCMC samples (one per row). 
#    - SSR: matrix of dimension (N_itr,P) containing the SSR (i.e. model-data misfit) MCMC samples.  
#           Only relevant for some algorithms. 
#    - cov_prop_scale: matrix of dimension (N_itr,D) containing the effective scale (see "adapt_cov_prop()")
#                      of the proposal covariance matrix at each iteration. This is the square root of the 
#                      diagonal of the proposal covariance times the square root of the global scaling factor. 


# -----------------------------------------------------------------------------
# Helper Functions. 
# -----------------------------------------------------------------------------

lpost_gp_integrated <- function(inputs_scaled, emulator_info_list, theta_prior_params, 
                                include_nugget=TRUE, inputs_unscaled=NULL) {
  # This function computes the approximate (unnormalized) posterior density defined by 
  # taking the expectation of the  GP-approximate posterior wrt the GP. This is equivalent to 
  # replacing the likelihood with  the expected likelihood. This is all conditional on fixed , 
  # variance parameters `sig2_eps` so we are really considering the conditional posterior here.  
  # If the log likelihood is represented by a GP(m(u), k(u)) then the unnormalized expected posterior 
  # is pi_0(u) * exp{m(u) + 0.5*k(u)} where pi_0(u) is the prior density. This function returns 
  # the log of this quantity, potentially for multiple values of `u`. Note that this is the 
  # log expected unnormalized posterior, not the expected log unnormalized posterior. 
  #
  # Args:
  #    inputs_scaled: matrix of dimension (M,D), with M the number of inputs at which to sample. Each row is an 
  #                   input at which to sample. .
  #    emulator_info_list: list, the emulator info list.
  #    theta_prior_params: the data.frame containing prior distribution info. 
  #    include_nugget: logical(1), whether or not to include the nugget variance in GP predictions. 
  #    inputs_unscaled: Unscaled version of `inputs_scaled`; if not passed will have to be computed in order 
  #                     to compute the prior. 
  #
  # Returns:
  #    numeric, the log density evaluations, as described above. 
  
  lpost_gp_pred_list <- predict_lpost_GP_approx(theta_vals_scaled=inputs_scaled, theta_vals_unscaled=inputs_unscaled, 
                                                emulator_info_list=emulator_info_list, sig2_eps=sig2_eps, 
                                                theta_prior_params=theta_prior_params, N_obs=N_obs, 
                                                return_vals=c("mean", "var"))
  
  return(lpost_gp_pred_list$mean + 0.5*lpost_gp_pred_list$var)
  
}


sample_emulator_cond <- function(inputs_scaled, emulator_info_list, cond_type, sig2_eps=NULL, gp_pred_list=NULL, 
                                 trunc_method="truncated", include_nugget=TRUE) {
  # Draws a sample either from p(phi|u) or p(phi|u,Sigma,Y). The former is the GP predictive 
  # distribution for the model-data misfit evaluated at input u. The latter is the conditional 
  # posterior of the model-data misfit value, given input u and variance parameters Sigma. 
  # Currently only supports drawing a single sample, and always includes the nugget variances.
  # Also, sets a floor of zero on the sampled values to prevent negative values. 
  #
  # Args:
  #    inputs_scaled: matrix of dimension (M,D), with M the number of inputs at which to sample. Each row is an 
  #                   input at which to sample. .
  #    emulator_info_list: list, the emulator info list. 
  #    cond_type: character(1), specifies which conditional to sample from. If "prior", 
  #               samples from the GP predictive distribution p(phi|u). If "post", samples
  #               from the GP conditional posterior distribution p(phi|u,Sigma,Y).
  #    sig2_eps: numeric(p), the vector of variance parameters. Only required if cond_type is "post". 
  #    gp_pred_list: list, as returned by the function `predict_independent_GPs()`. This can be optionally 
  #                  passed in if this list has already been pre-computed. This list is always assumed to 
  #                  contain information for the predictive distribution of the GP, NOT the conditional 
  #                  predictive distribution p(phi|u,Sigma,Y). 
  #    trunc_method: character, either "truncated" or "rectified". The method used to prevent negative SSR samples. 
  #
  # Returns:
  #    matrix, of dimension (N,P), with P being the number of output variables (i.e. number of GPs). 
  
  if(!(cond_type %in% c("prior", "post"))) stop("Invalid `cond_type`.")
  if(!(trunc_method %in% c("truncated", "rectified"))) stop("Invalid `trunc_method`")
  
  # Compute GP predictions. 
  if(is.null(gp_pred_list)) {
    gp_pred_list <- predict_independent_GPs(X_pred=inputs_scaled, gp_obj_list=emulator_info_list$gp_fits, 
                                            gp_lib=emulator_info_list$settings$gp_lib, include_cov_mat=FALSE, 
                                            denormalize_predictions=TRUE, output_stats=emulator_info_list$output_stats)
  }
  
  # Shift mean if sample from conditional posterior is desired. 
  var_sel <- ifelse(include_nugget, "var_comb", "var")
  if(cond_type=="post") {
    for(p in seq_along(gp_pred_list)) gp_pred_list[[p]]$mean <- gp_pred_list[[p]]$mean - 
                                                                0.5*(1/sig2_eps[p])*gp_pred_list[[p]][[var_sel]]
  }
  
  # Sample from GP predictive distribution. 
  SSR_samples <- sample_independent_GPs_pointwise(gp_pred_list, transformation_methods=trunc_method, 
                                                  include_nugget=include_nugget)
  
  return(SSR_samples)
  
}


# -----------------------------------------------------------------------------
# Proposal Covariance Adaptation 
# -----------------------------------------------------------------------------

adapt_cov_prop <- function(adapt_cov, adapt_scale, C, L, log_scale, sample_history, itr, accept_count,
                           MH_accept_prob, samp_mean, effective_log_scale, adapt_frequency=1000, 
                           accept_rate_target=0.24, scale_update_decay=0.7, scale_update_mult=1.0, init_threshold=3) {
  # Returns an adapted proposal covariance matrix. The covariance matrix is assumed to be of the form `scale * C`, where 
  # `scale` is a scaling factor. This function supports different algorithms to adapt C and to adapt `scale`, and the methods 
  # can be mixed and matched. Only a portion of the function arguments are needed for certain methods, so take care 
  # to ensure the correct arguments are being passed to the function. The function updates C and also computes the lower 
  # triangular Cholesky factor L of C. For certain methods (e.g. the AM algorithm) C will be updated every iteration, but 
  # L may only be updated intermittently depending on `adapt_cov_frequency`. The matrix `C_L` is defined as the covariance 
  # matrix that corresponds to the Cholesky factor `L`; i.e. C_L = LL^T. The covariance `C` will always be the most up-to-date, 
  # but one may instead choose to utilize `L` or `C_L` for proposals if the proposal covariance is not to be updated 
  # every MCMC iteration; this saves on computation (not having to compute a Cholesky decomposition every iteration) and 
  # can help prevent singular proposal covariances. 
  #
  # Args:
  #    adapt_cov: logical(1), if TRUE adapts the covariance matrix `C`. Otherwise, `C` is left unchanged. 
  #    adapt_scale: logical(1), if TRUE adapts the log scale `log_scale`; otherwise, it is left unchanged.
  #    C: matrix, the current d x d positive definite covariance matrix (not multiplied by the scaling factor). 
  #    L: matrix, the lower triangular Cholesky factor of C. This is what is actually used to generate proposals. In general, this 
  #       can be "out of alignment" with C; for example, C can be updated every iteration, but the Cholesky factor can be updated 
  #       less frequently to save computation time. 
  #    log_scale: numeric, the log of the scaling factor. 
  #    sample_history: matrix, itr_curr x p, where itr_curr is the current MCMC iteration. Note that this must be the entire 
  #                    sample history, even for methods like the AP algorithm, which only require a recent subset of the history. 
  #    itr: integer, the current MCMC iteration. 
  #    accept_count: integer(1), the number of accepted MH proposals over the subset of the sample history that 
  #                  will be used in the updates (not the 
  #                  acceptance count over the whole history!). This is only used for the "AP" and "pecan" methods.
  #    MH_accept_prob: numeric(1), the most recent Metropolis-Hastings acceptance probability.
  #    samp_mean: numeric(d), used only by "AM". This is the cumulative sample mean over the MCMC samples. 
  #               It is updated every iteration. 
  #    effective_log_scale: numeric(d), the (log) current standard deviations used in the proposal distributions for 
  #                              each parameter. These are obtained by multiplying the scale parameter (exp(log_scale))
  #                              by the diagonal of the proposal covariance matrix (then taking the log).
  #                              `effective_log_scale` should reflect  
  #                              the current proposal actually being used (i.e. it should align with `L`, and may not 
  #                              necessarily reflect the current value of `C`). 
  #    adapt_frequency: integer, number of iterations that specifies how often C or L will be updated. The AP algorithm 
  #                     updates both C and L together; the AM algorithm updates C every iteration, but only updates L according 
  #                     to `update_cov_frequency`. 
  #    accept_rate_target: numeric(1), the target acceptance rate.
  #    scale_update_decay: the log scale update formula used the term `scale_update_mult / itr^scale_update_decay`. 
  #    scale_update_mult: see `scale_update_decay`. 
  #    itr_threshold: integer, the proposal covariance matrix is not updated before this iteration. 
  #
  # Returns:
  #    list, with elements "C", "L", "log_scale", "samp_mean", "accept_count", and "effective_log_scale". 
  #    The first three are as described above. "samp_mean" is the mean of the MCMC samples up through the 
  #    current iteration. "accept_count" is the number of acceptences, which is reset at zero each time 
  #    `L` is updated. `effective_prop_scale` is described in "Args". 
  
  # Adapt proposal covariance. 
  if(adapt_cov) {
    samp_mean <- samp_mean + (1 / itr) * (sample_history[itr,] - samp_mean)
      
    if(itr==3) { # Sample covariance from first 3 samples. 
      samp_centered <- t(sample_history[1:3,]) - samp_mean
      C <- 0.5 * tcrossprod(samp_centered)
    } else if(itr > 3) { # Begin covariance updates every iteration. 
      samp_centered <- sample_history[itr,] - samp_mean
      w <- 1/(itr-1)
      C <- C + w * (itr*w*outer(samp_centered, samp_centered) - C)
    }
  }
  
  # Update Cholesky factor of covariance matrix (what is actually used for proposals).
  if((itr >= init_threshold) && (itr %% adapt_frequency==0)) {
    if(adapt_cov) L <- t(chol(C))
    # log_scale <- log_scale - log(sum(diag(C)))
    accept_count <- 0
    effective_log_scale <- log_scale + 0.5*log(diag(C))
  } else {
    # Adapt scaling factor for proposal covariance matrix.
    if(adapt_scale) {
      effective_log_scale <- effective_log_scale - log_scale
      log_scale <- log_scale + (scale_update_mult / itr^scale_update_decay) * (MH_accept_prob - accept_rate_target)
      effective_log_scale <- effective_log_scale + log_scale
    }
  }

  return(list(C=C, L=L, log_scale=log_scale, effective_log_scale=effective_log_scale, 
              samp_mean=samp_mean, accept_count=accept_count))
  
}


# -----------------------------------------------------------------------------
# Specific MCMC algorithms. 
# -----------------------------------------------------------------------------

mcmc_calibrate_ind_gp_gibbs <- function(computer_model_data, theta_prior_params, emulator_info_list,
                                        theta_init=NULL, sig2_eps_init=NULL, learn_sig_eps=FALSE, 
                                        sig_eps_prior_params=NULL, N_itr=50000, adapt_cov=TRUE, 
                                        adapt_scale=TRUE, adapt_frequency=1000,
                                        accept_rate_target=0.24, proposal_scale_decay=0.7,  
                                        proposal_scale_multiplier=1, cov_prop_init_=NULL, adapt_init_threshold=3, ...) {
  # Samples from the finite-dimensional GP emulator extended state space (theta, sig2_eps, phi), using 
  # a Metropolis-within-Gibbs approach with separate steps for `theta`, `sig2_eps`, and `phi`. 
  
  if(isTRUE(learn_sig_eps && sig_eps_prior_params$dist != "IG")) {
    stop("`mcmc_calibrate_ind_GP()` requires inverse gamma priors on observation variances.")
  }
  
  # Number observations in time series, number output variables, and dimension of parameter space.
  p <- length(computer_model_data$output_vars)
  d <- length(computer_model_data$pars_cal_names)
  
  # Ensure parameters in `theta_prior_params` are sorted correctly based on ordering in `computer_model_data`. 
  theta_prior_params <- theta_prior_params[computer_model_data$pars_cal_names,]
  
  # Objects to store samples.
  theta_samp <- matrix(nrow=N_itr, ncol=d)
  colnames(theta_samp) <- computer_model_data$pars_cal_names
  sig2_eps_samp <- matrix(nrow=N_itr, ncol=p)
  colnames(sig2_eps_samp) <- computer_model_data$output_vars
  SSR_samp <- matrix(nrow=N_itr, ncol=p)
  colnames(SSR_samp) <- computer_model_data$output_vars
  cov_prop_scales <- matrix(nrow=N_itr, ncol=d)
  colnames(cov_prop_scales) <- computer_model_data$pars_cal_names
  
  # Set initial conditions. 
  if(is.null(theta_init)) {
    theta_init <- sample_prior_theta(theta_prior_params)
  }
  
  if(learn_sig_eps) {
    if(is.null(sig2_eps_init)) {
      sig2_eps_init <- sample_prior_Sig_eps(sig_eps_prior_params)
    }
  } else {
    if(is.null(sig2_eps_init)) stop("Value for `sig2_eps_init` must be provided when `learn_sig_eps` is FALSE.")
  }
  
  theta_samp[1,] <- theta_init
  sig2_eps_samp[1,] <- sig2_eps_init
  
  theta_scaled_curr <- scale_input_data(theta_samp[1,,drop=FALSE], input_bounds=emulator_info_list$input_bounds)
  gp_pred_list_curr <- predict_independent_GPs(X_pred=theta_scaled_curr, gp_obj_list=emulator_info_list$gp_fits, 
                                               gp_lib=emulator_info_list$settings$gp_lib, include_cov_mat=FALSE, 
                                               denormalize_predictions=TRUE, output_stats=emulator_info_list$output_stats)
  gp_pred_means_curr <- sapply(gp_pred_list_curr, function(x) x$mean)
  gp_pred_vars_curr <- sapply(gp_pred_list_curr, function(x) x$var_comb)
  sig2_eps_curr <- sig2_eps_init
  lprior_theta_curr <- calc_lprior_theta(theta_init, theta_prior_params)
  SSR_curr <- sample_emulator_cond(theta_scaled_curr, emulator_info_list, cond_type="post", sig2_eps=sig2_eps_curr)
  SSR_samp[1,] <- SSR_curr
  
  # Proposal covariance.
  if(is.null(cov_prop_init)) Cov_prop <- diag(rep(1,d))
  else Cov_prop <- cov_prop_init
  Cov_prop <- diag(Cov_prop_init_diag, nrow=d)
  L_prop <- t(chol(Cov_prop))
  log_scale_prop <- log(2.38) - 0.5*log(d)
  effective_log_scale_prop <- log_scale_prop + 0.5*log(diag(Cov_prop))
  accept_count <- 0
  samp_mean <- theta_init
  cov_prop_scales[1,] <- exp(effective_log_scale_prop)
  
  for(itr in seq(2, N_itr)) {
    
    #
    # Metropolis step for theta.
    #
    
    # theta proposals.
    theta_prop <- theta_samp[itr-1,] + (exp(log_scale_prop) * L_prop %*% matrix(rnorm(d), ncol = 1))[,1]
    
    # Immediately reject if proposal is outside of prior bounds (i.e. prior density is 0).
    if(any(theta_prop < theta_prior_params[["bound_lower"]], na.rm=TRUE) ||
       any(theta_prop > theta_prior_params[["bound_upper"]], na.rm=TRUE)) {
      
      theta_samp[itr,] <- theta_samp[itr-1,]
      alpha <- 0
    } else {

      # Accept-Reject step. 
      lprior_theta_prop <- calc_lprior_theta(theta_prop, theta_prior_params)
      theta_scaled_prop <- scale_input_data(matrix(theta_prop, nrow=1), emulator_info_list$input_bounds)
      gp_pred_list_prop <- predict_independent_GPs(X_pred=theta_scaled_prop, gp_obj_list=emulator_info_list$gp_fits, 
                                                   gp_lib=emulator_info_list$settings$gp_lib, include_cov_mat=FALSE, 
                                                   denormalize_predictions=TRUE, output_stats=emulator_info_list$output_stats)
      gp_pred_means_prop <- sapply(gp_pred_list_prop, function(x) x$mean)
      gp_pred_vars_prop <- sapply(gp_pred_list_prop, function(x) x$var_comb)
      
      lpost_theta_curr <- lprior_theta_curr + sum(dnorm(drop(SSR_curr), gp_pred_means_curr, sqrt(gp_pred_vars_curr), log=TRUE))
      lpost_theta_prop <- lprior_theta_prop + sum(dnorm(drop(SSR_curr), gp_pred_means_prop, sqrt(gp_pred_vars_prop), log=TRUE))
      
      alpha <- min(1.0, exp(lpost_theta_prop - lpost_theta_curr))
      
      if(runif(1) <= alpha) {
        theta_samp[itr,] <- theta_prop
        gp_pred_list_curr <- gp_pred_list_prop
        theta_scaled_curr <- theta_scaled_prop
        lprior_theta_curr <- lprior_theta_prop
        accept_count <- accept_count + 1 
      } else {
        theta_samp[itr,] <- theta_samp[itr-1,]
      }
      
    }
    
    # Adapt proposal covariance matrix and scaling term.
    if(adapt_cov || adapt_scale) {
      adapt_list <- adapt_cov_prop(adapt_cov, adapt_scale, Cov_prop, L_prop, log_scale_prop, theta_samp, itr,
                                   accept_count, alpha, samp_mean, effective_log_scale_prop, adapt_frequency,
                                   accept_rate_target, proposal_scale_decay, proposal_scale_multiplier, adapt_init_threshold)
      Cov_prop <- adapt_list$C
      L_prop <- adapt_list$L
      log_scale_prop <- adapt_list$log_scale
      effective_log_scale_prop <- adapt_list$effective_log_scale
      samp_mean <- adapt_list$samp_mean
      accept_count <- adapt_list$accept_count
      cov_prop_scales[itr,] <- exp(effective_log_scale_prop)
    }
    
    #
    # Gibbs step for phi. 
    #
    
    SSR_curr <- sample_emulator_cond(theta_scaled_curr, emulator_info_list, cond_type="post", 
                                     sig2_eps=sig2_eps_curr, gp_pred_list=gp_pred_list_curr) 
    SSR_samp[itr,] <- SSR_curr
    
    #
    # Gibbs step for sig2_eps. 
    #
    
    if(learn_sig_eps) {
      sig2_eps_curr <- sample_cond_post_Sig_eps(SSR=SSR_curr, Sig_eps_prior_params=sig_eps_prior_params, n_obs=computer_model_data$n_obs)
      sig2_eps_samp[itr,] <- sig2_eps_curr
    } else {
      sig2_eps_samp[itr,] <- sig2_eps_init
    }
    
    #
    # Second Gibbs step for phi. 
    #
    
    # SSR_curr <- sample_emulator_cond(theta_scaled_curr, emulator_info_list, cond_type="post", 
    #                                  sig2_eps=sig2_eps_curr, gp_pred_list=gp_pred_list_curr) 
    
    
  }
  
  return(list(theta=theta_samp, sig_eps=sig2_eps_samp, Cov_prop=Cov_prop, 
              SSR=SSR_samp, cov_prop_scale=cov_prop_scales))
  
}


mcmc_calibrate_ind_gp_trajectory <- function(computer_model_data, theta_prior_params, emulator_info_list,
                                             theta_init=NULL, sig2_eps_init=NULL, learn_sig_eps=FALSE, 
                                             sig_eps_prior_params=NULL, N_itr=50000, adapt_cov=TRUE, 
                                             adapt_scale=TRUE, adapt_frequency=1000,
                                             accept_rate_target=0.24, proposal_scale_decay=0.7,  
                                             proposal_scale_multiplier=1, cov_prop_init=NULL, adapt_init_threshold=3,
                                             use_gp_cov=TRUE, second_gibbs_step=FALSE, ...) {
  # A GP-MCMC trajectory type algorithm, with the GP being viewed as an infinite-dimensional parameter. 
  # 
  # Function specific arguments: 
  #    - use_gp_cov: logical(1), if TRUE samples from the joint GP predictive distribution between 
  #                  the the current and proposed theta values. If FALSE, samples independently from 
  #                  the univariate predictive distributions (with the predictive covariance ignored). 
  #    - second_gibbs_step: logical(1), if TRUE performs a second Gibbs step for "phi" before the 
  #                         "sig2_eps" step. This second Gibbs step only requires sampling from the 
  #                         univariate predictive distribution at the current theta value. If FALSE, 
  #                         a second Gibbs step is not performed. 
  
  if(learn_sig_eps && sig_eps_prior_params$dist != "IG") {
    stop("`mcmc_calibrate_ind_GP()` requires inverse gamma priors on observation variances.")
  }
  
  # Number observations in time series, number output variables, and dimension of parameter space.
  p <- length(computer_model_data$output_vars)
  d <- length(computer_model_data$pars_cal_names)
  
  # Ensure parameters in `theta_prior_params` are sorted correctly based on ordering in `computer_model_data`. 
  theta_prior_params <- theta_prior_params[computer_model_data$pars_cal_names,]
  
  # Objects to store samples.
  theta_samp <- matrix(nrow=N_itr, ncol=d)
  colnames(theta_samp) <- computer_model_data$pars_cal_names
  sig2_eps_samp <- matrix(nrow=N_itr, ncol=p)
  colnames(sig2_eps_samp) <- computer_model_data$output_vars
  cov_prop_scales <- matrix(nrow=N_itr, ncol=d)
  colnames(cov_prop_scales) <- computer_model_data$pars_cal_names
  
  # Set initial conditions. 
  if(is.null(theta_init)) {
    theta_init <- sample_prior_theta(theta_prior_params)
  }
  
  if(learn_sig_eps) {
    if(is.null(sig2_eps_init)) {
      sig2_eps_init <- sample_prior_Sig_eps(sig_eps_prior_params)
    }
  } else {
    if(is.null(sig2_eps_init)) stop("Value for `sig2_eps_init` must be provided when `learn_sig_eps` is FALSE.")
  }
  
  theta_samp[1,] <- theta_init
  sig2_eps_samp[1,] <- sig2_eps_init
  
  theta_scaled_curr <- scale_input_data(theta_samp[1,,drop=FALSE], input_bounds=emulator_info_list$input_bounds)
  sig2_eps_curr <- sig2_eps_init
  lprior_theta_curr <- calc_lprior_theta(theta_init, theta_prior_params)
  
  # Proposal covariance.
  if(is.null(cov_prop_init)) Cov_prop <- diag(rep(1,d))
  else Cov_prop <- cov_prop_init
  L_prop <- t(chol(Cov_prop))
  log_scale_prop <- log(2.38) - 0.5*log(d)
  effective_log_scale_prop <- log_scale_prop + 0.5*log(diag(Cov_prop))
  accept_count <- 0
  samp_mean <- theta_init
  cov_prop_scales[1,] <- exp(effective_log_scale_prop)
  
  for(itr in seq(2, N_itr)) {
    
    #
    # Metropolis step for theta.
    #
    
    # theta proposals.
    theta_prop <- theta_samp[itr-1,] + (exp(log_scale_prop) * L_prop %*% matrix(rnorm(d), ncol=1))[,1]
    
    # Sample from GP at input locations corresponding to current and proposed parameter values. 
    theta_scaled_prop <- scale_input_data(matrix(theta_prop, nrow=1), emulator_info_list$input_bounds)
    gp_pred_list <- predict_independent_GPs(X_pred=rbind(theta_scaled_curr, theta_scaled_prop), 
                                            gp_obj_list=emulator_info_list$gp_fits, gp_lib=emulator_info_list$settings$gp_lib,
                                            include_cov_mat=use_gp_cov, denormalize_predictions=TRUE,
                                            output_stats=emulator_info_list$output_stats, cov_includes_nug=TRUE)
    
    if(use_gp_cov) {
      SSR_samples <- sample_independent_GPs_corr(gp_pred_list, transformation_methods="rectified", 
                                                 include_nugget=TRUE)
    } else {
      SSR_samples <- sample_independent_GPs_pointwise(gp_pred_list, 
                                                      transformation_methods=emulator_info_list$settings$transformation_method, 
                                                      include_nugget=TRUE)
    }
    
    # Immediately reject if proposal is outside of prior bounds (i.e. prior density is 0).
    if(any(theta_prop < theta_prior_params[["bound_lower"]], na.rm=TRUE) ||
       any(theta_prop > theta_prior_params[["bound_upper"]], na.rm=TRUE)) {
      
      theta_samp[itr,] <- theta_samp[itr-1,]
      SSR_curr <- SSR_samples[1,,drop=FALSE]
      alpha <- 0
      
    } else {

      # Accept-Reject step.
      lprior_theta_prop <- calc_lprior_theta(theta_prop, theta_prior_params)
      lpost_theta_vals <- calc_lpost_theta_product_lik(computer_model_data, lprior_vals=c(lprior_theta_curr,lprior_theta_prop), 
                                                       SSR=SSR_samples, vars_obs=sig2_eps_curr, 
                                                       normalize_lik=FALSE, na.rm=TRUE, return_list=FALSE)
      alpha <- min(1.0, exp(lpost_theta_vals[2] - lpost_theta_vals[1]))
      
      if(runif(1) <= alpha) {
        theta_samp[itr,] <- theta_prop
        theta_scaled_curr <- theta_scaled_prop
        lprior_theta_curr <- lprior_theta_prop
        SSR_curr <- SSR_samples[2,,drop=FALSE]
        pred_idx_curr <- 2
        accept_count <- accept_count + 1 
      } else {
        theta_samp[itr,] <- theta_samp[itr-1,]
        SSR_curr <- SSR_samples[1,,drop=FALSE]
        pred_idx_curr <- 1
      }
      
    }
    
    # Adapt proposal covariance matrix and scaling term.
    if(adapt_cov || adapt_scale) {
      adapt_list <- adapt_cov_prop(adapt_cov, adapt_scale, Cov_prop, L_prop, log_scale_prop, theta_samp, itr,  
                                   accept_count, alpha, samp_mean, effective_log_scale_prop, adapt_frequency,  
                                   accept_rate_target, proposal_scale_decay, proposal_scale_multiplier, adapt_init_threshold)
                                   
      Cov_prop <- adapt_list$C
      L_prop <- adapt_list$L
      log_scale_prop <- adapt_list$log_scale
      effective_log_scale_prop <- adapt_list$effective_log_scale
      samp_mean <- adapt_list$samp_mean
      accept_count <- adapt_list$accept_count
      cov_prop_scales[itr,] <- exp(effective_log_scale_prop)
    }
    
    #
    # Gibbs step for sig2_eps. 
    #
    
    if(learn_sig_eps) {
      if(second_gibbs_step) {
        SSR_curr <- sample_independent_GPs_pointwise(gp_pred_list, transformation_methods=emulator_info$settings$transformation_method,
                                                     idx_selector=pred_idx_curr, include_nugget=TRUE)
      }
      
      sig2_eps_curr <- sample_cond_post_Sig_eps(SSR=SSR_curr, Sig_eps_prior_params=sig_eps_prior_params, 
                                                n_obs=computer_model_data$n_obs)
      sig2_eps_samp[itr,] <- sig2_eps_curr
    } else {
      sig2_eps_samp[itr,] <- sig2_eps_init
    }

  }
  
  return(list(theta=theta_samp, sig_eps=sig2_eps_samp, Cov_prop=Cov_prop, cov_prop_scale=cov_prop_scales))
  
}


mcmc_calibrate_ind_gp_trajectory_trunc_prop <- function(computer_model_data, theta_prior_params, emulator_info_list,
                                                        theta_init=NULL, sig2_eps_init=NULL, learn_sig_eps=FALSE, 
                                                        sig_eps_prior_params=NULL, N_itr=50000, adapt_cov=TRUE, 
                                                        adapt_scale=TRUE, adapt_frequency=1000,
                                                        accept_rate_target=0.24, proposal_scale_decay=0.7,  
                                                        proposal_scale_multiplier=1, cov_prop_init=NULL, adapt_init_threshold=3,
                                                        use_gp_cov=TRUE, second_gibbs_step=FALSE, ...) {
  
  # `theta_prior_params` should already be truncated, if desired. 
  
  if(learn_sig_eps && sig_eps_prior_params$dist != "IG") {
    stop("`mcmc_calibrate_ind_GP()` requires inverse gamma priors on observation variances.")
  }
  
  # Number observations in time series, number output variables, and dimension of parameter space.
  p <- length(computer_model_data$output_vars)
  d <- length(computer_model_data$pars_cal_names)
  
  # Ensure parameters in `theta_prior_params` are sorted correctly based on ordering in `computer_model_data`. 
  theta_prior_params <- theta_prior_params[computer_model_data$pars_cal_names,]
  
  # Objects to store samples.
  theta_samp <- matrix(nrow=N_itr, ncol=d)
  colnames(theta_samp) <- computer_model_data$pars_cal_names
  sig2_eps_samp <- matrix(nrow=N_itr, ncol=p)
  colnames(sig2_eps_samp) <- computer_model_data$output_vars
  cov_prop_scales <- matrix(nrow=N_itr, ncol=d)
  colnames(cov_prop_scales) <- computer_model_data$pars_cal_names
  
  # Set initial conditions. 
  if(is.null(theta_init)) {
    theta_init <- sample_prior_theta(theta_prior_params)
  }
  
  if(learn_sig_eps) {
    if(is.null(sig2_eps_init)) {
      sig2_eps_init <- sample_prior_Sig_eps(sig_eps_prior_params)
    }
  } else {
    if(is.null(sig2_eps_init)) stop("Value for `sig2_eps_init` must be provided when `learn_sig_eps` is FALSE.")
  }
  
  theta_samp[1,] <- theta_init
  sig2_eps_samp[1,] <- sig2_eps_init
  
  theta_scaled_curr <- scale_input_data(theta_samp[1,,drop=FALSE], input_bounds=emulator_info_list$input_bounds)
  sig2_eps_curr <- sig2_eps_init
  lprior_theta_curr <- calc_lprior_theta(theta_init, theta_prior_params)
  
  # Proposal covariance.
  if(is.null(cov_prop_init)) Cov_prop <- diag(rep(1,d))
  else Cov_prop <- cov_prop_init
  L_prop <- t(chol(Cov_prop))
  log_scale_prop <- log(2.38) - 0.5*log(d)
  effective_log_scale_prop <- log_scale_prop + 0.5*log(diag(Cov_prop))
  accept_count <- 0
  samp_mean <- theta_init
  cov_prop_scales[1,] <- exp(effective_log_scale_prop)
  
  for(itr in seq(2, N_itr)) {
    
    #
    # Metropolis step for theta.
    #
    
    # theta proposals.
    theta_prop <- tmvtnorm::rtmvnorm(1, mean=theta_samp[itr-1,], sigma=exp(2*log_scale_prop)*Cov_prop, 
                                     lower=emulator_info_list$input_bounds["min",], 
                                     upper=emulator_info_list$input_bounds["max",])[1,]
    
    # Sample from GP at input locations corresponding to current and proposed parameter values. 
    theta_scaled_prop <- scale_input_data(matrix(theta_prop, nrow=1), emulator_info_list$input_bounds)
    gp_pred_list <- predict_independent_GPs(X_pred=rbind(theta_scaled_curr, theta_scaled_prop), 
                                            gp_obj_list=emulator_info_list$gp_fits, gp_lib=emulator_info_list$settings$gp_lib,
                                            include_cov_mat=use_gp_cov, denormalize_predictions=TRUE,
                                            output_stats=emulator_info_list$output_stats, cov_includes_nug=TRUE)
    
    if(use_gp_cov) {
      SSR_samples <- sample_independent_GPs_corr(gp_pred_list, transformation_methods="rectified", 
                                                 include_nugget=TRUE)
    } else {
      SSR_samples <- sample_independent_GPs_pointwise(gp_pred_list, 
                                                      transformation_methods=emulator_info_list$settings$transformation_method, 
                                                      include_nugget=TRUE)
    }
    
    # Accept-Reject step.
    lprior_theta_prop <- calc_lprior_theta(theta_prop, theta_prior_params)
    lpost_theta_vals <- calc_lpost_theta_product_lik(computer_model_data, lprior_vals=c(lprior_theta_curr,lprior_theta_prop), 
                                                     SSR=SSR_samples, vars_obs=sig2_eps_curr, 
                                                     normalize_lik=FALSE, na.rm=TRUE, return_list=FALSE)
    q_curr_prop <- tmvtnorm::dtmvnorm(theta_prop, mean=theta_samp[itr-1,], sigma=exp(2*log_scale_prop)*Cov_prop,
                                      lower=emulator_info_list$input_bounds[1,], upper=emulator_info_list$input_bounds[2,], log=TRUE)
    q_prop_curr <- tmvtnorm::dtmvnorm(theta_samp[itr-1,], mean=theta_prop, sigma=exp(2*log_scale_prop)*Cov_prop,
                                      lower=emulator_info_list$input_bounds[1,], upper=emulator_info_list$input_bounds[2,], log=TRUE)
    alpha <- min(1.0, exp(lpost_theta_vals[2] - lpost_theta_vals[1] + q_prop_curr - q_curr_prop))
      
    if(runif(1) <= alpha) {
      theta_samp[itr,] <- theta_prop
      theta_scaled_curr <- theta_scaled_prop
      lprior_theta_curr <- lprior_theta_prop
      SSR_curr <- SSR_samples[2,,drop=FALSE]
      pred_idx_curr <- 2
      accept_count <- accept_count + 1 
    } else {
      theta_samp[itr,] <- theta_samp[itr-1,]
      SSR_curr <- SSR_samples[1,,drop=FALSE]
      pred_idx_curr <- 1
    }

    # Adapt proposal covariance matrix and scaling term.
    if(adapt_cov || adapt_scale) {
      adapt_list <- adapt_cov_prop(adapt_cov, adapt_scale, Cov_prop, L_prop, log_scale_prop, theta_samp, itr,  
                                   accept_count, alpha, samp_mean, effective_log_scale_prop, adapt_frequency,  
                                   accept_rate_target, proposal_scale_decay, proposal_scale_multiplier, adapt_init_threshold)
      
      Cov_prop <- adapt_list$C
      L_prop <- adapt_list$L
      log_scale_prop <- adapt_list$log_scale
      effective_log_scale_prop <- adapt_list$effective_log_scale
      samp_mean <- adapt_list$samp_mean
      accept_count <- adapt_list$accept_count
      cov_prop_scales[itr,] <- exp(effective_log_scale_prop)
    }
    
    #
    # Gibbs step for sig2_eps. 
    #
    
    if(learn_sig_eps) {
      if(second_gibbs_step) {
        SSR_curr <- sample_independent_GPs_pointwise(gp_pred_list, transformation_methods=emulator_info$settings$transformation_method,
                                                     idx_selector=pred_idx_curr, include_nugget=TRUE)
      }
      
      sig2_eps_curr <- sample_cond_post_Sig_eps(SSR=SSR_curr, Sig_eps_prior_params=sig_eps_prior_params, 
                                                n_obs=computer_model_data$n_obs)
      sig2_eps_samp[itr,] <- sig2_eps_curr
    } else {
      sig2_eps_samp[itr,] <- sig2_eps_init
    }
    
  }
  
  return(list(theta=theta_samp, sig_eps=sig2_eps_samp, Cov_prop=Cov_prop, cov_prop_scale=cov_prop_scales))
  
}























