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

adapt_MH_proposal_cov <- function(cov_prop, log_scale_prop, times_adapted, adapt_cov, adapt_scale, 
                                  samp_interval, accept_rate, accept_rate_target=0.24,  
                                  adapt_factor_exponent=0.8, adapt_factor_numerator=10) {

  return_list <- list()
  adapt_factor <- 1 / (times_adapted + 3)^adapt_factor_exponent
  
  if(adapt_cov) {
    sample_cov_interval <- cov(samp_interval)
    cov_prop <- cov_prop + adapt_factor * (sample_cov_interval - cov_prop)
    L_cov_prop <- t(chol(cov_prop))
    return_list$L_cov <- L_cov_prop
  }
  
  if(adapt_scale) {
    log_scale_factor <- adapt_factor_numerator * adapt_factor * (accept_rate - accept_rate_target)
    log_scale_prop <- log_scale_prop + log_scale_factor
  }
  
  return_list$cov <- cov_prop
  return_list$log_scale <- log_scale_prop
  
  return(return_list)
}


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


#
# New functions (using new class structure). 
#

# Need to check that the llik_emulator is correct for this function. 
# Infer which sig2 need to be learned by the `lik_par_fixed` attribute.
# How to map `sig2_prior+_params` onto the proper parameters? 
# TODO: need to have checking in llikSumEmulator to make sure the same input parameters are 
# used for each term. 
mcmc_gp_noisy <- function(llik_emulator, par_prior_params, par_init=NULL, sig2_init=NULL, 
                          sig2_prior_params=NULL, N_itr=50000, cov_prop=NULL, 
                          log_scale_prop=NULL, mode="MCMH", use_gp_cov=FALSE,
                          adapt_cov_prop=TRUE, adapt_scale_prop=TRUE, 
                          adapt=adapt_cov_prop||adapt_scale_prop, accept_rate_target=0.24, 
                          adapt_factor_exponent=0.8, adapt_factor_numerator=10, adapt_interval=200, ...) {
  # TODO: Assuming `par_prior_params` is already truncated. Is this the best approach? 
  
  # Validation and setup for log-likelihood emulator. 
  # TODO: ensure that `sig2_init` is named vector, if non-NULL. And that the names include the 
  # names of the outputs where sig2 must be learned. Ensure `par_init` has names as well. 
  # validate_args_mcmc_gp_noisy(llik_emulator, par_prior_params, par_init, sig2_prior_params, N_itr,
  #                             cov_prop, adapt_cov_prop, adapt_scale, use_gp_cov)
  
  # Objects to store samples. 
  d <- llik_emulator$dim_input
  par_samp <- matrix(nrow=N_itr, ncol=d)
  colnames(par_samp) <- llik_emulator$input_names
  
  # Set initial conditions. 
  if(is.null(par_init)) par_init <- sample_prior_theta(par_prior_params)
  par_samp[1,] <- drop(par_init)
  par_curr <- par_samp[1,]
  lprior_par_curr <- calc_lprior_theta(par_curr, par_prior_params)
  
  # Setup for `sig2` (observation variances). Safe to assume that all of the 
  # non-fixed likelihood parameters are `sig2` since this is verified by 
  # `validate_args_mcmc_gp_noisy()` above. 
  learn_sig2 <- !unlist(llik_emulator$get_llik_term_attr("use_fixed_lik_par"))
  term_labels_learn_sig2 <- names(learn_sig2)[learn_sig2]
  include_sig2_Gibbs_step <- (length(term_labels_learn_sig2) > 0)
  sig2_curr <- sig2_init[term_labels_learn_sig2] # Only includes non-fixed variance params.

  if(include_sig2_Gibbs_step) {
    N_obs <- unlist(llik_emulator$get_llik_term_attr("N_obs", labels=term_labels_learn_sig2))
    sig2_samp <- matrix(nrow=N_itr, ncol=length(sig2_curr_learn))
    sig2_samp[1,] <- sig2_curr
  } else {
    sig2_curr <- NULL
    sig2_samp <- NULL
  }
  
  # Proposal covariance.
  if(is.null(cov_prop)) cov_prop <- diag(rep(1,d))
  if(is.null(log_scale_prop)) log_scale_prop <- log(2.38) - 0.5*log(d)
  L_cov_prop <- t(chol(cov_prop))
  accept_count <- 0
  times_adapted <- 0
  
  for(itr in 2:N_itr) {
    #
    # Metropolis step for calibration parameters.
    #
    
    # Random walk proposal. 
    par_prop <- par_curr + (exp(log_scale_prop) * L_cov_prop %*% matrix(rnorm(d), ncol=1))[,1]

    # Compute prior. 
    lprior_par_prop <- calc_lprior_theta(par_prop, par_prior_params)
    
    # Sample SSR. 
    emulator_samp_list <- llik_emulator$sample_emulator(rbind(par_curr,par_prop), use_cov=use_gp_cov, 
                                                        include_nugget=TRUE, ...)  

    # Immediately reject if proposal has prior density zero (which will often happen when the 
    # prior has been truncated to stay within the design bounds). 
    if(is.infinite(lprior_par_prop)) {
      par_samp[itr,] <- par_samp[itr-1,]
      SSR_idx <- 1
    } else {
      # Sample log-likelihood emulator. 
      llik_samp <- llik_emulator$assemble_llik(emulator_samp_list, lik_par_val=sig2_curr, conditional=TRUE,
                                               normalize=FALSE)

      # Accept-Reject step.
      lpost_par_curr <- lprior_par_curr + llik_samp[1]
      lpost_par_prop <- lprior_par_prop + llik_samp[2]
      alpha <- min(1.0, exp(lpost_par_prop - lpost_par_curr))

      if(runif(1) <= alpha) {
        par_samp[itr,] <- par_prop
        par_curr <- par_prop
        lprior_par_curr <- lprior_par_prop
        SSR_idx <- 2 
        accept_count <- accept_count + 1 
      } else {
        par_samp[itr,] <- par_curr
        SSR_idx <- 1
      }
      
    }
    
    # Adapt proposal covariance matrix and scaling term.
    if(adapt && (((itr-1) %% adapt_interval) == 0)) {
      times_adapted <- times_adapted + 1
      adapt_list <- adapt_MH_proposal_cov(cov_prop=cov_prop, log_scale_prop=log_scale_prop, 
                                          times_adapted=times_adapted, 
                                          adapt_cov=adapt_cov_prop, adapt_scale=adapt_scale_prop,
                                          samp_interval=par_samp[(itr-adapt_interval+1):itr,,drop=FALSE], 
                                          accept_rate=accept_count/adapt_interval, accept_rate_target, 
                                          adapt_factor_exponent, adapt_factor_numerator)
      cov_prop <- adapt_list$cov
      log_scale_prop <- adapt_list$log_scale
      if(adapt_cov_prop) L_cov_prop <- adapt_list$L_cov
      accept_count <- 0
    }
    
    #
    # Gibbs step for sig2. 
    #
    
    if(include_sig2_Gibbs_step) {
      # Update sum of squared residuals (SSR) sample. 
      SSR_curr <- sapply(term_labels_learn_sig2, function(lbl) emulator_samp_list[[lbl]][SSR_idx,])
      
      # TODO: how to deal with parameter ordering here? 
      # SSR_curr, sig2_prior_info, N_obs should all use llik labels. 
      sig2_curr <- sample_NIG_cond_post_sig2(SSR_curr, sig2_prior_info, N_obs)
      sig2_samp[itr,] <- sig2_curr
    }
    
  }
  
  return(list(par=par_samp, sig2=sig2_samp))
  
}


mcmc_gp_unn_post_dens_approx <- function(llik_emulator, par_prior_params, par_init=NULL, sig2_init=NULL, 
                                         sig2_prior_params=NULL, N_itr=50000, cov_prop=NULL, 
                                         log_scale_prop=NULL, approx_type="marginal",
                                         adapt_cov_prop=TRUE, adapt_scale_prop=TRUE, 
                                         adapt=adapt_cov_prop||adapt_scale_prop, accept_rate_target=0.24, 
                                         adapt_factor_exponent=0.8, adapt_factor_numerator=10, adapt_interval=200, ...) {
  # GP-accelerated MCMC algorithms which define an approximate posterior distribution by approximating 
  # the unnormalized posterior density. Examples of this include the "mean" approximation, in which 
  # the GP predictive mean is simply plugged into the unnormalized posterior density, and the "marginal" 
  # approximation, in which the expectation of the unnormalized posterior density is taken with 
  # respect to the GP. 
  # TODO: Assuming `par_prior_params` is already truncated. Is this the best approach? 
  # TODO: need to update the name of this function and the above one to reflect that these functions are valid
  # only for Gaussian likelihoods with sig2 params assigned independent IG priors, or likelihoods with 
  # fixed likelihood parameters. 
  # Supported values for `approx_type`: "marginal", "mean". 
  
  # Validation and setup for log-likelihood emulator. 
  # TODO: ensure that `sig2_init` is named vector, if non-NULL. And that the names include the 
  # names of the outputs where sig2 must be learned. Ensure `par_init` has names as well. 
  # validate_args_mcmc_gp_deterministic_approx(llik_emulator, par_prior_params, par_init, sig2_prior_params, N_itr,
  #                                           cov_prop, adapt_cov_prop, adapt_scale, use_gp_cov)
  
  # This should be moved to the argument validation function, once it is written. 
  if(approx_type=="marginal") assert_that(llik_emulator$llik_pred_dist == "Gaussian")
  
  # Objects to store samples. 
  d <- llik_emulator$dim_input
  par_samp <- matrix(nrow=N_itr, ncol=d)
  colnames(par_samp) <- llik_emulator$input_names
  
  # Setup for `sig2` (observation variances). Safe to assume that all of the 
  # non-fixed likelihood parameters are `sig2` since this is verified by 
  # `validate_args_mcmc_gp_noisy()` above. 
  learn_sig2 <- !unlist(llik_emulator$get_llik_term_attr("use_fixed_lik_par"))
  term_labels_learn_sig2 <- names(learn_sig2)[learn_sig2]
  include_sig2_Gibbs_step <- (length(term_labels_learn_sig2) > 0)
  sig2_curr <- sig2_init[term_labels_learn_sig2] # Only includes non-fixed variance params.
  
  if(include_sig2_Gibbs_step) {
    .NotYetImplemented()
    
    N_obs <- unlist(llik_emulator$get_llik_term_attr("N_obs", labels=term_labels_learn_sig2))
    sig2_samp <- matrix(nrow=N_itr, ncol=length(sig2_curr_learn))
    sig2_samp[1,] <- sig2_curr
  } else {
    sig2_curr <- NULL
    sig2_samp <- NULL
  }
  
  # Set initial conditions. 
  if(is.null(par_init)) par_init <- sample_prior_theta(par_prior_params)
  par_samp[1,] <- drop(par_init)
  par_curr <- par_samp[1,]
  lpost_pred_curr <- get_gp_lpost_approx(matrix(par_curr, nrow=1, dimnames=list(NULL, llik_emulator$input_names)), 
                                         approx_type, llik_emulator, par_prior_params,  
                                         lik_par_val=sig2_curr, conditional=TRUE, normalize=FALSE, ...)
  
  # Proposal covariance.
  if(is.null(cov_prop)) cov_prop <- diag(rep(1,d))
  if(is.null(log_scale_prop)) log_scale_prop <- log(2.38) - 0.5*log(d)
  L_cov_prop <- t(chol(cov_prop))
  accept_count <- 0
  times_adapted <- 0
  
  for(itr in 2:N_itr) {
    #
    # Metropolis step for calibration parameters.
    #
    
    # Random walk proposal. 
    par_prop <- par_curr + (exp(log_scale_prop) * L_cov_prop %*% matrix(rnorm(d), ncol=1))[,1]
    
    # Compute prior. 
    lprior_par_prop <- calc_lprior_theta(par_prop, par_prior_params)
    
    # Immediately reject if proposal has prior density zero (which will often happen when the 
    # prior has been truncated to stay within the design bounds). 
    if(is.infinite(lprior_par_prop)) {
      par_samp[itr,] <- par_samp[itr-1,]
      SSR_idx <- 1
    } else {
      # Compute log-posterior approximation. 
      lpost_pred_prop <- get_gp_lpost_approx(matrix(par_prop, nrow=1, dimnames=list(NULL, llik_emulator$input_names)), 
                                             approx_type, llik_emulator, par_prior_params,  
                                             lik_par_val=sig2_curr, conditional=TRUE, normalize=FALSE, ...)
      
      # Accept-Reject step.
      alpha <- min(1.0, exp(lpost_pred_prop - lpost_pred_curr))
      
      if(runif(1) <= alpha) {
        par_samp[itr,] <- par_prop
        par_curr <- par_prop
        lpost_pred_curr <- lpost_pred_prop
        SSR_idx <- 2 
        accept_count <- accept_count + 1 
      } else {
        par_samp[itr,] <- par_curr
        SSR_idx <- 1
      }
      
      # Adapt proposal covariance matrix and scaling term.
      if(adapt && (((itr-1) %% adapt_interval) == 0)) {
        times_adapted <- times_adapted + 1
        adapt_list <- adapt_MH_proposal_cov(cov_prop=cov_prop, log_scale_prop=log_scale_prop, 
                                            times_adapted=times_adapted, 
                                            adapt_cov=adapt_cov_prop, adapt_scale=adapt_scale_prop,
                                            samp_interval=par_samp[(itr-adapt_interval+1):itr,,drop=FALSE], 
                                            accept_rate=accept_count/adapt_interval, accept_rate_target, 
                                            adapt_factor_exponent, adapt_factor_numerator)
        cov_prop <- adapt_list$cov
        log_scale_prop <- adapt_list$log_scale
        if(adapt_cov_prop) L_cov_prop <- adapt_list$L_cov
        accept_count <- 0
      }
    }
    
    
    #
    # Gibbs step for sig2. 
    #
    
    if(include_sig2_Gibbs_step) {
      .NotYetImplemented()
    }
    
  }
  
  return(list(par=par_samp, sig2=sig2_samp))
  
}


mcmc_gp_acc_prob_approx <- function(llik_emulator, par_prior_params, par_init=NULL, sig2_init=NULL, 
                                    sig2_prior_params=NULL, N_itr=50000, cov_prop=NULL, 
                                    log_scale_prop=NULL, approx_type="marginal",
                                    adapt_cov_prop=TRUE, adapt_scale_prop=TRUE, 
                                    adapt=adapt_cov_prop||adapt_scale_prop, accept_rate_target=0.24, 
                                    adapt_factor_exponent=0.8, adapt_factor_numerator=10, adapt_interval=200, ...) {
  # TODO: I need to generalize this to allow for forward model emulators as well. 
  # TODO: I need to think about how to do the sig2 step with this approach. There seems to be no reason why 
  # there could not be different options for this. 
  #
  # GP-accelerated MCMC algorithms that approximate the acceptance probability of a Metropolis-Hastings MCMC
  # algorithm. Examples of this include the "mean" approximation, in which 
  # the GP predictive mean is simply plugged into acceptance ratio, the "joint-marginal", in which 
  # the expectation of the acceptance probability is taken with respect to the GP distribution, and the 
  # "marginal" approximation, which is like the "joint-marginal" but the GP covariance is ignored (i.e. set to 0). 
  # TODO: Assuming `par_prior_params` is already truncated. Is this the best approach? 
  # TODO: need to update the name of this function and the above one to reflect that these functions are valid
  # only for Gaussian likelihoods with sig2 params assigned independent IG priors, or likelihoods with 
  # fixed likelihood parameters. 
  # Supported values for `approx_type`: "mean", "marginal", "joint-marginal"
  
  # Validation and setup for log-likelihood emulator. 
  # TODO: ensure that `sig2_init` is named vector, if non-NULL. And that the names include the 
  # names of the outputs where sig2 must be learned. Ensure `par_init` has names as well. 
  # validate_args_mcmc_gp_deterministic_approx(llik_emulator, par_prior_params, par_init, sig2_prior_params, N_itr,
  #                                           cov_prop, adapt_cov_prop, adapt_scale, use_gp_cov)
  
  # This should be moved to the argument validation function, once it is written. 
  if(approx_type=="marginal") assert_that(llik_emulator$llik_pred_dist == "Gaussian")
  
  # Objects to store samples. 
  d <- llik_emulator$dim_input
  par_samp <- matrix(nrow=N_itr, ncol=d)
  colnames(par_samp) <- llik_emulator$input_names
  
  # Setup for `sig2` (observation variances). Safe to assume that all of the 
  # non-fixed likelihood parameters are `sig2` since this is verified by 
  # `validate_args_mcmc_gp_noisy()` above. 
  learn_sig2 <- !unlist(llik_emulator$get_llik_term_attr("use_fixed_lik_par"))
  term_labels_learn_sig2 <- names(learn_sig2)[learn_sig2]
  include_sig2_Gibbs_step <- (length(term_labels_learn_sig2) > 0)
  sig2_curr <- sig2_init[term_labels_learn_sig2] # Only includes non-fixed variance params.
  
  if(include_sig2_Gibbs_step) {
    .NotYetImplemented()
    
    N_obs <- unlist(llik_emulator$get_llik_term_attr("N_obs", labels=term_labels_learn_sig2))
    sig2_samp <- matrix(nrow=N_itr, ncol=length(sig2_curr_learn))
    sig2_samp[1,] <- sig2_curr
  } else {
    sig2_curr <- NULL
    sig2_samp <- NULL
  }
  
  # Set initial conditions. 
  if(is.null(par_init)) par_init <- sample_prior_theta(par_prior_params)
  par_samp[1,] <- drop(par_init)
  par_curr <- par_samp[1,]
  lprior_par_curr <- calc_lprior_theta(par_curr, par_prior_params)
  
  # Proposal covariance.
  if(is.null(cov_prop)) cov_prop <- diag(rep(1,d))
  if(is.null(log_scale_prop)) log_scale_prop <- log(2.38) - 0.5*log(d)
  L_cov_prop <- t(chol(cov_prop))
  accept_count <- 0
  times_adapted <- 0
  
  for(itr in 2:N_itr) {
    #
    # Metropolis step for calibration parameters.
    #
    
    # Random walk proposal. 
    par_prop <- par_curr + (exp(log_scale_prop) * L_cov_prop %*% matrix(rnorm(d), ncol=1))[,1]
    
    # Compute prior at proposed location. 
    lprior_par_prop <- calc_lprior_theta(par_prop, par_prior_params)
    
    # Immediately reject if proposal has prior density zero (which will often happen when the 
    # prior has been truncated to stay within the design bounds). 
    if(is.infinite(lprior_par_prop)) {
      par_samp[itr,] <- par_samp[itr-1,]
      SSR_idx <- 1
    } else {
      # Compute acceptance probability approximation. 
      alpha <- get_gp_mh_acc_prob_approx(par_curr, par_prop, llik_emulator, approx_type, 
                                         lik_par_val=sig2_curr, conditional=TRUE, normalize=FALSE, 
                                         lprior_vals=c(lprior_par_curr, lprior_par_prop), ...)
    
      # Accept-Reject step. 
      if(runif(1) <= alpha) {
        par_samp[itr,] <- par_prop
        par_curr <- par_prop
        lprior_par_curr <- lprior_par_prop
        accept_count <- accept_count + 1 
      } else {
        par_samp[itr,] <- par_curr
      }
      
      # Adapt proposal covariance matrix and scaling term.
      if(adapt && (((itr-1) %% adapt_interval) == 0)) {
        times_adapted <- times_adapted + 1
        adapt_list <- adapt_MH_proposal_cov(cov_prop=cov_prop, log_scale_prop=log_scale_prop, 
                                            times_adapted=times_adapted, 
                                            adapt_cov=adapt_cov_prop, adapt_scale=adapt_scale_prop,
                                            samp_interval=par_samp[(itr-adapt_interval+1):itr,,drop=FALSE], 
                                            accept_rate=accept_count/adapt_interval, accept_rate_target, 
                                            adapt_factor_exponent, adapt_factor_numerator)
        cov_prop <- adapt_list$cov
        log_scale_prop <- adapt_list$log_scale
        if(adapt_cov_prop) L_cov_prop <- adapt_list$L_cov
        accept_count <- 0
      }
    }
    
    #
    # Gibbs step for sig2. 
    #
    
    if(include_sig2_Gibbs_step) {
      .NotYetImplemented()
    }
    
  }
  
  return(list(par=par_samp, sig2=sig2_samp))
  
}


mcmc_gp_1d_discretization <- function(llik_emulator, par_prior_params, par_init=NULL, sig2_init=NULL, 
                                      sig2_prior_params=NULL, N_itr=50000, cov_prop=NULL, 
                                      log_scale_prop=NULL, N_grid_points=NULL,
                                      adapt_cov_prop=TRUE, adapt_scale_prop=TRUE, 
                                      adapt=adapt_cov_prop||adapt_scale_prop, accept_rate_target=0.24, 
                                      adapt_factor_exponent=0.8, adapt_factor_numerator=10, adapt_interval=200, ...) {
  # `N_grid_points` is the number of grid points in the discretization. 
  
  # This should be moved to the argument validation function, once it is written. 
  if(approx_type=="marginal") assert_that(llik_emulator$llik_pred_dist == "Gaussian")
  
  # Objects to store samples. 
  d <- llik_emulator$dim_input
  par_samp <- matrix(nrow=N_itr, ncol=d)
  colnames(par_samp) <- llik_emulator$input_names
  
  # Setup for `sig2` (observation variances). Safe to assume that all of the 
  # non-fixed likelihood parameters are `sig2` since this is verified by 
  # `validate_args_mcmc_gp_noisy()` above. 
  learn_sig2 <- !unlist(llik_emulator$get_llik_term_attr("use_fixed_lik_par"))
  term_labels_learn_sig2 <- names(learn_sig2)[learn_sig2]
  include_sig2_Gibbs_step <- (length(term_labels_learn_sig2) > 0)
  sig2_curr <- sig2_init[term_labels_learn_sig2] # Only includes non-fixed variance params.
  
  if(include_sig2_Gibbs_step) {
    .NotYetImplemented()
    
    N_obs <- unlist(llik_emulator$get_llik_term_attr("N_obs", labels=term_labels_learn_sig2))
    sig2_samp <- matrix(nrow=N_itr, ncol=length(sig2_curr_learn))
    sig2_samp[1,] <- sig2_curr
  } else {
    sig2_curr <- NULL
    sig2_samp <- NULL
  }
  
  # Set initial conditions. 
  if(is.null(par_init)) par_init <- sample_prior_theta(par_prior_params)
  par_samp[1,] <- drop(par_init)
  par_curr <- par_samp[1,]
  lprior_par_curr <- calc_lprior_theta(par_curr, par_prior_params)
  
  # Proposal covariance.
  if(is.null(cov_prop)) cov_prop <- diag(rep(1,d))
  if(is.null(log_scale_prop)) log_scale_prop <- log(2.38) - 0.5*log(d)
  L_cov_prop <- t(chol(cov_prop))
  accept_count <- 0
  times_adapted <- 0
  
  # Grid points for GP discretization. 
  grid <- get_tensor_product_grid(N_batch, prior_dist_info=par_prior_params)
                          
  # Gibbs step: GP. 
  
  
  # MH step: par. 
  
}



# ---------------------------------------------------------------------
# Unnormalized log posterior functions using log-likelihood emulators. 
# ---------------------------------------------------------------------

get_gp_lpost_approx <- function(input, approx_type, llik_emulator, par_prior_params, lik_par_val=NULL, 
                                conditional=FALSE, normalize=TRUE, llik_pred_list=NULL, ...) {
  if(approx_type == "mean") {
    return(gp_lpost_mean(input, llik_emulator, par_prior_params, lik_par_val=lik_par_val, 
                         conditional=conditional, normalize=normalize, llik_pred_list=llik_pred_list, ...))
  } else if (approx_type == "marginal") {
    return(gp_lpost_marginal(input, llik_emulator, par_prior_params, lik_par_val=lik_par_val, 
                             conditional=conditional, normalize=normalize, llik_pred_list=llik_pred_list, ...))
  } else {
    stop("Invalid `approx_type` ", approx_type)
  }
}

gp_lpost_mean <- function(input, llik_emulator, par_prior_params, lik_par_val=NULL, 
                          conditional=FALSE, normalize=TRUE, llik_pred_list=NULL, ...) {
  
  assert_that(llik_emulator$llik_pred_dist == "Gaussian")
  
  if(is.null(llik_pred_list)) {
    llik_vals <- drop(llik_emulator$predict(input, lik_par_val=lik_par_val, return_mean=TRUE, return_var=FALSE, 
                                            conditional=conditional, normalize=normalize, ...)$mean)
  } else {
    llik_vals <- drop(llik_pred_list$mean)
    assert_that(!is.null(llik_vals))
  }
  
  lprior_vals <- calc_lprior_theta(input, par_prior_params)
  
  return(llik_vals + lprior_vals)

}


gp_lpost_marginal <- function(input, llik_emulator, par_prior_params, lik_par_val=NULL, 
                              conditional=FALSE, normalize=TRUE, llik_pred_list=NULL, ...) {
  
  assert_that(llik_emulator$llik_pred_dist == "Gaussian")
  
  if(is.null(llik_pred_list)) {
    llik_pred_list <- llik_emulator$predict(input, lik_par_val=lik_par_val, return_mean=TRUE, return_var=TRUE, 
                                            conditional=conditional, normalize=normalize, ...)
  } else {
    assert_that(!is.null(llik_pred_list$mean) && !is.null(llik_pred_list$var))
  }
  
  llik_vals <- convert_Gaussian_to_LN(mean_Gaussian=llik_pred_list$mean, var_Gaussian=llik_pred_list$var,
                                      return_mean=TRUE, return_var=FALSE, log_scale=TRUE)$log_mean
  lprior_vals <- calc_lprior_theta(input, par_prior_params)
  
  return(llik_vals + lprior_vals)
  
}


# ---------------------------------------------------------------------
# Metropolis-Hastings acceptance probability approximations using 
# log-likelihood emulators. 
# ---------------------------------------------------------------------

get_gp_mh_acc_prob_approx <- function(input_curr, input_prop, llik_emulator, approx_type, 
                                      lik_par_val=sig2_curr, conditional=TRUE, normalize=FALSE, 
                                      par_prior_params=NULL, lprior_vals=NULL, llik_pred_list=NULL, ...) {
  # TODO: currently assumes symmetric proposal density, can generalize this when needed. 
  
  # Compute log-prior evaluations, if not already provided. 
  if(is.null(lprior_vals)) lprior_vals <- calc_lprior_theta(rbind(par_curr, par_prop), par_prior_params)
  
  if(approx_type == "mean") {
    .NotYetImplemented()
  } else if(approx_type %in% c("marginal", "joint-marginal")) {
    joint <- (approx_type == "joint-marginal")
    alpha <- gp_acc_prob_marginal(input_curr, input_prop, llik_emulator, joint, lik_par_val=lik_par_val, 
                                  conditional=conditional, normalize=normalize, llik_pred_list=llik_pred_list, 
                                  lprior_vals=lprior_vals, ...)
  } else {
    stop("Invalid `approx_type` ", approx_type)
  }
  
  return(alpha)
  
}


gp_acc_ratio_marginal <- function(input_curr, input_prop, llik_emulator, use_joint_dist, lik_par_val=NULL, 
                                  conditional=llik_emulator$default_conditional, 
                                  normalize=llik_emulator$default_normalize, llik_pred_list=NULL, 
                                  par_prior_params=NULL, lprior_vals=NULL, ...) {
  # For a Metropolis-Hastings acceptance probability of the form min{1, r(input_curr, input_prop)}, 
  # this function returns log E[r(input_curr, input_prop)], where the expectation is with respect to 
  # the `llik_emulator` distribution. Currently this function only works when the emulator distribution 
  # implies that r(input_curr, input_prop) ~ LN(m, s^2). 
  
  # Compute log-prior evaluations, if not already provided. 
  if(is.null(lprior_vals)) lprior_vals <- calc_lprior_theta(rbind(input_curr, input_prop), par_prior_params)
  
  # Compute predictive distribution at inputs.   
  if(is.null(llik_pred_list)) {
    input <- rbind(input_curr, input_prop)
    rownames(input) <- llik_emulator$input_names
    llik_pred_list <- llik_emulator$predict(input, lik_par_val=lik_par_val, return_mean=TRUE, return_var=TRUE, 
                                            return_cov=use_joint_dist, conditional=conditional, normalize=normalize, ...)
  } else {
    assert_that(!is.null(llik_pred_list$mean) && !is.null(llik_pred_list$var))
  }
  
  # Predictive mean and variance of the log acceptance ratio. This is the m and s^2 in 
  # r ~ LN(m, s^2). The variance depends on whether or not the joint dist is used. 
  # TODO: update the cov indexing below once this is changed. 
  m <- lprior_vals[2] - lprior_vals[1] + drop(llik_pred_list$mean)[2] - drop(llik_pred_list$mean)[1]
  s2 <- drop(llik_pred_list$var)[2] + drop(llik_pred_list$var)[1]
  if(use_joint_dist) s2 <- s2 - 2*llik_pred_list$cov[,,1][1,2]

  # Marginal acceptance ratio approximation. 
  return(exp(m + 0.5*s2))

}


gp_acc_prob_marginal <- function(input_curr, input_prop, llik_emulator, use_joint_dist, lik_par_val=NULL, 
                                 conditional=llik_emulator$default_conditional, 
                                 normalize=llik_emulator$default_normalize, llik_pred_list=NULL, 
                                 par_prior_params=NULL, lprior_vals=NULL, ...) {
  # For a Metropolis-Hastings acceptance probability of the form 
  # alpha(input_curr, input_prop) min{1, r(input_curr, input_prop)}, 
  # this function returns E[alpha(input_curr, input_prop)], where the expectation is with respect to 
  # the `llik_emulator` distribution.
  
  # Compute log-prior evaluations, if not already provided. 
  if(is.null(lprior_vals)) lprior_vals <- calc_lprior_theta(rbind(input_curr, input_prop), par_prior_params)
  
  # Compute predictive distribution at inputs.   
  if(is.null(llik_pred_list)) {
    input <- rbind(input_curr, input_prop)
    colnames(input) <- llik_emulator$input_names
    llik_pred_list <- llik_emulator$predict(input, lik_par_val=lik_par_val, return_mean=TRUE, return_var=TRUE, 
                                            return_cov=use_joint_dist, conditional=conditional, normalize=normalize, ...)
  } else {
    assert_that(!is.null(llik_pred_list$mean) && !is.null(llik_pred_list$var))
  }
  
  # Predictive mean and variance of the log acceptance ratio. This is the m and s^2 in 
  # r ~ LN(m, s^2). The variance depends on whether or not the joint dist is used. 
  # TODO: update the cov indexing below once this is changed. 
  m <- lprior_vals[2] - lprior_vals[1] + drop(llik_pred_list$mean)[2] - drop(llik_pred_list$mean)[1]
  s2 <- drop(llik_pred_list$var)[2] + drop(llik_pred_list$var)[1]
  if(use_joint_dist) s2 <- s2 - 2*llik_pred_list$cov[,,1][1,2]
  
  # Marginal acceptance probability is a linear combination of 1 and the marginal ratio approximation.
  log_acc_ratio_marginal <- m + 0.5*s2
  lw1 <- log(pnorm(m/sqrt(s2)))
  lw2 <- log(pnorm(-(m+s2)/sqrt(s2)))
  acc_ratio_marginal <- exp(matrixStats::logSumExp(c(lw1, lw2+log_acc_ratio_marginal)))
  
  return(acc_ratio_marginal)
  
}



# ---------------------------------------------------------------------
# Helper functions 
# ---------------------------------------------------------------------

run_llik_emulator_samplers <- function(sampler_settings_list, llik_emulator_default=NULL, 
                                       par_prior=NULL,  lik_par_prior_default=NULL, N_mcmc_default=NULL) {
  # Need to think about how to implement this. I want the ability to easily specify common settings that 
  # are shared across all algorithms, but also want the ability to have complete control over the 
  # settings for each algorithm, if desired. I'm currently thinking that there should be arguments 
  # to specify "defaults", which will be overwritten if they are also specified in "sampler_settings_list". 
  
  .NotYetImplemented()
}











