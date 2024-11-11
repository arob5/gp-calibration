#
# gp_mcmc_functions.r
# MCMC algorithms utilizing a Gaussian process (GP) emulator to accelerate 
# inference. 
#
# Andrew Roberts
#
# Depends:
#    gpWrapper, llikEmulator, general_helper_functions, mcmc_helper_functions

library(parallel)

# -----------------------------------------------------------------------------
# Code Standards and Conventions:
#
# This file contains MCMC algorithms designed to function with objects from the   
# `llikEmulator` class. The emphasis is on log-likelihood emulators constructed
# from Gaussian process (GP) approximations. The philosophy is for these 
# algorithms to be able to work with an exact (non-emulated) likelihood as well.
# That is, for llikEmulator objects with attribute `exact_llik = TRUE`, the 
# approximate MCMC algorithms here should sample from the exact posterior.
#
# Note that there are multiple types of parameters that may be samples in an 
# MCMC algorithm. The main parameters of interest here are those that 
# correspond to the inputs of the log-likelihood emulator, and are by 
# convention referred to as "par" in this file. The remaining set of parameters
# are typically those defining the likelihood; e.g., a variance parameter. 
# These are referred to as "likelihood parameters".
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Core MCMC algorithms. 
#
# At present, there are three primary approximate MCMC schemes, as well as one 
# wrapper function to allow access to the samplers in the `BayesianTools` R 
# package for exact MCMC.
#    1.) mcmc_noisy_llik: An MCMC algorithm defined with respect to a stochastic
#        approximation of the log-likelihood. This function implements methods 
#        that only require sampling from the log-likelihood surrogate, thus 
#        representing "noisy" versions of a traditional Metropolis-Hastings
#        scheme. The function arguments allow for a user to run a Monte Carlo 
#        within Metropolis or a pseudo marginal algorithm.
#    2.) mcmc_lik_approx: A Metropolis-Hastings algorithm defined with respect 
#        to a deterministic approximation of the true likelihood (alternatively,
#        can view this as an approximation to the unnormalized posterior 
#        density). This function relies of the `calc_lik_approx()` method of 
#        the llikEmulator object.
#    3.) mcmc_acc_prob_approx: An approximate Metropolis-Hastings like algorithm 
#        defined by a deterministic approximation to the true acceptance 
#        probability. This function relies on the helper function 
#        `get_gp_mh_acc_prob_approx()`.
#    4.) mcmc_bt_wrapper: A wrapper that provides access to the MCMC algorithms
#        implemented in the `BayesianTools` package. Hence, this function is 
#        only for inference with exact log-likelihood objects.
#
# At this point, these algorithms are only applicable when the true likelihood
# is of a multiplicative Gaussian form. This is due to the way that the 
# likelihood parameters are treated by the functions, and should be 
# generalized to more general likelihoods in the future.
# -----------------------------------------------------------------------------

# Need to check that the llik_emulator is correct for this function. 
# Infer which sig2 need to be learned by the `lik_par_fixed` attribute.
# How to map `sig2_prior+_params` onto the proper parameters? 
# TODO: need to have checking in llikSumEmulator to make sure the same input parameters are 
# used for each term. 

mcmc_noisy_llik <- function(llik_emulator, par_prior, par_init=NULL, sig2_init=NULL,
                            sig2_prior=NULL, mode="mcwmh", use_joint=TRUE,  
                            n_itr=50000L, cov_prop=NULL, log_scale_prop=NULL, 
                            adapt_cov_prop=TRUE, adapt_scale_prop=TRUE, 
                            adapt=adapt_cov_prop||adapt_scale_prop, 
                            return_prop_sd=FALSE, accept_rate_target=0.24, 
                            adapt_factor_exponent=0.8, adapt_factor_numerator=10, 
                            adapt_interval=200, ...) {
  # Args:
  #    llik_emulator: an object for which `is_llik_em(llik_emulator)` is TRUE.
  #    par_prior: data.frame defining the prior distribution for `par`. 
  #    par_init: numeric, d-dimensional vector, the initial condition. If NULL,
  #              samples from the prior to determine the initial condition.
  #    sig2_init: numeric, the initial condition for the variance parameters. 
  #               Not used if `llik_emulator$use_fixed_lik_par` is TRUE.
  #    sig2_prior: list defining the prior distribution on the variance 
  #                parameters.
  #    mode: character, either "mcwmh" or "pseudo-marg" to run either the Monte
  #          Carlo within Metropolis-Hastings or Pseudo-Marginal algorithm.
  #    use_joint: logical, whether to utilize the log-likelihood emulator
  #               joint distribution across different parameter values, or treat
  #               the emulator as independent on a parameter-by-parameter basis.
  #    n_itr: integer, the number of MCMC iterations.
  #    return_prop_sd: logical, if TRUE returns the trajectory of the standard 
  #                    deviations of the proposal distribution; that is, the 
  #                    square root of the diagonal of the proposal covariance.
  #    Remaining arguments define the proposal covariance adaptation, and are 
  #    passed to `adapt_MH_proposal_cov()`.
  #
  # Returns:
  # list, with the main element named "samp". This element is itself a list with
  # elements "par", "sig2", and "prop_sd_comb". Each of these are matrices with 
  # rows storing the MCMC samples. The "prop_sd_comb" matrix contains the 
  # square root of the diagonal of the proposal covariance at each iteration, 
  # and is only returned if `return_prop_sd = TRUE`. Note that this refers to 
  # the entire composite proposal covariance, which combines both `cov_prop`
  # and `log_scale_prop`. The elements of the outer list other than "samp" store 
  # other information regarding the MCMC run (initial conditions, proposal 
  # adaptation settings, etc.). 
  # TODO: ensure that `sig2_init` is named vector, if non-NULL. 
  # And that the names include the names of the outputs where sig2 must be 
  # learned. Ensure `par_init` has names as well. 

  # TODO: for some reason, `inherits(llik_emulator, "llikEmulator")` always 
  # evaluates to FALSE when this is called from parLapply(). Need to figure 
  # out why. 
  # assert_that(is_llik_em(llik_emulator))
  
  if(mode != "mcwmh") {
    .NotYetImplemented(msg="Only `mode=mcwmh` currently implemented.")
  }
  
  # Objects to store samples. 
  d <- llik_emulator$dim_input
  par_samp <- matrix(nrow=n_itr, ncol=d)
  colnames(par_samp) <- llik_emulator$input_names
  
  # Set initial conditions. 
  if(is.null(par_init)) par_init <- sample_prior(par_prior)
  par_samp[1,] <- drop(par_init)
  par_curr <- par_samp[1,]
  lprior <- get_lprior_dens(par_prior, check_bounds=TRUE)
  lprior_par_curr <- lprior(par_curr)
  
  # Setup for `sig2` (observation variances). Safe to assume that all of the 
  # non-fixed likelihood parameters are `sig2` since this is verified by 
  # `validate_args_mcmc_gp_noisy()` above. 
  learn_sig2 <- !unlist(llik_emulator$get_llik_term_attr("use_fixed_lik_par"))
  term_labels_learn_sig2 <- names(learn_sig2)[learn_sig2]
  include_sig2_Gibbs_step <- (length(term_labels_learn_sig2) > 0)
  sig2_curr <- sig2_init[term_labels_learn_sig2] # Only includes non-fixed variance params.
  
  if(include_sig2_Gibbs_step) {
    N_obs <- unlist(llik_emulator$get_llik_term_attr("N_obs", 
                                                     labels=term_labels_learn_sig2))
    sig2_samp <- matrix(nrow=n_itr, ncol=length(sig2_curr_learn))
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
  
  # If tracking the proposal standard deviations.
  if(return_prop_sd) {
    prop_sd_comb <- matrix(nrow=n_itr, ncol=d, 
                           dimnames=list(NULL, llik_emulator$input_names))
    prop_sd_comb[1,] <- exp(log_scale_prop) * sqrt(diag(cov_prop)) 
  } else {
    prop_sd_comb <- NULL
  }
  
  # Main MCMC loop. 
  tryCatch(
    {
      for(itr in 2L:n_itr) {
        #
        # Metropolis step for calibration parameters.
        #
        
        # Random walk proposal. 
        par_prop <- par_curr + (exp(log_scale_prop) * L_cov_prop %*% matrix(rnorm(d), ncol=1))[,1]
        
        # Compute prior. 
        lprior_par_prop <- lprior(par_prop)
        
        # Sample SSR. 
        emulator_samp_list <- llik_emulator$sample_emulator(rbind(par_curr,par_prop), 
                                                            use_cov=use_cov, 
                                                            include_nugget=TRUE, ...)  
        
        # Immediately reject if proposal has prior density zero (which will 
        # often happen when the prior has been truncated to stay within the 
        # design bounds).
        if(is.infinite(lprior_par_prop)) {
          par_samp[itr,] <- par_samp[itr-1,]
          curr_prop_idx <- 1L
        } else {
          # Sample log-likelihood emulator. 
          llik_samp <- llik_emulator$assemble_llik(emulator_samp_list, 
                                                   lik_par_val=sig2_curr)
          
          # Accept-Reject step.
          lpost_par_curr <- lprior_par_curr + llik_samp[1]
          lpost_par_prop <- lprior_par_prop + llik_samp[2]
          alpha <- min(1.0, exp(lpost_par_prop - lpost_par_curr))
          
          if(runif(1) <= alpha) {
            par_samp[itr,] <- par_prop
            par_curr <- par_prop
            lprior_par_curr <- lprior_par_prop
            curr_prop_idx <- 2L 
            accept_count <- accept_count + 1L
          } else {
            par_samp[itr,] <- par_curr
            curr_prop_idx <- 1L
          }
        }
        
        # Adapt proposal covariance matrix and scaling term.
        if(adapt && (((itr-1) %% adapt_interval) == 0)) {
          times_adapted <- times_adapted + 1L
          adapt_list <- adapt_MH_proposal_cov(cov_prop=cov_prop, 
                                              log_scale_prop=log_scale_prop, 
                                              times_adapted=times_adapted, 
                                              adapt_cov=adapt_cov_prop, 
                                              adapt_scale=adapt_scale_prop,
                                              samp_interval=par_samp[(itr-adapt_interval+1):itr,,drop=FALSE], 
                                              accept_rate=accept_count/adapt_interval, 
                                              accept_rate_target, adapt_factor_exponent, 
                                              adapt_factor_numerator)
          cov_prop <- adapt_list$cov
          log_scale_prop <- adapt_list$log_scale
          if(adapt_cov_prop) L_cov_prop <- adapt_list$L_cov
          accept_count <- 0L
        }
        
        if(return_prop_sd) {
          prop_sd_comb[itr,] <- exp(log_scale_prop) * sqrt(diag(crossprod(L_cov_prop)))
        }
        
        #
        # Gibbs step for sig2. 
        #
        
        if(include_sig2_Gibbs_step) {
          # Update sum of squared residuals (SSR) sample. 
          SSR_curr <- sapply(term_labels_learn_sig2, 
                             function(lbl) emulator_samp_list[[lbl]][curr_prop_idx,])
          
          # TODO: how to deal with parameter ordering here? 
          # SSR_curr, sig2_prior_info, N_obs should all use llik labels. 
          sig2_curr <- sample_NIG_cond_post_sig2(SSR_curr, sig2_prior_info, N_obs)
          sig2_samp[itr,] <- sig2_curr
        }
      }
    },  error = function(cond) {
      err <- cond
      message("mcmc_gp_noisy() MCMC error; iteration ", itr)
      message(conditionMessage(cond))
    }, finally = {
    }
  )
  
  return(list(samp=list(par=par_samp, sig2=sig2_samp, 
                        prop_sd_comb=prop_sd_comb), 
              log_scale_prop=log_scale_prop, L_cov_prop=L_cov_prop, 
              par_curr=par_curr, par_prop=par_prop, sig2_curr=sig2_curr,
              itr_curr=itr, condition=err))
  
}


mcmc_gp_unn_post_dens_approx <- function(llik_emulator, par_prior_params, par_init=NULL, sig2_init=NULL, 
                                         sig2_prior_params=NULL, N_itr=50000, cov_prop=NULL, 
                                         log_scale_prop=NULL, approx_type="marginal",
                                         adapt_cov_prop=TRUE, adapt_scale_prop=TRUE, 
                                         adapt=adapt_cov_prop||adapt_scale_prop, accept_rate_target=0.24, 
                                         adapt_factor_exponent=0.8, adapt_factor_numerator=10, adapt_interval=200, 
                                         alpha=0.9, ...) {
  # GP-accelerated MCMC algorithms which define an approximate posterior distribution by approximating 
  # the unnormalized posterior density. Examples of this include the "mean" approximation, in which 
  # the GP predictive mean is simply plugged into the unnormalized posterior density, and the "marginal" 
  # approximation, in which the expectation of the unnormalized posterior density is taken with 
  # respect to the GP. 
  # TODO: Assuming `par_prior_params` is already truncated. Is this the best approach? 
  # TODO: need to update the name of this function and the above one to reflect that these functions are valid
  # only for Gaussian likelihoods with sig2 params assigned independent IG priors, or likelihoods with 
  # fixed likelihood parameters. 
  # Supported values for `approx_type`: "marginal", "mean", "quantile". 
  # `alpha` is only relevant for the "quantile" approximation, and defines the 
  # pointwise quantile of the GP predictive distribution to use in the approximation. 
  
  # Validation and setup for log-likelihood emulator. 
  # TODO: ensure that `sig2_init` is named vector, if non-NULL. And that the names include the 
  # names of the outputs where sig2 must be learned. Ensure `par_init` has names as well. 
  # validate_args_mcmc_gp_deterministic_approx(llik_emulator, par_prior_params, par_init, sig2_prior_params, N_itr,
  #                                           cov_prop, adapt_cov_prop, adapt_scale, use_gp_cov)
  
  # This should be moved to the argument validation function, once it is written. 
  # TODO: need to update this check to include forward model emulator with Gaussian likelihood. 
  # if(approx_type=="marginal") assert_that(llik_emulator$llik_pred_dist == "Gaussian")
  
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
  if(is.null(par_init)) par_init <- sample_prior(par_prior_params)
  par_samp[1,] <- drop(par_init)
  par_curr <- par_samp[1,]
  par_curr_mat <- matrix(par_curr, nrow=1, dimnames=list(NULL, llik_emulator$input_names))
  llik_pred_curr <- llik_emulator$calc_lik_approx(par_curr_mat, approx_type, 
                                                  lik_par_val=sig2_curr, log_scale=TRUE, 
                                                  alpha=alpha, ...) 
  lpost_pred_curr <- llik_pred_curr + calc_lprior_theta(par_curr_mat, par_prior_params)
  
  # Proposal covariance.
  if(is.null(cov_prop)) cov_prop <- diag(rep(1,d))
  if(is.null(log_scale_prop)) log_scale_prop <- log(2.38) - 0.5*log(d)
  L_cov_prop <- t(chol(cov_prop))
  accept_count <- 0
  times_adapted <- 0
  
  tryCatch(
    {
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
          par_prop_mat <- matrix(par_prop, nrow=1, dimnames=list(NULL, llik_emulator$input_names))
          llik_pred_prop <- llik_emulator$calc_lik_approx(par_prop_mat, approx_type, 
                                                          lik_par_val=sig2_curr, log_scale=TRUE, alpha=alpha, ...) 
          lpost_pred_prop <- llik_pred_prop + calc_lprior_theta(par_prop_mat, par_prior_params)
          
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
      
    }, error = function(cond) {
      message("mcmc_gp_unn_post_dens_approx() MCMC error; iteration ", itr)
      message(conditionMessage(cond))
    }, finally = {
      return(list(samp=list(par=par_samp, sig2=sig2_samp), 
                  log_scale_prop=log_scale_prop, L_cov_prop=L_cov_prop, par_curr=par_curr,
                  par_prop=par_prop, sig2_curr=sig2_curr))
    }
  )
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
  if(is.null(par_init)) par_init <- sample_prior(par_prior_params)
  par_samp[1,] <- drop(par_init)
  par_curr <- par_samp[1,]
  lprior_par_curr <- calc_lprior_theta(par_curr, par_prior_params)
  
  # Proposal covariance.
  if(is.null(cov_prop)) cov_prop <- diag(rep(1,d))
  if(is.null(log_scale_prop)) log_scale_prop <- log(2.38) - 0.5*log(d)
  L_cov_prop <- t(chol(cov_prop))
  accept_count <- 0
  times_adapted <- 0
  
  tryCatch({
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
  }, error = function(cond) {
    message("mcmc_gp_acc_prob_approx() MCMC error; iteration ", itr)
    message(conditionMessage(cond))
  }, finally = {
    return(list(samp=list(par=par_samp, sig2=sig2_samp), 
                log_scale_prop=log_scale_prop, L_cov_prop=L_cov_prop, par_curr=par_curr,
                par_prop=par_prop, sig2_curr=sig2_curr))
  }
  )
}


mcmc_bt_wrapper <- function(llik_emulator, par_prior, par_init=NULL, 
                            n_itr=50000L, sampler="DEzs", settings_list=list(), ...) {
  # This function provides a wrapper around the BayesianTools package 
  # `runMCMC()` function, thus providing access to a variety of different 
  # MCMC algorithms. This wrapper is designed so that the argument interface
  # and returned value align with the other MCMC functions in this file.
  # In particular, `mcmc_bt_wrapper()` is set up to run a single chain of 
  # MCMC using a `llik_emulator` object. Parallelization for multichain runs 
  # should be accomplished using `run_mcmc_multichain()`, as for the other 
  # MCMC functions implemented here. The MCMC output from BayesianTools is 
  # transformed to have the same format as the other MCMC functions here (e.g.,
  # `mcmc_noisy_llik()`). Note that the `parallel` element of the BayesianTools
  # settings object refers to within-chain parallelization, for MCMC algorithms
  # for which parallel likelihood evaluations makes sense; it does not control 
  # the running of multiple chains in parallel.
  #
  # Note that BayesianTools appears to require `priorSampler` to be passed to 
  # the bayesianSetup object when the prior is a R function representing a 
  # log prior density. Certain algorithms (e.g., "DEzs") require generating 
  # an initial "population" of samples. In this case, it is recommended to 
  # set `par_init` to NULL and let the BayesianTools code setup the initial 
  # population. 
  #
  # Args:
  #    llik_emulator: an object for which `is_llik_em(llik_emulator)` is TRUE.
  #    par_prior: data.frame defining the prior distribution for `par`. 
  #    par_init: numeric, d-dimensional vector, the initial condition. If NULL,
  #              falls back on the default BayesianTools behavior for setting 
  #              the initial condition, which typically implies sampling from 
  #              the prior.
  #    n_itr: integer, the number of MCMC iterations.
  #    sampler: the MCMC algorithm used by `BayesianTools::runMCMC()`. See the
  #             BayesianTools documentation for options. The default "DEzs"
  #             aligns with the BayesianTools default, which refers to a 
  #             specific implementation of a differential evolution MCMC 
  #             algorithm.
  #   settings_list: There are many settings that can be altered in the MCMC 
  #                  algorithms implemented by `runMCMC()`, which are passed
  #                  via the `runMCMC()` "settings" argument. The list provided
  #                  in this argument is passed to this "settings" argument, 
  #                  with a couple caveats: the "n_itr" argument above is 
  #                  used to set `settings_list$iterations` which controls the 
  #                  number of MCMC iterations. Also, `settings_list$nrChains`
  #                  is always set to 1, since this function is designed to run 
  #                  a single MCMC chain. settings_list$startValue` is
  #                  set to `par_init`, which is sampled from the prior if 
  #                  `par_init` is NULL.
  #   ...: additional arguments passed to `llik_emulator$assemble_llik()`.
  #
  # Returns: 
  # Returns a list with elements "samp" and "bt_output". "samp" follows the 
  # formatting described in e.g. `mcmc_noisy_llik()`. "bt_output" is the object 
  # returned by the `runMCMC()` function. Note that some MCMC algorithms 
  # (e.g., the differential evolution algorithms) run multiple chains in order 
  # to design better proposal distributions. For these algorithms, the returned 
  # sample matrices incorporate samples combined from all of these chains. The 
  # original output from each chain separately can be found in
  # "bt_output$chain". 
  
  assert_that(is_llik_em(llik_emulator))
  assert_that(llik_emulator$exact_llik)
  
  # Log-likelihood function. Set up to accept numeric vector `par` as an 
  # argument.
  llik <- function(par) {
    par <- matrix(par, nrow=1)
    llik_emulator$assemble_llik(par, ...)
  }
  
  # Log prior density and sampling function.
  lprior <- get_lprior_dens(par_prior)
  prior_sampler <- get_prior_sampler(par_prior)
  
  # Create BayesianSetup object.
  bayesianSetup <- BayesianTools::createBayesianSetup(likelihood=llik, 
                                                      prior=lprior, 
                                                      parallel=FALSE,
                                                      priorSampler=prior_sampler, 
                                                      names=llik_emulator$input_names)
  
  # MCMC settings list. Ensure number of MCMC iterations and number of MCMC 
  # chains are in the list; other BayesianTools MCMC settings can optionally be 
  # included.
  settings_list$iterations <- n_itr
  settings_list$nrChains <- 1L
  
  # Set initial condition.
  if(!is.null(par_init) && is.numeric(par_init)) {
    settings_list$startValue <- matrix(par_init, nrow=1L)
  }
  
  # Run MCMC.
  bt_output <- BayesianTools::runMCMC(bayesianSetup=bayesianSetup, 
                                      sampler=sampler, settings=settings_list)
  
  # Convert output into a form that aligns with other MCMC functions. Note that 
  # for certain algorithms (e.g., "DEzs") that return multiple trajectories, 
  # this function will combine the trajectories in to a single chain.
  bt_samp <- BayesianTools::getSample(bt_output, parametersOnly=TRUE)

  # Return list with samples stored in matrices as well as the MCMC object 
  # directly returned by BayesianTools.
  return(list(samp=bt_samp, bt_output=bt_output))
}


# ------------------------------------------------------------------------------
# Proposal Covariance Adaptation 
# ------------------------------------------------------------------------------

adapt_MH_proposal_cov <- function(cov_prop, log_scale_prop, times_adapted, adapt_cov, adapt_scale, 
                                  samp_interval, accept_rate, accept_rate_target=0.24,  
                                  adapt_factor_exponent=0.8, adapt_factor_numerator=10) {
  # Follows the adaptive Metropolis scheme used in the Nimble probabilistic 
  # programming language. 
  # Described here: https://arob5.github.io/blog/2024/06/10/adapt-mh/

  return_list <- list()
  adapt_factor <- 1 / (times_adapted + 3)^adapt_factor_exponent
  
  if(adapt_cov) {
    if(accept_rate > 0) {
      sample_cov_interval <- cov(samp_interval)
      cov_prop <- cov_prop + adapt_factor * (sample_cov_interval - cov_prop)
    } 
    
    # Handle case that `cov_prop` may not be positive definite. 
    L_cov_prop <- try(t(chol(cov_prop)))
    if(inherits(L_cov_prop, "try-error")) {
      message("Problematic proposal cov: ", cov_prop)
      cov_prop <- Matrix::nearPD(cov_prop, ensureSymmetry=TRUE)
      L_cov_prop <- t(chol(cov_prop))
    }
    
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

# ------------------------------------------------------------------------------
# MCMC helper/wrapper functions.
#  For running chains in parallel and setting initial conditions.
# ------------------------------------------------------------------------------

run_mcmc_multichain <- function(mcmc_func_name, llik_emulator, n_chain=4L, 
                                par_init=NULL, par_prior=NULL, lik_par_init=NULL, 
                                lik_par_prior=NULL, ic_sample_method="LHS", 
                                try_parallel=TRUE, n_cores=NULL, package_list=NULL, 
                                dll_list=NULL, obj_list=NULL, 
                                test_label="default_lbl", ...) {
  # This wrapper function runs multiple MCMC chains in parallel. If initial 
  # conditions are not supplied by the user, then they will be sampled from 
  # their prior distribution. This function collects the output from all chains
  # into a single data.table satisfying the `samp_dt` requirements outlined 
  # in `mcmc_helper_functions.r`. Let d denote the dimension of the parameter 
  # space.
  #
  # Args:
  #    mcmc_func_name: character, the name of the MCMC function to call.
  #    llik_emulator: llikEmulator object encoding the exact or approximate 
  #                   log-likelihood.
  #    n_chain: integer, number of independent MCMC chains to run.
  #    par_init: matrix of dimension (n_chain,d), each row containing the 
  #              initial condition for one of the MCMC chains.
  #    par_prior: data.frame storing the prior distribution for the parameters.
  #    ic_sample_method: the sampling method used to sample the initial 
  #                      conditions for the chains. Default is Latin Hypercube 
  #                      sampling.
  #   ...: should contain all arguments to be passed to the MCMC function, other 
  #        than the initial conditions and priors.
  #
  # TODO: Need to generalize to handle `lik_par_init`/`lik_par_prior`.
  
  # Sample initial conditions from prior if not provided.
  if(is.null(par_init)) {
    par_init <- get_batch_design(ic_sample_method, N_batch=n_chain, 
                                 prior_params=par_prior)
  } else {
    assert_that(length(par_init) == n_chain)
  }
  par_init_list <- lapply(1:n_chain, function(i) par_init[i,])
  mcmc_func <- function(ic) get(mcmc_func_name)(llik_emulator=llik_emulator,
                                                par_prior=par_prior,
                                                par_init=ic, ...)
  
  if(try_parallel) {
    # Prepare for parallel run. 
    if(is.null(n_cores)) {
      n_cores <- parallel::detectCores()
    } else {
      n_cores <- min(n_cores, parallel::detectCores())
    }
    n_cores <- min(n_chain, n_cores)
    cl <- parallel::makeCluster(n_cores)
    
    # Ensure required packages, objects, and DLLs will be available during the 
    # parallel execution.
    export_list <- get_parallel_exports(package_list, dll_list, obj_list)
    parallel::clusterCall(cl, export_list$load_func, export_list$packages, 
                          export_list$dlls)
    parallel::clusterExport(cl, varlist=export_list$objects)
    mcmc_chain_list <- parLapply(cl, X=par_init_list, fun=mcmc_func)
    parallel::stopCluster(cl)
  } else {
    mcmc_chain_list <- lapply(par_init_list, mcmc_func)
  }

  # TODO: wrap above in trycatch and put stopCluster in finally block.
  
  # Convert to `samp_dt` data.table format
  samp_dt <- format_mcmc_output_multi_chain(lapply(mcmc_chain_list, function(l) l$samp), 
                                            test_label=test_label)
  
  # Non-sample output stored in its own list. 
  mcmc_outputs_list <- vector(mode="list", length=length(mcmc_chain_list)) 
  for(i in seq_along(mcmc_outputs_list)) {
    non_samp_elements <- setdiff(names(mcmc_chain_list[[i]]), "samp")
    mcmc_outputs_list[[i]] <- mcmc_chain_list[[i]][non_samp_elements]
  }
  
  return(list(samp=samp_dt, output_list=mcmc_outputs_list))
}


get_parallel_exports <- function(package_list=NULL, dll_list=NULL, obj_list=NULL) {
  # In using `parallel` apply functions, objects and packages must explicitly 
  # be made available in the cluster environment. This function returns lists 
  # containing R objects, packages, and dynamic-linked libraries (DLLs) that 
  # can then be passed to `parallel::clusterCall()` and 
  # `parallel::clusterExport()` to make all of these available in the cluster 
  # environment. This function is closely modeled after the function 
  # `generateParallelExecutor()` in the R `BayesianTools` package. By default, 
  # this function will return all loaded packages, DLLs, and objects in the 
  # global environment. This can be wasteful when the function being executed 
  # in parallel only relies on a small subset of these. The arguments of this 
  # function can thus explicitly pass the objects to include when the 
  # dependencies of the function are explicitly known. 
  #

  # Select all loaded packages if not explicitly passed.
  if(is.null(package_list)) package_list <- .packages()
  
  # Select all DLLs if not explicitly passed. Returns object of class 
  # `DLLInfoList`, which can be iterated and DLLs extracted as strings.
  if(is.null(dll_list)) {
    dll_list <- getLoadedDLLs()
    dll_list <- sapply(dll_list, function(x) x[[2]])
  }
  
  # Select all objects in global environment if not explicitly passed.
  if(is.null(obj_list)) obj_list <- ls(envir=.GlobalEnv)

  # Helper function to load packages and DLLs in the parallel cluster.
  load_package_dll <- function(package_list=NULL, dll_list=NULL) {
    # Load packages. 
    if(!is.null(package_list)) {
      for(package_name in package_list) library(package_name, character.only=TRUE)
    }
    
    # Load DLLs.
    if(!is.null(dll_list)) {
      for(dll_name in dll_list) try(dyn.load(dll_name), silent=TRUE)
    }
  }
  
  return_list <- list(packages=package_list, dlls=dll_list, objects=obj_list,
                      load_func=load_package_dll)
  return(return_list)
}


# ------------------------------------------------------------------------------
# Metropolis-Hastings acceptance probability approximations using 
# log-likelihood emulators. 
# ------------------------------------------------------------------------------

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


# -----------------------------------------------------------------------------
# Sample-Based Posterior Approximation: 
# The algorithms here attempt to approximately sample the so-called "sample-
# based" approximate posterior. Samples from this distribution are generated 
# by first sampling a random distribution (via sampling a GP trajectory) and 
# then sampling a parameter value from this distribution. In general, this is 
# impractical due to the requirement of sampling a random probability measure. 
# The functions here attempt to approximate this sampling procedure. 
# -----------------------------------------------------------------------------


sim_sample_based_approx_1d_grid <- function(llik_emulator, par_prior_params, bounds, 
                                            N_grid, N_samp, lik_par_val=NULL, 
                                            max_rejects=1000L, ...) {
  # This function (approximately) generates samples via: 
  #    (1) Sample a random log-likelihood L. 
  #    (2) Sample a value u from the posterior implied by L. 
  # To make this practical, L is sampled at a finite number of locations (an equally
  # spaced grid) and then linearly interpolated to produce an approximation to L^hat(u).
  # Then u is sampled from the distribution with unnormalized density 
  # prior(u)*exp{L^hat(u)} via rejection sampling. 
  
  assert_that(llik_emulator$dim_input==1L)
  
  # Construct grid. 
  par_grid <- seq(bounds[1], bounds[2], length.out=N_grid)
  par_grid_mat <- matrix(par_grid,ncol=1)
  colnames(par_grid_mat) <- llik_emulator$input_names
  
  # Iterate until desired number of samples is achieved. 
  par_samp <- vector(mode="numeric", length=N_samp)

  for(i in 1:N_samp) {
    
    # Sample log-likelihood values over dense grid. 
    llik_grid <- drop(llik_emulator$sample(par_grid_mat, lik_par_val=lik_par_val, N_samp=1, ...))
    lpost_grid <- llik_grid + calc_lprior_theta(par_grid_mat, par_prior_params)
    
    # Linearly interpolate to obtain approximate log-likelihood sample/trajectory. 
    # llik_interp <- get_linear_interpolation_1d(par_grid, llik_grid)
    # lpost_interp <- function(par) llik_interp(par) + calc_lprior_theta(par, par_prior_params)
    lik_interp <- get_linear_interpolation_1d(par_grid, exp(llik_grid))
    post_interp <- function(par) lik_interp(par) * exp(calc_lprior_theta(par, par_prior_params))
     
    # Upper bound for (unnormalized) density. 
    M <- max(exp(lpost_grid))
    assert_that(M > 0)
    assert_that(!is.infinite(M))
    
    # Rejection sampling until a sample is obtained. 
    N_rejects <- 0L 
    accepted <- FALSE
    while(!accepted && (N_rejects < max_rejects)) {
      par_prop <- runif(n=1, min=bounds[1], max=bounds[2])
      u <- runif(n=1, min=0, max=M)
      if(u <= post_interp(par_prop)) {
      # if(log(u) <= lpost_interp(par_prop)) {
        accepted <- TRUE
        par_samp[i] <- par_prop
      } else {
        N_rejects <- N_rejects+1L
        if(N_rejects==max_rejects) stop(paste0("Reached max rejects when attempting to simulate sample ", i))
      }
    }
  }
  
  return(par_samp)
}

estimate_sample_based_density_1d_grid <- function(llik_emulator, par_prior_params, input_grid, 
                                                  N_monte_carlo=1000L, lik_par_val=NULL, ...) {
  # This function is similar to `sim_sample_based_approx_1d_grid()` but its goal is not 
  # to draw parameter samples `u`, but instead to estimate the expectation of the random 
  # density of `u` implied by the algorithm 
  #    (1) Sample a random log-likelihood L. 
  #    (2) Sample a value u from the posterior implied by L.
  # This function estimates the density of the samples `u` returned by this algorithm. 
  # It does so by sampling density values on a grid, numerically normalizing each,
  # and then computing the average over the normalized values for each grid point. 
  # Note that this (approximately) computes the expected normalized density, whereas
  # the marginal approximation computes the expected unnormalized density and then 
  # normalizes after-the-fact. This approach has the benefit over 
  # `sim_sample_based_approx_1d_grid()` of not requiring rejection sampling, but
  # requires more care regarding numerical stability. I am completely ignoring this 
  # for now and working on the original (non-log) scale. If this function is to be 
  # used in cases where log-likelihood values might be very small, then this code 
  # should be updated. 
  
  assert_that(llik_emulator$dim_input==1L)
  
  # Simulate log-likelihood values. Return shape is (N_grid, N_monte_carlo).
  llik_samp <- llik_emulator$sample(input_grid, lik_par_val=lik_par_val, N_samp=N_monte_carlo, ...)

  # Assuming numbers are reasonable for now so exponentiating is fine. 
  lprior_dens <- calc_lprior_theta(input_grid, par_prior_params)
  post_dens_samp <- exp(add_vec_to_mat_cols(lprior_dens, llik_samp))
  
  # Normalize.
  dx <- abs(drop(input_grid)[2] - drop(input_grid)[1])
  for(j in 1:ncol(post_dens_samp)) {
    post_dens_samp[,j] <- post_dens_samp[,j] / int_trap(post_dens_samp[,j], dx)
  }
  
  # Average over Monte Carlo samples for each input.
  density_estimates <- rowMeans(post_dens_samp)
  
  return(density_estimates)
}


get_linear_interpolation_1d <- function(input_points, output_points) {
  
  point_order <- order(input_points)
  input_points <- input_points[point_order]
  output_points <- output_points[point_order]
  
  f <- function(x) {
    idx_l <- max(which(input_points <= x))
    idx_u <- idx_l+1L
    x_l <- input_points[idx_l]
    x_u <- input_points[idx_u]
    y_l <- output_points[idx_l]
    y_u <- output_points[idx_u]
    
    y_l + (y_u-y_l)/(x_u-x_l)*(x-x_l)
  }
  
  return(f)
  
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











