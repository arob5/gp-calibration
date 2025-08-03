# post_approx_grid.r
# 
# Various functions and helper functions for computing surrogate-based
# posterior approximations in one and two dimensions. These functions
# mostly operate by approximating normalizing constants over a dense
# grid and are thus not feasible in higher dimensions. The emphasis of
# these functions is on computing the normalized posterior density
# approximations.
#
# Andrew Roberts
#

normalize_density_grid <- function(ldens_vals, grid, log_scale=FALSE) {
  # A wrapper around `normalize_density_1d` and `normalize_density_2d`
  # that calls the proper function depending on the dimension of 
  # `grid`.
  
  # 1d input space.
  if(is.null(dim(grid)) || (ncol(grid) == 1L)) {
    normalize_density_1d(ldens_vals, grid, log_scale)
  } else if(ncol(grid) == 2L) {
    normalize_density_2d(ldens_vals, grid, log_scale)
  } else {
    stop("`normalize_density_grid` requires 1d or 2d `grid`.")
  }
  
}

normalize_density_1d <- function(ldens_vals, grid, log_scale=FALSE) {
  # Valid only for 1d input space.
  # Given unnormalized log-density values, returns the normalized
  # density, computed using LogSumExp. If `ldens_vals` is a matrix,
  # then each row is normalized. The normalizing constant is approximated
  # by the trapezoidal rule, and assumes that `grid` is ordered in ascending
  # order (it need not be equally spaced).
  
  if(is.null(dim(ldens_vals))) ldens_vals <- matrix(ldens_vals, nrow=1L) 
  n_grid <- ncol(ldens_vals)
  
  if(n_grid != length(drop(grid))) {
    stop("`ldens_vals` and `grid` dimension mismatch.")
  }
  
  l_du <- log(diff(drop(grid)))
  log_summands_1 <- add_vec_to_mat_rows(l_du, ldens_vals[, 1:(n_grid-1), drop=FALSE])
  log_summands_2 <- add_vec_to_mat_rows(l_du, ldens_vals[, 2:n_grid, drop=FALSE])
  
  log_term_1 <- matrixStats::rowLogSumExps(log_summands_1) - log(2)
  log_term_2 <- matrixStats::rowLogSumExps(log_summands_2) - log(2)
  log_norm_csts <- matrixStats::rowLogSumExps(cbind(log_term_1, log_term_2))
  ldens_norm <- add_vec_to_mat_cols(-log_norm_csts, ldens_vals)
  
  if(log_scale) return(ldens_norm)
  return(exp(ldens_norm))
}

normalize_density_2d <- function(ldens_vals, grid, log_scale=FALSE) {
  # Valid only for 12 input space.
  # Given unnormalized log-density values, returns the normalized
  # density, computed using LogSumExp. If `ldens_vals` is a matrix,
  # then each row is normalized. The normalizing constant is approximated
  # by Reimann sums, and assumes that `grid` is equally spaced (unlike
  # normalize_density_1d).
  
  if(!isTRUE(ncol(grid) == 2L)) stop("`grid` must have 2 columns.")
  if(is.null(dim(ldens_vals))) ldens_vals <- matrix(ldens_vals, nrow=1L) 
  
  # Determine grid spacing.
  x_unique <- sort(unique(grid[,1]))
  y_unique <- sort(unique(grid[,2]))
  dx <- unique(diff(x_unique))
  dy <- unique(diff(y_unique))
  
  # Allow for small differences in spacing consistent with numerical error.
  tol <- 1e-8
  if(abs(max(dx) - min(dx)) > tol) {
    stop("`grid` must be equally spaced in x-direction.")
  }
  if(abs(max(dy) - min(dy)) > tol) {
    stop("`grid` must be equally spaced in y-direction.")
  }
  
  dx <- dx[1]
  dy <- dy[1]
  
  # Riemann sum.
  dA <- dx * dy
  log_norm_csts <- matrixStats::rowLogSumExps(ldens_vals) + log(dA)
  ldens_norm <- add_vec_to_mat_cols(-log_norm_csts, ldens_vals)

  if(log_scale) return(ldens_norm)
  return(exp(ldens_norm))
}

get_EP_dens_grid <- function(lpost_em, input_grid, n_mc, log_scale=FALSE, ...) {
  # Approximates the density of the expected posterior approximation via
  # a grid-based approximation of the normalizing constant. 
  
  if(lpost_em$dim_input > 2L) stop("`lpost_em$dim_input` must be 1 or 2.")
  
  if(is.null(dim(input_grid))) input_grid <- matrix(input_grid, ncol=1L)
  n_grid <- nrow(input_grid)
  
  # Simulate log-likelihood values. Return shape is (n_grid, n_mc) before
  # transposing.
  lpost_samp <- t(lpost_em$sample(input_grid, N_samp=n_mc, ...))
  
  # Normalize each trajectory.
  lpost_samp_norm <- normalize_density_grid(lpost_samp, input_grid, log_scale=TRUE)
  
  # Average density over trajectories.
  lpost_ep_norm <- matrixStats::colLogSumExps(lpost_samp_norm) - log(n_mc)
  
  if(log_scale) return(lpost_ep_norm)
  return(exp(lpost_ep_norm))
}


# Exact posterior density.
get_post_exact <- function(grid_info, log_scale=FALSE) {
  drop(normalize_density_grid(grid_info$lpost, grid_info$input, log_scale=log_scale))
}

# Plug-In mean approximation.
get_post_mean <- function(grid_info, lpost_em, log_scale=FALSE, ...) {
  pred <- lpost_em$predict(grid_info$input, return_var=FALSE, ...)
  
  # Normalize density.
  drop(normalize_density_grid(pred$mean, grid_info$input, log_scale=log_scale))
}

# Expected likelihood approximation.
get_post_EL <- function(grid_info, lpost_em, log_scale=FALSE, ...) {
  pred <- lpost_em$predict_lik(grid_info$input, return_var=FALSE, log_scale=TRUE, ...)
  
  # Normalize density.
  drop(normalize_density_grid(pred$log_mean, grid_info$input, log_scale=log_scale))
}

get_post_noisy <- function(grid_info, lpost_em, par_prior, mode="mcwmh", use_joint=TRUE,
                           n_avg=1L, n_chains=4L, n_itr=10000L, itr_start=7000L, 
                           adjustment="rectified", lbl=NULL, cov_prop=NULL, 
                           ics=NULL, ...) {
  
  if(is.null(lbl)) lbl <- paste(mode, ifelse(use_joint, "joint", "ind"), sep="_")
  if(is.null(ics)) ics <- get_batch_design("simple", N_batch=n_chains, prior_params=par_prior)
  
  mcmc_settings <- list(mcmc_func_name="mcmc_noisy_llik", par_prior=par_prior, 
                        par_init=ics, mode=mode, use_joint=use_joint, n_avg=n_avg,
                        n_itr=n_itr, n_chain=n_chains, itr_start=itr_start, 
                        try_parallel=FALSE, cov_prop=cov_prop, 
                        log_scale_prop=0, adapt_cov_prop=FALSE, 
                        adapt_scale_prop=FALSE, llik_em=lpost_em, test_label=lbl)
  results <- do.call(run_mcmc_chains, mcmc_settings)
  
  max_rhat <- calc_R_hat(results$samp)$R_hat_vals$R_hat
  if(max_rhat > 1.05) message("Warning: R-hat above threshold.")
  
  return(results)
}




