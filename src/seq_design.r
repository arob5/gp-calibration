#
# seq_design.r
#
# General functions related to sequential design and optimization. While 
# `seq_design_gp.r` concerns sequential design for Gaussian processes (GPs)
# within the acquisition function minimization framework, and 
# `seq_design_for_post_approx.r` considers the same framework but specialized 
# to the used of log-likelihood emulators for posterior approximation, this 
# file contains functions related to experimental design which may not 
# fall into the acquisition function framework and may not involve GPs at all. 
# In particular, this includes "model free" space-filling methods such as 
# Latin Hypercube sampling (LHS) which are often used to generate an initial 
# design that will subsequently be refined by sequential methods. Also 
# included are functions tailored to solving Bayesian inverse problems 
# such as Ensemble Kalman Inversion (EKI). This also includes a class 
# of methods which construct a design by sub-sampling a set of candidate 
# points (e.g. support points and Stein points). 
# 
# Andrew Roberts
# 

get_batch_design <- function(method, N_batch, par_prior_params=NULL, design_candidates=NULL, 
                             design_candidate_weights=NULL, ...) {
  
  if(method == "LHS") return(get_LHS_sample(N_batch, par_prior_params, ...))
  else stop("Design method ", method, " not supported.")
  
}


get_LHS_sample <- function(N_batch, prior_params=NULL, bounds=NULL, order_1d=FALSE, ...) {
  # The dimension of the input space.
  d <- nrow(theta_prior_params)
  
  if(is.null(param_ranges)) {
    param_ranges <- matrix(NA, nrow = 2, ncol = d)
  }
  
  # Generate LHS design. 
  X_LHS <- randomLHS(N_points, d)
  
  # Apply inverse CDS transform using prior distributions.
  for(j in seq_len(d)) {
    if(theta_prior_params[j, "dist"] == "Uniform") {
      lower_bound <- ifelse(is.na(param_ranges[1, j]), theta_prior_params[j, "param1"], param_ranges[1, j])
      upper_bound <- ifelse(is.na(param_ranges[2, j]), theta_prior_params[j, "param2"], param_ranges[2, j])
      X_LHS[,j] <- qunif(X_LHS[,j], lower_bound, upper_bound) 
    } else if(theta_prior_params[j, "dist"] == "Gaussian") {
      lower_bound <- ifelse(is.na(param_ranges[1, j]), -Inf, param_ranges[1, j])
      upper_bound <- ifelse(is.na(param_ranges[2, j]), Inf, param_ranges[2, j])
      X_LHS[,j] <- qtruncnorm(X_LHS[,j], a = lower_bound, b = upper_bound, mean = theta_prior_params[j, "param1"], sd = theta_prior_params[j, "param2"])
    } else {
      stop("Unsupported prior distribution: ", theta_prior_params[j, "dist"])
    }
  }
  
  # For 1 dimensional data, order samples in increasing order. 
  if(order_1d && (d == 1)) {
    X_LHS <- X_LHS[order(X_LHS),,drop=FALSE]
  }
  
  colnames(X_LHS) <- rownames(theta_prior_params)
  
  return(X_LHS)
 
}















