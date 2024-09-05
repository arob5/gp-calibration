#
# general_helper_functions.r
#
# Andrew Roberts
#

add_vec_to_mat_rows <- function(v, M) {
  # `v`: numeric vector of length `n`. 
  # `M`: matrix with `n` columns. 
  
  assert_that(length(v) == ncol(M))
  M + rep(v, each=nrow(M))
}

add_vec_to_mat_cols <- function(v, M) {
  # `v`: numeric vector of length `n`. 
  # `M`: matrix with `n` rows 
  
  assert_that(length(v) == nrow(M))
  v + M
}

mult_vec_with_mat_rows <- function(v, M) {
  # `v`: numeric vector of length `n`. 
  # `M`: matrix with `n` rows 
  
  assert_that(length(v) == ncol(M))
  M * rep(v, each=nrow(M))
  
}

mult_vec_with_mat_cols <- function(v, M) {
  # `v`: numeric vector of length `n`. 
  # `M`: matrix with `n` rows 
  
  assert_that(length(v) == nrow(M))
  v * M
}

log_exp_minus_1 <- function(x, threshold=100) {
  # A numerically stable implementation of log(exp(x)-1). This computation 
  # shows up, for example, in computing the variance of a log-normal 
  # random variable. To prevent overflow, the computation uses the 
  # approximation log(exp(x)-1) ~ x for x greater than `threshold`. 
  # `x` may any object containing numeric values which can be 
  # indexed via `[`; e.g. a numeric vector or matrix. In these cases, 
  # the `log` and `exp` are interpreted as being applied elementwise. 
  #
  # Args:
  #    x: See above description. 
  #    threshold: numeric(1), the threshold used to define the cutoff
  #               point where the approximation is applied. 
  #
  # Returns: 
  #    The approximated value of log(exp(x)-1), which will be the same 
  #    class and shape as `x`. 
  
  idx_small_vals <- (x < threshold)
  x[idx_small_vals] <- log(exp(x[idx_small_vals]) - 1) 
  return(x)
}


log_diff_exp <- function(x, y, threshold=100) {
  # A numerically stable implementation of log{e^x - e^y}. The expression
  # is computed using:
  # log{e^x - e^y} = log{e^y * [e^{x-y} - 1]} = b + log{e^{x-y} - 1} = b + log_exp_minus_1(x-y).
  # When `x`, `y` are numeric vectors of length greater than 1, or matrices, then the computation 
  # is vectorized (computed elementwise). 
  #
  # Args:
  #    x,y: numeric vectors or matrices, of equal length/dimensions. 
  #    threshold: numeric(1), passed to `log_exp_minus_1`. 
  #
  # Returns: 
  #    Approximated value of log{e^x - e^y}, which will be the same class and 
  #    shape as `x` and `y`. 
  
  assert_that(class(x) == class(y))
  if(is.numeric(x) && is.null(dim(x))) assert_that(length(x) == length(y))
  else if(is.matrix(x)) assert_that(all(dim(x) == dim(y)))
  else stop("Unsupported class for `x` or `y`.")
  
  x + log_exp_minus_1(x-y, threshold=threshold)
}


# Basic trapezoidal rule for numerical integration. Mostly used in this code for 
# approximating normalizing constants for basic 1d integrals. 
int_trap <- function(x, dx) {
  x <- drop(x)
  n_pts <- length(x)
  
  return(0.5 * dx * (x[1] + x[n_pts]) +  dx * sum(x[2:(n_pts-1)]))
}



