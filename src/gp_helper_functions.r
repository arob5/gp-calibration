#
# gp_helper_functions.r
# Gaussian process (GP) helper functions, leveraged by the gpWrapper class and
# sub-classes. 
#
# Andrew Roberts
#


get_lengthscale_bounds <- function(X, include_one_half=FALSE, p=NULL) {
  # Returns empirically-determined bounds (and quantiles - see below)
  # for the lengthscale parameters of a Gaussian (exponentiated quadratic)
  # covariance function, parameterized for a d-dimensional input 
  # space as 
  #    k(x,y) = prod_{j=1}^{d} exp{-(x_j - y_j)^2 / ell_j^2}
  # if `include_one_half` is FALSE and 
  #    k(x,y) = prod_{j=1}^{d} exp{-0.5*(x_j - y_j)^2 / ell_j^2}
  # otherwise. The lengthscale parameters here refer to 
  # ell_1,...,ell_2 (not their squares). The matrix `X` should 
  # be of shape (n,d) and is interpretated as a design matrix/
  # training inputs. The bounds for ell_j are defined by 
  # computing all pairwise Euclidean distances between the scalars 
  # in the vector `X[,j]` (the jth dimension of the design points)
  # and setting the min/max bounds to the min/max of these 
  # pairwise distances. This is done independently for each 
  # dimension. Optionally, a vector of probabilities `p` can be 
  # provided, in which case the p-quantile of the pairwise 
  # distances will also be computed for each dimension, and stored
  # in additional rows of the matrix. 
  #
  # Returns:
  # matrix of dimension (2+length(p),d) [where d=ncol(X)] with the first 
  # row containing the lower bounds for ell_1,...,ell_d, 
  # respectively, and the second row containing the upper bounds. 
  # The rownames for these two rows is set to ("min","max")
  # The colnames of the matrix is set to `colnames(X)`. If `p` is non-NULL
  # then one row will be added corresponding to the respective quantiles 
  # specified by `p`. The rownames for these rows will be set to 
  # `q<100*pi>` where `pi` is the respective element of `p`; e.g., 
  # "q10" indicates the 10th percentile (pi = 0.1). 
  
  assert_that(is.matrix(X))
  if(!is.null(p)) {
    assert_that(is.vector(p))
    assert_that(is.numeric(p))
    assert_that(all(p <= 1.0))
    assert_that(all(p >= 0.0))
    quantile_names <- paste0("q", 100*p)
  } else {
    quantile_names <- NULL
  }
  
  d <- ncol(X)
  n_quantiles <- ifelse(is.null(p), 0L, length(p))
  bounds <- matrix(nrow=2+n_quantiles, ncol=d)
  
  for(j in 1:d) {
    # Pairwise Euclidean distance. 
    dists_j <- as.matrix(dist(X[,j], method="euclidean", diag=FALSE))
    dists_j <- dists_j[upper.tri(dists_j)]
    dists_j <- dists_j[dists_j > 0]
    
    # Bounds computed as min/max of pairwise distances. 
    bounds_j <- range(dists_j)
    
    # Compute quantiles of pairwise distances.
    if(!is.null(p)) {
      quantiles_j <- quantile(dists_j, p)
    } else {
      quantiles_j <- NULL
    }
    
    # Adjust if 1/2 factor is included in covariance parameterization.
    if(include_one_half) {
      bounds_j <- bounds_j / sqrt(2)
      quantiles_j <- quantiles_j / sqrt(2)
    }
    
    bounds[,j] <- c(bounds_j, quantiles_j)
  }
  
  # Row and column names. 
  rownames(bounds) <- c("min", "max", quantile_names)
  colnames(bounds) <- colnames(X)
  
  return(bounds)
}


get_marginal_variance_bounds <- function(M, p=0.9, return_variance=FALSE) {
  # A convenience function for determining bounds/default values for 
  # the marginal variance (i.e., scale) parameter of 
  # a Gaussian process (GP). For a correlation function C, the 
  # assumed parameterization is s^2 * C is `return_variance` is 
  # FALSE (the default), and otherwise s * C. This function 
  # assumes a mean-zero process. In particular, this function 
  # returns the scale `s` such that 
  #    P[-M <= N(0,s^2) <= M] = p. 
  # Multiple values of `M` can be passed (i.e., `M` can be a numeric
  # vector), in which case this process is repeated for each M and 
  # the vector of corresponding `s` values are returned. `s` is 
  # solved for by noting that the above probability can be rearranged
  # to obtain 
  #    s = M / Q(0.5*[1+p])
  # where Q is the standard normal inverse CDF/quantile function. 
  
  # Argument checking. 
  assert_that(is.vector(p))
  assert_that(length(p) == 1)
  assert_that((p>0) && (p<1))
  assert_that(is.vector(M))
  assert_that(all(M > 0))
  
  # Computes the standard deviation, not the variance. 
  calc_s <- function(m) {
    m / qnorm(0.5*(1+p))
  }
  
  s_vals <- sapply(M, calc_s)
  if(return_variance) return(s_vals^2)
  return(s_vals)
}








