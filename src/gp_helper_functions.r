#
# gp_helper_functions.r
# Gaussian process (GP) helper functions, leveraged by the gpWrapper class and
# sub-classes. 
#
# Andrew Roberts
#

get_lengthscale_bounds <- function(X, p_min=0.05, p_max=0.95, cor_min=0.01, 
                                   cor_max=0.5, p_extra=NULL, dim_by_dim=FALSE,
                                   include_one_half=FALSE, 
                                   convert_to_square=FALSE) {
  # Returns empirically-determined bounds (and quantiles - see below)
  # for the lengthscale parameters of a Gaussian (exponentiated quadratic)
  # covariance function, parameterized for a d-dimensional input 
  # space as 
  #    k(x,y) = prod_{j=1}^{d} exp{-(x_j - y_j)^2 / ell_j^2}
  # if `include_one_half` is FALSE and 
  #    k(x,y) = prod_{j=1}^{d} exp{-0.5*(x_j - y_j)^2 / ell_j^2}
  # otherwise. The lengthscale parameters here refer to 
  # ell_1,...,ell_2 (not their squares) by default. If `convert_to_square`,
  # then the returned bounds are defined with respect to the squares
  # of these lengthscales: ell_1^2,...,ell_2^2. The matrix `X` should 
  # be of shape (n,d) and is interpreted as a design matrix/
  # training inputs. The bounds for ell_j are defined by 
  # computing all pairwise Euclidean distances between 
  #  (1) the scalers in the vector `X[,j]` (the jth dimension of the design points),
  #      if `dim_by_dim = TRUE`; or 
  #  (2) the rows of `X` (i.e., using d-dimensional Euclidean distance)
  #      if `dim_by_dim = FALSE`. 
  # d_min[,j]/d_max[,j] are then defined as `p_min` and `p_max`
  # empirical quantiles of the pairwise distances in the jth
  # dimension. The minimum bound for ell_j is then defined as the 
  # value that sets the Gaussian correlation to `cor_min` when 
  # the inputs are distance `d_min` apart. Similarly, the maximum 
  # bound for ell_j is defined as the value that sets the Gaussian
  # correlation to `cor_max` when the inputs are distance `d_max`
  # apart. In addition to returning these bounds, this function 
  # can also return quantiles of the pairwise distance distribution 
  # (which can be useful for setting "starting" or "default" values
  # for the lengthscales). Note that the argument defaults are  
  # typically reasonable when `dim_by_dim = FALSE`; if setting  
  # `dim_by_dim = TRUE`, the defaults should probably be overwritten; for 
  # example, by raising the default `cor_min` and `cor_max` to the 1/d
  # power, where `d = ncol(X)`. 
  #
  # References:
  #  This function was inspired by the hetGP package function 
  # `auto_bounds()`. The behavior of this function agrees with that of 
  # `auto_bounds()` when `dim_by_dim = FALSE` and `convert_to_square = TRUE`. 
  #
  # Returns:
  # list, with elements "ell_bounds" and "dist_quantiles". The former is a 
  # matrix of dimension (2,d) [where d=ncol(X)] with rownames "lower" and 
  # "upper". The first contains the lower bounds for ell_1,...,ell_d, and 
  # the second row contains the upper bounds. "dist_quantiles" is 
  # the matrix returned by `get_pairwise_dist_quantiles_by_dim()` or 
  # `get_pairwise_dist_quantiles()`, but potentially adjusted to account 
  # for the "one half" scaling and converted to squared distances if 
  # `convert_to_square = TRUE`. This matrix contains the distance quantiles 
  # used in constructing the lengthscale bounds, as well as the min/max 
  # pairwise distances, and any other quantiles specified in `p_extra`.
  # For both matrices, column names are set to `colnames(X)`.

  # Compute min/max and quantiles of pairwise distances in each dimension.
  # Note that the first two rows of `dist_q` will be the min/max and the 
  # remaining rows will correspond to the quantiles. 
  if(dim_by_dim) {
    dist_q <- get_pairwise_dist_quantiles_by_dim(X, p=c(p_min, p_max, p_extra))
  } else {
    # First scale to [0,1]^d to reduce to isotropic setting.
    d <- ncol(X)
    X_bounds <- get_bounds(X)
    target_bounds <- rbind(rep(0,d), rep(1,d))
    X <- scale_inputs(X, target_bounds=target_bounds, source_bounds=X_bounds)
    dist_q <- get_pairwise_dist_quantiles(X, p=c(p_min, p_max, p_extra))
  }
  # Adjust if 1/2 factor is included in covariance parameterization.
  if(include_one_half) dist_q <- dist_q / sqrt(2)

  # Adjust min/max bounds to account for the correlation constraints.
  # The 3rd and 4th rows correspond to the distance quantiles to use in 
  # constructing the bounds.
  ell_bounds <- dist_q[c(3,4),,drop=FALSE]
  ell_bounds[1L,] <- ell_bounds[1L,] / sqrt(-log(cor_min))
  ell_bounds[2L,] <- ell_bounds[2L,] / sqrt(-log(cor_max))
  rownames(ell_bounds) <- c("lower", "upper")

  # Convert back to original scale. Note that the additive shift does not 
  # need to be done since these quantities are distances, so the shift 
  # cancels out.
  if(!dim_by_dim) {
    ell_bounds <- mult_vec_with_mat_rows(X_bounds[2,]-X_bounds[1,], ell_bounds)
    dist_q <- mult_vec_with_mat_rows(X_bounds[2,]-X_bounds[1,], dist_q)
  }
  
  if(convert_to_square) {
    ell_bounds <- ell_bounds^2
    dist_q <- dist_q^2
  }
  
  return(list(ell_bounds=ell_bounds, dist_quantiles=dist_q))
}


get_pairwise_dist_quantiles_by_dim <- function(X, p=NULL) {
  # Computes non-zero pairwise distances over the rows of `X` for each 
  # dimension (column) of `X` independently. Then returns the empirical 
  # min, max and quantiles of these pairwise distances. Note that 
  # each column of `X` is considered independently. The min/max of the 
  # pairwise distances is returned by default, and additional quantiles
  # are optionally included by specifying `p`.
  #
  # Args:
  #    X: matrix of shape (n,d), containing `n` points in `d` dimensions. 
  #    p: numeric vector of probabilities. The quantiles to compute in addition
  #       to the min and max. If NULL, only returns min and max. 
  #
  # Returns:
  # matrix of dimension (2+length(p),d) [where d=ncol(X)] with the first 
  # row containing the minimum pairwise distances in each respective dimension, 
  # and the second row containing the maximum pairwise distances. 
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
    
    bounds[,j] <- c(bounds_j, quantiles_j)
  }
  
  # Row and column names. 
  rownames(bounds) <- c("min", "max", quantile_names)
  colnames(bounds) <- colnames(X)
  
  return(bounds)
}

get_pairwise_dist_quantiles <- function(X, p=NULL) {
  # Like `get_pairwise_dist_quantiles()` but computes the pairwise distances 
  # using d-dimensional Euclidean distance (viewing each row of `X`) as a 
  # point in d-dimensional space, instead of computing the distances dimension
  # by dimension.
  #
  # Returns:
  #  For consistency with `get_pairwise_dist_quantiles()` returns a matrix 
  #  of the same dimensions, even though the values will be the same across
  #  the dimensions. See `get_pairwise_dist_quantiles()` for details on the 
  #  returned matrix. 
  
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
  
  # Compute pairwise Euclidean distances. 
  dists <- as.matrix(dist(X, method="euclidean", diag=FALSE))
  dists <- dists[upper.tri(dists)]
  dists <- dists[dists > 0]
  
  # Compute min, max, and quantiles of pairwise distances. 
  bounds[1,] <- min(dists)
  bounds[2,] <- max(dists)
  if(!is.null(p)) {
    quantiles <- quantile(dists, p)
    for(j in 1:d) bounds[3:nrow(bounds),j] <- quantiles
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


get_bounds <- function(X) {
  # For an nxd input matrix, returns a 2xd matrix with the first row storing
  # the minimum value of X in each dimension, and likewise with the maximum
  # in the second row.
  
  apply(X, 2, range)
}


scale_inputs <- function(X, source_bounds=NULL, target_bounds=NULL, invert=FALSE) {
  # Linearly scales data, dimension by dimension, from one hyper-rectangle to another 
  # hyper-rectangle.
  #
  # Args:
  #    X: matrix of dimension (n,d), `n` data points in `d` dimensions. 
  #    source_bounds: (2,d) matrix, defining the bounds in each dimension. 
  #                   Defaults to `get_bounds(X)`. 
  #    target_bounds: (2,d) matrix or vector of length 2, defaults to the unit 
  #                   hypercube [0,1]^d. If a vector of length 2, then treated 
  #                   as bounds that will be applied to all dimensions, hence 
  #                   defining a hypercube. 
  #    invert: if TRUE, the roles of source and targert bounds are reversed; i.e., 
  #            the inverse map is computed. 
  #
  # Returns:
  #    matrix of dimension `dim(X)` containing the transformed data `X`. 

  assert_that(is.matrix(X))
  d <- ncol(X)

  # Construct the source bounds.
  if(is.null(source_bounds)) source_bounds <- get_bounds(X) 

  # Construct the bounds in each dimension defining the hyperrectangle that
  # the data will be scaled to lie within.
  if(is.null(target_bounds)) {
    target_bounds <- rbind(rep(0,d), rep(1,d))
  } else if(is.vector(target_bounds)) {
    assert_that(length(target_bounds) == 2)
    target_bounds <- matrix(target_bounds,ncol=1)[,rep(1,d), drop=FALSE]
  } else {
    assert_that(is.matrix(target_bounds))
    assert_that(all(dim(target_bounds)==c(2,d)))
  }

  # Invert transformation. 
  if(invert) {
    temp <- source_bounds
    source_bounds <- target_bounds
    target_bounds <- temp
  }
  
  # Linearly map data from source bounds to target bounds.
  source_diff <- source_bounds[2,] - source_bounds[1,]
  target_diff <- target_bounds[2,] - target_bounds[1,]
  bound_ratio <- target_diff/source_diff
  add_vec_to_mat_rows(target_bounds[1,], mult_vec_with_mat_rows(bound_ratio, add_vec_to_mat_rows(-source_bounds[1,], X)))
}




