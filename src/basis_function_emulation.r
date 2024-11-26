# 
# basis_function_emulation.r
# Functions related to emulating functions with many outputs by approximating 
# the outputs as linear combinations of basis vectors.
#
# Xinyu Xu, Andrew Roberts
#
# Dependencies:
#    general_helper_functions.r
#

# ------------------------------------------------------------------------------
# Projection-based basis constructions.
#
# Functions assuming orthonormal bases. In particular, considers 
# decompositions of the form:
#
# G(u) = m + sum_{j=1}^{p} w_j(u) b_j, where 
# 
# * m is a p-dimensional vector.
# * The weights w_j(u) take the inner product form w_j(u) = <G(u)-m, b_j>.
# * b_1, ..., b_p are orthonormal p-dimensional vectors.
#
# An approximation can be achieved by truncating the above sum after r < p
# basis vectors. The canonical choice of basis vectors here are those provided 
# by principal components analysis (PCA).
# ------------------------------------------------------------------------------

eval_basis_weights <- function(U, fwd, B, m=NULL) {
  # Computes the weights w_j(u) = <G(u)-m, b_j> by first evaluating the forward 
  # model and then computing the inner product. Vectorized over multiple input 
  # points u and multible basis vectors b_j.
  #
  # Args:
  #    U: matrix of dimension (n,p), with rows containing the inputs.
  #    fwd: function, the forward model G(u). As of now, it must be vectorized
  #         such that `G(U)` returns a matrix of shape (n,p) containing the 
  #         model output for the ith input in its ith row.
  #    B: matrix of dimension (p,r), with columns containing orthonormal vectors
  #    m: optional numeric vector of length p, which will be subtracted from 
  #       forward model outputs if provided.
  #
  # Returns:
  #  matrix of dimension (n,r). The (i,j) entry contains w_j(u_i), with u_i
  #  being the ith row of `U`.
  
  if(is.null(m)) m <- rep(0, ncol(B))
  
  # Evaluate forward model at inputs and center the outputs.
  G <- add_vec_to_mat_rows(-m, fwd(U))
  
  # Compute the weights (magnitude of the projections).
  project_orthog_scalar(G, B)
}


project_orthog_scalar <- function(G, B) {
  # Compute the magnitude of the projections of vectors onto an orthonormal 
  # basis `B` - i.e., the "scalar projections". 
  #
  # Args:
  #    G: matrix of dimension (n,p), with rows containing the vectors to project.
  #    B: matrix of dimension (p,r), with columns containing orthonormal vectors
  #
  # Returns:
  #    matrix of dimension (n,r). The (i,j) entry contains <g_i,b_j>, the 
  #    inner product of the ith vector and the jth basis vector.
  
  G %*% B
}


pca <- function(G, r=ncol(G)) {
  # Given a (n,p) matrix `G` with observations given in the rows, computes 
  # a principal components decomposition, truncated after the first r 
  # dominant directions.
  #
  # Args:
  #    G: matrix, of dimension (n,p). Rows represent observations and columns 
  #       represent features/variables/dimensions.
  #    r: integer in [1,ncol(G)], the number of principal components to return. 
  #       Defaults is to not throw away any. 
  #
  # Returns:
  # list, with elements 
  #    mean: numeric vector of length p, the sample mean of the rows of 
  #          `G`. 
  #    vec: matrix of shape (p,r), the dominant r eigenvectors of the empirical 
  #         covariance matrix of the rows of `G`. The vectors are normalized to 
  #         have Euclidean norm 1, and are sorted in descending order according
  #         to the magnitude of their corresponding eigenvalues.
  #    val: numeric vector of length r, the dominant r eigenvalues.
  
  # Center the matrix.
  g_mean <- colMeans(G)
  G_centered <- add_vec_to_mat_rows(-g_mean, G)

  # Eigendecomposition of covariance.
  pca_result <- prcomp(G_centered, center=FALSE, scale.=FALSE, rank.=r)
  sqrt_eigvals <- pca_result$sdev
  eigvecs <- pca_result$rotation
   
  return(list(mean=g_mean, vec=eigvecs, sqrt_val=sqrt_eigvals))
}






