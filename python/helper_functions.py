#
# helper_functions.py
# Helper functions related to Gaussian process emulation.
#
# Andrew Roberts
#

import numpy as np
import scipy as sp

# For an nxd input array, returns a 2xd array with the first row storing
# the minimum value of X in each dimension, and likewise with the maximum
# in the second row.
def get_bounds(X):
    return np.vstack((np.min(X, axis=0), np.max(X, axis=0)))

# data_bounds: (2,d) array, defaults to `get_bounds(X)`.
# target_bounds: (2,d) or (2,) array, defaults to [0,1]^d.
# If `invert` is True, the roles of source and target bounds are reversed.
def scale_inputs(X, source_bounds=None, target_bounds=None, invert=False):
    d = X.shape[1]

    # Construct the source bounds.
    if source_bounds is None:
        source_bounds = get_bounds(X)

    # Construct the bounds in each dimension defining the hyperrectangle that
    # the data will be scaled to lie within.
    if target_bounds is None:
        target_bounds = np.tile(np.array([0,1]).reshape((2,1)), [1,d])
    elif target_bounds.shape == (2,):
        target_bounds = np.tile(target_bounds.reshape((2,1)), [1,d])
    else:
        assert target_bounds.shape == (2,d)

    if invert:
        source_bounds,target_bounds = target_bounds,source_bounds

    # Linearly map data from source bounds to target bounds.
    return target_bounds[0,:] + np.multiply(X - source_bounds[0,:],
                                           (target_bounds[1,:]-target_bounds[0,:]) / (source_bounds[1,:] - source_bounds[0,:]))


# By default, this function implicitly treats data as lying in the hypercube [-1,1]^d.
# One can think of [-1,1]^d as the hypercube associated with the bounds of
# the design points, so test points generated outside of this cube can be thought of
# as testing the extrapolation of the model.
# The output test points can be appropriately scaled with respect to a different
# hypercube by specifiying the `target_bounds` argument.
# Returns an array `X_test` of shape (d, N_points, d), where `X_test[j,:,:]` gives
# the `N_points` test inputs spread along the jth standard basis vector.
def gen_extrapolation_test_inputs(d, N_points, max_scaler=2.0, scale="linear", target_bounds=None):

    # These scalers are wrt the hypercube [-1,1]^d.
    reference_bounds = np.tile(np.array([-1,1]).reshape((2,1)), [1,d])
    if scale=="linear":
        scalers = np.linspace(start= -max_scaler, stop=max_scaler, num=N_points)
    elif scale=="log":
        raise NotImplementedError

    X_test = np.empty((d, N_points, d))
    for j in range(d):
        basis_vec = np.zeros(d)
        basis_vec[j] = 1.0
        X_test[j,:,:] = (np.tile(basis_vec.reshape((d,1)), [1,N_points]) * scalers).T

        if target_bounds is not None:
            X_test[j,:,:] = scale_inputs(X_test[j,:,:], source_bounds=reference_bounds, target_bounds=target_bounds)

    return X_test


# Given an array `X` of shape (n_point, n_dim), returns an array `pairwise_ranges`
# of shape (2, n_dim), where `pairwise_ranges[0,d]` is the minimum pairwise distance
# between distinct points in `X[:,d]`. Similarly, `pairwise_ranges[1,d]` contains
# the maximum pairwise distance in the dth dimension. The distance is in the
# Euclidean sense, but this could be generalized by passing additional arguments
# to `sp.spatial.distance.pdist`.
def get_pairwise_dist_range_per_dim(X):
    n_dim = X.shape[1]
    pairwise_ranges = np.empty((2,n_dim))

    for d in range(n_dim):
        pairwise_dists = sp.spatial.distance.pdist(X[:,d,np.newaxis])
        pairwise_ranges[0,d] = np.min(pairwise_dists)
        pairwise_ranges[1,d] = np.max(pairwise_dists)

    return pairwise_ranges

# Returns the inverse gamma parameters `alpha` and `beta` that achieve a given
# target mean `m` and standard deviation `s`. Returns a tuple containing the values
# of `alpha` and `beta`, respectively.
def get_IG_pars_from_mean_sd(m, s):
    v = np.power(s,2)
    alpha = 2 + np.power(m,2)/v
    beta = (alpha-1)*m

    # Condition for mean and variance to exist.
    assert alpha > 2

    return (alpha,beta)
