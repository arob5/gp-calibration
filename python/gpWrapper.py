#
# gpWrapper.r
# Class definitions for the `gpWrapper` class, and classes which inherit from
# this class.
#
# Andrew Roberts
#
# Depends: helper_functions.py
#

import numpy as np
import scipy as sp
import helper_functions as hf

# -----------------------------------------------------------------------------
# gpWrapper: the base Gaussian Process (GP) wrapper class.
#
# This class represents a "wrapper" around a GP implementation, potentially from
# a different Python package or a user-defined GP implementation.
# It serves as an interface between these other GP packages and PEcAn code.
# The idea is to define a modular object with
# a standardized interface providing access to common GP functionality
# (fit, predict, sample, plot, etc.) that can then be plugged into downstream
# tasks. This allows different GP packages to be swapped in without modifying
# the downstream code. `gpWrapper` defines a generic independent multi-output
# GP; i.e., `gpWrapper` objects may represent a collection of independent
# univariate GPs.
# -----------------------------------------------------------------------------

class gpWrapper:
    def __init__(self, X, Y, scale_input=False, normalize_output=False,
                 x_names=None, y_names=None,
                 default_jitter=np.sqrt(np.finfo(float).eps), *args):

        # Argument checking.
        assert isinstance(X, np.ndarray)
        assert isinstance(Y, np.ndarray)
        assert not np.isnan(X).any()
        assert not np.isnan(Y).any()
        assert X.ndim==2
        if Y.ndim==1:
            Y = Y.reshape((-1,1))
        assert Y.ndim==2

        # Initialize data attributes.
        self.X = X
        self.Y = Y
        self.X_dim = X.shape[1]
        self.Y_dim = Y.shape[1]

        # Default names if not provided.
        if x_names is None:
            x_names = ["x"+str(i+1) for i in range(self.X_dim)]
        if y_names is None:
            y_names = ["y"+str(i+1) for i in range(self.Y_dim)]

        # Initialize more data attributes. `X_names` and `Y_names` define
        # the official order of the X and Y variables, respectively.
        self.X_names = x_names
        self.Y_names = y_names
        self.scale_input = scale_input
        self.normalize_output = normalize_output
        self.X_bounds = hf.get_bounds(self.X)

        if self.normalize_output:
            self.Y_mean = np.mean(self.Y, axis=0)
            self.Y_std = np.std(self.Y, axis=0)
            self.Y_train = self.normalize(self.Y)
        else:
            self.Y_mean = None
            self.Y_std = None
            self.Y_train = self.Y # TODO: deep copy here instead?

        if self.scale_input:
            self.X_train = self.scale(self.X)
        else:
            self.X_train = self.X

        # The default jitter value; sometimes also called the "nugget" or
        # "nugget variance".
        self.default_jitter = default_jitter

    def normalize(self, Y_new, invert=False, transform_variance=False):
        if transform_variance:
            stop("not yet implemented.")
        return hf.normalize_matrix(Y_new, cst_add=-self.Y_mean,
                                   cst_mult=1/self.Y_std, invert=invert)

    def scale(self, X_new, invert=False):
        return hf.scale_inputs(X_new, source_bounds=self.X_bounds,
                               target_bounds=None, invert=invert)

    def fit_package(X_fit, Y_fit, kernel_name=None, mean_func_name=None, **kwargs):
        raise NotImplementedError

    def predict_package(X_new, output_name, return_mean=True, return_var=True,
                        return_cov=False, return_cross_cov=False, X_cross=None,
                        include_jitter=True):
        raise NotImplementedError
