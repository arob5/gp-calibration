#
# pymc_helper_functions.py
# Helper functions primary related to Gaussian processes in PyMC.
#
# Andrew Roberts
#

import pytensor.tensor as pt
import pymc as pm

class Quadratic(pm.gp.mean.Mean):
    R"""
    Quadratic mean function for Gaussian process. Currently
    uses basis functions {x,x^2} for each input dimension,
    plus an overall intercept. No interaction terms are included
    at present.

    Parameters
    ----------
    beta: variable, array or integer
          Coefficients on the basis functions, excluding intercept.
    beta0: variable, array (of length 1), or integer
           The intercept.
    """

    def __init__(self, coeffs, intercept=0):
        super().__init__()
        self.beta = coeffs
        self.beta0 = intercept

    def __call__(self, X):
        H =
        return pt.squeeze(pt.dot(H, self.beta) + self.beta0)


        return tt.alloc(1.0, X.shape[0]) * self.c
