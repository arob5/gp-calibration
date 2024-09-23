#
# init_design_setup.r
# A script that goes through some examples of setting up initial designs for 
# surrogate modeling applications. 
#

library(lhs)
library(ggplot2)
library(data.table)
library(docopt)

# File path setup. 
base_dir <- getwd()
src_dir <- file.path(base_dir, "src")

# Load functions for setting up tests with the VSEM model. 
source(file.path(src_dir, "inv_prob_test_functions.r"))

# Load functions for sampling from prior distribution.  
source(file.path(src_dir, "seq_design.r"))

# Load functions for plotting.   
source(file.path(src_dir, "plotting_helper_functions.r"))


# -----------------------------------------------------------------------------
# VSEM model setup.
#
# The Very Simple Ecosystem Model (VSEM) is a convenient toy model to use 
# for testing algorithms. It is implemented in the R package `BayesianTools`. 
# To make it easier to work with, I wrote a variety of wrapper functions which 
# are defined in `inv_prob_test_functions.r`. 
# -----------------------------------------------------------------------------

# Define the parameters that will be considered here. VSEM has 11 parameters, 
# but it is simpler to start with fewer and fix the rest at their default 
# values. Here we start by allowing 2 parameters to vary. The variable 
# `par_cal_names` is short for "calibration parameter names", the parameters you 
# will be varying.
par_cal_names <- c("KEXT", "tauV")
dim_par <- length(par_cal_names)
par_names <- get_vsem_par_names()
par_default <- get_vsem_default_pars()
print("Parameter defaults:")
print(data.frame(par_name=par_names, default_value=par_default))

# Define prior distribution on calibration parameters. I wrote a function that
# sets up default uniform priors, but we can consider alternatives in the 
# future. The data.frame `par_prior_params` that stores the information 
# storing the prior distributions is consistent with the formatting required 
# by the function `get_batch_design()` in `seq_design.r`
par_cal_idx <- which(par_names %in% par_cal_names)
par_prior_params <- get_vsem_default_priors()[par_cal_idx,,drop=FALSE]
rownames(par_prior_params) <- par_prior_params$par_name
print("Prior on calibration parameter:")
print(par_prior_params)


# -----------------------------------------------------------------------------
# Bayesian Inverse Problem Setup
#
# The Very Simple Ecosystem Model (VSEM) is a convenient toy model to use 
# for testing algorithms. It is implemented in the R package `BayesianTools`. 
# To make it easier to work with, I wrote a variety of wrapper functions which 
# are defined in `inv_prob_test_functions.r`. 
# -----------------------------------------------------------------------------










