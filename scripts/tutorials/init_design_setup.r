#
# init_design_setup.r
# A script that goes through some examples of setting up initial designs for 
# surrogate modeling applications. 
#

library(lhs)
library(ggplot2)
library(data.table)
library(docopt)

# For reproducibility.
set.seed(7543463)

# File path setup. 
base_dir <- getwd()
src_dir <- file.path(base_dir, "src")

# Load general helper functions.    
source(file.path(src_dir, "general_helper_functions.r"))

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

# The variable `n_time` defines the number of days that the model will be run 
# for. This will determine the dimension of the output space for the model G(u).
n_year <- 3
n_time <- 365*n_year
time_points <- 1:n_time

# The VSEM model is driven (i.e., "forced") by sunlight. This code simulates some 
# synthetic radiation data for this purpose. Don't worry about this for now, 
# but it needs to be set up to run the model. 
driver <- BayesianTools::VSEMcreatePAR(days=time_points)

# This sets up the map from model parameters to the time series of model outputs. 
# This model actually returns multiple time series of a few different output 
# variables, but we will start by looking at a single one: leaf area index (LAI). 
# This is a quantitative measure of the quantity of leaves in an area. 
output_var <- "LAI"
output_names <- get_vsem_output_names()
output_var_idx <- which(output_names==output_var) 
par_to_output_map <- get_vsem_fwd_model(driver, dim_par, par_cal_idx, par_default, simplify=FALSE)

# -----------------------------------------------------------------------------
# Bayesian Inverse Problem Setup
#
# We now set up the components defining a Bayesian inverse problem: the 
# forward model, data, likelihood, and prior.
# -----------------------------------------------------------------------------

# Forward model, a.k.a. the parameter-to-observable map. The function from 
# parameters to the model predictions of whatever quantity for which we 
# have data. As discussed below, we will be assuming we have daily LAI 
# observations, so the observation operator will simply extract the LAI 
# portion of the VSEM output. Note that this function is vectorized so 
# that you can evaluate at multiple values of `u`. To do so, the values 
# of `u` should be stacked as rows in a matrix. See examples later in this 
# script. 
fwd_model <- function(u) {
  # Run VSEM model with input parameters `u`. 
  model_outputs <- par_to_output_map(u)
  
  # Extract the output variable we are interested in. 
  single_run <- (dim(model_outputs)[1]==1L)
  var_trajectory <- model_outputs[,,output_var_idx]
  if(single_run) var_trajectory <- matrix(var_trajectory, nrow=1L)
  var_trajectory
}

# Simulating synthetic observations. Let's start simple and consider the 
# situation where we have daily noisy LAI observations. We will simulate 
# fake observations by assuming that the observations are generated 
# according to the VSEM model predictions at some "true" parameter value, 
# plus some Gaussian observation noise. For now, we will assume that we 
# know the variance of the observation noise. 

# Ground truth parameter and outputs. 
u_true <- get_vsem_default_pars()[par_cal_idx]
y_true <- drop(fwd_model(u_true))

# Observed data. `sig2_eps` is the observation variance. 
signal_to_noise_ratio <- 20
sig2_eps <- mean(y_true) / signal_to_noise_ratio
y <- y_true + sqrt(sig2_eps)*rnorm(n=length(y_true))

# Plot ground truth vs. observed data.
plt_data <- ggplot(data.frame(t=time_points, y_true=y_true, y=y)) + 
            geom_point(aes(x=t, y=y), color="black") + 
            geom_line(aes(x=t, y=y_true), color="red") + 
            xlab("days") + ylab("LAI") + ggtitle("True vs Observed LAI")
plot(plt_data)

# Likelihood: we will use the same likelihood that was used to generate the 
# data, which means we are assuming the statistical model is well-specified (i.e.,
# there is no model discrepancy). We can consider more realistic example later.
# The function defined here is actually the log-likelihood ("llik").
llik <- function(u) {
  # Run forward model (note that `u` can be a matrix of multiple parameter
  # values). 
  fwd <- fwd_model(u)
  
  # Compute Gaussian log-likelihood.
  squared_err <- add_vec_to_mat_rows(-y, fwd)^2
  -0.5*n_time*log(2*pi*sig2_eps) - 0.5*(1/sig2_eps)*rowSums(squared_err)
}


# Define prior distribution on calibration parameters. I wrote a function that
# sets up a default uniform prior. 
par_cal_idx <- which(par_names %in% par_cal_names)
par_prior_params <- get_vsem_default_priors()[par_cal_idx,,drop=FALSE]
rownames(par_prior_params) <- par_prior_params$par_name
print("Prior on calibration parameter:")
print(par_prior_params)





