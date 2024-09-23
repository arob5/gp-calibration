#
# basis_func_approx_setup.r
# A script that sets up the VSEM model, for Xinyu to begin testing representing
# the model outputs with respect to basis functions.
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
source(file.path(src_dir, "statistical_helper_functions.r"))

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
# values. Here we are starting with the simplest case of a single parameter
# to make visualization easy. The variable `par_cal_names` is short for 
# "calibration parameter names", the parameters you will be varying. The 
# varible `par_names` contains all 11 VSEM parameters.
par_cal_names <- "KEXT"
dim_par <- length(par_cal_names)
par_names <- get_vsem_par_names()
par_default <- get_vsem_default_pars()
print("Parameter defaults:")
print(data.frame(par_name=par_names, default_value=par_default))

# Define prior distribution on calibration parameters. I wrote a function that
# sets up a default uniform prior. 
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

# This sets up the map G(u); i.e., the function that takes input parameters 
# `u` and returns the time series of model outputs. This model actually returns 
# multiple time series of a few different output variables, but we will start 
# by looking at a single one: leaf area index (LAI). This is a quantitative 
# measure of the quantity of leaves in an area. 

# Map from calibration parameter to VSEM outputs.
output_var <- "LAI"
output_names <- get_vsem_output_names()
output_var_idx <- which(output_names==output_var) 
par_to_output_map <- get_vsem_fwd_model(driver, dim_par, par_cal_idx, par_default, simplify=FALSE)

G <- function(u) {
  # Run VSEM model with input parameters `u`. 
  model_outputs <- par_to_output_map(u)
  
  # Extract the output variable we are interested in. 
  single_run <- (dim(model_outputs)[1]==1L)
  var_trajectory <- model_outputs[,,output_var_idx]
  if(single_run) var_trajectory <- matrix(var_trajectory, nrow=1L)
  var_trajectory
}

# -----------------------------------------------------------------------------
# Exploring the forward model G(u)
# -----------------------------------------------------------------------------

# To start, try running the model at a single input and investigate the output. 
# We sample the input parameter from its prior distribution using the 
# `sample_prior_theta()` function. The function `G()` returns a matrix of shape 
# (1,n_time). We plot the output using ggplot. The periodic behavior of the 
# output is reflecting seasonal trends (e.g., there are few leaves in the winter 
# and many in the summer). The function `get_batch_design()` implements various 
# sampling methods: calling it with `method="simple"` and `N_batch=1` simply 
# draws a single value from the prior distribution. 
u_test <- get_batch_design("simple", N_batch=1, prior_params=par_prior_params)
g_test <- G(u_test)
plt_test <- ggplot(data=data.frame(t=time_points, y=drop(g_test))) + 
            geom_line(aes(x=t, y=y)) + 
            xlab("day") + ylab(output_var) + ggtitle("VSEM output")
plot(plt_test)

# Now let's try to run the model at a few different parameter values at the 
# same time. I set up the function `G()` so that this can be easily 
# accomplished. To do so, you must pass a matrix to `G()`, with each row 
# corresponding to a different input value. We now sample multiple values of the 
# parameter from its prior distribution using the `N_batch` argument. In this case, 
# the function `G()` will return a matrix of dimension `(n_samples, n_time)`. So 
# each row gives the model output for the corresponding input parameter value.
# To plot the output for all of these model runs, this time we call a helper 
# function `plot_curves_1d_helper()`, which is defined in `plotting_helper_functions.r`.
# This function plots one curve per column of the input matrix, so we must 
# transpose the model output `g_test` before passing it to the plotting function.
# Observe the effect that varying the parameter value has on the model outputs.
n_samples <- 5
u_test <- get_batch_design("simple", N_batch=n_samples, prior_params=par_prior_params)
g_test <- G(u_test)

plt_test <- plot_curves_1d_helper(time_points, t(g_test), 
                                  plot_title="VSEM output: varying parameter values",
                                  xlab="days", ylab=output_var)
plot(plt_test)














