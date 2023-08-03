# ---------------------------------------------------------------------------------------------------
# sequential_design_sim_study.r
# Simulation study for comparing different approaches to (batch) sequential design with GP emulator. 
#
# Andrew Roberts
# ---------------------------------------------------------------------------------------------------


library(lhs)
library(hetGP)
library(mlegp)
library(ggplot2)
library(viridis)
library(gridExtra)
library(data.table)
library(BayesianTools)

source("numerical_experiment_functions.r")
source("mcmc_calibration_functions.r")
source("gp_emulator_functions.r")
source("sequential_design_optimization.r")


# -----------------------------------------------------------------------------
# Setup:
#    - Specify forward model (i.e., computer model or simulator). 
#    - Set true values of calibration parameters. 
#    - Set true values of likelihood parameters. 
#    - Set prior distributions for calibration and likelihood parameters. 
# -----------------------------------------------------------------------------

# Computer model and synthetic data generation. 
computer_model_data <- generate_vsem_test_case(4)
print(computer_model_data$ref_pars[computer_model_data$pars_cal_sel,])

# Priors on calibration parameters. 
theta_prior_params <- computer_model_data$ref_pars[computer_model_data$pars_cal_sel,]
theta_prior_params[, "dist"] <- c("Uniform", "Uniform")
theta_prior_params[,"param1"] <- c(1.3, 0.4) 
theta_prior_params[,"param2"] <- c(1.7, 0.6)
theta_prior_params <- theta_prior_params[, c("dist", "param1", "param2")]
print(paste0(rep("-", 30), " Calibration parameters priors ", rep("-", 30)))
print(theta_prior_params)

# Priors on likelihood parameters. 
sig2_eps_prior_info <- get_IG_priors_numerical_test(sig2_true = diag(computer_model_data$Sig_eps), 
                                                    bias_frac = c(0.1, -0.15), bins = 50,
                                                    coef_var = c(0.3, 0.5), return_prior_plots = TRUE, 
                                                    output_variables = computer_model_data$output_vars)
sig2_eps_prior_params <- sig2_eps_prior_info$prior
print(paste0(rep("-", 30), " Likelihood parameters priors ", rep("-", 30)))
print(sig2_eps_prior_params)
















