#
# run_multidim_toy_example.r
# Bayesian inversion for synthetic examples with multidimensional 
# parameter spaces. This file is written to be able to run on a 
# computing cluster. 
#
# Andrew Roberts
#

# -----------------------------------------------------------------------------
# docopt string for parsing command line arguments.  
# -----------------------------------------------------------------------------

"Usage:
  test_docopt.r <run_id> <output_dir> [options]
  test_docopt.r (-h | --help)

Options:
  -h --help                                 Show this screen.
  --dim_par=<dim_par>                       Dimension of parameter space. 
  --dim_output=<dim_output>                 Dimension of output space. 
  --N_design=<N_design>                     Number of design points.
  --design_method=<design_method>           Algorithm to generate design points. 
  --N_design_test=<N_design_test>           Number of validation points for emulators. 
  --design_method_test=<design_method_test> Algorithm to generate validation points. 
  --mcmc_tags=<mcmc_tag,mcmc_tag>           GP-approx MCMC algorithms. 
  --N_mcmc=<N_mcmc>                         Number of MCMC iterations. 
" -> doc


# -----------------------------------------------------------------------------
# Setup 
# -----------------------------------------------------------------------------

print("---------------------------------- Setup ----------------------------------")

# Base directory (i.e., project directory). Required for loading the R project (renv). 
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")

# Load R project. 
renv::load(base_dir)
print(".libPaths()")
print(.libPaths())
renv::status()

library(lhs)
library(ggplot2)
library(viridis)
library(parallel)
library(gridExtra)
library(data.table)
library(docopt)

# Read command line arguments. 
arguments <- docopt(doc)
run_id <- arguments$run_id
output_dir <- arguments$output_dir
required_settings <- c("dim_par", "dim_output", "N_design", "design_method", 
                       "N_design_test", "design_method_test", "N_mcmc")
settings <- arguments[required_settings]

# Source and output paths.  
src_dir <- file.path(base_dir, "src")

# Source scripts. 
source(file.path(src_dir, "gpWrapper.r"))
source(file.path(src_dir, "general_helper_functions.r"))
source(file.path(src_dir, "statistical_helper_functions.r"))
source(file.path(src_dir, "plotting_helper_functions.r"))
source(file.path(src_dir, "mcmc_helper_functions.r"))
source(file.path(src_dir, "llikEmulator.r"))
source(file.path(src_dir, "gp_emulator_functions.r"))
source(file.path(src_dir, "sim_study_functions.r"))
source(file.path(src_dir, "mcmc_calibration_functions.r"))
source(file.path(src_dir, "seq_design_gp.r"))
source(file.path(src_dir, "seq_design.r"))
source(file.path(src_dir, "gp_mcmc_functions.r"))

# Seeds for random number generator for different portions of analysis. 
seed_data <- 82
seed_design <- 500

# Number of parameters to calibrate. Will determine the number of basis 
# functions in the forward model. 
if(!is.null(settings$dim_par)) settings$dim_par <- as.integer(settings$dim_par)

# Number of forward model outputs. For forward model emulation, this means 
# one GP must be fit per output. 
if(is.null(settings$dim_output)) {settings$dim_output <- 10L
} else {settings$dim_output <- as.integer(settings$dim_output)}


# Design settings for emulator construction. 
if(is.null(settings$N_design)) {
  settings$N_design <- 10L*settings$dim_par
} else {settings$N_design <- as.integer(settings$N_design)}

if(is.null(settings$design_method)) settings$design_method <- "LHS"

# Design settings for emulator validation. 
if(is.null(settings$N_design_test)) {
  settings$N_design_test <- min(100L*settings$dim_par, 500L)
} else {settings$N_design_test <- as.integer(settings$N_design_test)}

if(is.null(settings$design_method_test)) settings$design_method_test <- "LHS"

# Global likelihood normalization settings. 
default_conditional <- FALSE
default_normalize <- TRUE

# GP-accelerated MCMC algorithms to run. 
if(is.null(settings$mcmc_tags)) {
  settings$mcmc_tags <- c("gp-mean", "gp-marg", "gp-quantile", "mcwmh-joint", "mcwmh-ind")
} else {
  settings$mcmc_tags <- unlist(strsplit(settings$mcmc_tags, ","))
}

# Settings for MCMC runs.  
if(is.null(settings$N_mcmc)) {settings$N_mcmc <- 50000L
} else {settings$N_mcmc <- as.integer(settings$N_mcmc)}

# Create global variables with the setting values. 
missing_settings <- names(settings)[sapply(settings, is.null)]
if(length(missing_settings) > 0) stop("Missing settings: ", missing_settings)
dim_par <- settings$dim_par
dim_output <- settings$dim_output
N_design <- settings$N_design
design_method <- settings$design_method
N_design_test <- settings$N_design_test
design_method_test <- settings$design_method_test
mcmc_tags <- settings$mcmc_tags
N_mcmc <- settings$N_mcmc

#
# Create alphanumeric ID defining the run, and use to create directory. 
#

# MOVED THIS TO BASH FILE 
# # Run ID. 
# run_id <- paste(run_tag, paste0("d", dim_par), paste0("p", dim_output), 
#                 paste0("N", N_design), design_method, sep="_")
# 
# # Output directory. 
# output_dir <- file.path(base_output_dir, run_id)

# If `output_dir` already exists, append timestep and create new directory.
# Otherwise create the directory. 
if(dir.exists(output_dir)) {
  stop("`output_dir` already exists.")
  # timestamp <- as.character(Sys.time())
  # output_dir <- paste(output_dir, timestamp, sep="_")
  # message("`output_dir` already exists. Creating new `output_dir:`")
  # message(output_dir)
  # dir.create(output_dir)
} else {
  dir.create(output_dir, recursive=TRUE)
}
save(settings, file=file.path(output_dir, "settings.RData"))

# Print settings. 
print("--------------------- User specified settings -------------------------")

print("--> General settings")
print(paste0("run ID: ", run_id))
print(paste0("output directory: ", output_dir))
print(paste0("Data generation seed (seed_data): ", seed_data))
print(paste0("Design data seed (seed_design): ", seed_design))
print(paste0("Likelihood densities conditional: ", default_conditional))
print(paste0("Likelihood densities normalized: ", default_normalize))

print("--> Inverse problem")
print(paste0("Parameter space dimension: ", dim_par))
print(paste0("Forward model output dimension: ", dim_output))

print("--> Design")
print(paste0("Number design points: ", N_design))
print(paste0("Design method: ", design_method))
print(paste0("Number test points: ", N_design_test))
print(paste0("Test design method: ", design_method_test))

print("--> MCMC settings")
print(paste0("Number iterations: ", N_mcmc))
print(paste0("GP-accelerated MCMC algs: ", paste(mcmc_tags, collapse=", ")))


# -----------------------------------------------------------------------------
# Bayesian inverse problem setup. 
# -----------------------------------------------------------------------------

print("--------------------- Bayesian Inverse Problem Setup -------------------------")

set.seed(seed_data)

#
# Forward model: linear combination of basis functions. 
#

# Basis functions: restrict to the number specified by `dim_par`. 
basis_func_list <- list(phi0=function(t) rep(1, length(t)), phi1=function(t) t,
                        phi2=function(t) cos(2*pi*t), phi3=function(t) cos(2*pi*2*t),
                        phi4=function(t) sin(2*pi*t), phi5=function(t) sin(2*pi*t/3),
                        phi6=function(t) sin(2*pi*t/6), phi7=function(t) cos(2*pi*t/6),
                        phi8=function(t) cos(2*pi*t/5), phi9=function(t) cos(2*pi*t/8),
                        phi10=function(t) sin(2*pi*t/8))
basis_func_list <- basis_func_list[1:dim_par]

# Discretize the unit time interval to define the dimension of the output space. 
times <- seq(0, 1, length.out=dim_output)

# Forward model. 
eval_with_args <- function(FUN, ...) FUN(...)
Phi <- mapply(eval_with_args, basis_func_list, list(times))

fwd <- function(par) {
  Phi %*% matrix(par, ncol=1)
}

fwd_vect <- function(par_mat) {
  # `par_mat` assumed to have dimension (num inputs, dim_par). Returns 
  # output of dimension (num inouts, dim_output). 
  par_mat %*% t(Phi)
}

#
# Likelihood. 
#

# Observation covariance. 
cov_obs <- diag(0.5, nrow=dim_output)

#
# Prior. 
#

# Prior distribution. 
par_prior <- data.frame(dist="Gaussian", param1=rep(1,dim_par), param2=rep(1,dim_par), 
                        row.names=paste0("u", 1:dim_par))

#
# Synthetic data: ground truth and noisy observations. 
# 

# Ground truth. 
par_true <- sample_prior_theta(par_prior)
output_true <- fwd(par_true)

# Noisy observations. 
y <- output_true + t(chol(cov_obs)) %*% matrix(rnorm(dim_output), ncol=1)

#
# Collect inverse problem components in list, save plots summarizing data.  
#

# List containing inverse problem components. 
inv_prob <- list(y=y, fwd=fwd, fwd_vect=fwd_vect, par_prior=par_prior, 
                 cov_obs=cov_obs, par_true=par_true, output_true=output_true,
                 times=times, basis_func_list=basis_func_list)

# Plot ground truth and observations. 
df_inv_prob_data <- data.frame(t=times, y_true=output_true, y=y)
plt_outputs <- ggplot(data=df_inv_prob_data) +
               geom_point(aes(x=t, y=y), color="black") + 
               geom_point(aes(x=t, y=y_true), color="red") + 
               ylab("output") + ggtitle("True and Observed y")

plot(plt_outputs)
ggsave(file.path(output_dir, "output_data.png"), plt_outputs)


# -----------------------------------------------------------------------------
# Inversion using exact forward model. 
# -----------------------------------------------------------------------------

print("--------------------- Bayesian Inversion with True Forward Model -------------------------")

# Store exact log-likelihood object. 
llik_exact <- llikEmulatorExactGaussDiag(llik_lbl="exact", fwd_model=inv_prob$fwd, 
                                         fwd_model_vectorized=inv_prob$fwd_vect,
                                         y_obs=t(inv_prob$y), dim_par=dim_par,
                                         use_fixed_lik_par=TRUE, sig2=diag(cov_obs),
                                         par_names=rownames(inv_prob$par_prior), 
                                         default_conditional=default_conditional, 
                                         default_normalize=default_normalize)
inv_prob$llik_obj <- llik_exact
save(inv_prob, file=file.path(output_dir, "inv_prob_list.RData"))


# -----------------------------------------------------------------------------
# Initial design and test points.  
# -----------------------------------------------------------------------------

print("--------------------- Constructing design -------------------------")

set.seed(seed_design)

# Generate design. 
design_info <- list(design_method=design_method, N_design=N_design, seed=seed_design)
design_info$input <- get_batch_design("LHS", N_design, prior_params=inv_prob$par_prior)
design_info$fwd <- llik_exact$run_fwd_model(design_info$input)
design_info$llik <- llik_exact$assemble_llik(design_info$input)
design_info$lprior <- calc_lprior_theta(design_info$input, inv_prob$par_prior)
save(design_info, file=file.path(output_dir, "design_info.RData"))

# Test points. 
u_grid <- get_batch_design(design_method_test, N_design_test, prior_params=inv_prob$par_prior)
test_info <- list(input=u_grid, design_method=design_method_test, N_design=N_design_test)
test_info$fwd <- llik_exact$run_fwd_model(test_info$input) 
test_info$llik <- llik_exact$assemble_llik(test_info$input)
test_info$lprior <- calc_lprior_theta(test_info$input, inv_prob$par_prior)
save(test_info, file=file.path(output_dir, "test_info.RData"))


# -----------------------------------------------------------------------------
# Exact MCMC. 
# -----------------------------------------------------------------------------

# MCMC sampling using exact likelihood. 
mcmc_par_init <- design_info$input[which.max(design_info$llik + design_info$lprior),]
mcmc_exact_list <- mcmc_gp_noisy(inv_prob$llik_obj, inv_prob$par_prior, N_itr=N_mcmc, 
                                 mode="MCMH", par_init=mcmc_par_init)
samp_dt <- format_mcmc_output(mcmc_exact_list$samp, test_label="exact")
mcmc_list <- list(exact=mcmc_exact_list[setdiff(names(mcmc_exact_list), "samp")])


# -----------------------------------------------------------------------------
# Fitting Emulators.   
# -----------------------------------------------------------------------------

print("--------------------- Fitting Emulators -------------------------")

llik_em_list <- list()

# Log-likelihood emulator, constant mean. 
print("-> Fitting em_llik_const")
em_llik_gp <- gpWrapperHet(design_info$input, matrix(design_info$llik, ncol=1), 
                           normalize_output=TRUE, scale_input=TRUE)
em_llik_gp$fit("Gaussian", "constant", estimate_nugget=FALSE)
llik_em_list[["em_llik_const"]] <- llikEmulatorGP("em_llik_const", em_llik_gp, default_conditional=default_conditional, 
                                                  default_normalize=default_normalize, lik_par=diag(cov_obs), 
                                                  use_fixed_lik_par=TRUE)

# Log-likelihood emulator, quadratic mean.
print("-> Fitting em_llik_quad")
em_llik_gp2 <- gpWrapperKerGP(design_info$input, matrix(design_info$llik, ncol=1), 
                              normalize_output=TRUE, scale_input=TRUE)
em_llik_gp2$fit("Gaussian", "quadratic", estimate_nugget=FALSE, 
                optimFun="nloptr::nloptr", trace=TRUE, multistart=10)
llik_em_list[["em_llik_quad"]] <- llikEmulatorGP("em_llik_quad", em_llik_gp2, default_conditional=default_conditional, 
                                                 default_normalize=default_normalize, lik_par=diag(cov_obs), 
                                                 use_fixed_lik_par=TRUE)

# Log-likelihood emulator, constant mean, Gaussian plus quadratic kernel. 
print("-> Fitting em_llik_const_GaussQuad")
em_llik_gp3 <- gpWrapperKerGP(design_info$input, matrix(design_info$llik, ncol=1), 
                              normalize_output=TRUE, scale_input=TRUE)
em_llik_gp3$fit("Gaussian_plus_Quadratic", "constant", estimate_nugget=FALSE, 
                optimFun="nloptr::nloptr", trace=TRUE, multistart=10)
llik_em_list[["em_llik_const_GaussQuad"]] <- llikEmulatorGP("em_llik_const_GaussQuad", em_llik_gp3, 
                                                            default_conditional=default_conditional, 
                                                            default_normalize=default_normalize, 
                                                            lik_par=diag(cov_obs), use_fixed_lik_par=TRUE)
                                                            
# Forward model emulator.
print("-> Fitting em_fwd")
em_fwd_gp <- gpWrapperHet(design_info$input, design_info$fwd, normalize_output=TRUE, scale_input=TRUE)
em_fwd_gp$fit("Gaussian", "constant", estimate_nugget=FALSE)
llik_em_list[["em_fwd"]] <- llikEmulatorFwdGaussDiag("em_fwd", em_fwd_gp, t(inv_prob$y), sig2=diag(cov_obs), 
                                                     default_conditional=default_conditional, 
                                                     default_normalize=default_normalize, use_fixed_lik_par=TRUE, 
                                                     par_names=rownames(inv_prob$par_prior))

# Save log-likelihood emulators. 
save(llik_em_list, file=file.path(output_dir, "llik_em_list.RData"))


# -----------------------------------------------------------------------------
# Log-Likelihood predictions at test points.   
# -----------------------------------------------------------------------------

# Emulator predictions. 
print("--------------------- Emulator Predictions at Test Points -------------------------")
emulator_pred_list <- list()
for(em_name in names(llik_em_list)) {
  print(paste0("Predicting with emulator: ", em_name))
  emulator_pred_list[[em_name]] <- llik_em_list[[em_name]]$predict_emulator(test_info$input)
}
save(emulator_pred_list, file=file.path(output_dir, "emulator_pred_list.RData"))


# Log-likelihood predictions. 
print("--------------------- Log-Likelihood Predictions at Test Points -------------------------")
for(em_name in names(llik_em_list)) {
  plt <- llik_em_list[[em_name]]$plot_pred_validation(test_info$input, true_llik=test_info$llik, 
                                                      emulator_pred_list=emulator_pred_list[[em_name]], 
                                                      plot_type="llik", conditional=default_conditional,
                                                      normalize=default_normalize, include_interval=TRUE,
                                                      interval_method="pm_std_dev", N_std_dev=1.64, 
                                                      plot_title=paste0("llik predictions: ", em_name))
  ggsave(file.path(output_dir, paste0("llik_pred_", em_name, ".png")), plt)
}


# -----------------------------------------------------------------------------
# Run MCMC Samplers.  
# -----------------------------------------------------------------------------

print("--------------------- GP-Accelerated MCMC -------------------------")

for(em_name in names(llik_em_list)) {
  print(paste0("----- ", em_name, " -----"))
  mcmc_info_list <- run_approx_mcmc_comparison(inv_prob_list=inv_prob, 
                                               llik_em_obj=llik_em_list[[em_name]], 
                                               mcmc_tags=mcmc_tags, save_dir=output_dir, 
                                               samp_dt=samp_dt, mcmc_list=mcmc_list, 
                                               test_label_suffix=em_name, N_itr=N_mcmc, 
                                               par_init=mcmc_par_init, overwrite=TRUE)
  
  samp_dt <- mcmc_info_list$samp
  mcmc_list <- mcmc_info_list$mcmc_list
}


print("-------------------------- End of Script ------------------------------")




