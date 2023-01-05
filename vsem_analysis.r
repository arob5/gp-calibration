#
# vsem_analysis.r
# Testing parameter calibration methods for the "Very Simple Ecosystem Model" (VSEM)
# implemented in the BayesianTools R package. 
#
# Andrew Roberts
#

library(data.table)
library(BayesianTools)

set.seed(5)

#
# Setup: settings control synthetic data generation and which parameters to calibrate
#

# Number days in time series 
N_days <- 1000

# Standard deviation of observation error (for now assuming Gaussian noise) for each 
# output. Also store minimum and lower bounds on the standard deviation to align with the 
# reference parameter formatting in VSEMgetDefaults(). Note that the NEE output of the
# model VSEM() is scales by 1000, so the observation error magnitude should correspond 
# to this scale. 
ref_pars_lik <- data.frame(best = c(2.0, 2.0, 2.0, 2.0), 
                           lower = c(0.1, 0.1, 0.1, 0.1), 
                           upper = c(4.0, 4.0, 4.0, 4.0))
rownames(ref_pars_lik) <- c("NEE", "Cv", "Cs", "CR")

# Names of parameters to calibrate. This can include the likelihood parameters
# to calibrate the observation noise variances, etc. 
pars_cal <- c("KEXT", "sig_eps")

# Identify which outputs to constrain in the model. Each constrained output factors 
# into the likelihood. Choices are "NEE" and the three carbon pools "Cv", "Cs", and "CR".
output_vars <- c("NEE", "Cv")

#
# Synthetic data generation
#

# Create time series of Phosynthetically Active Radiation (PAR) which is the forcing variable
# for the VSEM model. 
PAR <- VSEMcreatePAR(seq_len(N_days))

# The "best" column of the reference parameters are used to generate the "true" 
# output data. We will add noise to this data to simulate the field observations.
ref_pars <- VSEMgetDefaults()
data_ref <- as.data.table(VSEM(ref_pars$best, PAR))

# Add observational noise; NEE is scaled by 1000. Fow now considering independent 
# outputs, but can consider models that fill in the off-diagonal Sig_eps values
# in the future. 
data_ref[, NEE := NEE*1000]
Sig_eps <- diag(ref_pars_lik[colnames(data_ref), "best"])
Lt <- sqrt(Sig_eps)
Z <- matrix(rnorm(N_days*N_outputs), N_days, N_outputs) 
data_obs <- data_ref + Z * Lt

# Index selector for calibration parameters
pars_cal_sel <- which(rownames(ref_pars) %in% pars_cal)

#
# Likelihood
#

llik_Gaussian <- function(par, par_lik, par_ref, par_cal_sel, output_vars, PAR, data_obs) {
  # This is more or less the simplest possible likelihood to use in this setting. 
  # It assumes iid Gaussian likelihoods for each output, and assumes independence 
  # across outputs. 
  
  # Parameters not calibrated are fixed at default values
  theta <- par_ref$best
  theta[par_cal_sel] <- par
  
  # Run forward model, re-scale NEE.
  pred_model <- as.data.table(VSEM(theta[1:11], PAR))[, ..output_vars]
  if("NEE" %in% output_vars) {
    pred_model[, NEE := NEE*1000]
  }
  
  # Evaluate log-likelihood for each output
  model_errs <- pred_model - data_obs[, ..output_vals]
  
  # Sum log-likelihoods, per independence assumption. 
}


# here is the likelihood 
likelihood <- function(par, sum = TRUE){
  # set parameters that are not calibrated on default values 
  x = refPars$best
  x[parSel] = par
  predicted <- VSEM(x[1:11], PAR) # replace here VSEM with your model 
  predicted[,1] = 1000 * predicted[,1] # this is just rescaling
  diff <- c(predicted[,1:4] - obs[,1:4]) # difference betweeno observed and predicted
  # univariate normal likelihood. Note that there is a parameter involved here that is fit
  llValues <- dnorm(diff, sd = x[12], log = TRUE)  
  if (sum == FALSE) return(llValues)
  else return(sum(llValues))
}


#
# Plots of synthetic data
# 

# Plot PAR, the forcing term in the model
plot(PAR, main = "PAR (driving the model)", xlab = "Day", ylab = "PAR")

# Ground-truth NEE
plot(data_ref[,NEE], main = "Ground Truth NEE", xlab = "Day", ylab = "NEE", col = "red")

# Normalize NEE and PAR and plot them together
# Note: NEE > 0 implies positive flux to atmosphere
NEE_ref_norm <- data_ref[, (NEE - mean(NEE))/sd(NEE)]
PAR_norm <- (PAR - mean(PAR)) / sd(PAR)
plot(NEE_ref_norm, main = "Ground Truth NEE vs. PAR (Normalized)", 
     xlab = "Day", ylab = "Z-scaore NEE/PAR", col = "red")
points(PAR_norm, col = "blue")

legend(x = "bottomright",     
       legend = c("NEE", "PAR"), 
       pch = c(1, 1),
       col = c("red", "blue"))   

# Plot the three ground-truth carbon pools 
plot(data_ref[,Cs], main = "Ground Truth C pools", xlab = "Day", type = "l", lwd = 3,
     ylab = "kg C/m^2", ylim = data_ref[, range(Cs, CR, Cv)], col = "black")
lines(data_ref[,Cv], col = "green")
lines(data_ref[,CR], col = "brown", lty = 3, lwd = 3)
legend(x = "right",     
       legend = c("SOM", "Above-Ground", "Below-Ground"), 
       lty = c(1, 1, 3),
       lwd = 3,
       col = c("black", "green", "brown"))  

# De-mean Carbon pool time series and plot 
# Peaks in Cv (above) and CR (below) correspond to peaks in NEE and dips in PAR
plot(data_ref[,Cs-mean(Cs)], main = "Ground Truth C pools (Normalized)", xlab = "Day", type = "l", lwd = 3,
     ylab = "kg C/m^2", col = "black")
lines(data_ref[,Cv-mean(Cv)], col = "green")
lines(data_ref[,CR-mean(CR)], col = "brown", lty = 3, lwd = 3)
legend(x = "bottomright",     
       legend = c("SOM", "Above-Ground", "Below-Ground"), 
       lty = c(1, 1, 3),
       lwd = 3,
       col = c("black", "green", "brown"))  

# Plot above-ground vegetation by itself
plot(data_ref[,Cv], main = "Ground Truth Above-Ground Carbon", xlab = "Day", type = "l", lwd = 3,
     ylab = "kg C/m^2", col = "green")

# BayesianTools has a nice function to plot observed vs. ground-truth time series
plotTimeSeries(observed = data_obs[[1]], 
               predicted = data_ref[[1]], main = colnames(data_ref)[1])
plotTimeSeries(observed = data_obs[[2]], 
               predicted = data_ref[[2]], main = colnames(data_ref)[2])
plotTimeSeries(observed = data_obs[[3]], 
               predicted = data_ref[[3]], main = colnames(data_ref)[3])
plotTimeSeries(observed = data_obs[[4]], 
               predicted = data_ref[[4]], main = colnames(data_ref)[4])







