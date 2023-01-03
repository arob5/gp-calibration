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
# Setup: settings control synthetic data generation
#

# Number days in time series 
N_days <- 1000

# Standard deviation of observation error (for now assuming Gaussian noise). Also 
# store minimum and lower bounds on the standard deviation to align with the 
# reference parameter formatting in VSEMgetDefaults()
sig_eps <- 2.0
sig_eps_rng <- c(0.1, 4.0)

#
# Synthetic data generation
#

# Create time series of Phosynthetically Active Radiation (PAR) which is the forcing variable
# for the VSEM model. 
PAR <- VSEMcreatePAR(seq_len(N_days))

# The "best" column of the reference parameters are used to generate the "true" 
# output data. We will add noise to this data to simulate the field observations.
ref_pars <- VSEMgetDefaults()
data_ref <- as.data.table(VSEM(refPars$best, PAR))

# Add observational noise, assumes same observational noise on each output so re-scale 
# NEE to be on similar scale as the C pools so the noise level makes sense.
ref_pars["sig_eps",] <- c(sig_eps, sig_eps_rng)
data_ref[, NEE := NEE*1000]
data_obs <- data_ref + rnorm(nrow(data_ref), sd = ref_pars["sig_eps", "best"])


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







