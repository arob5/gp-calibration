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


#
# Synthetic data generation
#

# Create time series of Phosynthetically Active Radiation (PAR) which is the forcing variable
# for the VSEM model. 
PAR <- VSEMcreatePAR(seq_len(N_days))
plot(PAR, main = "PAR (driving the model)", xlab = "Day", ylab = "PAR")

# The "best" column of the reference parameters are used to generate the "true" 
# output data. We will add noise to this data to simulate the field observations.
ref_pars <- VSEMgetDefaults()
ref_data <- as.data.table(VSEM(refPars$best, PAR))

# Ground-truth NEE
plot(ref_data[,"NEE"], main = "Ground Truth NEE", xlab = "Day", ylab = "NEE", col = "red")

# Normalize NEE and PAR and plot them together
# Note: NEE > 0 implies positive flux to atmosphere
NEE_ref_norm <- ref_data[, (NEE - mean(NEE))/sd(NEE)]
PAR_norm <- (PAR - mean(PAR)) / sd(PAR)
plot(NEE_ref_norm, main = "Ground Truth NEE vs. PAR (Normalized)", 
     xlab = "Day", ylab = "Z-scaore NEE/PAR", col = "red")
points(PAR_norm, col = "blue")

legend(x = "bottomright",     
       legend = c("NEE", "PAR"), 
       pch = c(1, 1),
       col = c("red", "blue"))   

# Plot the three ground-truth carbon pools 
plot(ref_data[,Cs], main = "Ground Truth C pools", xlab = "Day", type = "l", lwd = 3,
     ylab = "kg C/m^2", ylim = ref_data[, range(Cs, CR, Cv)], col = "black")
lines(ref_data[,Cv], col = "green")
lines(ref_data[,CR], col = "brown", lty = 3, lwd = 3)
legend(x = "right",     
       legend = c("SOM", "Above-Ground", "Below-Ground"), 
       lty = c(1, 1, 3),
       lwd = 3,
       col = c("black", "green", "brown"))  

# De-mean Carbon pool time series and plot 
# Peaks in Cv (above) and CR (below) correspond to peaks in NEE and dips in PAR
plot(ref_data[,Cs-mean(Cs)], main = "Ground Truth C pools (Normalized)", xlab = "Day", type = "l", lwd = 3,
     ylab = "kg C/m^2", col = "black")
lines(ref_data[,Cv-mean(Cv)], col = "green")
lines(ref_data[,CR-mean(CR)], col = "brown", lty = 3, lwd = 3)
legend(x = "bottomright",     
       legend = c("SOM", "Above-Ground", "Below-Ground"), 
       lty = c(1, 1, 3),
       lwd = 3,
       col = c("black", "green", "brown"))  







