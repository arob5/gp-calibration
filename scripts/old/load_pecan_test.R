# --------------------------------------------------------------------------------
# load_pecan_test.R
#
# Loads data file saved from 'pda.emulator.R' and transforms to be run on my test 
# algorithms. 
#
# Andrew Roberts
# Working Directory: /projectnb/dietzelab/arober/pecan_personal
# --------------------------------------------------------------------------------

library(PEcAn.assim.batch)
library(PEcAn.utils)
library(PEcAn.emulator)

setwd("/projectnb/dietzelab/arober/pecan_personal")
base_dir <- getwd()
rdata_path <- file.path(base_dir, "..", "test_Rdata", "history.pda1000031535.Rdata")

# Load R data file and check the current step
rdata <- load(rdata_path)
print(current.step)

#
# Different datasets to be calibrated against
#

# Number of datasets
cat("Number of datasets: ", length(inputs))

# First dataset. What is the UST column? Looks like time stamp is every 
# 30 minutes. 
inputs1.data <- inputs[[1]]$data
dim(inputs1.data)
head(inputs1.data)

# Variables associated with this dataset are "FC" and "UST". Presumably these 
# are two output variables to be calibrated against. 
length(inputs[[1]]$variable.name)
inputs[[1]]$variable.name[1]$variable.name$variable.drv
inputs[[1]]$variable.name[2]$variable.name$variable.drv

# Only a portion of the "FC" column is not NA, and same for the UST column. 
sum(!is.na(inputs1.data$FC))
sum(!is.na(inputs1.data$UST))
head(inputs1.data[!is.na(inputs1.data$FC),])
head(inputs1.data[!is.na(inputs1.data$UST),])

# Length of "obs" equal to number of rows of inputs[[1]]$data
length(inputs[[1]]$obs)
sum(!is.na(inputs[[1]]$obs))
inputs[[1]]$obs[!is.na(inputs1.data$FC)][1:10]

# Number observations when both variables are not NA
sum(!is.na(inputs1.data$FC) & !is.na(inputs1.data$UST))

# "obs" is equal to the column FC of the data
all.equal(inputs1.data$FC, inputs[[1]]$obs)

# What is inputs[[1]]$par?
inputs[[1]]$par

# inputs[[1]]$n is just the number of observations of FC to be used in the dataset
# (number non-NA values of obs)
inputs[[1]]$n == sum(!is.na(inputs1.data$FC))

# inputs[[1]]$n_eff is probably some sort of effective sample size calculation 
# taking into account the correlation in the FC time series; need to check how 
# they calculate this. 
inputs[[1]]$n_eff

# Q: why do "obs" and "n" only correspond to FC, and there are no analogous 
# variables for UST? 
# Need to figure out what UST is. UST is also listed as the second variable 
# under inputs[[2]], so the variable listed first seems to be the main output 
# of interest. This output is "LE" for inputs[[2]].
inputs[[2]]$variable.name

# Why is the first column of inputs[[2]]$data called "Qle" and not "LE"?
# This is partially addressed in a comment in pda.get.model.output.R, which calls
# "LE" a data variable name, and "Qle" a model output name. Seems like there are 
# just a few different naming schemes in use. 
head(inputs[[2]]$data)

# Plot some of the FC data (Note: only printing non-NA observations, so messing up the
# time dependence in the plot, so take with grain of salt)
sel.FC.non.NA <- !is.na(inputs1.data$FC)
n.plot <- 60
plot(1:n.plot, inputs1.data$FC[sel.FC.non.NA][1:n.plot])

# Plot some of the UST data 
sel.UST.non.NA <- !is.na(inputs1.data$UST)
plot(1:n.plot, inputs1.data$UST[sel.UST.non.NA][1:n.plot])

# Plot some of the LE data 
sel.Qle.non.NA <- !is.na(inputs[[2]]$data$Qle)
plot(1:n.plot, inputs[[2]]$data$Qle[sel.Qle.non.NA][1:n.plot])


#
# Likelihood associated with different datasets to be calibrated against
#

# "inputs" element of the "settings" object stores information that pertains to 
# the "inputs" object explored above.
names(settings$assim.batch$inputs[[1]])
settings$assim.batch$inputs[[1]]$path # Path to data I guess? 
settings$assim.batch$inputs[[1]]$likelihood # Laplace likelihood associated with this dataset
settings$assim.batch$inputs[[1]]$variable.id




# TODO: plot the data above (obs and data)


#
# Design points/knots
#

# Investigating the form of the "knots" (i.e. design points) which are fed into 
# the computer model. The variable "run.ids" contains one run ID per knot. 
head(run.ids)

# Investigating the form of the output of the computer model, as returned by 
# 'pda.get.model.output()'. In pda.emulator.R the outputs from the computer model
# are saved in a list called 'model.out', which has length equal to the number of 
# design points. 
names(model.out)

# Investigating data used to fit GP 
# TODO: look where variable "x" in pda.emulator.R comes from.






