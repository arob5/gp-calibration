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






