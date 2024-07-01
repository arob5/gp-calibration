#!/bin/bash -l
#$ -N multidim_toy_example
#$ -P dietzelab                 # Specify project. 
#$ -l buyin                     # Request buyin node. 
#$ -j y                         # Merge the error and output streams into a single file.
#$ -m beas                      # Email when job begins/ends/aborted/suspended

# TODO: need to redirect the output/error files to the run directory. 

# Load modules. 
module load gcc/8.3.0
module load R/4.3.1

# Ensure the R project library for project `gp-calibration` is being used. 
# export R_LIBS="/projectnb/dietzelab/arober/gp-calibration/renv/library/R-4.3/x86_64-pc-linux-gnu"

# Required for the `Rscript` command to be found. 
export PATH="$PATH:/share/pkg.8/r/4.3.1/install/lib64/R"

Rscript ../run_multidim_toy_example.r linGauss --dim_par=3 --dim_output=10 --N_design=30 --N_design_test=400

