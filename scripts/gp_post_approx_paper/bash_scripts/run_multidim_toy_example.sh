#!/bin/bash -l
#$ -N multidim_toy_example
#$ -P dietzelab                 # Specify project. 
#$ -l buyin                     # Request buyin node. 
#$ -j y                         # Merge the error and output streams into a single file.
#$ -m beas                      # Email when job begins/ends/aborted/suspended

# NOTE: 
# The R library path is not set explicitly in this file. Therefore, if using an 
# R project, one should ensure that renv::load() is called within the R script 
# being run so that the correct library path is set. 

# TODO: need to redirect the output/error files to the run directory. 

# Settings to pass to Rscript. 
RUN_TAG="linGauss"
DIM_PAR=11
DIM_OUTPUT=10
N_DESIGN=100
N_DESIGN_TEST=600
DESIGN_METHOD="LHS"
DESIGN_METHOD_TEST="LHS"
N_MCMC=50000
MCMC_TAGS="gp-mean,gp-marg,mcwmh-joint,mcwmh-ind"

# Create output directory. Creation of this directory is handled (if necessary)
# within the R script. 
SIM_RUN_ID="d${DIM_PAR}_p${DIM_OUTPUT}_N${N_DESIGN}_${DESIGN_METHOD}"
LOG_FILENAME="${JOB_NAME}.o${JOB_ID}"
OUTPUT_DIR="/projectnb/dietzelab/arober/gp-calibration/output/gp_post_approx_paper/${RUN_TAG}/${SIM_RUN_ID}"

# Load modules, including specification of the R version. 
module load gcc/8.3.0
module load R/4.3.1

# Required for the `Rscript` command to be found. 
export PATH="$PATH:/share/pkg.8/r/4.3.1/install/lib64/R"

# Run Rscript. 
Rscript ../run_multidim_toy_example.r ${SIM_RUN_ID} ${OUTPUT_DIR} \
        --dim_par=${DIM_PAR} --dim_output=${DIM_OUTPUT} \
        --N_design=${N_DESIGN} --N_design_test=${N_DESIGN_TEST} \
        --design_method=${DESIGN_METHOD} \
        --design_method_test=${DESIGN_METHOD_TEST} \
        --N_mcmc=${N_MCMC} --mcmc_tags=${MCMC_TAGS}

# Move output file to the output directory (if it exists). 
if [ -d ${OUTPUT_DIR} ]; then
  mv "${LOG_FILENAME}" "${OUTPUT_DIR}/${LOG_FILENAME}"
fi



