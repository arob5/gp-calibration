#!/bin/bash -l
#$ -N run_rmarkdown
#$ -P dietzelab                 # Specify project. 
#$ -l buyin                     # Request buyin node. 
#$ -j y                         # Merge the error and output streams into a single file.
#$ -m beas                      # Email when job begins/ends/aborted/suspended

# NOTE: 
# The R library path is not set explicitly in this file. Therefore, if using an 
# R project, one should ensure that renv::load() is called within the R script 
# being run so that the correct library path is set. 

# Settings: used to identify the correct filepaths. 
# `SCRIPT_NAME` should not include the file extension. This is used both 
# to identify the r script to run (with assumed ".r" extension), as well 
# as the directory for the RMarkdown output.
# At present, this is set up to only allow passing a single argument to the 
# Rmarkdown file; specifically, the experiment tag.
EXPERIMENT_TAG="vsem"
SCRIPT_NAME="inv_prob_setup"

echo "${EXPERIMENT_TAG}"
echo "${SCRIPT_NAME}"

# For logging.
LOG_FILENAME="${JOB_NAME}.o${JOB_ID}"

# Directory and filepath to which RMarkdown output will be directed.
OUT_DIR="/projectnb/dietzelab/arober/gp-calibration/output/gp_inv_prob/${EXPERIMENT_TAG}/${SCRIPT_NAME}"
OUT_PATH="${OUT_DIR}/${SCRIPT_NAME}.html"

# Path to RMarkdown file to execute.
RMD_IN_PATH="/projectnb/dietzelab/arober/gp-calibration/scripts/gp_post_approx_paper/sim_study/notebooks/${SCRIPT_NAME}.Rmd"

# Define R code used to run Rmarkdown.
R_RUN_CMD="rmarkdown::render(input='${RMD_IN_PATH}', output_file='${OUT_PATH}', output_format='html_document', param=list(args='${EXPERIMENT_TAG}'))"

# Load modules, including specification of the R version. 
module load gcc/8.3.0
module load R/4.3.1

# Required for the `Rscript` command to be found. 
export PATH="$PATH:/share/pkg.8/r/4.3.1/install/lib64/R"

# Render Rmarkdown.
echo "${R_RUN_CMD}"
Rscript -e "${R_RUN_CMD}"

# Move output file to the output directory (if it exists). 
if [ -d ${OUT_DIR} ]; then
mv "${LOG_FILENAME}" "${OUT_DIR}/${LOG_FILENAME}"
fi
