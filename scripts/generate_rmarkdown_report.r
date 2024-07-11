#
# generate_rmarkdown_report.r
# A convenience script to render an RMarkdown document and save it to a 
# specific file path. The main purpose of this script is to name the file 
# something other than the name of the RMarkdown file itself, which is 
# not directly possible when using the "knit" button in RStudio. 
#
# Andrew Roberts
#

# TODO: pass settings to Rmarkdown file 
# https://stackoverflow.com/questions/18929782/how-can-i-pass-variables-into-an-r-markdown-rmd-file

#
# Settings 
#

base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")
scripts_dir <- file.path(base_dir, "scripts")
output_dir <- file.path(base_dir, "output")
tag <- "linGauss"
run_id <- "d11_p10_N100_LHS"
run_dir_relative <- file.path("gp_post_approx_paper", tag, run_id)

rmd_path <- file.path(scripts_dir, "gp_post_approx_paper", "multidim_toy_examples.Rmd")
save_dir <- file.path(output_dir, run_dir_relative, "rmarkdown")
                      
base_filename <- run_id
file_extension <- ".html"

#
# Set up directory, avoid overwriting. 
#

# If output filename already already exists, append timestamp to filename to 
# avoid overwriting. 
if(dir.exists(save_dir)) {
  
  filename <- paste0(base_filename, file_extension)
  if(file.exists(file.path(save_dir, filename))) {
    timestamp <- as.character(Sys.time())
    base_filename <- paste(base_filename, timestamp, sep="_")
  }
  
} else {
  dir.create(save_dir, recursive=TRUE)
}
filename <- paste0(base_filename, file_extension)


#
# Render and save file. 
#

rmarkdown::render(input=rmd_path, output_file=filename, output_dir=save_dir)
                  



