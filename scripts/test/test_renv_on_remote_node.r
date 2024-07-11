#
# test_renv_on_remote_node.r
# Testing loading an renv project on a compute node on the shared computing 
# cluster (SCC). Assuming R version 4.3. 
#
# Andrew Roberts
# 

print("R.version:")
R.version

print(".libPaths() before loading project:")
.libPaths()

# Load project. 
project_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")
print(paste0("`project_dir`: ", project_dir))
print("Calling renv::load() from project `project_dir`")
setwd(project_dir)
renv::load(getwd())

print(".libPaths() after loading project:")
.libPaths()

print("Try loading package: kergp")
library(kergp)




