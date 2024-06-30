"Usage:
  test_docopt.r [options]
  test_docopt.r (-h | --help)

Options:
  -h --help                                 Show this screen.
  --dim_par=<dim_par>                       Dimension of parameter space. 
  --dim_output=<dim_output>                 Dimension of output space. 
  --N_design=<N_design>                     Number of design points.
  --design_method=<design_method>           Algorithm to generate design points. 
  --N_design_test=<N_design_test>           Number of validation points for emulators. 
  --design_method_test=<design_method_test> Algorithm to generate validation points. 
  --mcmc_tags=<mcmc_tag,mcmc_tag>           GP-approx MCMC algorithms. 
  --N_mcmc=<N_mcmc>                         Number of MCMC iterations. 
" -> doc

# print(.libPaths())
# renv::load("/projectnb/dietzelab/arober/gp-calibration")
# print(.libPaths())
# 
# .libPaths(new="/usr3/graduate/arober/.cache/R/renv/sandbox/R-4.3/x86_64-pc-linux-gnu/b4907fb1")
# print(.libPaths())

.libPaths()

library(docopt)
arguments <- docopt(doc)
print(arguments)

unlist(strsplit(arguments$mcmc_tags, ","))
as.integer(arguments$N_design)











