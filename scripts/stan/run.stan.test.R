# Working directory: /projectnb2/dietzelab/arober/pecan

# ---------------------
# Setup
# ---------------------

library(PEcAn.assim.batch)
library(PEcAn.emulator)

# Set working directory for cluster runs
setwd('/projectnb2/dietzelab/arober/pecan')

# Execution gates
use.gp.approx <- TRUE
use.stan <- FALSE

# File path to .RData file 
Rdata.dir <- file.path('..', 'test_Rdata')
base.Rdata.filename <- 'history.pda1000031535'
pre.mcmc.data.path <- file.path(Rdata.dir, paste0(base.Rdata.filename, '_manualtest.Rdata'))

# Path to outputs
out.dir <- 'test_outputs'

# Paths to helper functions
test.code.dir <- file.path('..', 'test_code')
helper.function.filenames <- c('mcmc.GP.stan.R', 'stan.helper.functions.R', 'pda.emulator.stan.test.R')
for(f in file.path(test.code.dir, helper.function.filenames)) source(f)


# ---------------------
# Load Data
# ---------------------

# Should load .RData file saved from pda.emulator.R with `current.step` equal to "pre-MCMC"
load(pre.mcmc.data.path)


if(use.gp.approx && !use.stan) {
  pda.emulator.run.mcmc(settings, out.dir, Rdata.dir, base.Rdata.filename,  
                        gp, init.list, rng, mix, jmp.list, prior.fn.all,  
                        run.normal, run.round, n.of.obs, llik.fn, 
                        hyper.pars, resume.list)
}







# Test code:

stan.string <- generate.stan.script(prior.table = prior.all[prior.ind.all, ], 
                                    write.to.file = TRUE)
cat(stan.string)

gp.params <- create.gp.params.list(gp[[1]], 'mlegp')
print(gp.params)

