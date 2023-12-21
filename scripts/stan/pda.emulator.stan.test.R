# pda.emulator.stan.test.R
# This file contains helper functions used by run.stan.test.R. Some of the 
# functions are just portions of the function pda.emulator() that have been 
# slightly modified to be convenient to run the stan tests. 
#
# Andrew Roberts


pda.emulator.run.mcmc <- function(settings, out.dir, Rdata.dir, base.Rdata.filename,  
                                  gp, init.list, rng, mix, jmp.list, prior.fn.all,  
                                  run.normal, run.round, n.of.obs, llik.fn, 
                                  hyper.pars, resume.list) {
  # Code taken from pda.emulator immediately after step "pre-MCMC". Minor 
  # modifications were made to the original code for the convenience of 
  # running tests.

  # start the clock
  ptm.start <- proc.time()

  # prepare for parallelization
  dcores <- parallel::detectCores() - 1
  ncores <- min(max(dcores, 1), settings$assim.batch$chain)


  logfile_path <- file.path(out.dir, "pda.log")
  PEcAn.logger::logger.setOutputFile(logfile_path)

  # Can't get this to work with outfile
  cl <- parallel::makeCluster(ncores, type="FORK") # , outfile = logfile_path) 

  ## Sample posterior from emulator
  mcmc.out <- parallel::parLapply(cl, 1:settings$assim.batch$chain, function(chain) {
    PEcAn.emulator::mcmc.GP(gp          = gp, ## Emulator(s)
                            x0          = init.list[[chain]],     ## Initial conditions
                            nmcmc       = settings$assim.batch$iter,       ## Number of reps
                            rng         = rng,       ## range
                            format      = "lin",      ## "lin"ear vs "log" of LogLikelihood
                            mix         = mix,     ## Jump "each" dimension independently or update them "joint"ly
                            jmp0        = jmp.list[[chain]],  ## Initial jump size
                            ar.target   = settings$assim.batch$jump$ar.target,   ## Target acceptance rate
                            priors      = prior.fn.all$dprior[prior.ind.all], ## priors
                            settings    = settings,
                            run.block   = (run.normal | run.round),
                            n.of.obs    = n.of.obs,
                            llik.fn     = llik.fn,
                            hyper.pars  = hyper.pars,
                            resume.list = resume.list[[chain]]
  )
})

parallel::stopCluster(cl)

# Stop the clock
ptm.finish <- proc.time() - ptm.start
PEcAn.logger::logger.info(paste0("Emulator MCMC took ", paste0(round(ptm.finish[3])), 
                                 " seconds for ", paste0(settings$assim.batch$iter), " iterations."))

current.step <- "post-MCMC"
save(list = ls(all.names = TRUE),envir=environment(), 
     file=file.path(Rdata.dir, base.Rdata.filename, '_manual_postMCMC.Rdata'))

}


