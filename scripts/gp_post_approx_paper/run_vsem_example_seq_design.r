#
# run_vsem_example_seq_design.r
# Sequential acqusition of design points following the definition of an initial 
# surrogate model in `run_vsem_example.r`. The main function defined below 
# takes two inputs: a fit `llikEmulator` object and (optionally) samples from
# a current approximate posterior distribution. These are then used to acquire
# a new batch of design points. 
#
# Andrew Roberts
#

