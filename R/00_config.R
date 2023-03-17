## Setting up paths and global variables
## -----------------------------------------------------------------------------

library(future)

# Prepare future:: parallel and increase memory max size for Seurat
plan("multisession", workers = 8)  
options(future.globals.maxSize = 2097152000)

# Cores for parallel::
ncore <- 8

plot_dir <- "/home/amorin/Plots/TR_singlecell/"