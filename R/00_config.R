## Setting up paths and global variables
## -----------------------------------------------------------------------------

library(future)

# Prepare future:: parallel and increase memory max size for Seurat
plan("multisession", workers = 8)
options(future.globals.maxSize = 2097152000)

# Cores for parallel::
ncore <- 8

plot_dir <- "/home/amorin/Plots/TR_singlecell/"


# Human Protein Atlas

expr_mat_l_path <- "/home/amorin/scratch/R_objects/HPA_expression_mat_list.RDS"
cor_mat_l_path <- "/home/amorin/scratch/R_objects/HPA_cor_mat_list.RDS"