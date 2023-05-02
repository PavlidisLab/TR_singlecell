## Setting up paths and global variables
## -----------------------------------------------------------------------------

# library(future)

# Prepare future:: parallel and increase memory max size for Seurat
# plan("multisession", workers = 8)
# options(future.globals.maxSize = 2097152000)

# Cores for parallel::
ncore <- 8

plot_dir <- "/home/amorin/Plots/TR_singlecell/"


# Location of aggregate coexpression matrices
amat_dir <- "~/scratch/R_objects/"


# 1:1 orthologous protein coding genes
pc_ortho_path <- "/space/grp/amorin/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv"

# TR-target rankings from genomics evidence
evidence_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"


# Human Protein Atlas

expr_mat_l_path <- "/home/amorin/scratch/R_objects/HPA_expression_mat_list.RDS"
cor_mat_l_path <- "/home/amorin/scratch/R_objects/HPA_cor_mat_list.RDS"