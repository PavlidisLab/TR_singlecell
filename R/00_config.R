## Setting up paths and global variables
## TODO: pcoding download and common filter needs to be formalized
## TODO: mito gene was manually downloaded from biomart... doc this https://www.biostars.org/p/310641/
## -----------------------------------------------------------------------------


# Cores for parallel::
ncore <- 8

plot_dir <- "/home/amorin/Plots/TR_singlecell/"


# Location of aggregate coexpression matrices
amat_dir <- "~/scratch/R_objects/"


# 1:1 orthologous protein coding genes
pc_ortho_path <- "/space/grp/amorin/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv"

# Pcoding paths
ref_hg_path <- "/space/grp/amorin/Metadata/refseq_select_hg38.tsv"
ref_mm_path <- "/space/grp/amorin/Metadata/refseq_select_mm10.tsv"
ens_hg_path <- "/space/grp/amorin/Metadata/ensembl_human_protein_coding_105.tsv"
ens_mm_path <- "/space/grp/amorin/Metadata/ensembl_mouse_protein_coding_105.tsv"


# TR-target rankings from genomics evidence
evidence_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"


# Human Protein Atlas

expr_mat_l_path <- "/home/amorin/scratch/R_objects/HPA_expression_mat_list.RDS"
cor_mat_l_path <- "/home/amorin/scratch/R_objects/HPA_cor_mat_list.RDS"