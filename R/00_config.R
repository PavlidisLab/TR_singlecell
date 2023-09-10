## Setting up paths and global variables
## TODO: pcoding download and common filter needs to be formalized
## TODO: mito gene was manually downloaded from biomart... doc this https://www.biostars.org/p/310641/
## TODO: ribo human ribo S https://www.genenames.org/data/genegroup/#!/group/728  https://www.genenames.org/cgi-bin/genegroup/download?id=728&type=node
## TODO: ribo human ribo L https://www.genenames.org/data/genegroup/#!/group/729  https://www.genenames.org/cgi-bin/genegroup/download?id=729&type=node
## TODO: TF data acquisition
## -----------------------------------------------------------------------------


# Cores for parallel::
ncore <- 8


# Top level of plot directory
plot_dir <- "/home/amorin/Plots/TR_singlecell/"


# Location of (most) single cell data sets + metadata
sc_dir <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datasets"


# Location of aggregate coexpression matrices
amat_dir <- "/space/scratch/amorin/TR_singlecell/"


# 1:1 orthologous protein coding genes
pc_ortho_path <- "/space/grp/amorin/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv"


# ENSEMBL and Refseq Select protein coding paths
ref_hg_path <- "/space/grp/amorin/Metadata/refseq_select_hg38.tsv"
ref_mm_path <- "/space/grp/amorin/Metadata/refseq_select_mm10.tsv"
ens_hg_path <- "/space/grp/amorin/Metadata/ensembl_human_protein_coding_105.tsv"
ens_mm_path <- "/space/grp/amorin/Metadata/ensembl_mouse_protein_coding_105.tsv"


# TF gene paths
tfs_hg_path <- "/space/grp/amorin/Metadata/human_tfs.tsv"
tfs_mm_path <- "/space/grp/amorin/Metadata/mouse_tfs.tsv"


# Human L/S ribo genes
sribo_hg_path <- "/space/grp/amorin/Metadata/HGNC_human_Sribosomal_genes.csv"
lribo_hg_path <- "/space/grp/amorin/Metadata/HGNC_human_Lribosomal_genes.csv"


# TR-target rankings from genomics evidence
evidence_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"


# Gsheets ID for TRsc datasets
gsheets_id <- "1SQx_rFdBNBwOEdChaHkbQY1vStcbuA3fF0CZk6YECEc"


# Local copy of meta
sc_meta_path <- "/space/grp/amorin/Metadata/single_cell_dataset_meta.tsv"


# List of cell types per dataset
celltype_list_path <- "/space/scratch/amorin/TR_singlecell/celltype_list.RDS"


# Binary matrices tracking if a gene was measured in at least once cell type
msr_mat_hg_path <- "/space/scratch/amorin/TR_singlecell/binary_measurement_matrix_human.RDS"
msr_mat_mm_path <- "/space/scratch/amorin/TR_singlecell/binary_measurement_matrix_mouse.RDS"


# List of the most correlated gene pair per experiment
avg_coexpr_hg_path <- "/space/scratch/amorin/R_objects/avg_coexpr_hg.RDS"
avg_coexpr_mm_path <- "/space/scratch/amorin/R_objects/avg_coexpr_mm.RDS"


# List of TF and L/S ribo aggregate matrices
agg_tf_hg_path <- "/space/scratch/amorin/R_objects/TF_agg_mat_list_human.RDS"
agg_tf_mm_path <- "/space/scratch/amorin/R_objects/TF_agg_mat_list_mouse.RDS"
agg_ribo_hg_path <- "/space/scratch/amorin/R_objects/ribo_agg_mat_list_human.RDS"
agg_ribo_mm_path <- "/space/scratch/amorin/R_objects/ribo_agg_mat_list_mouse.RDS"


# List of TF and ribo vector similarities
sim_tf_hg_path <- "/space/scratch/amorin/R_objects/similarity_TF_human.RDS"
sim_tf_mm_path <- "/space/scratch/amorin/R_objects/similarity_TF_mouse.RDS"
sim_ribo_hg_path <- "/space/scratch/amorin/R_objects/similarity_ribo_human.RDS"
sim_ribo_mm_path <- "/space/scratch/amorin/R_objects/similarity_ribo_mouse.RDS"


# List of null sampled gene vector topk similarities
null_topk_hg_path <- "/space/scratch/amorin/R_objects/sampled_null_topk_intersect_human.RDS"
null_topk_mm_path <- "/space/scratch/amorin/R_objects/sampled_null_topk_intersect_mouse.RDS"


# List of summarized TF rankings
rank_tf_hg_path <- "/space/scratch/amorin/R_objects/ranking_TF_human.RDS"
rank_tf_mm_path <- "/space/scratch/amorin/R_objects/ranking_TF_mouse.RDS"


# Curated targets from on going curation and from Eric's 2021 paper
gsheets_curated <- "1PB2P-9Xk2zV0RSZnkY5BdnV6E4KkDpEKvo68Sw_Rnx8"
pavlab_curation_path <- "/home/amorin/Data/Metadata/pavlab_curation_sep2023.tsv"
chu2021_records_path <- "/home/amorin/Data/Metadata/Chu2021_records_DTRI.tsv"
chu2021_all_path <- "/home/amorin/Data/Metadata/Chu2021_all_DTRI.tsv"

# Output of formatted curated targets
curated_all_path <- "/home/amorin/Data/Metadata/Curated_targets_all_Sept2023.tsv"


# List of ROCR performances
coexpr_recover_curated_hg_path <- "/space/scratch/amorin/R_objects/coexpr_recover_curated_hg.RDS"
coexpr_recover_curated_mm_path <- "/space/scratch/amorin/R_objects/coexpr_recover_curated_mm.RDS"
unibind_recover_curated_hg_path <- "/space/scratch/amorin/R_objects/unibind_recover_curated_hg.RDS"
unibind_recover_curated_mm_path <- "/space/scratch/amorin/R_objects/unibind_recover_curated_mm.RDS"




# Human Protein Atlas
# expr_mat_l_path <- "/home/amorin/scratch/R_objects/HPA_expression_mat_list.RDS"
# cor_mat_l_path <- "/home/amorin/scratch/R_objects/HPA_cor_mat_list.RDS"