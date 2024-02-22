## Setting up paths and global variables
## -----------------------------------------------------------------------------


# Cores for parallel::
ncore <- 8


# Top level of plot directory
plot_dir <- "/home/amorin/Plots/TR_singlecell/"


# Location of (most) single cell data sets + metadata
sc_dir <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datasets"


# Location of aggregate coexpression matrices
amat_dir <- "/space/scratch/amorin/TR_singlecell/"


# DIOPT raw input and 1:1 orthologous protein coding genes
diopt_path <- "/space/grp/DIOPT/DIOPTvs8_export_Sanja Rogic.txt" 
pc_ortho_path <- "/space/grp/amorin/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv"


# ENSEMBL and Refseq Select protein coding paths
ref_hg_path <- "/space/grp/amorin/Metadata/refseq_select_hg38.tsv"
ref_mm_path <- "/space/grp/amorin/Metadata/refseq_select_mm10.tsv"
ens_hg_path <- "/space/grp/amorin/Metadata/ensembl_human_protein_coding_105.tsv"
ens_mm_path <- "/space/grp/amorin/Metadata/ensembl_mouse_protein_coding_105.tsv"


# TF gene paths
tfs_hg_path <- "/space/grp/amorin/Metadata/TFs_human.tsv"
tfs_mm_path <- "/space/grp/amorin/Metadata/TFs_mouse.tsv"


# Ortho L/S ribo genes
ribo_path <- "/space/grp/amorin/Metadata/L_S_ribosomal_ortho_genes.tsv"


# TR-target rankings from genomics evidence
evidence_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"


# Gsheets ID for TRsc datasets, local copy of raw meta, and final
gsheets_id <- "1SQx_rFdBNBwOEdChaHkbQY1vStcbuA3fF0CZk6YECEc"
gsheets_meta_raw_path <- "/space/grp/amorin/Metadata/gsheets_single_cell_dataset_meta_raw.tsv"
sc_meta_path <- "/space/grp/amorin/Metadata/single_cell_dataset_meta.tsv"


# List of cell types per dataset
celltype_list_path <- "/space/scratch/amorin/R_objects/celltype_list.RDS"


# Binary matrices tracking if a gene was measured in at least once cell type
msr_mat_hg_path <- "/space/scratch/amorin/R_objects/binary_measurement_matrix_hg.RDS"
msr_mat_mm_path <- "/space/scratch/amorin/R_objects/binary_measurement_matrix_mm.RDS"


# List of the most correlated gene pair per experiment
avg_coexpr_hg_path <- "/space/scratch/amorin/R_objects/avg_coexpr_hg.RDS"
avg_coexpr_mm_path <- "/space/scratch/amorin/R_objects/avg_coexpr_mm.RDS"


# List of TF and L/S ribo aggregate matrices
agg_tf_hg_path <- "/space/scratch/amorin/R_objects/agg_mat_TF_list_hg.RDS"
agg_tf_mm_path <- "/space/scratch/amorin/R_objects/agg_mat_TF_list_mm.RDS"
agg_ribo_hg_path <- "/space/scratch/amorin/R_objects/agg_mat_ribo_list_hg.RDS"
agg_ribo_mm_path <- "/space/scratch/amorin/R_objects/agg_mat_ribo_list_mm.RDS"


# List of TF and ribo vector similarities
sim_tf_hg_path <- "/space/scratch/amorin/R_objects/similarity_TF_hg.RDS"
sim_tf_mm_path <- "/space/scratch/amorin/R_objects/similarity_TF_mm.RDS"
sim_ribo_hg_path <- "/space/scratch/amorin/R_objects/similarity_ribo_hg.RDS"
sim_ribo_mm_path <- "/space/scratch/amorin/R_objects/similarity_ribo_mm.RDS"


# List of null sampled gene vector topk similarities
null_topk_hg_path <- "/space/scratch/amorin/R_objects/sampled_null_topk_intersect_hg.RDS"
null_topk_mm_path <- "/space/scratch/amorin/R_objects/sampled_null_topk_intersect_mm.RDS"


# List of summarized aggregate coexpression rankings
rank_tf_hg_path <- "/space/scratch/amorin/R_objects/ranking_agg_TF_hg.RDS"
rank_tf_mm_path <- "/space/scratch/amorin/R_objects/ranking_agg_TF_mm.RDS"
rank_ribo_hg_path <- "/space/scratch/amorin/R_objects/ranking_agg_ribo_hg.RDS"
rank_ribo_mm_path <- "/space/scratch/amorin/R_objects/ranking_agg_ribo_mm.RDS"


# Curated targets from on going curation and from Eric's 2021 paper
gsheets_curated <- "1PB2P-9Xk2zV0RSZnkY5BdnV6E4KkDpEKvo68Sw_Rnx8"
pavlab_curation_path <- "/home/amorin/Data/Metadata/pavlab_curation_sep2023.tsv"
chu2021_records_path <- "/home/amorin/Data/Metadata/Chu2021_records_DTRI.tsv"
chu2021_all_path <- "/home/amorin/Data/Metadata/Chu2021_all_DTRI.tsv"


# Output of formatted curated targets
curated_all_path <- "/home/amorin/Data/Metadata/Curated_targets_all_Sept2023.tsv"


# Unibind summarized ChIP-seq data
bind_dat_path <- "/space/scratch/amorin/R_objects/processed_unibind_data.RDS"
bind_summary_path <- "/space/scratch/amorin/R_objects/unibind_bindscore_summary.RDS"
bind_model_path <- "/space/scratch/amorin/R_objects/unibind_bindscore_modelfit.RDS"


# List of ROCR performances
coexpr_auc_hg_path <- "/space/scratch/amorin/R_objects/coexpr_recover_curated_hg.RDS"
coexpr_auc_mm_path <- "/space/scratch/amorin/R_objects/coexpr_recover_curated_mm.RDS"
unibind_auc_hg_path <- "/space/scratch/amorin/R_objects/unibind_recover_curated_hg.RDS"
unibind_auc_mm_path <- "/space/scratch/amorin/R_objects/unibind_recover_curated_mm.RDS"
avg_vs_ind_auc_hg_path <- "/space/scratch/amorin/R_objects/individual_vs_average_recover_curated_hg.RDS"
avg_vs_ind_auc_mm_path <- "/space/scratch/amorin/R_objects/individual_vs_average_recover_curated_mm.RDS"
rev_coexpr_auc_hg_path <- "/space/scratch/amorin/R_objects/reverse_coexpr_recover_curated_hg.RDS"
rev_coexpr_auc_mm_path <- "/space/scratch/amorin/R_objects/reverse_coexpr_recover_curated_mm.RDS"

