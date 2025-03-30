## Setting up paths and global variables
## -----------------------------------------------------------------------------


# Cores for parallel::
ncore <- 8


# Main directories
# ---

# Meta table dir (protein coding, TFs, ribo, curated targets, final experiments)
meta_dir <- "/space/grp/amorin/Metadata/"


# Top level of plot directory
plot_dir <- "/home/amorin/Plots/TR_singlecell/"


# Location of unprocessed single cell data sets + metadata
sc_dir <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datasets"


# Location of processed single cell data and aggregate coexpression matrices
amat_dir <- "/space/scratch/amorin/TR_singlecell/"


# Various outputs like .RDS files
output_dir <- "/space/scratch/amorin/TRsc_output/"

# ---



# DIOPT raw input table for orthology mapping and 1:1 human to mouse mapping
# NOTE: The DIOPT input is a hard path from Sanja
# https://www.flyrnai.org/cgi-bin/DRSC_orthologs.pl
diopt_path <- "/space/grp/DIOPT/DIOPTvs9_export_human2other_Sanja_20220901.txt" 
pc_ortho_path <- paste0(meta_dir, "hg_mm_1to1_ortho_genes_DIOPT_V9.tsv")


# ENSEMBL and Refseq Select protein coding paths
ref_hg_path <- paste0(meta_dir, "refseq_select_hg38_jan2024.tsv")
ref_mm_path <- paste0(meta_dir, "refseq_select_mm10_jan2024.tsv")
ens_hg_path <- paste0(meta_dir, "ensembl_human_protein_coding_105.tsv")
ens_mm_path <- paste0(meta_dir, "ensembl_mouse_protein_coding_105.tsv")


# AnimalTFDB paths for TF genes
tfs_hg_path <- paste0(meta_dir, "AnimalTFDB_human_V4.tsv")
tfs_mm_path <- paste0(meta_dir, "AnimalTFDB_mouse_V4.tsv")


# Ortho L/S ribo genes
ribo_path <- paste0(meta_dir, "L_S_ribosomal_ortho_genes.tsv")


# Gsheets ID for TRsc datasets, local copy of raw meta, and final
gsheets_id <- "1SQx_rFdBNBwOEdChaHkbQY1vStcbuA3fF0CZk6YECEc"
gsheets_meta_raw_path <- paste0(meta_dir, "gsheets_single_cell_dataset_meta_raw.tsv")
sc_meta_path <- paste0(meta_dir, "single_cell_dataset_meta.tsv")


# List of author-annotated cell types per dataset
celltype_list_path <- paste0(output_dir, "celltype_list.RDS")


# Binary matrices tracking if genes were measured in experiments
msr_mat_hg_path <- paste0(output_dir, "binary_measurement_matrix_hg.RDS")
msr_mat_mm_path <- paste0(output_dir, "binary_measurement_matrix_mm.RDS")
comsr_mat_hg_path <- paste0(output_dir, "comeasurement_matrix_hg.RDS")
comsr_mat_mm_path <- paste0(output_dir, "comeasurement_matrix_mm.RDS")


# List of the most correlated gene pairs across experiments
avg_coexpr_hg_path <- paste0(output_dir, "avg_coexpr_hg.RDS")
avg_coexpr_mm_path <- paste0(output_dir, "avg_coexpr_mm.RDS")


# List of TF and L/S ribo aggregate coexpression matrices
agg_tf_hg_path <- paste0(output_dir, "agg_mat_TF_list_hg.RDS")
agg_tf_mm_path <- paste0(output_dir, "agg_mat_TF_list_mm.RDS")
agg_ribo_hg_path <- paste0(output_dir, "agg_mat_ribo_list_hg.RDS")
agg_ribo_mm_path <- paste0(output_dir, "agg_mat_ribo_list_mm.RDS")


# List of summarized aggregate coexpression rankings
rank_tf_hg_path <- paste0(output_dir, "ranking_agg_coexpr_TF_hg.RDS")
rank_tf_mm_path <- paste0(output_dir, "ranking_agg_coexpr_TF_mm.RDS")
rank_tf_ortho_path <- paste0(output_dir, "ranking_agg_coexpr_TF_ortho.RDS")
rank_ribo_hg_path <- paste0(output_dir, "ranking_agg_coexpr_ribo_hg.RDS")
rank_ribo_mm_path <- paste0(output_dir, "ranking_agg_coexpr_ribo_mm.RDS")


# List of integrated rankings
rank_int_hg_path <- paste0(output_dir, "ranking_agg_integrated_TF_hg.RDS")
rank_int_mm_path <- paste0(output_dir, "ranking_agg_integrated_TF_mm.RDS")
rank_int_ortho_path <- paste0(output_dir, "ranking_agg_integrated_TF_ortho.RDS")


# Curated targets from on going curation and from Eric's 2021 paper
gsheets_curated <- "1rKu0inmJt67q50cz18ccmMARXd9E_LQNElsWW7tn5j8"
pavlab_curation_path <- paste0(meta_dir, "pavlab_curation_sep2023.tsv")
chu2021_records_path <- paste0(meta_dir, "Chu2021_records_DTRI.tsv")
chu2021_all_path <- paste0(meta_dir, "Chu2021_all_DTRI.tsv")
curated_all_path <- paste0(meta_dir, "Curated_targets_all_Sept2023.tsv")


# Unibind summarized ChIP-seq data
# This lives on the Pavlab servers, and was generated via: 
# https://github.com/PavlidisLab/Unibind_analysis/
bind_summary_path <- "/space/scratch/amorin/Unibind/unibind_Permissive_bindscore_summary.RDS"


# List of literature curation benchmark performances
coexpr_auc_hg_path <- paste0(output_dir, "coexpr_recover_curated_hg.RDS")
coexpr_auc_mm_path <- paste0(output_dir, "coexpr_recover_curated_mm.RDS")
unibind_auc_hg_path <- paste0(output_dir, "unibind_recover_curated_hg.RDS")
unibind_auc_mm_path <- paste0(output_dir, "unibind_recover_curated_mm.RDS")
avg_vs_ind_auc_hg_path <- paste0(output_dir, "individual_vs_average_recover_curated_hg.RDS")
avg_vs_ind_auc_mm_path <- paste0(output_dir, "individual_vs_average_recover_curated_mm.RDS")
rev_coexpr_auc_hg_path <- paste0(output_dir, "reverse_coexpr_recover_curated_hg.RDS")
rev_coexpr_auc_mm_path <- paste0(output_dir, "reverse_coexpr_recover_curated_mm.RDS")
int_auc_hg_path <- paste0(output_dir, "integrated_recover_curated_hg.RDS")
int_auc_mm_path <- paste0(output_dir, "integrated_recover_curated_mm.RDS")


# Summary tables of literature curation benchmark performance
auc_table_hg_path <- paste0(output_dir, "curation_benchmark_summary_hg.tsv")
auc_table_mm_path <- paste0(output_dir, "curation_benchmark_summary_mm.tsv")


# Tiered evidence based on K cutoff
tiered_evidence_path <- paste0(output_dir, "tiered_evidence_list.RDS")
tiered_evidence_flat_path <- paste0(output_dir, "tiered_evidence_flat_list.RDS")


# GO annotations
go_date <- "2024-09-08"
go_path <- paste0(output_dir, "go_terms_", go_date, ".xml")
anno_hg_path <- paste0(output_dir, "gemma_generic_human_anno")
anno_mm_path <- paste0(output_dir, "gemma_generic_mouse_anno")

# GO output results
erminer_coexpr_hg_path <- paste0(output_dir, "coexpr_erminer_pr_hg.RDS")
erminer_coexpr_mm_path <- paste0(output_dir, "coexpr_erminer_pr_mm.RDS")
erminer_coexpr_ortho_path <- paste0(output_dir, "coexpr_erminer_pr_ortho.RDS")
