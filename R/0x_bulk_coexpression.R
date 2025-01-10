## 
## -----------------------------------------------------------------------------

library(tidyverse)
library(reshape2)
library(WGCNA)
library(parallel)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")

n_samps <- 1000
set.seed(5)
force_resave <- FALSE


# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)



gtex_expr_url <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz"
gtex_expr_path <- "/home/amorin/Data/Expression_files/GTEx/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz"


gtex_meta_url <- "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
gtex_meta_path <- "/home/amorin/Data/Expression_files/GTEx/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"


if (!file.exists(gtex_expr_path)) {
  download.file(url = gtex_expr_url, destfile = gtex_expr_path)
}


if (!file.exists(gtex_meta_path)) {
  download.file(url = gtex_meta_url, destfile = gtex_meta_path)
}


gtex <- read.delim(gtex_expr_path, skip = 2, check.names = FALSE)
gtex_meta <- read.delim(gtex_meta_path)



# Subset meta to IDs in expression table

gtex_meta <- filter(gtex_meta, SAMPID %in% colnames(gtex))


# Common protein coding genes

common_genes <- intersect(pc_hg$Symbol, gtex$Description)

pc_hg <- pc_hg %>% 
  filter(Symbol %in% common_genes) %>% 
  arrange(match(Symbol, common_genes))





# Format GTEX expression table into a sample x gene matrix

gtex_mat <- gtex %>% 
  filter(Description %in% common_genes) %>% 
  distinct(Description, .keep_all = TRUE) %>%   # note 26 duplicated symbols
  arrange(match(Description, common_genes))

gtex_mat <- t(as.matrix(gtex_mat[, gtex_meta$SAMPID]))

colnames(gtex_mat) <- common_genes


identical(
  gtex_mat["GTEX-1117F-0011-R4b-SM-GI4VM", "ASCL1"],
  filter(gtex, Description == "ASCL1")$`GTEX-1117F-0011-R4b-SM-GI4VM`
)

identical(rownames(gtex_mat), gtex_meta$SAMPID)
identical(pc_hg$Symbol, colnames(gtex_mat))





aggregate_gtex_coexpr <- function(gtex_mat, gtex_meta, pc_df) {
  
  tissues <- unique(gtex_meta$SMTS)
  
  amat <- aggtools::init_agg_mat(pc_df)
  
  for (tissue in tissues) {
  
    message(paste(tissue, Sys.time()))
  
    gtex_meta_sub <- filter(gtex_meta, SMTS == tissue)
    gtex_sub <- gtex_mat[gtex_meta_sub$SAMPID, ]
  
    cmat <- WGCNA::cor(gtex_sub)
    cmat <- na_to_zero(cmat)
    rmat <- aggtools::transform_correlation_mat(cmat, "allrank")
    amat <- amat + rmat
    rm(cmat, rmat, gtex_sub)
    gc(verbose = FALSE)
    
  }
  
  amat <- allrank_mat(amat) 
  max_rank <- max(amat, na.rm = TRUE)
  amat <- amat / max_rank
  amat <- diag_to_one(amat)
  amat <- lowertri_to_symm(amat)
  
  return(amat)
}



agg_path <- "/space/scratch/amorin/TRsc_output/GTEX_aggregate_coexpr.RDS"


if (!file.exists(agg_path)) {
  
  agg <- aggregate_gtex_coexpr(gtex_mat, gtex_meta, pc_hg)
  saveRDS(agg, agg_path)
  
 } else {
   
   agg <- readRDS(agg_path)
  
}









# Get ortho-matched symbols of TFs with available data, as well as all ortho
# targets, which are used for null
# ------------------------------------------------------------------------------


# Human

ortho_tf_hg <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$TF_Symbol | Symbol_hg %in% curated$TF_Symbol) %>% 
  filter(Symbol_hg %in% names(rank_tf_hg)) %>% 
  pull(Symbol_hg)

ortho_target_hg <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$Target_Symbol | Symbol_hg %in% curated$Target_Symbol) %>% 
  filter(Symbol_hg %in% rownames(rank_tf_hg[[1]])) %>% 
  pull(Symbol_hg)

tf_hg <- union(
  intersect(names(rank_tf_hg), str_to_upper(curated$TF_Symbol)),
  ortho_tf_hg)

target_hg <- union(
  intersect(rank_tf_hg[[1]]$Symbol, str_to_upper(curated$Target_Symbol)),
  ortho_target_hg)


# Mouse

ortho_tf_mm <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$TF_Symbol | Symbol_hg %in% curated$TF_Symbol) %>% 
  filter(Symbol_mm %in% names(rank_tf_mm)) %>% 
  pull(Symbol_mm)

ortho_target_mm <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$Target_Symbol | Symbol_hg %in% curated$Target_Symbol) %>% 
  filter(Symbol_mm %in% rownames(rank_tf_mm[[1]])) %>% 
  pull(Symbol_mm)

tf_mm <- union(
  intersect(names(rank_tf_mm), str_to_title(curated$TF_Symbol)),
  ortho_tf_mm)

target_mm <- union(
  intersect(rank_tf_mm[[1]]$Symbol, str_to_title(curated$Target_Symbol)),
  ortho_target_mm)



# Aggregate coexpr rankings relative to single dataset profiles
# ------------------------------------------------------------------------------


# Human
save_function_results(
  path = avg_vs_ind_auc_hg_path,
  fun = get_colwise_curated_auc_list,
  args = list(
    tfs = tf_hg,
    agg_l = agg_tf_hg,
    msr_mat = msr_hg,
    curated_df = curated,
    ortho_df = pc_ortho,
    pc_df = pc_hg,
    species = "Human",
    ncores = ncore,
    verbose = TRUE
  ),
  force_resave = force_resave
)



# Mouse
save_function_results(
  path = avg_vs_ind_auc_mm_path,
  fun = get_colwise_curated_auc_list,
  args = list(
    tfs = tf_mm,
    agg_l = agg_tf_mm,
    msr_mat = msr_mm,
    curated_df = curated,
    ortho_df = pc_ortho,
    pc_df = pc_mm,
    species = "Mouse",
    ncores = ncore,
    verbose = TRUE
  ),
  force_resave = force_resave
)



# Positive coexpr (top of aggr coexpr rankings)
# ------------------------------------------------------------------------------


# Human
save_function_results(
  path = coexpr_auc_hg_path,
  fun = curated_obs_and_null_auc_list,
  args = list(
    tfs = tf_hg,
    rank_l = rank_tf_hg,
    score_col = "Avg_aggr_coexpr",
    curated_df = curated,
    label_all = target_hg,
    ortho_df = pc_ortho,
    pc_df = pc_hg,
    species = "Human",
    n_samps = n_samps,
    ncores = ncore,
    verbose = TRUE
  ),
  force_resave = force_resave
)





# tissue <- "Kidney"
# 
# 
# gtex_meta_sub <- filter(gtex_meta, SMTS == tissue)
# gtex_sub <- gtex_mat[gtex_meta_sub$SAMPID, ]
# 
# cmat <- WGCNA::cor(gtex_sub)
# rmat <- aggtools::transform_correlation_mat(cmat, "allrank")
# amat <- amat + rmat
# rm(cmat, rmat, gtex_sub)
# gc(verbose = FALSE)


