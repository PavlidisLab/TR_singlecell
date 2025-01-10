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


# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)

# Aggregate GTEX coexpression matrix path
agg_path <- "/space/scratch/amorin/TRsc_output/GTEX_aggregate_coexpr.RDS"



# Download and load GTEX V10
# ------------------------------------------------------------------------------

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



# Subset GTEX meta to IDs in expression table

gtex_meta <- filter(gtex_meta, SAMPID %in% colnames(gtex))


# Common protein coding genes and TFs

common_genes <- intersect(pc_hg$Symbol, gtex$Description)
common_tfs <- intersect(tfs_hg$Symbol, gtex$Description)

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



# Aggregate coexpression matrices per tissue
# ------------------------------------------------------------------------------


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




if (!file.exists(agg_path)) {
  
  agg <- aggregate_gtex_coexpr(gtex_mat, gtex_meta, pc_hg)
  saveRDS(agg, agg_path)
  
 } else {
   
   agg <- readRDS(agg_path)
  
}



# Literature curation benchmark
# ------------------------------------------------------------------------------


n_samps <- 1000
set.seed(5)



# Get ortho-matched symbols of TFs with available data, as well as all ortho
# targets, which are used for null


ortho_tf_hg <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$TF_Symbol | Symbol_hg %in% curated$TF_Symbol) %>% 
  filter(Symbol_hg %in% common_tfs) %>% 
  pull(Symbol_hg)

ortho_target_hg <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$Target_Symbol | Symbol_hg %in% curated$Target_Symbol) %>% 
  filter(Symbol_hg %in% common_genes) %>% 
  pull(Symbol_hg)

tf_hg <- union(
  intersect(common_tfs, str_to_upper(curated$TF_Symbol)),
  ortho_tf_hg)

target_hg <- union(
  intersect(common_genes, str_to_upper(curated$Target_Symbol)),
  ortho_target_hg)



# Split TF agg coexpr matrix into a list of dataframes for function call

split_mat_to_list <- function(mat) {
  
  l <- asplit(mat, 2)
  
  l <- lapply(l, function(x) {
    data.frame(Aggr_coexpr = x) %>% 
      rownames_to_column(var = "Symbol") %>% 
      arrange(desc(Aggr_coexpr))
  })
}


agg_tf <- agg[, common_tfs]
agg_l <- split_mat_to_list(agg_tf)


gtex_auc_path <- "/space/scratch/amorin/TRsc_output/gtex_recover_curated_hg.RDS"


save_function_results(
  path = gtex_auc_path,
  fun = curated_obs_and_null_auc_list,
  args = list(
    tfs = tf_hg,
    rank_l = agg_l,
    score_col = "Aggr_coexpr",
    curated_df = curated,
    label_all = target_hg,
    ortho_df = pc_ortho,
    pc_df = pc_hg,
    species = "Human",
    n_samps = n_samps,
    ncores = ncore,
    verbose = TRUE
  ),
  force_resave = FALSE
)
