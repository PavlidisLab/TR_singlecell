library(tidyverse)
library(Seurat)
library(future)
library(parallel)
library(WGCNA)
source("R/utils/functions.R")
source("R/00_config.R")

sc_dir_hg <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/alex_sc_requests/human/has_celltype_metadata/"
sc_dir_mm <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/alex_sc_requests/mice/has_celltype_metadata/"

files_hg <- list.files(sc_dir_hg, recursive = TRUE)
files_mm <- list.files(sc_dir_mm, recursive = TRUE)

"GSE180928/GSE180928_filtered_cell_counts.csv" 
"GSE180928/GSE180928_metadata.csv" 

"GSE188516/GSE188516_GEO_RNA_barcodes.tsv"                              
"GSE188516/GSE188516_GEO_RNA_counts.mtx"                                
"GSE188516/GSE188516_GEO_RNA_features.tsv"

"GSE195445/GSE195445_Human_obj.rds"                                     
"GSE195445/GSE195445_Mouse_metadata.csv"                                
"GSE195445/GSE195445_Mouse_normalized_count.csv" 

"GSE217511/GSE217511_Combined_log-normalized_AveExp_snRNAseq_Human.xlsx"

"GSE139600/GSE139600_Processed_SS2_Neuron_TPM.txt" 

"GSE160797/GSE160797_normalized_counts_eye_macs.csv"                                       

"GSE160797/GSE160797_normalized_counts_retina_baseline_cnv.csv" 

"GSE163252/GSE163252_DRG.processed.umi.raw.counts.txt" 

"GSE168031/GSE168031_expression_mats_counts.csv"                                           
"GSE168031/GSE168031_expression_mats_groups.csv"                                           
"GSE168031/GSE168031_expression_mats_TPM.csv"    

"GSE93374/GSE93374_cell_metadata.txt"                                                      
"GSE93374/GSE93374_Merged_all_020816_BatchCorrected_LNtransformed_doubletsremoved_Data.txt"



dat <- read.delim(file.path(sc_dir_mm, "GSE93374/GSE93374_Merged_all_020816_BatchCorrected_LNtransformed_doubletsremoved_Data.txt"))
meta <- read.delim(file.path(sc_dir_mm, "GSE93374/GSE93374_cell_metadata.txt"))

meta_ids <- str_replace(meta$X1.ID, ".*_", "")
dat_ids <- str_replace(colnames(dat), ".*_", "")

meta_ids_dup <- duplicated(meta_ids)
dat_ids_dup <- duplicated(dat_ids)

dat2 <- dat[, !dat_ids_dup]
meta2 <- meta[!meta_ids_dup, ]

colnames(dat2) <- str_replace(colnames(dat2), ".*_", "")

meta2 <- meta2 %>% 
  mutate(ID = str_replace(X1.ID, ".*_", "")) %>% 
  filter(ID %in% intersect(colnames(dat2), ID)) %>% 
  dplyr::rename(Cell_type = X7.clust_all)

dat2 <- dat2[, meta2$ID]

stopifnot(identical(colnames(dat2), meta2$ID))

cts <- unique(meta2$Cell_type)


# Min count of non-zero expressing cells to keep gene for correlating
min_count <- 20

# Init aggregated matrix
amat <- matrix(0, nrow = nrow(dat2), ncol = nrow(dat2))
rownames(amat) <- colnames(amat) <- rownames(dat2)


# Loop over TFs and cor with every other gene

thresh <- floor(nrow(dat2) * 0.005)
thresh <- thresh + 1  # +1 for the moment to account for keeping self cor


for (i in cts) {
  
  message(paste(i, Sys.time()))
  
  # Isolate cells from cell type. If not enough cells meet the minimum non-zero
  # threshold, set all counts to NA (producing an NA cor)
  
  cells <- filter(meta2, Cell_type == i)$ID
  mat <- t(as.matrix(dat2[, cells]))
  na_genes <- apply(mat, 2, function(x) sum(x != 0)) <= min_count
  mat[, na_genes] <- NA
  
  # Pearson cor mat
  cmat <- WGCNA::corAndPvalue(x = mat, use = "pairwise.complete.obs", alternative = "greater")
  
  # Threshold by pval
  
  binmat <- mclapply(1:ncol(cmat$p), function(x) {
    
    vec <- cmat$p[, x]
    passes <- names(sort(vec)[1:thresh])
    vec <- ifelse(names(vec) %in% passes, 1, 0)
    return(vec)
    
  }, mc.cores = ncore)
  
  binmat <- do.call(cbind, binmat)
  
  # Running addition to aggregate matrix
  amat <- amat + binmat
  
  # Clean up
  gc()
  
}


# Rank of sum of ranks
# amat <- colrank_mat(-amat)

# Save out
saveRDS(amat, "~/scratch/R_objects/GSE93374_ranksumrank_celltype_cor.RDS")
