library(tidyverse)
library(Seurat)
library(future)
library(parallel)
library(WGCNA)
source("R/utils/functions.R")
source("R/00_config.R")

sc_dir_hg <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/alex_sc_requests/human/has_celltype_metadata/"
sc_dir_mm <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/alex_sc_requests/mice/has_celltype_metadata/"

dat <- read.delim(file.path(sc_dir_mm, "GSE93374/GSE93374_Merged_all_020816_BatchCorrected_LNtransformed_doubletsremoved_Data.txt"))
meta <- read.delim(file.path(sc_dir_mm, "GSE93374/GSE93374_cell_metadata.txt"))

dat <- as.matrix(dat)

# Get common IDs (different between data and metadata...)
meta_ids <- str_replace(meta$X1.ID, ".*_", "")
dat_ids <- str_replace(colnames(dat), ".*_", "")

# Remove handful of duplicate IDs
meta_ids_dup <- duplicated(meta_ids)
dat_ids_dup <- duplicated(dat_ids)

dat <- dat[, !dat_ids_dup]
meta <- meta[!meta_ids_dup, ]

colnames(dat) <- str_replace(colnames(dat), ".*_", "")

meta <- meta %>% 
  mutate(ID = str_replace(X1.ID, ".*_", "")) %>% 
  filter(ID %in% intersect(colnames(dat), ID)) %>% 
  dplyr::rename(Cell_type = X7.clust_all)

dat <- dat[, meta$ID]

stopifnot(identical(colnames(dat), meta$ID))

# Cell types to iterate over
cts <- unique(meta$Cell_type)

tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

# Min count of non-zero expressing cells to keep gene for correlating
min_count <- 20

# Init aggregated matrix
# amat_binthresh <- matrix(0, nrow = nrow(dat), ncol = nrow(dat))
# rownames(amat_binthresh) <- colnames(amat_binthresh) <- rownames(dat)
# amat_rsr <- amat_binthresh

amat_binthresh <- matrix(0, nrow = nrow(dat), ncol = length(tfs))
rownames(amat_binthresh) <- rownames(dat)
colnames(amat_binthresh) <- tfs
amat_rsr <- amat_binthresh




# Threshold of top 0.5% of correlation per vector (selected by one-sided pval)
thresh <- floor(nrow(dat) * 0.005)
thresh <- thresh + 1  # +1 for the moment to account for keeping self cor



# Loop over TFs and cor with every other gene

for (i in cts) {
  
  message(paste(i, Sys.time()))
  
  # Isolate cells from cell type. If not enough cells meet the minimum non-zero
  # threshold, set all counts to NA (producing an NA cor)
  
  cells <- filter(meta, Cell_type == i)$ID
  mat <- t(dat[, cells])
  na_genes <- apply(mat, 2, function(x) sum(x != 0)) <= min_count
  mat[, na_genes] <- NA
  
  
  if (sum(is.na(mat)) == length(mat)) {
    message(paste(i, "all NAs"))
    next()
  }
  
  # Pearson cor mat
  # cmat <- WGCNA::corAndPvalue(x = mat, use = "pairwise.complete.obs", alternative = "greater")
  
  cmat <- WGCNA::corAndPvalue(x = mat, 
                              y = mat[, tfs],
                              use = "pairwise.complete.obs", 
                              alternative = "greater")
  
  # Threshold by pval
  
  binmat <- mclapply(1:ncol(cmat$p), function(x) {
    
    vec <- cmat$p[, x]
    passes <- names(sort(vec)[1:thresh])
    vec <- ifelse(names(vec) %in% passes, 1, 0)
    return(vec)
    
  }, mc.cores = ncore)
  
  binmat <- do.call(cbind, binmat)
  
  
  # Convert to ranks and set NA ranks to mean of rank matrix
  rmat <- colrank_mat(cmat$cor)
  rmat <- na_to_mean(rmat)
  
  
  # Running addition to aggregate matrix
  amat_binthresh <- amat_binthresh + binmat
  amat_rsr <- amat_rsr + rmat
  
  # Clean up
  gc()
  
}


# Rank of sum of ranks
amat_rsr <- colrank_mat(-amat_rsr)
amat_binthresh2 <- colrank_mat(amat_binthresh)


## CHECKING discrenpency in ranking

genes <- rownames(dat)


df <- data.frame(
  Symbol = genes,
  RSR = amat_rsr[genes, "Ascl1"],
  BT_count = amat_binthresh[genes, "Ascl1"],
  BT_rank = amat_binthresh2[genes, "Ascl1"])


filter(df, Symbol == "Prdx6")


cor_l <- lapply(cts, function(i) {
  
  cells <- filter(meta, Cell_type == i)$ID
  mat <- t(dat[, cells])
  na_genes <- apply(mat, 2, function(x) sum(x != 0)) <= min_count
  mat[, na_genes] <- NA
  
  if (sum(is.na(mat)) == length(mat)) {
    message(paste(i, "all NAs"))
    next()
  }
  
  cmat <- WGCNA::corAndPvalue(x = mat, 
                              y = mat[, "Ascl1"],
                              use = "pairwise.complete.obs", 
                              alternative = "greater")
  
})
names(cor_l) <- cts
  

cor_mat <- do.call(cbind, lapply(cor_l, function(x) x$cor))
colnames(cor_mat) <- cts
cor_mat[1:5, 1:5]


binmat <- mclapply(1:ncol(cor_mat), function(x) {
  vec <- cor_mat[, x]
  passes <- names(sort(vec)[1:thresh])
  vec <- ifelse(names(vec) %in% passes, 1, 0)
  return(vec)
}, mc.cores = ncore)


binmat <- do.call(cbind, binmat)
colnames(binmat) <- cts
rownames(binmat) <- genes
binmat[1:5, 1:5]
head(sort(binmat[, 1], decreasing = TRUE), 100)
binmat["Ascl1", ]


# Convert to ranks and set NA ranks to mean of rank matrix
rmat <- colrank_mat(cor_mat)
rmat <- na_to_mean(rmat)


##


# Save out
saveRDS(amat_rsr, "~/scratch/R_objects/GSE93374_ranksumrank_celltype_cor.RDS")
saveRDS(amat_binthresh, "~/scratch/R_objects/GSE93374_binthresh_celltype_cor.RDS")
