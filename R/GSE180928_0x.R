library(tidyverse)
library(Seurat)
library(future)
library(parallel)
library(WGCNA)
source("R/utils/functions.R")
source("R/00_config.R")


sc_dir_hg <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/alex_sc_requests/human/has_celltype_metadata/"
dat <- read.delim(file.path(sc_dir_hg, "GSE180928/GSE180928_filtered_cell_counts.csv"), sep = ",")
meta <- read.delim(file.path(sc_dir_hg, "GSE180928/GSE180928_metadata.csv"), sep = ",")

dat2 <- dat
rownames(dat2) <- dat$X
dat2$X <- NULL
colnames(dat2) <- str_replace_all(colnames(dat2), "\\.", "_")
dat2 <- as.matrix(dat2)

meta <- meta %>% 
  dplyr::rename(ID = X, Cell_type = Cluster) %>% 
  mutate(ID = str_replace_all(ID, "-", "_"))

stopifnot(all(colnames(dat2) %in% meta$ID))

dat2 <- dat2[, meta$ID]

stopifnot(identical(colnames(dat2), meta$ID))

cts <- unique(meta$Cell_type)



# Min count of non-zero expressing cells to keep gene for correlating
min_count <- 20

# Init aggregated matrix
amat_binthresh <- matrix(0, nrow = nrow(dat2), ncol = nrow(dat2))
rownames(amat_binthresh) <- colnames(amat_binthresh) <- rownames(dat2)
amat_rsr <- amat_binthresh

# Loop over TFs and cor with every other gene

thresh <- floor(nrow(dat2) * 0.005)
thresh <- thresh + 1  # +1 for the moment to account for keeping self cor


for (i in cts) {
  
  message(paste(i, Sys.time()))
  
  # Isolate cells from cell type. If not enough cells meet the minimum non-zero
  # threshold, set all counts to NA (producing an NA cor)
  
  cells <- filter(meta, Cell_type == i)$ID
  mat <- t(dat2[, cells])
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

# Save out
saveRDS(amat_rsr, "~/scratch/R_objects/GSE180928_ranksumrank_celltype_cor.RDS")
saveRDS(amat_binthresh, "~/scratch/R_objects/GSE180928_binthresh_celltype_cor.RDS")
