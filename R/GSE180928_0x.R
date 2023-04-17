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

rownames(dat) <- dat$X
dat$X <- NULL
colnames(dat) <- str_replace_all(colnames(dat), "\\.", "_")
dat <- as.matrix(dat)

meta <- meta %>% 
  dplyr::rename(ID = X, Cell_type = Cluster) %>% 
  mutate(ID = str_replace_all(ID, "-", "_"))

stopifnot(all(colnames(dat) %in% meta$ID))

dat <- dat[, meta$ID]

stopifnot(identical(colnames(dat), meta$ID))

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

# Loop over TFs and cor with every other gene

thresh <- floor(nrow(dat) * 0.005)
thresh <- thresh + 1  # +1 for the moment to account for keeping self cor


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


## CHECKING discrenpency in ranking


genes <- rownames(dat)



cor_l <- lapply(cts, function(i) {
  
  cells <- filter(meta, Cell_type == i)$ID
  mat <- t(dat[, cells])
  na_genes <- apply(mat, 2, function(x) sum(x != 0)) <= min_count
  mat[, na_genes] <- NA
  
  if (sum(is.na(mat)) == length(mat)) {
    message(paste(i, "all NAs"))
    return(NA)
  }
  
  cmat <- WGCNA::corAndPvalue(x = mat, 
                              y = mat[, "ASCL1"],
                              use = "pairwise.complete.obs", 
                              alternative = "greater")
  
})
names(cor_l) <- cts
cor_l <- cor_l[!is.na(cor_l)]


cor_mat <- do.call(cbind, lapply(cor_l, function(x) x$cor))
pval_mat <- do.call(cbind, lapply(cor_l, function(x) x$p))
colnames(cor_mat) <- colnames(pval_mat) <- names(cor_l)


cor_binmat <- mclapply(1:ncol(cor_mat), function(x) {
  vec <- cor_mat[, x]
  passes <- names(sort(vec, decreasing = TRUE)[1:thresh])
  vec <- ifelse(names(vec) %in% passes, 1, 0)
  return(vec)
}, mc.cores = ncore)

cor_binmat <- do.call(cbind, cor_binmat)



pval_binmat <- mclapply(1:ncol(pval_mat), function(x) {
  vec <- pval_mat[, x]
  passes <- names(sort(vec)[1:thresh])
  vec <- ifelse(names(vec) %in% passes, 1, 0)
  return(vec)
}, mc.cores = ncore)
pval_binmat <- do.call(cbind, pval_binmat)


colnames(cor_binmat) <- colnames(pval_binmat) <- names(cor_l)
rownames(cor_binmat) <- rownames(pval_binmat) <- genes


sum_cor_binmat <- rowSums(cor_binmat, na.rm = TRUE)
rank_cor_binmat <- rank(-sum_cor_binmat, ties.method = "min", na.last = "keep")


sum_pval_binmat <- rowSums(pval_binmat, na.rm = TRUE)
rank_pval_binmat <- rank(-sum_pval_binmat, ties.method = "min", na.last = "keep")



# Convert to ranks and set NA ranks to mean of rank matrix
rmat <- colrank_mat(cor_mat)
rmat <- na_to_mean(rmat)
sum_rmat <- rowSums(rmat)
rank_sum_rmat <- rank(sum_rmat)


df <- data.frame(
  Symbol = genes,
  Bin_cor_count = sum_cor_binmat[genes],
  Bin_pval_count = sum_pval_binmat[genes],
  Rank_sum = sum_rmat[genes],
  Rank_bin_cor = rank_cor_binmat[genes],
  Rank_bin_pval = rank_pval_binmat[genes],
  Rank_sum_rank = rank_sum_rmat[genes]
)


gene <- "DLL1"


df2 <- data.frame(
  Cor = cor_mat[gene, ],
  Binary_threshold = pval_binmat[gene, ],
  Rank_cor = rmat[gene, ]
) %>% 
  arrange(desc(Cor))



##


# Save out
saveRDS(amat_rsr, "~/scratch/R_objects/GSE180928_ranksumrank_celltype_cor.RDS")
saveRDS(amat_binthresh, "~/scratch/R_objects/GSE180928_binthresh_celltype_cor.RDS")
