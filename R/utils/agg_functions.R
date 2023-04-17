library(WGCNA)
library(tidyverse)


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
genes <- rownames(dat)

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


mat <- dat[, filter(meta, Cell_type == cts[1])$ID]




# TODO: ranksumrank function
# TODO: need to distinguish core agg function from apply over all
# TODO: binthresh function
# TODO: minrank function (should return best cell type and value as well)




# Return a gene x gene matrix of 0s used to track coexpression aggregation
# across cell types. Rows are the full set of provided genes, while columns
# are an optional subset of genes.


init_agg_mat <- function(row_genes, col_genes = NULL) {
  
  if (!is.null(col_genes)) {
    
    stopifnot(all(col_genes %in% row_genes))
    
    amat <- matrix(0, nrow = length(row_genes), ncol = length(col_genes))
    rownames(amat) <- row_genes
    colnames(amat) <- col_genes
    
  } else {
      
    amat <- matrix(0, nrow = length(row_genes), ncol = length(row_genes))
    rownames(amat) <- colnames(amat) <- row_genes
  
  }
  
  return(amat)
}



# Rank sum rank from Harris et al., 2021 (Jesse Gillis) 
# https://pubmed.ncbi.nlm.nih.gov/34015329/
# Rank coexpression (1=best) across cell types. Set NAs to network mean. Sum 
# the cell-type ranks, and then rank order these sums (1=best).


aggregate_cor <- function(cmat_list, impute_na = TRUE, ncores = 1) {
  
  # Convert cors to ranks (1=best)
  rank_list <- mclapply(cmat_list, colrank_mat, mc.cores = ncores)
  
  # Set NAs
  if (impute_na) {
    rank_list <- lapply(rank_list, na_to_mean)
  }
  
  # Sum list of rank matrices into a single matrix
  # https://stackoverflow.com/questions/42628385/sum-list-of-matrices-with-nas
  sum_rank <- apply(simplify2array(rank_list), 1:2, sum, na.rm = TRUE)
  
  # Convert sum of ranks into a final rank (1=best)
  final_rank <- colrank_mat(-sum_rank)
  
  return(final_rank)
}






# Check init agg mat
amat1 <- init_agg_mat(row_genes = genes)
amat2 <- init_agg_mat(row_genes = genes, col_genes = tfs)


# Check row rank
rowrank1 <- rowrank_mat(mat) # default is min
rowrank2 <- rowrank_mat(mat, ties_arg = "random")
rowrank3 <- rowrank_mat(mat, ties_arg = "max")
assertthat::are_equal(dim(mat), dim(rowrank1))

