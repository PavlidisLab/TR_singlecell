library(WGCNA)
library(tidyverse)
source("R/utils/functions.R")


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


thresh <- floor(nrow(dat) * 0.005)
thresh <- thresh + 1  # +1 for the moment to account for keeping self cor



# Subset for speed of testing

# mat <- dat[1:1000, filter(meta, Cell_type == cts[1])$ID]
# meta_sub <- filter(meta, Cell_type == meta$Cell_type[1])
# meta_sub2 <- filter(meta, Cell_type %in% cts[1:2])
# stopifnot(identical(ncol(mat), nrow(meta_sub)))
# saveRDS(list(mat, meta_sub, meta_sub2), file = "~/scratch/R_objects/17_04_2023_GSE180928_subset_for_test.RDS")


dat2 <- readRDS("~/scratch/R_objects/17_04_2023_GSE180928_subset_for_test.RDS")
mat <- dat2[[1]]
meta_sub <- dat2[[2]]
meta_sub2 <- dat2[[3]]



# TODO: ranksumrank function
# TODO: need to distinguish core agg function from apply over all
# TODO: binthresh function
# TODO: minrank function (should return best cell type and value as well)
# TODO: does `...` with col_subset declared work?
# TODO: aggregate mat symmetric?
# TODO: col subset is clumsy, actually needed?
# TODO: be explicit about shape of input mat (is t() burried bad?)
# TODO: col/rowrank sign confusing?



# Rank matrix columns such that 1=best. Return as same dimension as input.

colrank_mat <- function(mat, ties_arg = "min", na_arg = "keep") {
  rank_mat <- apply(-mat, 2, rank, ties.method = ties_arg, na.last = na_arg)
  return(rank_mat)
}


# Rank matrix rows such that 1=best. Return as same dimension as input.

rowrank_mat <- function(mat, ties_arg = "min", na_arg = "keep") {
  rank_mat <- apply(-mat, 1, rank, ties.method = ties_arg, na.last = na_arg)
  return(t(rank_mat))
}


# If a col of mat has fewer than min_count elements that are non-zero, set that
# col to NAs. This is done to produce an NA during correlation, instead of 
# allowing values resulting from fewer observations.

under_min_count_to_na <- function(mat, min_count = 20) {
  
  na_genes <- apply(mat, 2, function(x) sum(x != 0)) <= min_count
  mat[, na_genes] <- NA
  return(mat)
}


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



# 

cor_and_rank <- function(mat, 
                         cor_method = "pearson",
                         na_arg = "last",
                         ties_arg = "min") {
  
  stopifnot(na_arg %in% c("last", "mean"))
  
  cmat <- WGCNA::cor(
    x = mat, method = cor_method, use = "pairwise.complete.obs")
  
  if (na_arg == "last") {
    rmat <- colrank_mat(cmat, na_arg = TRUE, ties_arg = ties_arg)
  } else {
    rmat <- colrank_mat(cmat, na_arg = "keep", ties_arg = ties_arg)
    rmat <- na_to_mean(rmat)
  }
  
  return(rmat)
}



# Rank sum rank from Harris et al., 2021 (Jesse Gillis) 
# https://pubmed.ncbi.nlm.nih.gov/34015329/
# Rank coexpression (1=best) across cell types. Set NAs to network mean. Sum 
# the cell-type ranks, and then rank order these sums (1=best).


all_RSR_aggregate <- function(mat,
                              meta, 
                              min_cell = 20, 
                              # cor_method = "pearson",
                              # na_action,
                              # ties_arg = "min",
                              ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    rmat <- cor_and_rank(ct_mat)
    amat <- amat + rmat
  
  }
  
  amat <- colrank_mat(-amat, na_arg = "last", ties_arg = "min")
  
  return(amat)
}



# 


# Check init agg mat
genes <- rownames(mat)
amat1 <- init_agg_mat(row_genes = genes)
amat2 <- init_agg_mat(row_genes = genes, col_genes = str_to_upper(tfs))


# Check min to NA
na_mat <- under_min_count_to_na(t(mat))
na_gene <- genes[which(is.na(na_mat), arr.ind = TRUE)[, "col"][1]]
sum(mat[na_gene, ] != 0) < min_count
sum(is.na(mat))
sum(is.na(na_mat))


# Check row rank
# TODO: necessary? colrank
rowrank1 <- rowrank_mat(mat) # default is min
rowrank2 <- rowrank_mat(mat, ties_arg = "random")
rowrank3 <- rowrank_mat(mat, ties_arg = "max")
assertthat::are_equal(dim(mat), dim(rowrank1))


# Check cor and rank
cr1 <- cor_and_rank(mat, na_arg = "last")
cr2 <- cor_and_rank(na_mat, na_arg = "last")

cr3 <- cor_and_rank(mat, na_arg = "mean")
cr4 <- cor_and_rank(na_mat, na_arg = "mean")

cr5 <- cor_and_rank(mat, ties_arg = "min")
cr6 <- cor_and_rank(na_mat, ties_arg = "random")


cr1[1:5, 1:5]
cr2[1:5, 1:5]
cr3[1:5, 1:5]
cr4[1:5, 1:5]
cr5[1:5, 1:5]
cr6[1:5, 1:5]

identical(cr1, cr2)

# Check RSR
rsr1 <- all_RSR_aggregate(mat, meta_sub)
rsr2 <- all_RSR_aggregate(mat, meta_sub, na_action = "mean")
rsr3 <- all_RSR_aggregate(mat, meta_sub, na_acion = "last")
rsr4 <- all_RSR_aggregate(mat, meta_sub, ties_arg = "random")

rsr1[1:5, 1:5]
rsr2[1:5, 1:5]
rsr3[1:5, 1:5]
rsr4[1:5, 1:5]

identical(rsr1, rsr3)


# na_action,

head(sort(cmat[, "BEX2"], decreasing = TRUE))
head(sort(rmat[, "BEX2"], decreasing = FALSE))
head(sort(amat[, "BEX2"], decreasing = FALSE))


check_na <- cmat[is.na(cmat[, 1]), 1]
rmat[names(check_na), 1]
sum(is.na(cmat[,1]))

mat[1:5, 1:5]
cmat[1:5, 1:5]
rmat[1:5, 1:5]
amat[1:5, 1:5]



# Speed of cor. ncore8 actually slightly slows down, likely NA related. Pval
res <- microbenchmark::microbenchmark(
  F1 = WGCNA::cor(mat, nThreads = 1),
  F2 = WGCNA::cor(mat, nThreads = 8),
  F3 = WGCNA::corAndPvalue(mat),
  F4 = stats::cor(mat),
  times = 10
)




binary_threshold <- function(mat, threshold) {
  
  # Threshold by pval
  
  binmat <- mclapply(1:ncol(cmat$p), function(x) {
    
    vec <- cmat$p[, x]
    passes <- names(sort(vec)[1:thresh])
    vec <- ifelse(names(vec) %in% passes, 1, 0)
    return(vec)
    
  }, mc.cores = ncore)
  
  binmat <- do.call(cbind, binmat)
  
}



all_threshold_aggregate <- function(mat,
                                    meta, 
                                    min_cell = 20, 
                                    cor_method = "pearson",
                                    threshold = 0,
                                    ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    cells <- filter(meta, Cell_type == ct)$ID
    ct_mat <- t(mat[, cells])
    na_genes <- apply(ct_mat, 2, function(x) sum(x != 0)) <= min_cell
    ct_mat[, na_genes] <- NA
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "all NAs"))
      next()
    }
    
    cmat <- WGCNA::corAndPvalue(
      x = mat, 
      method = cor_method, 
      use = "pairwise.complete.obs", 
      alternative = "greater")
    
    
    
    
    
    rmat <- colrank_mat(cmat, na_arg = TRUE)
    
    amat <- amat + rmat
    
  }
  
  amat <- colrank_mat(-amat)
  
  return(amat)
}
