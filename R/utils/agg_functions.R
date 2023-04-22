library(WGCNA)
library(tidyverse)
source("R/utils/functions.R")


sc_dir_hg <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/alex_sc_requests/human/has_celltype_metadata/"
dat_path <- file.path(sc_dir_hg, "GSE180928/GSE180928_filtered_cell_counts.csv")
meta_path <- file.path(sc_dir_hg, "GSE180928/GSE180928_metadata.csv")
out_path <- "/space/scratch/amorin/R_objects/20_04_2023_GSE180928.RDS"


if (!file.exists(out_path)) {
  
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
  
  saveRDS(list(dat, meta), file = out_path)
  
} else {
  
  dat <- readRDS(out_path)
  meta <- dat[[2]]
  dat <- dat[[1]]
  
}



stopifnot(identical(colnames(dat), meta$ID))

cts <- unique(meta$Cell_type)
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")
genes <- rownames(dat)


# Min count of non-zero expressing cells to keep gene for correlating
min_cell <- 20


thresh <- floor(nrow(dat) * 0.005)
thresh <- thresh + 1  # +1 for the moment to account for keeping self cor



# Subset for speed of testing

mat <- dat[1:1000, ]
meta_sub <- filter(meta, Cell_type %in% cts[1:2])


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


# TODO: 

allrank_mat <- function(mat, ties_arg = "min", na_arg = "keep") {
  
  rmat <- rank(-mat, ties.method = ties_arg, na.last = na_arg)
  rmat <- matrix(rmat, nrow = nrow(mat), ncol = ncol(mat))
  rownames(rmat) <- colnames(rmat) <- rownames(mat)
  return(rmat)
}


# TODO:

get_cor_mat <- function(mat, 
                        cor_method = "pearson",
                        lower_tri = TRUE,
                        ncores = 1) {
  
  cmat <- WGCNA::cor(
    mat, method = cor_method, use = "pairwise.complete.obs", nThreads = ncores)
  
  if (lower_tri) {
    diag(cmat) <- NA
    cmat[upper.tri(cmat)] <- NA
  }
  
  return(cmat)
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


# 1: Column rank where NAs cors are set to 0.
all_RSR_aggregate1 <- function(mat,
                              meta, 
                              min_cell = 20, 
                              ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    cmat[is.na(cmat)] <- 0
    # diag(cmat) <- NA
    rmat <- colrank_mat(cmat)
    amat <- amat + rmat
    
  }
  
  amat <- colrank_mat(-amat, na_arg = "last", ties_arg = "min")
  
  return(amat)
}




# 2: All rank where NAs cors are set to 0.
all_RSR_aggregate2 <- function(mat,
                               meta,
                               min_cell = 20,
                               ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    # Get full symmetric corr mat. Then replace NAs with 0 and make triangular
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    cmat[is.na(cmat)] <- 0
    # diag(cmat) <- NA
    cmat[upper.tri(cmat)] <- NA
    
    rmat <- allrank_mat(cmat)
    amat <- amat + rmat
    
  }
  
  amat <- allrank_mat(-amat)
  
  return(amat)
}



# 3: Col rank where NA rank sent to mean rank
all_RSR_aggregate3 <- function(mat,
                               meta,
                               min_cell = 20,
                               ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    # diag(cmat) <- NA
    rmat <- colrank_mat(cmat, na_arg = "keep")
    rmat <- na_to_mean(rmat)
    amat <- amat + rmat
    
  }
  
  amat <- colrank_mat(-amat, na_arg = "last", ties_arg = "min")
  
  return(amat)
}



# 4: All rank where NA rank sent to mean rank
all_RSR_aggregate4 <- function(mat,
                               meta,
                               min_cell = 20,
                               ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    #
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    cmat[upper.tri(cmat)] <- NA
    
    rmat <- allrank_mat(cmat, na_arg = "keep")
    rmat[lower.tri(rmat) & is.na(rmat)] <- mean(rmat[lower.tri(rmat)], na.rm = TRUE)
  
    amat <- amat + rmat
    
  }
  
  amat <- allrank_mat(-amat)
  
  return(amat)
}



# 5: Col rank where NA rank sent to last
all_RSR_aggregate5 <- function(mat,
                               meta,
                               min_cell = 20,
                               ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    #
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    # diag(cmat) <- NA
    
    # na_arg = TRUE results in placing NAs in random order after last non-NA rank
    # rmat <- colrank_mat(cmat, na_arg = TRUE)
    
    # This sets all NA ranks to length of genes
    rmat <- colrank_mat(cmat, na_arg = "keep")
    rmat[is.na(rmat)] <- length(genes)
    
    amat <- amat + rmat
    
  }
  
  amat <- colrank_mat(-amat, na_arg = "last", ties_arg = "min")
  
  return(amat)
}



# 6: All rank where NA rank sent to last
all_RSR_aggregate6 <- function(mat,
                               meta,
                               min_cell = 20,
                               ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    cmat[upper.tri(cmat)] <- NA
    
    rmat <- allrank_mat(cmat, na_arg = "keep")
    
    n_comparisons <- length(genes) * (length(genes) - 1) / 2
    max_nonNA_rank <- max(rmat, na.rm = TRUE)
    n_na <- length(rmat[lower.tri(rmat) & is.na(rmat)])
    
    # This sets NA ranks to +1 increasing values after last non-NA rank
    # rmat[lower.tri(rmat) & is.na(rmat)] <- (max_nonNA_rank + 1) : (max_nonNA_rank + n_na)
    
    # This sets all NA ranks to the fixed last non-NA rank+1
    rmat[lower.tri(rmat) & is.na(rmat)] <- max_nonNA_rank + 1
    
    amat <- amat + rmat
    
  }
  
  amat <- allrank_mat(-amat)
  
  return(amat)
}


# Average Z-score
all_zscore_aggregate <- function(mat,
                                  meta,
                                  min_cell = 20,
                                  ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  na_mat <- amat
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    
    na_ix <- which(is.na(cmat), arr.ind = TRUE)
    na_mat[na_ix] <- na_mat[na_ix] + 1
    
    cmat[is.na(cmat)] <- 0
    # diag(cmat) <- NA
  
    zmat <- scale(cmat)
    zmat[is.na(zmat)] <- 0
    
    amat <- amat + zmat

  }
  
  agg_list <- list(
    Amat = amat,
    Avg_all = (amat / length(cts)),
    Avg_nonNA = (amat / (length(cts) - na_mat)),
    NA_mat = na_mat)
  
  return(agg_list)
}



# Check init agg mat
genes <- rownames(mat)
amat1 <- init_agg_mat(row_genes = genes)
# amat2 <- init_agg_mat(row_genes = genes, col_genes = str_to_upper(tfs))


# Check min to NA. Transpose as intended that genes are columns
na_mat <- under_min_count_to_na(t(mat))
na_gene <- genes[which(is.na(na_mat), arr.ind = TRUE)[, "col"][1]]
sum(mat[na_gene, ] != 0) < min_count
sum(is.na(mat))
sum(is.na(na_mat))
na_mat[1:5, 1:5]


# Check cor. Transpose mat so that genes are cols.
cmat1 <- get_cor_mat(t(mat))
cmat2 <- get_cor_mat(t(mat), lower_tri = FALSE)
cmat3 <- get_cor_mat(na_mat)
cmat4 <- get_cor_mat(na_mat, lower_tri = FALSE)

all(is.na(cmat3[na_gene, ]))
table(!is.na(cmat2[na_gene, ]))

cmat1[1:5, 1:5]
cmat2[1:5, 1:5]
cmat3[1:5, 1:5]
cmat4[1:5, 1:5]


# Check ranking of each element of matrix
allrank1 <- allrank_mat(cmat1)
allrank2 <- allrank_mat(cmat2)
allrank3 <- allrank_mat(cmat3)
allrank4 <- allrank_mat(cmat4)

allrank1[1:5, 1:5]
allrank2[1:5, 1:5]

which(cmat1 == max(cmat1, na.rm = TRUE), arr.ind = TRUE)
which(allrank1 == min(allrank1, na.rm = TRUE), arr.ind = TRUE)

cmat1[630, 498, drop = FALSE]
allrank1[630, 498, drop = FALSE]

df <- mat_to_df(allrank1) %>% 
  dplyr::rename(Rank_cor = Value) %>% 
  mutate(Cor = mat_to_df(cmat1)[, "Value"])



# Check row rank. Default is ties as min and keep NA. mat will be ordered by
# cells, na_mat will be ordered by genes

rowrank1 <- rowrank_mat(mat)
rowrank2 <- rowrank_mat(na_mat)
rowrank3 <- rowrank_mat(mat, ties_arg = "random")
rowrank4 <- rowrank_mat(na_mat, ties_arg = "random")
rowrank5 <- rowrank_mat(mat, na_arg = "last")
rowrank6 <- rowrank_mat(na_mat, na_arg = "last")

assertthat::are_equal(dim(mat), dim(rowrank1))


head(sort(mat[1, ], decreasing = TRUE))
head(sort(rowrank1[1, ], decreasing = FALSE))
head(sort(rowrank2[1, ], decreasing = FALSE))

rowrank1[1:5, 1:5]
rowrank2[1:5, 1:5]
rowrank3[1:5, 1:5]
rowrank4[1:5, 1:5]
rowrank5[1:5, 1:5]
rowrank6[1:5, 1:5]



# Check col rank
colrank1 <- colrank_mat(t(mat)) # default is min

colrank1[1:5, 1:5]
t(mat)[1:5, 1:5]


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
rsr1 <- all_RSR_aggregate1(mat, meta)
rsr2 <- all_RSR_aggregate2(mat, meta)
rsr3 <- all_RSR_aggregate3(mat, meta)
rsr4 <- all_RSR_aggregate4(mat, meta)
rsr5 <- all_RSR_aggregate5(mat, meta)
rsr6 <- all_RSR_aggregate6(mat, meta)

summary(sapply(1:ncol(rsr1), function(x) cor(rsr1[,x], rsr3[,x])))
summary(sapply(1:ncol(rsr1), function(x) cor(rsr1[,x], rsr5[,x])))
summary(sapply(1:ncol(rsr4), function(x) cor(rsr6[,x], rsr5[,x], use = "pairwise.complete.obs")))







## Need to keep track of NAs -> 0s without ranking symmetric NAs

ct_mat <- t(mat[, filter(meta, Cell_type == cts[1])$ID])
ct_mat <- under_min_count_to_na(ct_mat, min_cell)

cmat_tri <- get_cor_mat(ct_mat)
cmat_full <- get_cor_mat(ct_mat, lower_tri = FALSE)

cmat_full_nato0 <- cmat_full
cmat_full_nato0[is.na(cmat_full_nato0)] <- 0

cmat_full_nato0_tri <- cmat_full_nato0
diag(cmat_full_nato0_tri) <- NA
cmat_full_nato0_tri[upper.tri(cmat_full_nato0_tri)] <- NA

rmat_full_nato0_tri <- allrank_mat(cmat_full_nato0_tri)



ct_mat[1:5, 1:5]
cmat_tri[1:5, 1:5]
cmat_full[1:5, 1:5]
cmat_full_nato0[1:5, 1:5]
cmat_full_nato0_tri[1:5, 1:5]
rmat_full_nato0_tri[1:5, 1:5]


sum(is.na(cmat_full))
which(is.na(cmat_full), arr.ind = TRUE)[1,]
cmat_full[24, 1, drop = FALSE]

sum(is.na(cmat_full_nato0))
sum(is.na(cmat_tri))
sum(is.na(cmat_full_nato0_tri))

which(cmat_full_nato0_tri == 0, arr.ind = TRUE)[1,]
cmat_full_nato0_tri[24, 1, drop = FALSE]


which(rmat_full_nato0_tri == min(rmat_full_nato0_tri, na.rm = TRUE), arr.ind = TRUE)
cmat_full_nato0_tri[630, 498, drop = FALSE]
rmat_full_nato0_tri[630, 498, drop = FALSE]

##


cmat[1:5, 1:5]
rmat[1:5, 1:5]
amat[1:5, 1:5]
which(cmat == max(cmat, na.rm = TRUE), arr.ind = TRUE)
which(rmat == min(rmat, na.rm = TRUE), arr.ind = TRUE)
which(amat == min(amat, na.rm = TRUE), arr.ind = TRUE)# Speed of cor. ncore8 actually slightly slows down, likely NA related.
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
