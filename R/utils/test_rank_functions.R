## Testing functions used to generate the aggregate correlation matrices
## -----------------------------------------------------------------------------

library(assertthat)
library(parallel)
source("R/utils/functions.R")
source("R/00_config.R")

id <- "GSE180928"
processed_path <- file.path(amat_dir, id, paste0(id, "_clean_mat_and_meta.RDS"))
# allrank_path <- file.path(out_dir, paste0(id, "_RSR_allrank.RDS"))
# corlist_path <- file.path(out_dir, paste0(id, "_celltype_corlist.RDS"))
# zcor_path <- file.path(out_dir, paste0(id, "_fishersZ.RDS"))


dat <- readRDS(processed_path)
meta <- dat[[2]]
mat <- as.matrix(dat[[1]])
# rsr <- readRDS(allrank_path)
# zcor <- readRDS(zcor_path)
# na_mat <- zcor$NA_mat


# Subset for speed of testing
test_mat <- mat[1:1000, ]



# Check init agg mat: matrix of 0s used to track iterative sum of ranks
# ------------------------------------------------------------------------------


genes <- rownames(test_mat)
amat1 <- init_agg_mat(row_genes = genes)
assert_that(all(amat1 == 0))
assert_that(nrow(amat1) == length(genes))
assert_that(ncol(amat1) == length(genes))



# Check subset and filter: used to filter the full count matrix to just the
# current cell type, and to set genes that are below the min count filter to NA
# ------------------------------------------------------------------------------


min_count <- 20
ct <- meta$Cell_type[1]
ct_ids <- filter(meta, Cell_type == ct)$ID

assert_that(length(ct_ids) > 0)


# Subset with no filtering
ct_mat1 <- t(test_mat[, ct_ids])


# Subset with filter
ct_mat2 <- subset_and_filter(test_mat, 
                             meta, 
                             cell_type = ct, 
                             min_count = min_count)


# Ensure genes are now columns

assert_that(all(colnames(ct_mat1) %in% rownames(test_mat)))
assert_that(all(colnames(ct_mat2) %in% rownames(test_mat)))


# Coerce a highly expressed gene within the cell type to mostly 0s so as to get
# caught by the filter

expr_gene <- names(which.max(colMeans(ct_mat1)))

assert_that(length(expr_gene) == 1)

# Applying to a copy of the original full matrix (note transposed dimensions)

test_mat_expr_gene_to_0 <- test_mat
test_mat_expr_gene_to_0[expr_gene, ] <- 0

ct_mat3 <- subset_and_filter(test_mat_expr_gene_to_0, 
                             meta,
                             cell_type = ct, 
                             min_count = min_count)


assert_that(all(!is.na(ct_mat1[, expr_gene])))
assert_that(all(!is.na(ct_mat2[, expr_gene])))
assert_that(all(is.na(ct_mat3[, expr_gene])))



# Testing generation of correlation matrix for a single cell type
# ------------------------------------------------------------------------------







# Default is random
colrank1 <- colrank_mat(test_mat) 
colrank2 <- colrank_mat(test_mat) 
are_equal(colrank1, colrank2)


# Self as best rank=1
all(diag(colrank1) == 1)



# Default is random
allrank1 <- allrank_mat(test_mat)
allrank2 <- allrank_mat(test_mat)
are_equal(allrank1, allrank2)


