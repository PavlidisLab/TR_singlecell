## Note that ChatGPT4 was used to initiate most tests, all of which I reviewed
## and further added to.
## -----------------------------------------------------------------------------

library(testthat)
source("R/utils/dev_functions.R")
source("R/00_config.R")


# TODO: more targeted tests for NA counting in aggregate


# Mock data
# ------------------------------------------------------------------------------



set.seed(123)

# 100 cells split 2 cell types, 100 genes, gene 1 mostly 0s

mat_dense <- matrix(sample(1:100, 10000, replace = TRUE), nrow = 100)
mat_dense[1, 1:95] <- 0  
rownames(mat_dense) <- paste0("Gene", 1:100)
colnames(mat_dense) <- paste0("Cell", 1:100)
mat_sparse <- Matrix(mat_dense, sparse = TRUE)


meta <- data.frame(
  ID = paste0("Cell", 1:100),
  Cell_type = rep(c("Type1", "Type2"), each = 50),
  stringsAsFactors = FALSE
)




# colrank_mat(): focus on default ties/NA arguments
# ------------------------------------------------------------------------------



# Test for basic functionality
test_that("colrank_mat ranks columns correctly", {
  
  mock_data <- matrix(c(3, 1, 4, 2, 5, 9, 8, 6, 7), nrow = 3, ncol = 3)
  expected_output <- matrix(c(2, 3, 1, 3, 2, 1, 1, 3, 2), nrow = 3, ncol = 3)
  result <- colrank_mat(mock_data)
  
  expect_equal(result, expected_output)
})



# Test for handling NA values
test_that("colrank_mat handles NA values correctly", {
  
  mock_data <- matrix(c(3, NA, 4, 2, 5, NA, 8, 6, 7), nrow = 3, ncol = 3)
  expected_output <- matrix(c(2, NA, 1, 2, 1, NA, 1, 3, 2), nrow = 3, ncol = 3)
  result <- colrank_mat(mock_data)
  
  expect_equal(result, expected_output)
})



# Test for assigning ties
test_that("colrank_mat assigns ties to minimum values", {
  
  mock_data <- matrix(c(2, 0, 0, -1, 3, 1, 4, 3, 3, 3, 4, 4), nrow = 4, ncol = 3)
  expected_output <- matrix(c(1, 2, 2, 4, 2, 4, 1, 2, 3, 3, 1, 1), nrow = 4, ncol = 3)
  result <- colrank_mat(mock_data)

  expect_equal(result, expected_output)
})




# allrank_mat(): focus on default ties/arg arguments
# ------------------------------------------------------------------------------



# Test for basic functionality
test_that("allrank_mat ranks columns correctly", {
  
  mock_data <- matrix(c(3, 1, 4, 2, 5, 9, 8, 6, 7), nrow = 3, ncol = 3)
  expected_output <- matrix(c(7, 9, 6, 8, 5, 1, 2, 4, 3), nrow = 3, ncol = 3)
  result <- allrank_mat(mock_data)
  
  expect_equal(result, expected_output)
})



# Test for handling NA values
test_that("allrank_mat handles NA values correctly", {
  
  mock_data <- matrix(c(3, NA, 4, 2, 5, NA, 8, 6, 7), nrow = 3, ncol = 3)
  expected_output <- matrix(c(6, NA, 5, 7, 4, NA, 1, 3, 2), nrow = 3, ncol = 3)
  result <- allrank_mat(mock_data)
  
  expect_equal(result, expected_output)
})



# Test for assigning ties
test_that("allrank_mat assigns ties to minimum values", {
  
  mock_data <- matrix(c(-1, 0, 6, 3, 0, 4, -2, 0, 4), nrow = 3, ncol = 3)
  expected_output <- matrix(c(8, 5, 1, 4, 5, 2, 9, 5, 2), nrow = 3, ncol = 3)
  result <- allrank_mat(mock_data)
  
  expect_equal(result, expected_output)
})




# sparse_pcor()
# ------------------------------------------------------------------------------



test_that("sparse_pcor works correctly with valid input", {
  
  set.seed(1)
  
  mat_sparse <- Matrix::rsparsematrix(10, 10, density = 0.5)
  mat_sparse[, 1] <- 0
  colnames(mat_sparse) <- paste0("gene", 1:10)
  
  mat_dense <- as.matrix(mat_sparse)
  
  result_sparse <- sparse_pcor(mat_sparse)
  
  result_dense <- 
    suppressWarnings(cor(mat_dense, method = "pearson", use = "everything"))
  
  expect_equal(result_sparse, result_dense, tolerance = 1e-6)
  
})




test_that("sparse_pcor throws an error with invalid input", {
  
  mat_dense <- matrix(c(1, 0, 3, 0, 5, 6, 0, 8, 9), nrow = 3)
  
  expect_error(sparse_pcor(mat_dense))
  
})




# sparse_scor()
# ------------------------------------------------------------------------------



test_that("sparse_scor works correctly with valid input", {
  
  set.seed(1)
  
  mat_sparse <- Matrix::rsparsematrix(10, 10, density = 0.5)
  mat_sparse[, 1] <- 0
  colnames(mat_sparse) <- paste0("gene", 1:10)
  
  mat_dense <- as.matrix(mat_sparse)
  
  result_sparse <- sparse_scor(mat_sparse)
  
  result_dense <- 
    suppressWarnings(cor(mat_dense, method = "spearman", use = "everything"))
  
  expect_equal(result_sparse, result_dense, tolerance = 1e-6)
  
})




test_that("sparse_pcor throws an error with invalid input", {
  
  mat_dense <- matrix(c(1, 0, 3, 0, 5, 6, 0, 8, 9), nrow = 3)
  
  expect_error(sparse_pcor(mat_dense))
  
})






# calc_correlation()
# ------------------------------------------------------------------------------



test_that("calc_correlation works for sparse method", {
  
  set.seed(1)
  mat <- Matrix::rsparsematrix(10, 10, density = 0.1)
  result <- calc_correlation(mat, cor_method = "sparse")
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), ncol(mat))
  expect_equal(nrow(result), nrow(mat))
  
})



test_that("calc_correlation works for pearson method", {
  
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  result <- calc_correlation(mat, cor_method = "pearson")
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), ncol(mat))
  expect_equal(nrow(result), nrow(mat))
  
})



test_that("calc_correlation works for spearman method", {
  
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  result <- calc_correlation(mat, cor_method = "spearman")
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), ncol(mat))
  expect_equal(nrow(result), nrow(mat))
  
})



test_that("calc_correlation throws an error for invalid cor_method", {
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  expect_error(calc_correlation(mat, cor_method = "invalid"))
})



test_that("calc_correlation handles NAs correctly for sparse method", {
  
  set.seed(1)
  mat <- Matrix::rsparsematrix(10, 10, density = 0.1)
  mat[, 1] <- 0
  
  result <- calc_correlation(mat, cor_method = "sparse")
  expect_true(is.matrix(result))
  expect_true(any(is.na(result)))
  
})



test_that("calc_correlation throws an error for sparse matrix with dense method", {
  mat <- Matrix::rsparsematrix(10, 10, density = 0.1)
  expect_error(calc_correlation(mat, cor_method = "pearson", ncores = 1))
})



test_that("calc_correlation throws an error for dense matrix with sparse method", {
  mat <- matrix(rnorm(100), nrow = 10)
  expect_error(calc_correlation(mat, cor_method = "sparse", ncores = 1))
})



# set_under_min_count()
# ------------------------------------------------------------------------------



test_that("set_under_min_count sets columns to 0 correctly for sparse matrix", {
  
  mat <- Matrix(c(1, 0, 3, 0, 5, 0, 7, 0, 0, 0), nrow = 5, ncol = 2, sparse = TRUE)
  result <- set_under_min_count(mat, min_count = 2)
  
  expect_true(all(result[, 2] == 0))
  expect_equal(mat[, 1], result[, 1])
  
})



test_that("set_under_min_count handles edge cases for sparse matrix", {
  
  mat <- Matrix(c(1, 2, 3, 4, 5, 1, 2, 3, 4, 0), nrow = 5, ncol = 2, sparse = TRUE)
  result <- set_under_min_count(mat, min_count = 5)

  expect_true(all(result[, 2] == 0))
  expect_error(set_under_min_count(mat, min_count = 6))
  
})



test_that("set_under_min_count handles minimum count of 0 for sparse matrix", {
  
  mat <- Matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1, sparse = TRUE)
  result <- set_under_min_count(mat, min_count = 0)
  expect_false(any(result[, 1] == 0))
  
})



test_that("set_under_min_count checks arguments", {
  
  mat_dense <- matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1)
  mat_sparse <- Matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1, sparse = TRUE)
  
  expect_error(set_under_min_count(mat_dense, min_count = 1))
  expect_error(set_under_min_count(mat_sparse, min_count = -1))
  expect_error(set_under_min_count(mat_sparse, min_count = 10))

})




# subset_celltype_and_filter()
# ------------------------------------------------------------------------------



test_that("subset_celltype_and_filter works with sparse matrix and default settings", {
  
  result <- subset_celltype_and_filter(mat_sparse, meta, "Type1")
 
  expect_true(all(result[, 1] == 0))
  expect_equal(nrow(result), 50)
  expect_equal(ncol(result), 100)
  expect_true(all(rownames(result) %in% meta$ID))
  
})



test_that("subset_celltype_and_filter handles case where no genes meet the min_count", {
  
  mat_dense <- mat_dense
  mat_dense[2:100, ] <- mat_dense[1, ]
  mat_sparse <- Matrix(mat_dense, sparse = TRUE)
  
  result_sparse <- subset_celltype_and_filter(mat_sparse, meta, "Type1")
  
  expect_true(all(result_sparse == 0))
  
})



test_that("subset_celltype_and_filter handles invalid cell_type", {
  expect_error(subset_celltype_and_filter(mat_sparse, meta, cell_type = "InvalidType"))
})



test_that("subset_celltype_and_filter handles incorrect metadata format", {
  invalid_meta <- data.frame(ID = paste0("Cell", 1:10), stringsAsFactors = FALSE)
  expect_error(subset_celltype_and_filter(mat_sparse, invalid_meta, "Type1"))
})



test_that("subset_celltype_and_filter handles arguments", {
  expect_error(subset_celltype_and_filter(mat_dense, meta, "Type1"))
  expect_error(subset_celltype_and_filter(mat_sparse, meta, "Type1", min_count = -1))
  expect_error(subset_celltype_and_filter(mat_sparse, meta, "Type1", min_count = nrow(mat_dense) + 1))
})



test_that("subset_celltype_and_filter checks all cell IDs are in matrix", {
  
  add_mock <- data.frame(ID = paste0("Cell", 101:110),
                         Cell_type = rep("Type1", 10))
  
  meta <- rbind(meta, add_mock)
  
  expect_error(subset_celltype_and_filter(mat_sparse, meta, "Type1"))

})



# Matrix helpers
# ------------------------------------------------------------------------------



# Set NAs in matrix to 0
test_that("na_to_zero works correctly", {
  
  mat <- matrix(c(1, NA, 3, 4), 2, 2)
  result <- na_to_zero(mat)
  expected <- matrix(c(1, 0, 3, 4), 2, 2)
  
  expect_equal(result, expected)
  
})



# Set diag in matrix to 1
test_that("diag_to_one works correctly", {
  
  mat <- matrix(c(1, 2, 3, 4), 2, 2)
  result <- diag_to_one(mat)
  expected <- matrix(c(1, 2, 3, 1), 2, 2)
  
  expect_equal(result, expected)
  
})



# Set upper triangle to NA
test_that("uppertri_to_na works correctly", {
  
  mat <- matrix(c(1, 2, 3, 4), 2, 2)
  result <- uppertri_to_na(mat)
  expected <- matrix(c(1, 2, NA, 4), 2, 2)
  
  expect_equal(result, expected)
  
})



# Replace upper tri of mat with lower tri
test_that("lowertri_to_symm works correctly", {
  
  mat <- matrix(c(1, 2, 3, 4), 2, 2)
  result <- lowertri_to_symm(mat)
  expected <- matrix(c(1, 2, 2, 4), 2, 2)
  
  expect_equal(result, expected)
  
})



# Initiate a matrix of 0s for holding aggregate correlations
test_that("init_agg_mat works correctly", {
  
  mat <- matrix(c(1, 2, 3, 4), 2, 2)
  rownames(mat) <- c("gene1", "gene2")
  result <- init_agg_mat(mat)
  expected <- matrix(0, 2, 2)
  rownames(expected) <- colnames(expected) <- c("gene1", "gene2")
  
  expect_equal(result, expected)
  
})



# Count the NAs in cmat and increment the corresponding indices of na_mat
test_that("increment_na_mat works correctly", {
  
  cmat <- matrix(c(NA, 2, 3, NA), 2, 2)
  na_mat <- matrix(0, 2, 2)
  result <- increment_na_mat(cmat, na_mat)
  expected <- matrix(c(1, 0, 0, 1), 2, 2)
  
  expect_equal(result, expected)
  
})



# transform_correlation_mat()
# ------------------------------------------------------------------------------



test_that("transform_correlation_mat works correctly with allrank", {
  
  mat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  rownames(mat) <- colnames(mat) <- c("gene1", "gene2")
  result <- transform_correlation_mat(mat, "allrank")
  expected <- matrix(c(1, 3, NA, 1), 2, 2)
  rownames(expected) <- colnames(expected) <- c("gene1", "gene2")
  
  expect_equal(result, expected)
  
})



test_that("transform_correlation_mat works correctly with colrank", {
  
  mat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  rownames(mat) <- colnames(mat) <- c("gene1", "gene2")
  result <- transform_correlation_mat(mat, "colrank")
  expected <- matrix(c(1, 2, 2, 1), 2, 2)
  rownames(expected) <- colnames(expected) <- c("gene1", "gene2")
  
  expect_equal(result, expected)
  
})



test_that("transform_correlation_mat works correctly with FZ", {
  
  mat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  rownames(mat) <- colnames(mat) <- c("gene1", "gene2")
  result <- transform_correlation_mat(mat, "FZ")
  expected <- matrix(c(Inf, 0.5493061, 0.5493061, Inf), 2, 2)
  rownames(expected) <- colnames(expected) <- c("gene1", "gene2")
  
  expect_equal(result, expected, tolerance = 1e-6)
  
})



test_that("transform_correlation_mat fails with incorrect arguments", {
  
  mat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  rownames(mat) <- colnames(mat) <- c("gene1", "gene2")
  mat_df <- data.frame(mat)
  mat_badname <- mat
  rownames(mat_badname) <- rev(rownames(mat))
  
  expect_error(transform_correlation_mat(mat, "invalid_method"))
  expect_error(transform_correlation_mat(mat_df, "FZ"))
  expect_error(transform_correlation_mat(mat_badname, "FZ"))
  
})



test_that("transform_correlation_mat handles NA values correctly", {
  
  mat <- matrix(c(1, NA, NA, 1), 2, 2)
  
  rownames(mat) <- colnames(mat) <- c("gene1", "gene2")
  
  result_allrank <- transform_correlation_mat(mat, "allrank")
  expected_allrank <- matrix(c(1, 3, NA, 1), 2, 2)
  rownames(expected_allrank) <- colnames(expected_allrank) <- c("gene1", "gene2")
  
  result_colrank <- transform_correlation_mat(mat, "colrank")
  expected_colrank <- matrix(c(1, 2, 2, 1), 2, 2)
  rownames(expected_colrank) <- colnames(expected_colrank) <- c("gene1", "gene2")
  
  result_fz <- transform_correlation_mat(mat, "FZ")
  expected_fz <- matrix(c(Inf, 0, 0, Inf), 2, 2)
  rownames(expected_fz) <- colnames(expected_fz) <- c("gene1", "gene2")
  
  expect_equal(result_allrank, expected_allrank)
  expect_equal(result_colrank, expected_colrank)
  expect_equal(result_fz, expected_fz)
  
})



# finalize_agg_mat()
# ------------------------------------------------------------------------------



test_that("finalize_agg_mat works correctly with allrank", {
  
  mat <- matrix(c(2, 6, NA, 2), 2, 2)
  rownames(mat) <- colnames(mat) <- c("gene1", "gene2")
  na_mat <- matrix(0, 2, 2)
  result <- finalize_agg_mat(mat, agg_method = "allrank", 2, na_mat = na_mat)
  expected <- matrix(c(1, 0.333333, 0.333333, 1), 2, 2)
  rownames(expected) <- colnames(expected) <- c("gene1", "gene2")
  
  expect_equal(result, expected, tolerance = 1e-6)
  
})



test_that("finalize_agg_mat works correctly with colrank", {
  
  mat <- matrix(c(2, 4, 4, 2), 2, 2)
  rownames(mat) <- colnames(mat) <- c("gene1", "gene2")
  na_mat <- matrix(0, 2, 2)
  result <- finalize_agg_mat(mat, "colrank", 2, na_mat)
  expected <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  rownames(expected) <- colnames(expected) <- c("gene1", "gene2")
  
  expect_equal(result, expected, tolerance = 1e-6)
})



test_that("finalize_agg_mat works correctly with FZ", {
  
  mat <- matrix(c(Inf, 1.098612, 1.098612, Inf), 2, 2)
  rownames(mat) <- colnames(mat) <- c("gene1", "gene2")
  na_mat <- matrix(0, 2, 2)
  result <- finalize_agg_mat(mat, "FZ", 2, na_mat)
  expected <- matrix(c(Inf, 0.5493061, 0.5493061, Inf), 2, 2)
  rownames(expected) <- colnames(expected) <- c("gene1", "gene2")
  
  expect_equal(result, expected, tolerance = 1e-6)
})



# aggregate_celltype_correlation()
# ------------------------------------------------------------------------------



test_that("aggregate_celltype_correlation works for dense matrix with Pearson", {
  
  result_allrank <- aggregate_celltype_correlation(mat_dense, meta, "pearson", "allrank")
  result_colrank <- aggregate_celltype_correlation(mat_dense, meta, "pearson", "colrank")
  result_fz <- aggregate_celltype_correlation(mat_dense, meta, "pearson", "FZ")
  
  expect_true(is.matrix(result_allrank$Agg_mat))
  expect_true(is.matrix(result_colrank$Agg_mat))
  expect_true(is.matrix(result_fz$Agg_mat))
  
  expect_true(is.matrix(result_allrank$NA_mat))
  expect_equal(dim(result_allrank$Agg_mat), dim(result_allrank$NA_mat))
  expect_equal(result_allrank$NA_mat, result_colrank$NA_mat)
  expect_equal(result_allrank$NA_mat, result_fz$NA_mat)
  
  expect_equal(dim(result_allrank$Agg_mat), c(nrow(mat_dense), nrow(mat_dense)))
  expect_equal(dim(result_colrank$Agg_mat), c(nrow(mat_dense), nrow(mat_dense)))
  expect_equal(dim(result_fz$Agg_mat), c(nrow(mat_dense), nrow(mat_dense)))
  
  expect_equal(diag(result_allrank$Agg_mat), diag(result_colrank$Agg_mat))
  

})



test_that("aggregate_celltype_correlation works for sparse matrix with Pearson", {
  
  result_allrank <- aggregate_celltype_correlation(mat_sparse, meta, "sparse", "allrank")
  result_colrank <- aggregate_celltype_correlation(mat_sparse, meta, "sparse", "colrank")
  result_fz <- aggregate_celltype_correlation(mat_sparse, meta, "sparse", "FZ")
  
  expect_true(is.matrix(result_allrank$Agg_mat))
  expect_true(is.matrix(result_colrank$Agg_mat))
  expect_true(is.matrix(result_fz$Agg_mat))
  
  expect_true(is.matrix(result_allrank$NA_mat))
  expect_equal(dim(result_allrank$Agg_mat), dim(result_allrank$NA_mat))
  expect_equal(result_allrank$NA_mat, result_colrank$NA_mat)
  expect_equal(result_allrank$NA_mat, result_fz$NA_mat)
  
  expect_equal(dim(result_allrank$Agg_mat), c(nrow(mat_dense), nrow(mat_dense)))
  expect_equal(dim(result_colrank$Agg_mat), c(nrow(mat_dense), nrow(mat_dense)))
  expect_equal(dim(result_fz$Agg_mat), c(nrow(mat_dense), nrow(mat_dense)))
  
  expect_equal(diag(result_allrank$Agg_mat), diag(result_colrank$Agg_mat))
  
  
})




test_that("aggregate_celltype_correlation handles insufficient counts", {
  
  # Set min_cell to a high number to force insufficient counts
  result <- aggregate_celltype_correlation(mat_dense, meta, "pearson", "allrank", min_cell = 51)
  
  expect_true(is.matrix(result$Agg_mat))
  expect_true(is.matrix(result$NA_mat))
  
  expect_equal(dim(result$NA_mat), c(10, 10))
  expect_true(all(result$NA_mat == 1))
  
})

