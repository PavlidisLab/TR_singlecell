## Tests and sanity checks for functions used to compare ranked vectors
## -----------------------------------------------------------------------------

library(assertthat)
source("R/utils/vector_comparison_functions.R")

k <- 1000 # topk cutoff
genes <- rownames(agg_hg[[1]])

# Define a test gene to focus on, as well as two datasets that will be used for
# copying vectors to check for matches

gene <- "RPL3"
ds1 <- "GSE180928"
ds2 <- "GSE216019"

mat1 <- agg_hg[[ds1]]
mat2 <- agg_hg[[ds2]]

vec1 <- mat1[, gene]
vec2 <- mat2[, gene]

assert_that(!are_equal(vec1, vec2))


# Here, create a copy of a test mat, copying a test gene vector into another 
# dataset to look for perfect matches

mat_copy <- mat1
mat_copy[, gene] <- vec2

agg_hg_copy <- agg_hg
agg_hg_copy[[ds1]] <- mat_copy

assert_that(are_equal(agg_hg_copy[[ds1]][, gene], vec2))
assert_that(!are_equal(mat1[, gene], vec2))


# Now, create a test mat where the order of the test gene is reversed to ensure
# minimal similarity

vec_sort <- sort(vec1, decreasing = TRUE)
vec_rev <- sort(vec1, decreasing = FALSE)
names(vec_rev) <- names(vec_sort)
vec_rev <- vec_rev[genes]

mat_rev <- mat2
mat_rev[, gene] <- vec_rev

agg_hg_rev <- agg_hg
agg_hg_rev[[ds2]] <- mat_rev

assert_that(are_equal(agg_hg[[ds1]][, gene], vec1))
assert_that(!are_equal(agg_hg_rev[[ds2]][, gene], vec2))



# Gene similarity matrix to generate raw similarity of given gene across all
# experiment pairs. Ensure that the copy results in identical results
# ------------------------------------------------------------------------------


gsm_topk_copy <- gene_similarity_matrix(
  agg_l = agg_hg_copy,
  gene = gene,
  msr = "Topk",
  k = k,
  ncores = ncore
)


gsm_auprc_copy <- gene_similarity_matrix(
  agg_l = agg_hg_copy,
  gene = gene,
  msr = "AUPRC",
  k = k,
  ncores = ncore
)


gsm_cor_copy <- gene_similarity_matrix(
  agg_l = agg_hg_copy,
  gene = gene,
  msr = "Cor",
  k = k,
  ncores = ncore
)



assert_that(are_equal(gsm_topk_copy[ds1, ds2], gsm_topk_copy[ds1, ds1], k))
assert_that(are_equal(gsm_cor_copy[ds1, ds2], gsm_cor_copy[ds1, ds1], 1))
assert_that(are_equal(gsm_auprc_copy[ds1, ds2], gsm_auprc_copy[ds1, ds1], 1))


# Now reverse, to ensure minimal similarity. Not sure what this entails for 
# AUPRC, so instead check that the value is lower than the baseline (# pos/total)


gsm_topk_rev <- gene_similarity_matrix(
  agg_l = agg_hg_rev,
  gene = gene,
  msr = "Topk",
  k = k,
  ncores = ncore
)


gsm_auprc_rev <- gene_similarity_matrix(
  agg_l = agg_hg_rev,
  gene = gene,
  msr = "AUPRC",
  k = k,
  ncores = ncore
)


gsm_cor_rev <- gene_similarity_matrix(
  agg_l = agg_hg_rev,
  gene = gene,
  msr = "Cor",
  k = k,
  ncores = ncore
)


assert_that(are_equal(gsm_topk_rev[ds1, ds2], 0))
assert_that(are_equal(gsm_cor_rev[ds1, ds2], -1))
assert_that(gsm_auprc_rev[ds1, ds2] < (k / length(genes)))

