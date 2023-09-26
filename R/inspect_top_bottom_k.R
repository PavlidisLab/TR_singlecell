source("R/utils/plot_functions.R")



k <- 200
check_k_arg <- TRUE
agg_l <- agg_tf_mm 
msr_mat <- msr_mm 
genes <- tfs_mm$Symbol
x <- "Mecp2"

# For the given gene, create a matrix from the experiments that it is measured
# in. Remove the gene itself to prevent inflated overlap
gene_mat <- gene_vec_to_mat(agg_l, x)
gene_mat <- gene_mat[setdiff(rownames(gene_mat), x), ]
gene_mat <- subset_to_measured(gene_mat, msr_mat = msr_mat, gene = x)
gene_mat_rev <- -gene_mat

# Experiment x experiment similarity matrices
gene_cor <- colwise_cor(gene_mat, cor_method = "spearman")
gene_topk <- colwise_topk_intersect(gene_mat, k = k, check_k_arg = check_k_arg)
gene_btmk <- colwise_topk_intersect(gene_mat_rev, k = k, check_k_arg = check_k_arg)
gene_jacc <- colwise_jaccard(binarize_topk_btmk(gene_mat, k = k))




get_similarity_pair_df <- function(cor_mat, topk_mat, jacc_mat, btmk_mat) {
  
  df <-
    mat_to_df(cor_mat, symmetric = TRUE, value_name = "Scor") %>%
    cbind(
      Topk = mat_to_df(topk_mat, symmetric = TRUE)$Value,
      Bottomk = mat_to_df(btmk_mat, symmetric = TRUE)$Value,
      Jaccard = mat_to_df(jacc_mat, symmetric = TRUE)$Value
    )
  
  return(df)
}




# Data frame of unique dataset pairs and their similarities
sim_df <- get_similarity_pair_df(gene_cor, gene_topk, gene_jacc, gene_btmk)
sim_df <- arrange(sim_df, desc(Bottomk))

# Inspecting the max bottom k pair - ensure negative cors

id1 <- sim_df$Row[1]
id2 <- sim_df$Col[1]

dat <- load_dat_list(c(id1, id2))  # CPM
# dat <- load_dat_list(c(id1, id2), pattern = "_clean_mat_and_meta.RDS")  # lognorm

sort1 <- sort(gene_mat[, id1], decreasing = TRUE)
sort2 <- sort(gene_mat[, id2], decreasing = TRUE)

sort_rev1 <- sort(gene_mat[, id1])
sort_rev2 <- sort(gene_mat[, id2])


# head(sort1, k)
# head(sort2, k)
# head(sort_rev1, k)
# head(sort_rev2, k)


topk1 <- topk_sort(gene_mat[, id1], k = k)
topk2 <- topk_sort(gene_mat[, id2], k = k)
topk_common <- intersect(topk1, topk2)
topk_max <- names(sort(rowMeans(gene_mat[topk_common, c(id1, id2)]), decreasing = TRUE))[1]


btmk1 <- topk_sort(gene_mat_rev[, id1], k = k)
btmk2 <- topk_sort(gene_mat_rev[, id2], k = k)
btmk_common <- intersect(btmk1, btmk2)
btmk_max <- names(sort(rowMeans(gene_mat[btmk_common, c(id1, id2)])))[1]


# This looks at the highest rank (most pos cor) for each and averaged across the paired datasets
cor_max1 <- all_celltype_cor(dat[[id1]]$Mat, dat[[id1]]$Meta, gene1 = x, gene2 = topk1[1])
cor_max2 <- all_celltype_cor(dat[[id2]]$Mat, dat[[id2]]$Meta, gene1 = x, gene2 = topk2[1])
cor_max <- all_celltype_cor(dat[[id2]]$Mat, dat[[id2]]$Meta, gene1 = x, gene2 = topk_max)
p_max1 <- all_celltype_scatter(dat[[id1]]$Mat, dat[[id1]]$Meta, gene1 = x, gene2 = topk_max)[names(cor_max1)]
p_max2 <- all_celltype_scatter(dat[[id2]]$Mat, dat[[id2]]$Meta, gene1 = x, gene2 = topk_max)[names(cor_max2)]


# This looks at the lowest rank (most negative cor) for the paired datasets
cor_min1 <- all_celltype_cor(dat[[id1]]$Mat, dat[[id1]]$Meta, gene1 = x, gene2 = btmk1[1])
cor_min2 <- all_celltype_cor(dat[[id2]]$Mat, dat[[id2]]$Meta, gene1 = x, gene2 = btmk2[1])
cor_min <- all_celltype_cor(dat[[id2]]$Mat, dat[[id2]]$Meta, gene1 = x, gene2 = btmk_max)
p_min1 <- all_celltype_scatter(dat[[id1]]$Mat, dat[[id1]]$Meta, gene1 = x, gene2 = btmk_max)[names(cor_min1)]
p_min2 <- all_celltype_scatter(dat[[id2]]$Mat, dat[[id2]]$Meta, gene1 = x, gene2 = btmk_max)[names(cor_min2)]


# Inspect the first element before/after ties - "weakest" cors that are not 0

k1 <- check_k(sort1, k = 5e3)
k2 <- check_k(sort2, k = 10e3)

k_rev1 <- check_k(sort_rev1, k = 17e3)
k_rev2 <- check_k(sort_rev2, k = 10e3)

sort1[(k1 - 2):(k1 + 2)]
sort2[(k2 - 2):(k2 + 2)]

sort_rev1[(k_rev1 - 2):(k_rev1 + 2)]
sort_rev2[(k_rev2 - 2):(k_rev2 + 2)]


# Last pos cor before ties - weakest pos cor
last_pos1 <- all_celltype_cor(dat[[id1]]$Mat, dat[[id1]]$Meta, gene1 = x, gene2 = names(sort1)[k1])
last_pos2 <- all_celltype_cor(dat[[id2]]$Mat, dat[[id2]]$Meta, gene1 = x, gene2 = names(sort2)[k2])

p_lp1 <- all_celltype_scatter(dat[[id1]]$Mat, dat[[id1]]$Meta, gene1 = x, gene2 = names(sort1)[k1])[names(last_pos1)]
p_lp2 <- all_celltype_scatter(dat[[id2]]$Mat, dat[[id2]]$Meta, gene1 = x, gene2 = names(sort2)[k2])[names(last_pos2)]

# First neg cor after ties - weakest neg cor
first_neg1 <- all_celltype_cor(dat[[id1]]$Mat, dat[[id1]]$Meta, gene1 = x, gene2 = names(sort_rev1)[k_rev1])
first_neg2 <- all_celltype_cor(dat[[id2]]$Mat, dat[[id2]]$Meta, gene1 = x, gene2 = names(sort_rev2)[k_rev2])

p_fn1 <- all_celltype_scatter(dat[[id1]]$Mat, dat[[id1]]$Meta, gene1 = x, gene2 = names(sort_rev1)[k_rev1])[names(first_neg1)]
p_fn2 <- all_celltype_scatter(dat[[id2]]$Mat, dat[[id2]]$Meta, gene1 = x, gene2 = names(sort_rev2)[k_rev2])[names(first_neg2)]



plot(sort1, 
     main = paste(x, id1),
     xlab = "Sorted genes",
     ylab = "Standardized RSR",
     cex.axis = 1.5,
     cex.lab = 1.5)

abline(v = k1, col = "red")
abline(v = which(names(sort1) == names(sort_rev1[k_rev1])), col = "red")



plot(sort2, 
     main = paste(x, id2),
     xlab = "Sorted genes",
     ylab = "Standardized RSR",
     cex.axis = 1.5,
     cex.lab = 1.5)

abline(v = k2, col = "red")
abline(v = which(names(sort2) == names(sort_rev2[k_rev2])), col = "red")

