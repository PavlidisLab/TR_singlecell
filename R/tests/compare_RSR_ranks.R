## This script looks to compare the rankings from coexpression aggregation
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "GSE180928"
dat <- readRDS(file.path("/space/scratch/amorin/TR_singlecell", id, paste0(id, "_clean_mat_and_meta.RDS")))
rsr_all <- readRDS(file.path("/space/scratch/amorin/TR_singlecell", id, paste0(id, "_RSR_allrank.RDS")))
rsr_all <- lowertri_to_symm(rsr_all)  # RSR_all lower triangle of ranks, make symmetric for col subsetting
rsr_col <- readRDS(file.path("/space/scratch/amorin/TR_singlecell", id, paste0(id, "_RSR_colrank.RDS")))
zcor <- readRDS(file.path("/space/scratch/amorin/TR_singlecell", id, paste0(id, "_fishersZ.RDS")))
na_mat <- zcor$NA_mat
mat <- dat[[1]]
meta <- dat[[2]]
genes <- rownames(rsr_all)
cts <- unique(meta$Cell_type)

stopifnot(identical(rownames(rsr_all), rownames(rsr_col)))
stopifnot(identical(rownames(rsr_all), rownames(zcor$Avg_all)))
stopifnot(identical(rownames(rsr_all), rownames(zcor$NA_mat)))


med_na <- apply(na_mat, 1, median)
min_na <- min(na_mat)
max_na <- max(na_mat)
which_min_na <- which(na_mat == min_na, arr.ind = TRUE)


cor_l <- vapply(genes, function(x) {
  cor(rsr_all[, x], rsr_col[, x], use = "pairwise.complete.obs", method = "spearman")
}, numeric(1))

cor_l <- sort(unlist(cor_l))


cor_summ <- summary(cor_l)
hist(cor_l, breaks = 100, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, main = "Agreement between column and all ranking")
best_cor <- cor_l[length(cor_l)] 
worst_cor <- cor_l[1] # comparing noise


# This is awful and slow, but the idea is to get a tally of cell types that
# the given gene has a minimum counts

expr_ct_counts <- parallel::mclapply(genes, function(x) {
  
  counts <- vapply(cts, function(y) {
    sum(mat[x, meta$Cell_type == y] != 0) >= 20
  }, logical(1))
  
  return(sum(counts))
  
}, mc.cores = ncore)
names(expr_ct_counts) <- genes


# And here, averaging the count of NA pairs for each gene

na_mean <- rowMeans(na_mat)

na_df <- data.frame(
  Symbol = genes,
  NA_mean = na_mean[genes],
  RSR_cor = cor_l[genes],
  Expr_counts = unlist(expr_ct_counts)
)


px1 <- ggplot(na_df, aes(x = Expr_counts, y = RSR_cor)) +
  geom_point(shape = 21, size = 2.5) +
  xlab("Count of cell types with at least 20 expressing cells") +
  ylab("Spearman's correlation between column and all ranking") +
  ggtitle("Agreement between column and all ranking") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))


plot(na_df$NA_mean, na_df$RSR_cor)
plot(na_df$Expr_counts, na_df$RSR_cor)

# Lowest cor between RSR all/col among maximally expressed genes

diff_cor <- na_df %>% 
  filter(Expr_counts == max(Expr_counts)) %>% 
  slice_min(RSR_cor, n = 1)


# Organize the ranks for a given gene

# gene <- names(worst_cor)
# gene <- names(best_cor)
# gene <- diff_cor$Symbol
gene <- "ASCL1"


rank_df <- data.frame(
  Symbol = rownames(rsr_all),
  RSR_all = rsr_all[, gene],
  RSR_col = rsr_col[, gene],
  Zcor_all = zcor$Avg_all[, gene],
  Zcor_nonNA = zcor$Avg_nonNA[, gene],
  N_na = na_mat[, gene]
)


cor(select_if(rank_df, is.numeric), use = "pairwise.complete.obs", method = "spearman")


px2 <- ggplot(rank_df, aes(x = RSR_all, y = RSR_col)) +
  geom_point(shape = 21, size = 2.5) +
  xlab("All ranking") +
  ylab("Column ranking") +
  ggtitle(gene) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))
  


# Here, splitting the ranks into 10 bins, decreasing by column or all rank. 
# Then correlating the ranks in each bin to show how the random ranking
# of NA->0 cor appears

n_bins <- 50
a1 <- arrange(rank_df, desc(RSR_all)) %>% filter(Symbol != gene)
a2 <- split(a1, cut(1:nrow(a1), n_bins))
a3 <- unlist(lapply(a2, function(x) cor(x$RSR_all, x$RSR_col, method = "spearman")))


# Assume that dip is noise from random ranking of NA cor -> 0, while raise in
# last bin is for preserved negative correlation. Or is it ~all NAs resulting
# in blocky similarity?


px3 <- data.frame(Scor = a3, k = 1:n_bins) %>% 
  ggplot(aes(x = k, y = Scor)) +
  geom_point(shape = 21, size = 2.5) +
  xlab("K bins in decreasing order of all ranking") +
  ylab("Spearman's correlation between column and all ranking") +
  ggtitle(gene) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))




# looking at tf


gene1 <- "HES1"
na_mean[gene1]
gene1_rank_all <- sort(rsr_all[gene1, ])
gene1_rank_all <- gene1_rank_all[names(gene1_rank_all) != gene1]
gene2 <- names(gene1_rank_all[1])
# gene2 <- names(gene1_rank_all[length(gene1_rank_all)])
ct_cors <- all_celltype_cor(mat, meta, gene1 = gene1, gene2 = gene2)
ct_maxcor <- names(ct_cors[1])
ct_mincor <- names(ct_cors[length(ct_cors)])
plot(mat[gene1, meta$Cell_type == ct_maxcor], mat[gene2, meta$Cell_type == ct_maxcor])
plot(mat[gene1, meta$Cell_type == ct], mat[gene2, meta$Cell_type == ct])


# Explicit rank cutoff and cor

k <- nrow(a2[[1]])
sub_all <- slice_max(rank_df, RSR_all, n = k)
sub_col <- slice_max(rank_df, RSR_col, n = k)
cor(select_if(sub_all, is.numeric), method = "spearman")
cor(select_if(sub_col, is.numeric), method = "spearman")
plot(sub_all$RSR_all, sub_all$RSR_col)
plot(sub_col$RSR_all, sub_col$RSR_col)


# Looking at rank correlation of a given gene relative to the rest of gene
# vectors within all rank. This is distinct than the best aggregate cor pair.
# It considers which gene's rankings most align with the given gene.

top_cor <- lapply(genes, function(x) cor(rsr_all[, gene], rsr_all[, x], method = "spearman"))
names(top_cor) <- genes
top_cor <- sort(unlist(top_cor), decreasing = TRUE)
top_cor <- top_cor[names(top_cor) != gene]
rsr_all[gene, names(top_cor[1])]
plot(rsr_all[, gene], rsr_all[, names(top_cor[1])])

tt <- data.frame(Symbol = names(top_cor),
                 Gene_vec_cor = top_cor, 
                 Aggregate_rank = rsr_all[names(top_cor), gene])
