## This script looks to compare the rankings from coexpression aggregation


# TODO: replace with proper output

saveRDS(
  list(
    Col_std = amat_colrank_std,
    All_std = amat_allrank_std,
    All_nostd = amat_allrank_nostd
  ),
  file = "/space/scratch/amorin/R_objects/Posner2022_subset_RSR_comparison.RDS"
)

gene <- "Rpl37a"

cor(amat_colrank_std[, gene], amat_allrank_std[, gene], method = "spearman", use = "pairwise.complete.obs")
plot(amat_colrank_std[, gene], amat_allrank_std[, gene])


cor(amat_allrank_std[, gene], amat_allrank_nostd[, gene], method = "spearman", use = "pairwise.complete.obs")
plot(amat_allrank_std[, gene], amat_allrank_nostd[, gene])


cor(amat_colrank_std[, gene], amat_allrank_nostd[, gene], method = "spearman", use = "pairwise.complete.obs")
plot(amat_colrank_std[, gene], amat_allrank_nostd[, gene])




#
cmat_col1
rmat_col1
amat_col1

cmat_col2
rmat_col2
amat_col2

amat_col_std
amat_col_nostd

plot(amat_col_std[, 1], amat_col_nostd[, 1])


cmat_all1
rmat_all1
amat_all1

cmat_all2
rmat_all2
amat_all2

amat_all_std
amat_all_nostd

amat_all_std_full <- lowertri_to_symm(amat_all_std)
amat_all_nostd_full <- lowertri_to_symm(amat_all_nostd)

plot(amat_all_std_full[, gene], amat_all_nostd_full[, gene])
plot(amat_all_std_full[, gene], amat_col_std[, gene])

amat_all_std[1:5, 1:5]

cmat_col1[1:5, 1:5]
cmat_all1[1:5, 1:5]
