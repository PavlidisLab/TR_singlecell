library(WGCNA)
library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "GSE180928"
out_dir <- file.path("/space/scratch/amorin/TR_singlecell", id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta.RDS"))

allrank_rand_std_path <- file.path(out_dir, paste0(id, "_RSR_allrank.RDS"))
allrank_min_std_path <- file.path(out_dir, paste0(id, "_RSR_allrank_minrank_std.RDS"))
allrank_min_nostd_path <- file.path(out_dir, paste0(id, "_RSR_allrank_minrank_nostd.RDS"))
colrank_rand_std_path <- file.path(out_dir, paste0(id, "_RSR_colrank.RDS"))
colrank_min_std_path <- file.path(out_dir, paste0(id, "_RSR_colrank_minrank_std.RDS"))
colrank_min_nostd_path <- file.path(out_dir, paste0(id, "_RSR_colrank_minrank_nostd.RDS"))

load_df <- data.frame(
  Paths = c(
    allrank_rand_std_path,
    allrank_min_std_path,
    allrank_min_nostd_path,
    colrank_rand_std_path,
    colrank_min_std_path,
    colrank_min_nostd_path
  ),
  ID = c(
    "All_random_std",
    "All_min_std",
    "All_min_nostd",
    "Col_random_std",
    "Col_min_std",
    "Col_min_nostd"
  )
)


mat_l <- load_agg_mat_list(ids = load_df$ID, paths = load_df$Paths, make_symmetric = FALSE)
mat_l[1:3] <- lapply(mat_l[1:3], function(x) lowertri_to_symm(x))


# Check for range of standardized ranks within a column


gene <- "MEF2A"  # ASCL1  RPL3

gene_mat <- gene_vec_to_mat(mat_l, gene)
rank_mat <- colrank_mat(gene_mat, ties_arg = "min")
gene_df <- data.frame(gene_mat) %>% rownames_to_column(var = "Symbol")
rank_df <- data.frame(rank_mat) %>% rownames_to_column(var = "Symbol")
cor(gene_mat, method = "spearman")


# Check for build ups of min ranked values (is this present in aggregate form?)


hist_l <- lapply(load_df$ID, function(x) {
  
  list(
    Top = head(sort(table(gene_df[[x]]), decreasing = TRUE)), 
    Hist = hist(gene_df[[x]], breaks = 100, main = x)
  )
  
})
names(hist_l) <- load_df$ID


plot(gene_df$All_min_std, gene_df$All_min_nostd)
