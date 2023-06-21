## Generate cell type correlation into a list for inspection of intermediates of
## aggregation.
## -----------------------------------------------------------------------------

source("R/00_config.R")
source("R/utils/functions.R")

id <- "GSE180928"
species <- "Human" 

sc_dir <- file.path("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datasets", id)
dat_path <- file.path(sc_dir, paste0(id, "_cellxgene_seurat.RDS"))
out_dir <- file.path("/space/scratch/amorin/TR_singlecell", id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta.RDS"))
colrank_path <- file.path(out_dir, paste0(id, "_RSR_colrank_minrank_std.RDS"))
corlist_path <- file.path(out_dir, paste0(id, "_celltype_corlist.RDS"))
zcor_path <- file.path(out_dir, paste0(id, "_fishersZ.RDS"))


dat <- readRDS(processed_path)
meta <- dat[[2]]
mat <- as.matrix(dat[[1]])
rsr <- readRDS(colrank_path)
zcor <- readRDS(zcor_path)
na_mat <- zcor$NA_mat


# List of raw correlation matrices

get_cor_list <- function(mat, meta, ncores = 1) {
  
  cts <- unique(meta$Cell_type)
  
  # Split mat into list of cell type matrices, setting low count genes to NA,
  # calculating gene-cor cor for each cell type, and setting NA cors to 0
  
  ct_l <- mclapply(cts, function(x) {
    
    ct_mat <- subset_and_filter(mat = mat, 
                                meta = meta, 
                                cell_type = x, 
                                min_count = 20)
    
    cor_mat <- ct_mat %>%
      get_cor_mat(lower_tri = FALSE) %>%
      na_to_zero() %>%
      diag_to_one()
    
    return(cor_mat)
    
  }, mc.cores = ncores)
  
  names(ct_l) <- cts
  
  return(ct_l)
}



if (!file.exists(corlist_path)) {
  cor_l <- get_cor_list(mat = mat, meta = meta, ncores = ncore)
  saveRDS(cor_l, corlist_path)
} else {
  cor_l <- readRDS(corlist_path)
}


# Get the average cor across cell types
cor_avg <- Reduce("+", cor_l) / length(cor_l)


# Rank cor list + aggregate: all min rank
rank_l_all <- lapply(cor_l, function(x) allrank_mat(upper_to_na(x), ties_arg = "min"))
agg_mat_all <- Reduce("+", rank_l_all)
rank_mat_all <- allrank_mat(agg_mat_all, ties_arg = "min")
srank_mat_all <- rank_mat_all / sum(!is.na(rank_mat_all))


# Rank cor list + aggregate: col min rank
rank_l_col <- lapply(cor_l, function(x) colrank_mat(x, ties_arg = "min"))
agg_mat_col <- Reduce("+", rank_l_col)
rank_mat_col <- colrank_mat(agg_mat_col, ties_arg = "min")
srank_mat_col <- rank_mat_col / nrow(rank_mat_col)


# Ensure this matches what was previously generated
identical(srank_mat_col, rsr)



# Compare average cor to standardized rank/RSR, including count of NA pairs

set.seed(3)
rsr_df_all <- srank_mat_all %>% 
  mat_to_df(symmetric = TRUE) %>% 
  dplyr::rename(RSR_all = Value) %>% 
  mutate(Rank_RSR_all = rank(-RSR_all, ties.method = "min"))


# Ranking doesn't make sense in this context as RSRs are independent per column
set.seed(3)
rsr_df_col <- srank_mat_col %>% 
  mat_to_df(symmetric = TRUE) %>% 
  dplyr::rename(RSR_col = Value)
  # mutate(Rank_RSR_col = rank(-RSR_col, ties.method = "min"))


set.seed(3)
cor_df <- cor_avg %>% 
  mat_to_df(symmetric = TRUE) %>% 
  dplyr::rename(Avg_cor = Value) %>% 
  mutate(Rank_cor = rank(-Avg_cor, ties.method = "min"))


# cor(rsr_df$Avg_cor, cor_df$RSR_rand, method = "spearman")


all_df <- cbind(rsr_df_all,
                rsr_df_col[, c("RSR_col", "Rank_RSR_col")],
                cor_df[, c("Avg_cor", "Rank_cor")],
                Count_NA = mat_to_df(na_mat, symmetric = TRUE)$Value)


rm(rsr_df_all, rsr_df_col, cor_df)
gc()


# Comparing the top values for a given gene

gene <- "ZSWIM6"

gene_df <- data.frame(
  Symbol = rownames(cor_avg),
  Avg_cor = cor_avg[, gene],
  RSR_all = lowertri_to_symm(srank_mat_all)[, gene],
  RSR_col = srank_mat_col[, gene],
  Count_NA = na_mat[, gene]
)


gene_df$Rank_cor <- rank(-gene_df$Avg_cor, ties.method = "min")
gene_df$Rank_RSR_all <- rank(-gene_df$RSR_all, ties.method = "min")
gene_df$Rank_RSR_col <- rank(-gene_df$RSR_col, ties.method = "min")


gene_df %>% 
  filter(Rank_RSR_all < 2e3 & Rank_RSR_col < 2e3) %>% 
  ggplot(aes(x = Rank_RSR_all, y = Rank_RSR_col)) +
  geom_point() + 
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))



gene1 <- gene   
gene2 <- "BRK1"   #  PRDX5  BRK1
 
unlist(lapply(cor_l, function(x) x[gene1, gene2]))
unlist(lapply(rank_l_all, function(x) x[gene1, gene2]))
unlist(lapply(rank_l_col, function(x) x[gene1, gene2]))

