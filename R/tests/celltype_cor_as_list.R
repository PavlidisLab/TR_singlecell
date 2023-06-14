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
allrank_path <- file.path(out_dir, paste0(id, "_RSR_allrank.RDS"))
corlist_path <- file.path(out_dir, paste0(id, "_celltype_corlist.RDS"))
zcor_path <- file.path(out_dir, paste0(id, "_fishersZ.RDS"))


dat <- readRDS(processed_path)
meta <- dat[[2]]
mat <- as.matrix(dat[[1]])
rsr <- readRDS(allrank_path)
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



# Rank cor list
rank_l <- lapply(cor_l, function(x) allrank_mat(upper_to_na(x)))
# rank_l2 <- lapply(cor_l, function(x) allrank_mat(upper_to_na(x)))
# rank_l_cp <- rank_l
# rank_l <- rank_l2




# Aggregate rank + standardize
agg_mat <- Reduce("+", rank_l)
rank_mat <- allrank_mat(agg_mat)
srank_mat <- rank_mat / sum(!is.na(rank_mat))


# Get the average cor across cell types
cor_avg <- Reduce("+", cor_l) / length(cor_l)


# Compare average cor to standardized rank generated here and the previously
# generated RSR mat


srank_df <- srank_mat %>% 
  mat_to_df(symmetric = TRUE) %>% 
  dplyr::rename(Srank = Value) %>% 
  mutate(Rank_srank = rank(-Srank, ties.method = "random"))


rsr_df <- rsr %>% 
  mat_to_df(symmetric = TRUE) %>% 
  dplyr::rename(RSR = Value) %>% 
  mutate(Rank_RSR = rank(-RSR, ties.method = "random"))


cor_df <- cor_avg %>% 
  mat_to_df(symmetric = TRUE) %>% 
  dplyr::rename(Avg_cor = Value) %>% 
  mutate(Rank_cor = rank(-Avg_cor, ties.method = "random"))


# cor(rsr_df$Avg_cor, cor_df$RSR, method = "spearman")


all_df <- cbind(rsr_df,
                cor_df[, c("Avg_cor", "Rank_cor")],
                srank_df[, c("Srank", "Rank_srank")],
                Count_NA = mat_to_df(na_mat, symmetric = TRUE)$Value)


top_df <- all_df %>% 
  filter(Rank_cor < 10e3 | Rank_RSR < 10e3 | Rank_srank < 10e3) %>% 
  mutate(Diff_RSR_Srank = abs(Rank_srank - Rank_RSR),
         Diff_RSR_cor = abs(Rank_srank - Rank_cor))


# Find relationship between the difference in ranks and the count of NAs

boxplot(top_df$Diff_RSR_cor ~ top_df$Count_NA)


# Inspect most different rank with all measurements (no NAs)

most_diff <- top_df %>% filter(Count_NA == 0) %>% slice_max(Diff_RSR_cor)


gene1 <- most_diff$Row  # "RHOXF2"  "FAM72A"
gene2 <- most_diff$Col   # "RHOXF2B"  "FAM72B"
all_celltype_cor(mat, meta, gene1, gene2)
mean(all_celltype_cor(mat, meta, gene1, gene2))
mean(unlist(lapply(cor_l, function(x) x[gene1, gene2])))
unlist(lapply(rank_l, function(x) x[gene1, gene2]))



# Comparing the top values for a given TF

tf <- "ASCL1"

tf_df <- data.frame(
  Symbol = rownames(cor_avg),
  Cor = cor_avg[, tf],
  RSR = lowertri_to_symm(rsr)[, tf],
  Srank = lowertri_to_symm(srank_mat)[, tf],
  Count_NA = na_mat[, tf]
)


tf_df$Rank_cor <- rank(-tf_df$Cor, ties.method = "random")
tf_df$Rank_RSR <- rank(-tf_df$RSR, ties.method = "random")
tf_df$Rank_Srank <- rank(-tf_df$Srank, ties.method = "random")


tf_df %>% 
  filter(Rank_RSR < 2e3 & Rank_Srank < 2e3) %>% 
  ggplot(aes(x = Rank_RSR, y = Rank_Srank)) +
  geom_point() + 
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))
