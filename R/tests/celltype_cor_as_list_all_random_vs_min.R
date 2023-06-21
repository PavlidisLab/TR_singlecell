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


# Get the average cor across cell types
cor_avg <- Reduce("+", cor_l) / length(cor_l)


# Rank cor list + aggregate: random rank default seed
rank_l_rand <- lapply(cor_l, function(x) allrank_mat(upper_to_na(x)))
agg_mat_rand <- Reduce("+", rank_l_rand)
rank_mat_rand <- allrank_mat(agg_mat_rand)
srank_mat_rand <- rank_mat_rand / sum(!is.na(rank_mat_rand))


# Ensure this matches what was previously generated
identical(srank_mat_rand, rsr)


# Rank cor list + aggregate: min rank
rank_l_min <- lapply(cor_l, function(x) allrank_mat(upper_to_na(x), ties_arg = "min"))
agg_mat_min <- Reduce("+", rank_l_min)
rank_mat_min <- allrank_mat(agg_mat_min, ties_arg = "min")
srank_mat_min <- rank_mat_min / sum(!is.na(rank_mat_min))


# Compare average cor to standardized rank/RSR, including count of NA pairs

set.seed(3)
rsr_df_rand <- srank_mat_rand %>% 
  mat_to_df(symmetric = TRUE) %>% 
  dplyr::rename(RSR_rand = Value) %>% 
  mutate(Rank_RSR_rand = rank(-RSR_rand, ties.method = "min"))


set.seed(3)
rsr_df_min <- srank_mat_min %>% 
  mat_to_df(symmetric = TRUE) %>% 
  dplyr::rename(RSR_min = Value) %>% 
  mutate(Rank_RSR_min = rank(-RSR_min, ties.method = "min"))


set.seed(3)
cor_df <- cor_avg %>% 
  mat_to_df(symmetric = TRUE) %>% 
  dplyr::rename(Avg_cor = Value) %>% 
  mutate(Rank_cor = rank(-Avg_cor, ties.method = "min"))


# cor(rsr_df$Avg_cor, cor_df$RSR_rand, method = "spearman")


all_df <- cbind(rsr_df_rand,
                rsr_df_min[, c("RSR_min", "Rank_RSR_min")],
                cor_df[, c("Avg_cor", "Rank_cor")],
                Count_NA = mat_to_df(na_mat, symmetric = TRUE)$Value)


rm(rsr_df_rand, rsr_df_min, cor_df)
gc()



top_df <- all_df %>% 
  filter(Rank_cor < 10e3 | Rank_RSR_rand < 10e3 | Rank_RSR_min < 10e3) %>% 
  mutate(Diff_rand_min = abs(Rank_RSR_min - Rank_RSR_rand),
         Diff_rand_cor = abs(Rank_RSR_rand - Rank_cor),
         Diff_min_cor = abs(Rank_RSR_min - Rank_cor))



# Find relationship between the difference in ranks and the count of NAs

ggplot(top_df, aes(x = as.factor(Count_NA), y = Diff_rand_cor)) +
  geom_boxplot(width = 0.3) +
  xlab("Count of NAs") +
  ylab("Difference in ranks") +
  ggtitle("Effect of NAs on difference in average cor and RSR ranks") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20))


top_df %>% 
  filter(Count_NA == 0) %>% 
  ggplot(., aes(x = RSR_rand, y = RSR_min)) +
  geom_point()

  

# Inspect most different rank with all measurements (no NAs)

most_diff <- top_df %>% filter(Count_NA == 0) %>% slice_max(Diff_rand_cor)


gene1 <- most_diff$Row  # "RHOXF2"  "FAM72A"
gene2 <- most_diff$Col   # "RHOXF2B"  "FAM72B"

unlist(lapply(cor_l, function(x) x[gene1, gene2]))
unlist(lapply(rank_l_rand, function(x) x[gene1, gene2]))
unlist(lapply(rank_l_min, function(x) x[gene1, gene2]))

rank_mat_min[gene1, gene2]
rank_mat_rand[gene1, gene2]


agg_mat_rand[gene1, gene2]
agg_mat_min[gene1, gene2]



head(sort(table(agg_mat_rand)))


# Comparing the top values for a given TF

tf <- "ASCL1"

tf_df <- data.frame(
  Symbol = rownames(cor_avg),
  Avg_cor = cor_avg[, tf],
  RSR = lowertri_to_symm(rsr)[, tf],
  RSR2 = lowertri_to_symm(srank_mat_min)[, tf],
  Count_NA = na_mat[, tf]
)


tf_df$Rank_cor <- rank(-tf_df$Avg_cor, ties.method = "min")
tf_df$Rank_RSR <- rank(-tf_df$RSR, ties.method = "min")
tf_df$Rank_RSR2 <- rank(-tf_df$RSR2, ties.method = "min")


tf_df %>% 
  filter(Rank_RSR < 2e3 & Rank_RSR2 < 2e3) %>% 
  ggplot(aes(x = Rank_RSR, y = Rank_RSR2)) +
  geom_point() + 
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))
