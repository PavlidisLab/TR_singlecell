## Generate cell type correlation into a list for inspection of intermediates of
## aggregation.
## NOTE: This is extremely resource intensive. Main factor is number of cell 
## types, as a gene x gene matrix is generated for each into a list, and then
## another list of matched size is created for the correlation ranks.
## -----------------------------------------------------------------------------

source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "GSE165003"
species <- "Mouse" 
corlist_path <- paste0("/space/scratch/amorin/R_objects/", id, "_celltype_cor_list.RDS")


pc <- if (str_to_lower(species) %in% c("human", "hg")) {
  read.delim(ens_hg_path, stringsAsFactors = FALSE)
} else if (str_to_lower(species) %in% c("mouse", "mm")) {
  read.delim(ens_mm_path, stringsAsFactors = FALSE)
} else {
  stop("Species not recognized")
}


dat <- load_dat_list(id)
meta <- dat[[id]]$Meta
mat <- suppressWarnings(as.matrix(dat[[id]]$Mat))


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
    
    message(paste(x, Sys.time()))
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      return(NA)
    }
    
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



# All rank cor list - coerce numeric as was running into integer overflow
rank_l <- lapply(cor_l, function(x) {
  rmat <- allrank_mat(upper_to_na(x), ties_arg = "min")
  # rmat <- allrank_mat(upper_to_na(x), ties_arg = "random")
  rmat <- matrix(as.numeric(rmat), ncol = ncol(rmat), nrow = nrow(rmat))
  rownames(rmat) <- colnames(rmat) <- rownames(x)
  return(rmat)
})


# Sum ranks
agg_mat1 <- Reduce("+", rank_l)

# Final rank of summed ranks                      
agg_mat2 <- allrank_mat(agg_mat1, ties_arg = "min")
# agg_mat2 <- allrank_mat(agg_mat1, ties_arg = "random")

# Standardize ranks and convert from lower tri to symmetric
agg_mat3 <- agg_mat2 / sum(!is.na(agg_mat2))
agg_mat3 <- lowertri_to_symm(agg_mat3)


# Ensure this matches what was previously generated
# agg <- load_agg_mat_list(id, genes = pc$Symbol)[[1]]
# stopifnot(all(abs(agg - agg_mat3) < 1e-5))
# rm(agg)
# gc()


# Inspecting gene pair

gene1 <- "Mef2c"  #  "Runx1"
gene2 <- "Dgcr8"  #  "Ankrd11"   "Rbm5"

agg_df <- data.frame(Symbol = rownames(agg_mat3),
                     RSR = agg_mat3[, gene1],
                     Avg_cor = cor_avg[, gene1]
                     ) %>%
  filter(Symbol != gene1) %>%
  arrange(desc(Avg_cor))


qplot(agg_df, yvar = "Avg_cor", xvar = "RSR")


# The rank of gene 2 in the vector of RSRs for gene 1
rsr_rank <- which(arrange(agg_df, desc(RSR))$Symbol == gene2)

# The raw correlations for gene1-gene2 for each cell type
raw_cor <- vapply(cor_l, function(x) x[gene2, gene1], numeric(1))

# The average raw correlation across cell types
raw_cor_avg <- cor_avg[gene2, gene1]

# The rank of gene2 in gene1's column vector across cell types
# (This ranks just within the gene1 vector, while RSR ranks the entire matrix)
col_rank <- vapply(cor_l, function(x) {
  vec <- x[, gene1]
  vec <- vec[names(vec) != gene1]
  rank_vec <- rank(-vec, ties.method = "min")
  rank_vec[gene2]
}, numeric(1))

# The average column rank
col_rank_avg <- mean(col_rank)
col_rank_med <- median(col_rank)


# The ranks of gene 1 and gene 2 for each cell type
# (This rank is relative to the entire matrix)
all_rank <- vapply(rank_l, function(x) x[gene1, gene2], numeric(1))


cor_df <- data.frame(
  Cell_type = names(col_rank),
  Raw_pcor = round(raw_cor, 4),
  Col_rank = col_rank,
  All_rank = all_rank
)


check_ct <- slice_min(cor_df, Raw_pcor)$Cell_type
all_rank[check_ct] / sum(all_rank[names(all_rank) != check_ct])
