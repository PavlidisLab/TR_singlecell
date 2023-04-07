## 
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(parallel)
library(pheatmap)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

out_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_oddsratio_cellcounts.RDS"

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

# Ranked targets from paper
rank_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"
rank_l <- readRDS(rank_path)

tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")


# Convert expression matrix to binary expression for quicker look up

count_mat <- GetAssayData(sdat, slot = "counts")
bin_mat <- count_mat > 0


# Calls fisher.test() for a 2x2 table of the count of cells that are
# +/- non-zero counts for the TF and +/- non-zero counts for a gene.

get_cellcount_or <- function(bin_mat, 
                             tf, 
                             gene,
                             tf_label = c("- TF", "+ TF"),
                             target_label = c("- target", "+ target")) {
  
  n_tf_only <- sum(bin_mat[tf, ] & !bin_mat[gene, ]) 
  n_target_only <- sum(!bin_mat[tf, ] & bin_mat[gene, ])
  n_target_and_tf <- sum(bin_mat[tf, ] & bin_mat[gene, ])
  n_none <- sum(!bin_mat[tf, ] & !bin_mat[gene, ])
  
  stopifnot(identical(
    n_tf_only + n_target_only + n_target_and_tf + n_none, ncol(bin_mat)))
  
  cont_mat <- matrix(
    c(n_none, n_tf_only, n_target_only, n_target_and_tf), 
    nrow = 2, ncol = 2,
    dimnames = list("TF expression" = tf_label, "Target" = target_label))
  
  or <- fisher.test(cont_mat)
  
  df <- data.frame(Symbol = gene,
                   Odds_ratio = or$estimate,
                   Lower = or$conf.int[1],
                   Upper = or$conf.int[2],
                   Pval = or$p.value,
                   row.names = NULL)
  
  
  return(df)
}



get_or_list <- function(bin_mat, tf, gene_vec, ...) {
  
  or_list <- mclapply(gene_vec, function(x) {
    
    message(paste(x, Sys.time()))
    
    res <- tryCatch({
      get_cellcount_or(bin_mat, tf, x)
    }, error = function(e) {
      return(NA)
    })
    
    return(res)
    
  }, mc.cores = ncore)
  
  names(or_list) <- input_genes
  
  return(or_list)
}



tf <- "Ascl1"
n_genes <- 100


# input_genes <- rank_l$Mouse[[tf]] %>%
#   filter(Symbol %in% rownames(bin_mat)) %>%
#   filter(Symbol != tf) %>% 
#   slice_head(n = n_genes) %>%
#   pull(Symbol)


input_genes <- rank_l$Mouse[[tf]] %>%
  filter(Symbol %in% rownames(bin_mat)) %>%
  pull(Symbol)



if (!file.exists(out_path)) {
  or_list <- get_or_list(bin_mat, tf, input_genes)
  saveRDS(or_list, out_path)
} else {
  or_list <- readRDS(out_path)
}





or_df <- do.call(rbind, or_list) %>% 
  mutate(Symbol = factor(Symbol, levels = unique(Symbol)))


log_or_df <- data.frame(cbind(
  Symbol = or_df$Symbol, 
  log(dplyr::select(or_df, Odds_ratio, Upper, Lower)),
  Pval = or_df$Pval
))



# ggplot(or_df, aes(x = Symbol, y = Odds_ratio)) +
#   geom_point(size = 2.4, colour = "firebrick4") +
#   geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.4) +
#   geom_hline(yintercept = 1, linetype = "dashed") +
#   ylab("Odds ratio") +
#   xlab("Genes ordered by integrated evidence") +
#   theme_classic() +
#   theme(
#     axis.title = element_text(size = 25),
#     axis.text.x = element_blank(),
#     axis.text.y = element_text(size = 20)
#   )
# 
# 
# 
# ggplot(log_or_df, aes(x = Symbol, y = Odds_ratio)) +
#   geom_point(size = 2.4, colour = "firebrick4") +
#   geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.4) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   ylab("Log odds ratio") +
#   xlab("Genes ordered by integrated evidence") +
#   theme_classic() +
#   theme(
#     axis.title = element_text(size = 25),
#     axis.text.x = element_blank(),
#     axis.text.y = element_text(size = 20)
#   )
# 
# 
# heatmap_df <- log_or_df %>% 
#   mutate(Odds_ratio = ifelse(Pval < 0.05, Odds_ratio, NA)) %>% 
#   arrange(desc(abs(Odds_ratio)))
# 
# 
# pheatmap(heatmap_df$Odds_ratio,
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          na_col = "black",
#          cellwidth = 5,
#          cellheight = 5)
