## TODO: fix save
## TODO: need a min count filter
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(parallel)
library(pheatmap)
library(cowplot)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

out_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_oddsratio_cellcounts"

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



get_or_list <- function(bin_mat, tf, gene_vec, ncores = 1) {
  
  or_list <- mclapply(gene_vec, function(x) {
    
    res <- tryCatch({
      get_cellcount_or(bin_mat, tf, x)
    }, error = function(e) {
      return(NA)
    })
    
    return(res)
    
  }, mc.cores = ncores)
  
  names(or_list) <- gene_vec
  
  return(or_list)
}



if (!file.exists(out_path)) {

  all_or_list <- mclapply(tfs, function(x) {

    message(paste(x, Sys.time()))
    
    input_genes <- rank_l$Mouse[[x]] %>%
      filter(Symbol %in% rownames(bin_mat)) %>%
      pull(Symbol)

    or_list <- get_or_list(bin_mat, x, input_genes, ncores = 4)
    
    saveRDS(or_list, paste0(out_path, "_", x, ".RDS"))
    
    return(or_list)

  }, mc.cores = 4)

  names(all_or_list) <- tfs

  saveRDS(all_or_list, paste0(out_path, "_all.RDS"))

} else {

  or_list <- readRDS(paste0(out_path, "_all.RDS"))

}



# Post hoc count filter (should be done earlier)

keep_gene <- apply(bin_mat, 1, function(x) sum(x != 0) >= 20)
keep_gene <- names(keep_gene[keep_gene])


or_df_list <- lapply(tfs, function(x) {
  
  data.frame(do.call(rbind, or_list[[x]])) %>% 
    mutate(Log_OR = log(Odds_ratio),
           Log_upper = log(Upper),
           Log_lower = log(Lower)) %>% 
    left_join(rank_l$Mouse[[x]], by = "Symbol") %>% 
    filter(Symbol %in% keep_gene & Symbol != x) %>% 
    mutate(Top500 = Rank_integrated <= 500)

})
names(or_df_list) <- tfs



# Cherry pick the top OR that also has annotated evidence


or_boxplot <- function(df, xvar, xlab) {

  ggplot(df, aes(x = !!sym(xvar), y = Log_OR)) +
    geom_point(size = 2.4, colour = "firebrick4") +
    geom_errorbar(aes(ymin = Log_lower, ymax = Log_upper), width = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Log odds ratio") +
    xlab(xlab) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 25),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 20))
}


# Density plot of log OR by group status


or_densplot <- function(df, group) {
  
  ggplot(df, aes(x = Log_OR, fill = !!sym(group))) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ylab("Density") +
    xlab("Log odds ratio") +
    scale_fill_manual(values = c("lightgrey", "orchid4")) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 25),
      axis.text = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20))
}




plot_l <- lapply(or_df_list, function(x) {
  
  p1_df <- x %>%
    filter(Curated_target & Top500) %>%
    mutate(Symbol = factor(Symbol, levels = unique(Symbol)))
  
  p1 <- or_boxplot(p1_df, xvar = "Symbol", xlab = "Symbol")
  p2 <- or_densplot(x, group = "Curated_target")
  p3 <- or_densplot(x, group = "Top500")

  return(list(P1 = p1, P2 = p2, P3 = p3))  
})



left <- plot_l$Ascl1$P1 + ylim(-5, 2)
right <- plot_grid(plot_l$Ascl1$P2, plot_l$Ascl1$P3, nrow = 2)
plot_grid(left, right, nrow = 1, rel_widths = c(1, 0.75))


wilx_l <- lapply(or_df_list, function(x) {
  c(Curated_target = wilcox.test(x$Log_OR ~ x$Curated_target)$p.value,
    Top500 = wilcox.test(x$Log_OR ~ x$Top500)$p.value)
})



or_df_list$Ascl1 %>% 
  filter(Top500 | Curated_target) %>% 
  slice_max(Odds_ratio, n = 5)
  

or_df_list$Ascl1 %>% 
  filter((Top500 | Curated_target) & Log_OR != -Inf) %>% 
  slice_min(Odds_ratio, n = 5)


FeaturePlot(sdat, c("Ascl1", "Traf4", "Igf2"), ncol = 3)
