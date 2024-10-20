## TODO
## -----------------------------------------------------------------------------

library(clusterProfiler)
library(parallel)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
source("R/00_config.R")
source("R/utils/plot_functions.R")


# go_coexpr_hg_path <- "/space/scratch/amorin/R_objects/coexpr_goenrich_hg.RDS"
# go_coexpr_mm_path <- "/space/scratch/amorin/R_objects/coexpr_goenrich_mm.RDS"
# go_coexpr_ortho_path <- "/space/scratch/amorin/R_objects/coexpr_goenrich_ortho.RDS"
# score_col <- "p.adjust"
# term_col <- "Description"

# go_coexpr_hg_path <- "/space/scratch/amorin/R_objects/coexpr_erminer_hg.RDS"
# go_coexpr_mm_path <- "/space/scratch/amorin/R_objects/coexpr_erminer_mm.RDS"
# go_coexpr_ortho_path <- "/space/scratch/amorin/R_objects/coexpr_erminer_ortho.RDS"

go_coexpr_hg_path <- "/space/scratch/amorin/R_objects/coexpr_erminer_pr_hg.RDS"
go_coexpr_mm_path <- "/space/scratch/amorin/R_objects/coexpr_erminer_pr_mm.RDS"
go_coexpr_ortho_path <- "/space/scratch/amorin/R_objects/coexpr_erminer_pr_ortho.RDS"

score_col <- "CorrectedMFPvalue"
term_col <- "Name"

pc_ortho <- read.delim(pc_ortho_path)
go_coexpr_hg <- readRDS(go_coexpr_hg_path)
go_coexpr_mm <- readRDS(go_coexpr_mm_path)
go_coexpr_ortho <- readRDS(go_coexpr_ortho_path)


# Example of one table - ASCL1 used in paper
# example_hg <- "ASCL1"
# example_mm <- "Ascl1"
# go_coexpr_hg[[example_hg]]
# go_coexpr_mm[[example_mm]]
# go_coexpr_ortho[[example_hg]]


# ErmineR sporadic bug that shifts columns over by one for a subset of terms 
# such that Name/term column is all NAs. This bug results in loss of GeneMembers
# data for affected elements. Check and shift all columns over by one.

check_and_fix_shifted_cols <- function(go_l, term_col) {
  
  tfs <- names(go_l)
  
  go_l <- lapply(tfs, function(x) {
    
    df <- go_l[[x]]
    na_terms <- which(is.na(df[[term_col]]))
    
    if (length(na_terms) > 0) {
      fix_df <- as.data.frame(df)  # tibble throws error if mismatch types
      fix_df[na_terms, 1:(ncol(fix_df) - 1)] <- fix_df[na_terms, 2:ncol(fix_df)]
      df <- as_tibble(fix_df)
    }
    
   return(df)
  })
  
  names(go_l) <- tfs
  return(go_l)
}



go_coexpr_hg <- check_and_fix_shifted_cols(go_coexpr_hg, term_col)
go_coexpr_mm <- check_and_fix_shifted_cols(go_coexpr_mm, term_col)
go_coexpr_ortho <- check_and_fix_shifted_cols(go_coexpr_ortho, term_col)


# Only consider results reaching sig. cut-off

filter_results <- function(go_l, cutoff_col, cutoff = 0.05) {
  lapply(go_l, function(x) filter(x, !!sym(cutoff_col) < cutoff))
}


go_coexpr_filt_hg <- filter_results(go_coexpr_hg, score_col)
go_coexpr_filt_mm <- filter_results(go_coexpr_mm, score_col)
go_coexpr_filt_ortho <- filter_results(go_coexpr_ortho, score_col)


# Count of significant terms per TF

n_terms_per_tf <- function(go_l) {
  n_l <- lapply(go_l, nrow)
  n_df <- data.frame(Symbol = names(go_l), N = unlist(n_l))
  n_df <- arrange(n_df, desc(N))
  return(n_df)
}

tf_n_hg <- n_terms_per_tf(go_coexpr_filt_hg)
tf_n_mm <- n_terms_per_tf(go_coexpr_filt_mm)
tf_n_ortho <- n_terms_per_tf(go_coexpr_filt_ortho)


# TFs with no associated terms

zero_hg <- tf_n_hg[tf_n_hg$N == 0, "Symbol"]
zero_mm <- tf_n_hg[tf_n_mm$N == 0, "Symbol"]
zero_ortho <- tf_n_hg[tf_n_ortho$N == 0, "Symbol"]


# Are ortho TFs correlated in their count of terms?

tf_n_all <- pc_ortho %>% 
  left_join(tf_n_hg, by = c("Symbol_hg" = "Symbol")) %>% 
  dplyr::rename(N_hg = N) %>% 
  left_join(tf_n_mm, by = c("Symbol_mm" = "Symbol")) %>%
  dplyr::rename(N_mm = N) %>% 
  filter(!is.na(N_hg), !is.na(N_mm)) %>% 
  left_join(tf_n_ortho, by = c("Symbol_hg" = "Symbol")) %>% 
  dplyr::rename(N_ortho = N) 


cor.test(tf_n_all$N_hg, tf_n_all$N_mm)
human0 <- tf_n_all$N_hg == 0
mouse0 <- tf_n_all$N_mm == 0
fisher.test(table(human0, mouse0))


# How well overlapping are the terms for ortho TFs?

# jacc_terms <- function(all_df, term_col, hg_l, mm_l) {
#   
#   jacc <- NA
#   tf_hg <- all_df$Symbol_hg[x]
#   tf_mm <- all_$Symbol_mm[x]
#   terms_hg <- hg_l[[tf_hg]][[term_col]]
#   terms_mm <- mm_l[[tf_mm]][[term_col]]
#   
#   
#   length(intersect(terms_hg, terms_mm)) / length(union(terms_hg, terms_mm))
#   
# }


tf_n_all$Jaccard <- unlist(lapply(1:nrow(tf_n_all), function(x) {
  tf_hg <- tf_n_all$Symbol_hg[x]
  tf_mm <- tf_n_all$Symbol_mm[x]
  terms_hg <- go_coexpr_filt_hg[[tf_hg]]$Name
  terms_mm <- go_coexpr_filt_mm[[tf_mm]]$Name
  length(intersect(terms_hg, terms_mm)) / length(union(terms_hg, terms_mm))
}))


summary(filter(tf_n_all, N_hg > 1 | N_mm > 1)$Jaccard)
summary(filter(tf_n_all, N_hg > 1 & N_mm > 1)$Jaccard)



# What are the most frequent/rare terms? Include the associated TFs

n_tfs_per_term <- function(go_l, term_col) {
  
  tf_terms <- lapply(names(go_l), function(tf) {
    term <- go_l[[tf]][[term_col]]
    if (length(term) == 0) return(NA)
    data.frame(Term = term, TF = tf, stringsAsFactors = FALSE)
  })
  
  tf_terms <- tf_terms[!is.na(tf_terms)]
  
  all_terms_df <- do.call(rbind, tf_terms)
  
  # Group by Term and tally occurrences
  result <- all_terms_df %>%
    group_by(Term) %>%
    summarise(
      N = n(),  # Count occurrences of each term
      TF = paste(unique(TF), collapse = ", ")  # Concat unique TFs
    ) %>% 
    arrange(desc(N))
  
  return(result)
}




# TODO Note: running on subset of list with at least 1 sig. term (error if include 0)
term_n_hg <- n_tfs_per_term(go_coexpr_filt_hg[tf_n_hg[tf_n_hg$N > 0, "Symbol"]], term_col)
term_n_mm <- n_tfs_per_term(go_coexpr_filt_mm[tf_n_mm[tf_n_mm$N > 0, "Symbol"]], term_col)
term_n_ortho <- n_tfs_per_term(go_coexpr_filt_ortho[tf_n_ortho[tf_n_ortho$N > 0, "Symbol"]], term_col)



# Are there terms that show species enrichment?

term_n_all <- data.frame(Term = union(term_n_hg$Term, term_n_mm$Term)) %>% 
  left_join(term_n_hg, by = "Term") %>% 
  dplyr::rename(N_hg = N, TF_hg = TF) %>% 
  left_join(term_n_mm, by = "Term") %>% 
  dplyr::rename(N_mm = N, TF_mm = TF) %>% 
  left_join(term_n_ortho, by = "Term") %>% 
  dplyr::rename(N_ortho = N, TF_ortho = TF) %>% 
  replace(is.na(.), 0)


cor.test(term_n_all$N_hg, term_n_all$N_mm)



# Plots
# ------------------------------------------------------------------------------


plot_bar <- function(df, count_col, group_col, xlab, ylab, title) {
  
  ggplot(df, aes(x = !!sym(count_col), 
                 y = reorder(!!sym(group_col), !!sym(count_col)))) +
    geom_bar(stat = "identity", colour = "slategrey", fill = "slategrey") +
    ylab(ylab) +
    xlab(xlab) +
    ggtitle(title) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text.y = element_blank(),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
  
}



plot_bar(tf_n_hg, 
         count_col = "N", 
         group_col = "Symbol", 
         xlab = "Count of GO terms (biological process)", 
         ylab = "TR", 
         title = "Human")


plot_bar(tf_n_mm, 
         count_col = "N", 
         group_col = "Symbol", 
         xlab = "Count of GO terms (biological process)", 
         ylab = "TR", 
         title = "Mouse")



plot_bar(term_n_hg, 
         count_col = "N", 
         group_col = "Term", 
         xlab = "Count of TRs belonging to term", 
         ylab = "GO terms (biological process)", 
         title = "Human")


plot_bar(term_n_mm, 
         count_col = "N", 
         group_col = "Term", 
         xlab = "Count of TRs belonging to term", 
         ylab = "GO terms (biological process)", 
         title = "Mouse")




ggplot(tf_n_all, aes(x = N_hg, y = N_mm)) +
  geom_point(shape = 21, size = 2.4) +
  geom_smooth(method = "lm") +
  xlab("Count of significant terms (human)") +
  ylab("Count of significant terms (mouse)") +
  ggtitle("N=1,241 orthologous TRs") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))



ggplot(term_n_all, aes(x = N_hg, y = N_mm)) +
  geom_point(shape = 21, size = 2.4) +
  geom_smooth(method = "lm") +
  xlab("Count of associated TRs (human)") +
  ylab("Count of associated TRs (mouse)") +
  ggtitle("N=5,692 unique terms between species") +
  


plot_hist(tf_n_hg,
          stat_col = "N",
          xlab = "Count of GO terms (biological process)",
          title = "Human") +
  ylab("Count of TRs")




plot_hist(term_n_hg,
          stat_col = "N",
          xlab = "Count of TRs",
          title = "Human") +
  ylab("Count of GO terms (biological process)")




# Demo simple barchart of GO terms and sig

go_coexpr_filt_ortho$ASCL1 %>% 
  mutate(Score = -log10(CorrectedMFPvalue)) %>% 
  arrange(Score) %>% 
  mutate(Name = factor(Name, levels = unique(Name))) %>% 
  ggplot(., aes(x = Score, y = Name)) +
  geom_bar(stat = "identity") +
  xlab("-log10 adjusted p-value") +
  ggtitle("ASCL1 orthologous rankings") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))
  



## Indidvidual


# dat_l <- readRDS(rank_tf_hg_path)
# dat_l <- readRDS(rank_tf_mm_path)
dat_l <- readRDS(rank_tf_ortho_path)

species <- "org.Hs.eg.db"
# species <- "org.Mm.eg.db"

tf <- "ASCL1"
# tf <- "Ascl1"

# universe <- dat_l[[1]]$Symbol
universe <- dat_l[[1]]$Symbol_hg

gene_entrez <- bitr(universe, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = species)

# top_genes <- slice_min(dat_l[[tf]], Rank_aggr_coexpr, n = topn)[["Symbol"]]
top_genes <- slice_min(dat_l[[tf]], Rank_aggr_coexpr_ortho, n = topn)[["Symbol_hg"]]
top_genes <- filter(gene_entrez, SYMBOL %in% top_genes)
  
go <- enrichGO(gene = top_genes$ENTREZID,
               universe = gene_entrez$ENTREZID,
               OrgDb = species,
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE)


dotplot(go, font.size = 17, showCategory = 20)
