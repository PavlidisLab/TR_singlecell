## Combine the coexpression and binding rankings to nominate TF-gene interactions
## with reproducible evidence across species and methods
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(cowplot)
library(ggrepel)
library(pheatmap)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 500

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Saved list of the aggregate coexpression ranks
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)

# Average bind scores
collection <- "Permissive"
stopifnot(collection %in% c("Robust", "Permissive"))
bind_summary_path <- paste0("/space/scratch/amorin/R_objects/unibind_", collection, "_bindscore_summary.RDS")
bind_summary <- readRDS(bind_summary_path)
bind_dat <- readRDS("/space/scratch/amorin/R_objects/processed_unibind_data.RDS")

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)



# Functions
# ------------------------------------------------------------------------------


# For TFs with available ChIP-seq data, add the aggregated binding profile to
# the aggregate coexpression ranking, and calc the rank product of the binding
# and coexpression profiles

join_coexpr_bind <- function(coexpr_l, bind_mat) {
  
  common_genes <- intersect(names(coexpr_l), colnames(bind_mat))
  
  join_l <- lapply(common_genes, function(x) {
    
    bind_df <- data.frame(Symbol = rownames(bind_mat), 
                          Bind_score = bind_mat[, x])
    
    join_df <- left_join(coexpr_l[[x]], bind_df, by = "Symbol") %>%
      filter(Symbol != x) %>% 
      mutate(
        Rank_bind = rank(-Bind_score, ties.method = "min"),
        RP = rank(Rank_RSR * Rank_bind)
      ) %>% 
      relocate(Rank_RSR, .before = RP) %>% 
      arrange(RP)
  })
  
  names(join_l) <- common_genes
  return(join_l)
}


# Take the list of ranked dfs for each species, and for the ortho TFs with 
# binding data, subset their ranked dfs to ortho genes and then join the two
# species. Re-rank the coexpression and binding scores with the reduced gene 
# list, and generate the rank product across all four ranks

join_and_rank_ortho <- function(rank_hg, rank_mm, pc_ortho, ncores = 1) {
  
  # Common ortho TFs
  tf_ortho <- filter(pc_ortho, 
                     Symbol_hg %in% names(rank_hg) & 
                     Symbol_mm %in% names(rank_mm))
  
  # Only keeping key columns for each ranking df
  sub_cols <- c("Symbol", "Avg_RSR", "Bind_score")
  
  # Iteratively subset to ortho genes, join human and mouse, and re-rank
  ortho_l <- mclapply(1:nrow(tf_ortho), function(x) {
    
    df_hg <- left_join(rank_hg[[tf_ortho$Symbol_hg[x]]][, sub_cols],
                       pc_ortho[, c("Symbol_hg", "ID")],
                       by = c("Symbol" = "Symbol_hg")) %>%
      filter(!is.na(ID))
    
    
    df_mm <- left_join(rank_mm[[tf_ortho$Symbol_mm[x]]][, sub_cols], 
                       pc_ortho[, c("Symbol_mm", "ID")], 
                       by = c("Symbol" = "Symbol_mm")) %>% 
      filter(!is.na(ID))
    
    
    df_ortho <- left_join(df_hg, df_mm,
                          by = "ID",
                          suffix = c("_hg", "_mm")) %>%
      filter(!is.na(Avg_RSR_hg) & !is.na(Avg_RSR_mm)) %>%
      mutate(
        Rank_RSR_hg = rank(-Avg_RSR_hg, ties.method = "min"),
        Rank_RSR_mm = rank(-Avg_RSR_mm, ties.method = "min"),
        Rank_bind_hg = rank(-Bind_score_hg, ties.method = "min"),
        Rank_bind_mm = rank(-Bind_score_mm, ties.method = "min"),
        RP = rank(log(Rank_RSR_hg + Rank_RSR_mm + Rank_bind_hg + Rank_bind_mm))
      )
  }, mc.cores = ncores)
  
  names(ortho_l) <- tf_ortho$Symbol_hg
  return(ortho_l)
}


# Taking the list of ortho rankings, subset each rank df into 3 tiers of
# evidence: stringent (top K in all 4 rankings), middle (top K in 3/4), and
# relaxed, which allows genes with top K in one data type in the opposite species
# (but not both data types in just one species - want to emphasize ortho aspect)


# TODO:
# TODO: DUX4
# TODO: k from species-specific versus ortho
# TODO: Note species rank versus ortho rank (minimal difference)

subset_tiered_evidence <- function(rank_ortho, 
                                   rank_hg,
                                   rank_mm,
                                   pc_ortho,
                                   k, 
                                   ncores = 1) {
  # Ortho symbol map
  tfs <- filter(pc_ortho, 
                Symbol_hg %in% names(rank_hg) |
                Symbol_mm %in% names(rank_mm))
  
  
  tiered_l <- mclapply(1:nrow(tfs), function(x) {
    
    # NULL if no data across species
    specific_hg <- specific_mm <- stringent <- elevated <- mixed <- NULL
    
    # Species and ortho rankings for the given TF
    df_ortho <- rank_ortho[[tfs$Symbol_hg[x]]]
    df_hg <- rank_hg[[tfs$Symbol_hg[x]]]
    df_mm <- rank_mm[[tfs$Symbol_mm[x]]]
    
    # Convert ortho ranks to logical status for cut-off. Filter species to genes
    # meeting cut-off in both data types. May need to decrease k if the current 
    # k includes tied values
    
    if (!is.null(df_hg)) {
      
      k_rsr_hg <- check_k(sort(df_hg$Avg_RSR, decreasing = TRUE), k = k)
      k_bind_hg <- check_k(sort(df_hg$Bind_score, decreasing = TRUE), k = k)
      df_hg <- filter(df_hg, Rank_RSR <= k_rsr_hg & Rank_bind <= k_bind_hg)
      
    }
    
    if (!is.null(df_mm)) {
      
      k_rsr_mm <- check_k(sort(df_mm$Avg_RSR, decreasing = TRUE), k = k)
      k_bind_mm <- check_k(sort(df_mm$Bind_score, decreasing = TRUE), k = k)
      df_mm <- filter(df_mm, Rank_RSR <= k_rsr_mm & Rank_bind <= k_bind_mm)
      
    }
    
    if (!is.null(df_ortho)) {
      
      k_rsr_hg <- check_k(sort(df_ortho$Avg_RSR_hg, decreasing = TRUE), k = k)
      k_bind_hg <- check_k(sort(df_ortho$Bind_score_hg, decreasing = TRUE), k = k)
      k_rsr_mm <- check_k(sort(df_ortho$Avg_RSR_mm, decreasing = TRUE), k = k)
      k_bind_mm <- check_k(sort(df_ortho$Bind_score_mm, decreasing = TRUE), k = k)
      
      df_ortho <- mutate(df_ortho,
                         Cutoff_RSR_hg = Rank_RSR_hg <= k_rsr_hg,
                         Cutoff_RSR_mm = Rank_RSR_mm <= k_rsr_mm,
                         Cutoff_bind_hg = Rank_bind_hg <= k_bind_hg,
                         Cutoff_bind_mm = Rank_bind_mm <= k_bind_mm)
      
      # Stringent requires top k in all four of human and mouse coexpr and binding 
      stringent <- df_ortho %>%  
        filter(Cutoff_RSR_hg & Cutoff_RSR_mm & Cutoff_bind_hg & Cutoff_bind_mm) %>%
        select(-contains("Cutoff")) %>%
        arrange(RP)
    
      # Elevated requires evidence in 3/4 rankings
      elevated <- df_ortho %>%
        mutate(n = rowSums(.[grep("Cutoff", names(.))])) %>%
        filter(n == 3) %>%
        select(-c(contains("Cutoff"), n)) %>%
        arrange(RP)
    
      # Species specific
      specific_hg <- filter(df_hg, Symbol %!in% c(stringent$Symbol_hg, elevated$Symbol_hg))
      specific_mm <- filter(df_mm, Symbol %!in% c(stringent$Symbol_mm, elevated$Symbol_mm))
    
      # Mixed species requires evidence in one data type for each species
      mixed <- df_ortho %>%
        filter(
        (Cutoff_RSR_hg & Cutoff_bind_mm) | (Cutoff_RSR_mm & Cutoff_bind_hg)) %>% 
        filter(
        Symbol_hg %!in% c(stringent$Symbol_hg, elevated$Symbol_hg, specific_hg$Symbol) &
        Symbol_mm %!in% specific_mm$Symbol) %>% 
        select(-contains("Cutoff")) %>% 
        arrange(RP)
    
    } 
    
    list(Stringent = stringent, 
         Elevated = elevated,
         Human = df_hg,
         Human_specific = specific_hg,
         Mouse = df_mm,
         Mouse_specific = specific_mm,
         Mixed = mixed)
    
  }, mc.cores = ncores)
  
  names(tiered_l) <- tfs$Symbol_hg
  return(tiered_l)
}


# Join aggregate coexpression and binding and generate tiered evidence
# ------------------------------------------------------------------------------


# Produces lists of TF gene rank dfs for coexpression and binding data

rank_hg <- join_coexpr_bind(rank_tf_hg, bind_summary$Human_TF)
rank_mm <- join_coexpr_bind(rank_tf_mm, bind_summary$Mouse_TF)
rank_ortho <- join_and_rank_ortho(rank_hg, rank_mm, pc_ortho, ncores = ncore)


# Produces a list for each TF within binding evidence in at least one species,
# grouping genes by their status in the top K across species/data type rankings

tiered_l <- subset_tiered_evidence(rank_ortho = rank_ortho,
                                   rank_hg = rank_hg,
                                   rank_mm = rank_mm,
                                   pc_ortho = pc_ortho,
                                   k = k, 
                                   ncores = ncore)


# Inspecting the collections
# ------------------------------------------------------------------------------


# Tally the number of interactions at each tier for each TF

n_interactions <- lapply(tiered_l, function(tf) {
  unlist(lapply(tf, function(x) ifelse(is.null(x), NA, nrow(x))))
})

n_interactions <- do.call(rbind, n_interactions)


# Summary of total interactions at each tier

summ_tier <- lapply(colnames(n_interactions), function(x) {
  list(
    N_pairs = sum(n_interactions[, x], na.rm = TRUE),
    N_TFs = sum(!is.na(n_interactions[, x]) & n_interactions[, x] > 0),
    Summ = summary(n_interactions[n_interactions[, x] > 0, x]))
})
names(summ_tier) <- colnames(n_interactions)



# Genes can reoccur across different TFs reproducible interactions. Tally these
# reoccurences, and group which TFs are affiliated with each unique gene
# NOTE: Pretty clunky, as generating different versions of list with 3 levels
# (tier, TF, and gene symbols in that tier) to make summaries.
# ------------------------------------------------------------------------------


tier_names <- colnames(n_interactions)


# Extract the gene symbols for each TF at each tier
# list: TF -> tier -> gene symbols
tf_symbols <- lapply(tiered_l, function(x) {
  lapply(x, function(y) {
    y[, intersect(colnames(y), c("Symbol", "Symbol_hg"))]
  })
})


# Invert this list and remove empty so list is now tier with TF elements
# list: tier -> TF -> gene symbols
tier_symbols <- lapply(tier_names, function(x) {
  l <- lapply(tf_symbols, pluck, x)
  l <- l[lapply(l, length) > 0]
})
names(tier_symbols) <- tier_names


# Tally genes at each tier
tally_symbols <- lapply(tier_symbols, function(x) sort(table(unlist(x))))


# Get the count of genes affiliated with more than 1 TF at each tier
n_multi <- lapply(tally_symbols, function(x) sum(x > 1))


# Inverted list: tier -> unique gene symbols -> TFs affiliated with each gene
gene_tf  <- lapply(tier_names, function(tier) {
  
  gene_names <- names(tally_symbols[[tier]])
  tf_l <- tier_symbols[[tier]]
  
  gene_l <- lapply(gene_names, function(gene) {
      names(unlist(lapply(tf_l, function(tf) which(gene %in% tf))))
  })
  names(gene_l) <- gene_names
  
  # order by size of 'TF sets' each gene belongs to
  gene_l <- gene_l[order(vapply(gene_l, length, integer(1)))]
  
  return(gene_l)
})
names(gene_tf) <- tier_names



# Get the corresponding TFs from the genes with the largest TF sets
max_genes <- lapply(tally_symbols, function(x) names(x[x == max(x)]))
max_tfs <- lapply(tier_names, function(x) gene_tf[[x]][max_genes[[x]]])
names(max_tfs) <- tier_names



# Focus on stringent interactions in paper
# tail(tally_symbols$Stringent, 10)
# n_multi$Stringent
# max_genes$Stringent
# max_tfs$Stringent
# unique(unlist(max_tfs$Stringent))


# Which data type is typically absent from the elevated collection?
# ------------------------------------------------------------------------------


missing_counts <- lapply(tiered_l, function(x) {
  
  if (is.null(x$Elevated)) {
    return(NA)
  }

  df <- data.frame(
    N = nrow(x$Elevated),
    N_coexpr_hg = sum(x$Elevated$Rank_RSR_hg > k),
    N_coexpr_mm = sum(x$Elevated$Rank_RSR_mm > k),
    N_bind_hg = sum(x$Elevated$Rank_bind_hg > k),
    N_bind_mm = sum(x$Elevated$Rank_bind_mm > k))
  
  df$N_human <- sum(df$N_bind_hg, df$N_coexpr_hg)
  df$N_mouse <- sum(df$N_bind_mm, df$N_coexpr_mm)
  df$N_coexpr <- sum(df$N_coexpr_hg, df$N_coexpr_mm)
  df$N_bind <- sum(df$N_bind_hg, df$N_bind_mm)
  
  return(df)
})


missing_counts <- do.call(rbind, missing_counts[!is.na(missing_counts)])


# Proportions
missing_prop <- missing_counts
colnames(missing_prop) <- str_replace(colnames(missing_prop), "N_", "Prop_")
missing_prop[, 2:ncol(missing_prop)] <- t(apply(missing_counts, 1, function(x) x[2:ncol(missing_counts)] / x["N"]))
missing_prop <- round(missing_prop, 3)

# missing_counts["MAFB", ]  # MAFB largest elevated collection with 0 stringent
# summary(missing_counts)
# summary(missing_prop)
# summary(filter(missing_prop, N >= 5))


# How many species-specific interactions are not in the ortho set?


common_hg <- intersect(names(rank_hg), names(rank_ortho))


specific_hg <- lapply(names(rank_hg), function(x) {

  topk <- topk_hg[[x]]$Symbol
  
  not_ortho <- length(setdiff(topk, pc_ortho$Symbol_hg))
    
  if (x %!in% names(rank_ortho)) {
    
    ortho <- FALSE
    specific = length(topk)
    
  } else {
    
    ortho <- TRUE
    specific <- length(setdiff(
    topk,
    # unlist(lapply(tiered_l[[x]], pluck, "Symbol_hg")))  # Including Mixed
    unlist(lapply(tiered_l[[x]][c("Stringent", "Elevated")], pluck, "Symbol_hg"))
    ))
  }
  
  prop_not_ortho <- not_ortho / specific
  
  data.frame(Symbol = x, 
             Ortho = ortho, 
             N_specific = specific, 
             N_not_ortho = not_ortho,
             Prop_not_ortho = prop_not_ortho)
  
})


specific_hg <- do.call(rbind, specific_hg)




specific_mm <- lapply(names(rank_mm), function(x) {

  topk <- topk_mm[[x]]$Symbol
  
  tf_ortho <- filter(pc_ortho, Symbol_mm == x)
  
  not_ortho <- length(setdiff(topk, pc_ortho$Symbol_mm))
    
  if (tf_ortho$Symbol_hg %!in% names(rank_ortho)) {
    
    ortho <- FALSE
    specific = length(topk)
    
  } else {
    
    ortho <- TRUE
    specific <- length(setdiff(
    topk,
    # unlist(lapply(tiered_l[[x]], pluck, "Symbol_mm")))  # Including Mixed
    unlist(lapply(tiered_l[[tf_ortho$Symbol_hg]][c("Stringent", "Elevated")], pluck, "Symbol_mm"))
    ))
  }
  
  prop_not_ortho <- not_ortho / specific
  
  data.frame(Symbol = x, 
             Ortho = ortho, 
             N_specific = specific, 
             N_not_ortho = not_ortho,
             Prop_not_ortho = prop_not_ortho)
  
})


specific_mm <- do.call(rbind, specific_mm)




n_topk_hg <- left_join(specific_hg, n_topk_ortho, by = "Symbol") %>% 
  dplyr::rename("Human_specific" = "N_specific")
n_topk_hg[is.na(n_topk_hg)] <- 0


n_topk_mm <- specific_mm %>%
  left_join(pc_ortho, by = c("Symbol" = "Symbol_mm")) %>% 
  left_join(., n_topk_ortho, by = c("Symbol_hg" = "Symbol")) %>% 
  dplyr::rename("Mouse_specific" = "N_specific")
n_topk_mm[is.na(n_topk_mm)] <- 0


# TODO: make histogram plots of species-specific gains

summary(filter(n_topk_hg, Ortho))
summary(filter(n_topk_mm, Ortho))

summary(filter(n_topk_hg, !Ortho))
summary(filter(n_topk_mm, !Ortho))

sum(filter(n_topk_hg, Ortho)$N_not_ortho == 0)
sum(filter(n_topk_mm, Ortho)$N_not_ortho == 0)


# Genes in species-specific that are not in ortho set

setdiff(topk_hg$STAT1$Symbol, pc_ortho$Symbol_hg)



# Candidate evolutionary divergent interactions

mid_ortho <- nrow(rank_ortho[[1]]) / 2

evo_div <- lapply(rank_ortho, function(rank_df) {

  # May need to decrease k if the current k includes tied values
  k_bind_hg <- check_k(sort(rank_df$Bind_score_hg, decreasing = TRUE), k = k)
  k_bind_mm <- check_k(sort(rank_df$Bind_score_mm, decreasing = TRUE), k = k)
  k_rsr_hg <- check_k(sort(rank_df$Avg_RSR_hg, decreasing = TRUE), k = k)
  k_rsr_mm <- check_k(sort(rank_df$Avg_RSR_mm, decreasing = TRUE), k = k)
  
  human = filter(rank_df,
                 Rank_RSR_hg <= k_rsr_hg & Rank_bind_hg <= k_bind_hg &
                 Rank_RSR_mm >= mid_ortho & Rank_bind_mm >= mid_ortho)
  
  mouse = filter(rank_df,
                 Rank_RSR_mm <= k_rsr_mm & Rank_bind_mm <= k_bind_mm &
                 Rank_RSR_hg >= mid_ortho & Rank_bind_hg >= mid_ortho)
  
  list(Human = human, Mouse = mouse)
  
})



n_evo_div <- data.frame(
  Symbol = names(rank_ortho),
  do.call(rbind, lapply(evo_div, function(x) data.frame(Human = nrow(x$Human), Mouse = nrow(x$Mouse))))
)



n_evo_div %>% 
  pivot_longer(cols = c("Human", "Mouse"),
               names_to = "Species",
               values_to = "Count") %>% 
  ggplot(
  # aes(x = reorder(Symbol, Count, FUN = median), y = Count, fill = Species, colour = Species)) +
  aes(x = reorder(Symbol, Count, FUN = median), y = Count, fill = Species)) +
  # aes(x = Symbol, y = Count, fill = Scheme, colour = Scheme)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  ylab("Count of candidate interactions") +
  xlab("Transcription factor") +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  # scale_colour_manual(values = c("royalblue", "goldenrod")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))



# Plotting interactions grouped by ortho collection


p_dfx <- pivot_longer(n_topk_hg,
                      cols = c("Human_specific", "Stringent", "Elevated", "Mixed"),
                      names_to = "Scheme",
                      values_to = "Count") %>% 
  mutate(Scheme = factor(Scheme, levels = c("Mixed", "Human_specific", "Elevated", "Stringent")))
  


p_dfx$Symbol <- factor(p_dfx$Symbol, levels = arrange(n_topk_hg, Stringent)$Symbol)
# p_dfx$Symbol <- factor(p_dfx$Symbol, levels = arrange(n_topk_hg, Human_specific)$Symbol)
  

# https://github.com/dgrtwo/drlib/blob/master/R/reorder_within.R
ggplot(
  p_dfx, 
  # aes(x = reorder(Symbol, Count, FUN = median), y = Count, fill = Scheme, colour = Scheme)) +
  aes(x = Symbol, y = Count, fill = Scheme, colour = Scheme)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~ Ortho, scales = "free") +
  ylab("Count of reproducible interactions") +
  xlab("Transcription factor") +
  ggtitle(paste0("N=", n_distinct(p_dfx$Symbol))) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))



p_dfx2a <- n_topk_hg %>% 
  # mutate(N_ortho = rowSums(n_topk_hg[, c("Stringent", "Elevated", "Mixed")])) %>% 
  mutate(N_ortho = rowSums(n_topk_hg[, c("Stringent", "Elevated")])) %>% 
  select(-c(Stringent, Elevated, Mixed)) %>% 
  pivot_longer(cols = c("Human_specific", "N_ortho"),
               names_to = "Scheme",
               values_to = "Count")



p_dfx2b <- n_topk_mm %>% 
  # mutate(N_ortho = rowSums(n_topk_mm[, c("Stringent", "Elevated", "Mixed")])) %>% 
  mutate(N_ortho = rowSums(n_topk_mm[, c("Stringent", "Elevated")])) %>% 
  select(-c(Stringent, Elevated, Mixed)) %>% 
  pivot_longer(cols = c("Mouse_specific", "N_ortho"),
               names_to = "Scheme",
               values_to = "Count")


p2a <- ggplot(
  p_dfx2a, 
  aes(x = reorder(Symbol, Count, FUN = median), y = Count, fill = Scheme, colour = Scheme)) +
  # aes(x = Symbol, y = Count, fill = Scheme, colour = Scheme)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~Ortho, scales = "free") +
  ylab("Count of reproducible interactions") +
  xlab("Transcription factor") +
  # ggtitle(paste0("N=", nrow(p_dfx))) +
  scale_fill_manual(values = c("royalblue", "darkgrey")) +
  scale_colour_manual(values = c("royalblue", "darkgrey")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        # axis.title = element_text(size = 20),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        legend.position = c(0.75, 0.75),
        plot.margin = margin(c(10, 20, 10, 10)))


p2b <- ggplot(
  p_dfx2b, 
  aes(x = reorder(Symbol, Count, FUN = median), y = Count, fill = Scheme, colour = Scheme)) +
  # aes(x = Symbol, y = Count, fill = Scheme, colour = Scheme)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~Ortho, scales = "free") +
  ylab("Count of reproducible interactions") +
  xlab("Transcription factor") +
  # ggtitle(paste0("N=", nrow(p_dfx))) +
  scale_fill_manual(values = c("goldenrod", "darkgrey")) +
  scale_colour_manual(values = c("goldenrod", "darkgrey")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        # axis.title = element_text(size = 20),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        legend.position = c(0.75, 0.75),
        plot.margin = margin(c(10, 20, 10, 10)))


# https://stackoverflow.com/questions/33114380/centered-x-axis-label-for-muliplot-using-cowplot-package
p2 <- plot_grid(p2a, p2b, ncol = 1) 

library(grid)
library(gridExtra)

y.grob <- textGrob("Count of reproducible interactions", 
                   gp = gpar(fontsize = 20), rot = 90)

x.grob <- textGrob("Transcription factor", 
                   gp = gpar(fontsize = 20))

#add to plot

p2 <- grid.arrange(arrangeGrob(p2, left = y.grob, bottom = x.grob))
  





# Demo overlap of AP1 members


ap1_genes <- c("FOS", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND")

ap1_hg <- lapply(ap1_genes, function(x) {
  rank_hg[[x]] %>% 
    arrange(match(Symbol, pc_hg$Symbol)) %>%
    select(Symbol, Rank_bind, Rank_RSR, RP) %>% 
    rename_with(~paste0(., "_", x), -c("Symbol")) 
})


ap1_hg <- plyr::join_all(ap1_hg, by = "Symbol")



# Plots
# ------------------------------------------------------------------------------


# Stacked barchart of tiered evidence counts

plot_df1 <- pivot_longer(n_topk_ortho,
                         cols = c("Stringent", "Elevated", "Mixed"),
                         names_to = "Scheme",
                         values_to = "Count") %>% 
  mutate(Scheme = factor(Scheme, levels = c("Mixed", "Elevated", "Stringent")))
  

plot_df1$Symbol <- factor(plot_df1$Symbol, levels = arrange(n_topk_ortho, Stringent)$Symbol)


p1 <- ggplot(
  plot_df1, 
  aes(x = reorder(Symbol, Count, FUN = sum), y = Count, fill = Scheme, colour = Scheme)) +
  # aes(x = Symbol, y = Count, fill = Scheme, colour = Scheme)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Count of reproducible interactions") +
  xlab("Transcription factor") +
  # ggtitle(paste0("N=", nrow(n_topk_ortho))) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        legend.position = c(0.75, 0.75),
        plot.margin = margin(c(10, 20, 10, 10)))


plot_grid(p1, p2, rel_widths = c(1, 0.5))



# Heatmap of tiered evidence for a given TF

plot_tf <- "ASCL1"

scale01 <- function(x) (x - min(x)) / (max(x) - min(x))


genes <- c(tiered_l[[plot_tf]]$Stringent$Symbol_hg,
           tiered_l[[plot_tf]]$Elevated$Symbol_hg,
           tiered_l[[plot_tf]]$Mixed$Symbol_hg)


plot_df2 <- rank_ortho[[plot_tf]] %>% 
  mutate(Avg_RSR_hg = scale01(Avg_RSR_hg),
         Avg_RSR_mm = scale01(Avg_RSR_mm),
         Bind_score_hg = scale01(Bind_score_hg),
         Bind_score_mm = scale01(Bind_score_mm),
         Rank_RSR_hg = as.integer(Rank_RSR_hg <= k), 
         Rank_RSR_mm = as.integer(Rank_RSR_mm <= k), 
         Rank_bind_hg = as.integer(Rank_bind_hg <= k), 
         Rank_bind_mm = as.integer(Rank_bind_mm <= k)) %>%  
  filter(Symbol_hg %in% genes) %>% 
  arrange(match(Symbol_hg, genes))

rownames(plot_df2) <- plot_df2$Symbol_hg


# Padding between species (columns) and genes in evidence tiers (rows)

ns <- nrow(tiered_l[[plot_tf]]$Stringent)
nm <- nrow(tiered_l[[plot_tf]]$Elevated)

gaps_row <- rep(c(ns, nm + ns), each = 4)
gaps_col <- rep(1, 4)


pheatmap(plot_df2[, c("Avg_RSR_hg", "Avg_RSR_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000'),
         breaks = seq(0, 1, length.out = 9),
         border_color = "black",
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_stringent_coexpr_heatmap.png")
)


pheatmap(plot_df2[, c("Bind_score_hg", "Bind_score_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858'),
         breaks = seq(0, 1, length.out = 9),
         border_color = "black",
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_stringent_binding_heatmap.png")
)


# Binarizing status


pheatmap(plot_df2[, c("Rank_RSR_hg", "Rank_RSR_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c('white', '#7f0000'),
         border_color = "black",
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_stringent_coexpr_binary_heatmap.png")
)


pheatmap(plot_df2[, c("Rank_bind_hg", "Rank_bind_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c('white', '#045a8d'),
         border_color = "black",
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_stringent_binding_binary_heatmap.png")
)


# Bin status


plot_df2 <- rank_ortho[[plot_tf]] %>% 
  
  mutate(
    
    Avg_RSR_hg = scale01(Avg_RSR_hg),
    Avg_RSR_mm = scale01(Avg_RSR_mm),
    Bind_score_hg = scale01(Bind_score_hg),
    Bind_score_mm = scale01(Bind_score_mm), 
    
    Rank_RSR_hg = case_when(Rank_RSR_hg <= k ~ 0,
                            Rank_RSR_hg > k & Rank_RSR_hg <= k + 500 ~ 1,
                            TRUE ~ 2), 
    
    Rank_RSR_mm = case_when(Rank_RSR_mm <= k ~ 0,
                            Rank_RSR_mm > k & Rank_RSR_mm <= k + 500 ~ 1,
                            TRUE ~ 2),
    
    Rank_bind_hg = case_when(Rank_bind_hg <= k ~ 0,
                             Rank_bind_hg > k & Rank_bind_hg <= k + 500 ~ 1,
                             TRUE ~ 2),
    
    Rank_bind_mm = case_when(Rank_bind_mm <= k ~ 0,
                             Rank_bind_mm > k & Rank_bind_mm <= k + 500 ~ 1,
                             TRUE ~ 2)
    ) %>% 
 
  filter(Symbol_hg %in% genes) %>% 
  arrange(match(Symbol_hg, genes))

rownames(plot_df2) <- plot_df2$Symbol_hg



pheatmap(plot_df2[, c("Rank_RSR_hg", "Rank_RSR_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rev(c('white', "#fc8d59", '#7f0000')),
         border_color = "black",
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_stringent_coexpr_step_heatmap.png")
)


pheatmap(plot_df2[, c("Rank_bind_hg", "Rank_bind_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rev(c('white', "#a6bddb",'#045a8d')),
         border_color = "black",
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_stringent_binding_step_heatmap.png")
)


# Binary curated status


labels_curated <- get_curated_labels(tf = plot_tf, 
                                     curated_df = curated, 
                                     ortho_df = pc_ortho,
                                     pc_df = pc_hg, 
                                     species = "Human", 
                                     remove_self = TRUE)

curated_vec <- setNames(as.integer(genes %in% labels_curated), genes)


pheatmap(curated_vec,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "black"),
         border_color = "black",
         gaps_row = gaps_row,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_curated_heatmap.png")
)


# Using RP as order, take the top whatever, and show the evidence of the 
# other groups as blocked colours



plot_df3 <- rank_ortho[[plot_tf]] %>% 
  
  slice_min(RP, n = 40) %>% 
  
  mutate(
    
    Rank_RSR_hg = case_when(Rank_RSR_hg <= k ~ 0,
                            Rank_RSR_hg > k & Rank_RSR_hg <= k + 500 ~ 1,
                            TRUE ~ 2), 
    
    Rank_RSR_mm = case_when(Rank_RSR_mm <= k ~ 0,
                            Rank_RSR_mm > k & Rank_RSR_mm <= k + 500 ~ 1,
                            TRUE ~ 2),
    
    Rank_bind_hg = case_when(Rank_bind_hg <= k ~ 0,
                             Rank_bind_hg > k & Rank_bind_hg <= k + 500 ~ 1,
                             TRUE ~ 2),
    
    Rank_bind_mm = case_when(Rank_bind_mm <= k ~ 0,
                             Rank_bind_mm > k & Rank_bind_mm <= k + 500 ~ 1,
                             TRUE ~ 2)
    )

rownames(plot_df3) <- plot_df3$Symbol_hg


pheatmap(plot_df3[, c("Rank_RSR_hg", "Rank_RSR_mm", "Rank_bind_hg", "Rank_bind_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rev(c('white', "#a6bddb",'#045a8d')),
         border_color = "black",
         gaps_col = rep(1:4, each = 4),
         cellwidth = 20,
         cellheight = 20,
         legend = FALSE,
         filename = file.path(plot_dir, "demo_rp_order_step_heatmap.png")
)



curated_vec <- setNames(as.integer(plot_df3$Symbol_hg %in% labels_curated), plot_df3$Symbol_hg)


pheatmap(curated_vec,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "black"),
         border_color = "black",
         cellwidth = 20,
         cellheight = 20,
         legend = FALSE,
         filename = file.path(plot_dir, "demo_curated_heatmap.png")
)



# Join species specific and tiered


plot_tf <- "PAX6"

rank_l <- tiered_l[[plot_tf]]


rm_genes <- filter(pc_ortho, Symbol_hg %in% c(rank_l$Stringent$Symbol_hg, rank_l$Elevated$Symbol_hg))


gain_hg <- topk_hg[[plot_tf]] %>%
  filter(Symbol %!in% rm_genes$Symbol_hg) %>%
  dplyr::rename(Rank_RSR_hg = Rank_RSR, Rank_bind_hg = Rank_bind) %>%
  mutate(Symbol_hg = Symbol) %>% 
  left_join(pc_ortho, by = "Symbol_hg")


gain_mm <- topk_mm[[str_to_title(plot_tf)]] %>%
  filter(Symbol %!in% rm_genes$Symbol_mm) %>%
  dplyr::rename(Rank_RSR_mm = Rank_RSR, Rank_bind_mm = Rank_bind) %>%
  mutate(Symbol_mm = Symbol) %>% 
  left_join(pc_ortho, by = "Symbol_mm")



gain_hg <- left_join(
  gain_hg,
  rank_mm[[str_to_title(plot_tf)]][, c("Symbol", "Rank_bind", "Rank_RSR", "RP")],
  by = c("Symbol_mm" = "Symbol")) %>%
  dplyr::rename(Rank_RSR_mm = Rank_RSR, Rank_bind_mm = Rank_bind)


gain_mm <- left_join(
  gain_mm,
  rank_hg[[plot_tf]][, c("Symbol", "Rank_bind", "Rank_RSR", "RP")],
  by = c("Symbol_hg" = "Symbol")) %>%
  dplyr::rename(Rank_RSR_hg = Rank_RSR, Rank_bind_hg = Rank_bind) %>% 
  mutate(Symbol_hg = ifelse(is.na(Symbol_hg), Symbol_mm, Symbol_hg))



demo_tiered <- list(
  Stringent = tiered_l[[plot_tf]]$Stringent,
  Elevated = tiered_l[[plot_tf]]$Elevated,
  Human = gain_hg,
  Mouse = gain_mm,
  Mixed = filter(rank_l$Mixed, Symbol_hg %!in% c(gain_hg$Symbol_hg, gain_mm$Symbol_hg))
)



plot_df4 <- lapply(demo_tiered, `[`, c("Rank_RSR_hg", "Rank_RSR_mm", "Rank_bind_hg", "Rank_bind_mm")) %>%
  do.call(rbind, .) %>% 

  mutate(

    Rank_RSR_hg = case_when(Rank_RSR_hg <= k ~ 0,
                            Rank_RSR_hg > k & Rank_RSR_hg <= k + 500 ~ 1,
                            is.na(Rank_RSR_hg) ~ NA_real_,
                            TRUE ~ 2),

    Rank_RSR_mm = case_when(Rank_RSR_mm <= k ~ 0,
                            Rank_RSR_mm > k & Rank_RSR_mm <= k + 500 ~ 1,
                            is.na(Rank_RSR_mm) ~ NA_real_,
                            TRUE ~ 2),

    Rank_bind_hg = case_when(Rank_bind_hg <= k ~ 0,
                             Rank_bind_hg > k & Rank_bind_hg <= k + 500 ~ 1,
                             is.na(Rank_bind_hg) ~ NA_real_,
                             TRUE ~ 2),

    Rank_bind_mm = case_when(Rank_bind_mm <= k ~ 0,
                             Rank_bind_mm > k & Rank_bind_mm <= k + 500 ~ 1,
                             is.na(Rank_bind_mm) ~ NA_real_,
                             TRUE ~ 2)
  )




rownames(plot_df4) <- unlist(lapply(demo_tiered, pluck, "Symbol_hg"))


# Padding between species and genes in evidence tiers 


gap_genes <- rep(
  head(cumsum(unlist(lapply(demo_tiered, nrow))), -1),
  each = 4)


gap_evidence <- rep(1:4, each = 4)


pheatmap(t(plot_df4),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rev(c('white', "#a6bddb",'#045a8d')),
         border_color = "black",
         na_col = "black",
         gaps_col = gap_genes,
         gaps_row = gap_evidence,
         # cellwidth = 10,
         cellheight = 10,
         legend = FALSE,
         angle_col = 90
         # filename = file.path(plot_dir, "demo_rp_order_step_heatmap.png")
)


# Versus just using the ortho list

plot_tf <- "PAX6"

rank_l <- tiered_l[[plot_tf]]

rm_genes <- filter(pc_ortho, Symbol_hg %in% c(rank_l$Stringent$Symbol_hg, rank_l$Elevated$Symbol_hg))


gain_hg2 <- rank_ortho[[plot_tf]] %>%
  filter(Rank_RSR_hg <= k & Rank_bind_hg <= k) %>% 
  filter(Symbol_hg %!in% rm_genes$Symbol_hg)


gain_mm2 <- rank_ortho[[plot_tf]] %>%
  filter(Rank_RSR_mm <= k & Rank_bind_mm <= k) %>% 
  filter(Symbol_hg %!in% rm_genes$Symbol_hg)


demo_tiered2 <- list(
  Stringent = tiered_l[[plot_tf]]$Stringent,
  Elevated = tiered_l[[plot_tf]]$Elevated,
  Human = gain_hg2,
  Mouse = gain_mm2,
  Mixed = filter(rank_l$Mixed, Symbol_hg %!in% c(gain_hg2$Symbol_hg, gain_mm2$Symbol_hg))
)


plot_df5 <- lapply(demo_tiered2, `[`, c("Rank_RSR_hg", "Rank_RSR_mm", "Rank_bind_hg", "Rank_bind_mm")) %>%
  do.call(rbind, .) %>% 

  mutate(

    Rank_RSR_hg = case_when(Rank_RSR_hg <= k ~ 0,
                            Rank_RSR_hg > k & Rank_RSR_hg <= k + 500 ~ 1,
                            is.na(Rank_RSR_hg) ~ NA_real_,
                            TRUE ~ 2),

    Rank_RSR_mm = case_when(Rank_RSR_mm <= k ~ 0,
                            Rank_RSR_mm > k & Rank_RSR_mm <= k + 500 ~ 1,
                            is.na(Rank_RSR_mm) ~ NA_real_,
                            TRUE ~ 2),

    Rank_bind_hg = case_when(Rank_bind_hg <= k ~ 0,
                             Rank_bind_hg > k & Rank_bind_hg <= k + 500 ~ 1,
                             is.na(Rank_bind_hg) ~ NA_real_,
                             TRUE ~ 2),

    Rank_bind_mm = case_when(Rank_bind_mm <= k ~ 0,
                             Rank_bind_mm > k & Rank_bind_mm <= k + 500 ~ 1,
                             is.na(Rank_bind_mm) ~ NA_real_,
                             TRUE ~ 2)
  )


rownames(plot_df5) <- unlist(lapply(demo_tiered2, pluck, "Symbol_hg"))


# Padding between species and genes in evidence tiers 


gap_genes <- rep(
  head(cumsum(unlist(lapply(demo_tiered2, nrow))), -1),
  each = 4)


gap_evidence <- rep(1:4, each = 4)


pheatmap(t(plot_df5),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rev(c('white', "#a6bddb",'#045a8d')),
         border_color = "black",
         na_col = "black",
         gaps_col = gap_genes,
         gaps_row = gap_evidence,
         # cellwidth = 10,
         cellheight = 10,
         legend = FALSE,
         angle_col = 90
)



labels_curated <- get_curated_labels(tf = plot_tf, 
                                     curated_df = curated, 
                                     ortho_df = pc_ortho,
                                     pc_df = pc_hg, 
                                     species = "Human", 
                                     remove_self = TRUE)


curated_vec <- setNames(as.integer(rownames(plot_df5) %in% labels_curated), rownames(plot_df5))


pheatmap(t(curated_vec),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "black"),
         border_color = "black",
         gaps_col = gap_genes,
         # cellwidth = 20,
         cellheight = 10,
         legend = FALSE,
         angle_col = 90
)

