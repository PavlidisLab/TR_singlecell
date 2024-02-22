## Combine the coexpression and binding rankings to nominate TF-gene interactions
## with reproducible evidence across species and methods
# TODO: name anonymous functions; consider splitting joining/re-ranking
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(cowplot)
library(ggrepel)
library(pheatmap)
library(grid)
library(gridExtra)
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


tier_names <- names(tiered_l[[1]])


# Tally the number of interactions at each tier for each TF

n_interactions <- lapply(tiered_l, function(tf) {
  unlist(lapply(tf, function(x) ifelse(is.null(x), NA, nrow(x))))
})

n_interactions <- do.call(rbind, n_interactions) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Symbol")


# Summary of total interactions at each tier

summ_tier <- lapply(tier_names, function(x) {
  list(
    N_pairs = sum(n_interactions[, x], na.rm = TRUE),
    N_TFs = sum(!is.na(n_interactions[, x]) & n_interactions[, x] > 0),
    Summ = summary(n_interactions[n_interactions[, x] > 0, x]))
})
names(summ_tier) <- tier_names



# Genes can reoccur across different TFs reproducible interactions. Tally these
# reoccurences, and group which TFs are affiliated with each unique gene
# NOTE: Pretty clunky, as generating different versions of list with 3 levels
# (tier, TF, and gene symbols in that tier) to make summaries.
# ------------------------------------------------------------------------------


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


# Reporting
# missing_counts["MAFB", ]  # MAFB largest elevated collection with 0 stringent
# summary(missing_counts)
# summary(missing_prop)
# summary(filter(missing_prop, N >= 5))


# How many species-specific interactions are not in the ortho set?
# ------------------------------------------------------------------------------


not_ortho <- lapply(names(tiered_l), function(tf) {
  
  tier_hg <- tiered_l[[tf]]$Human_specific
  tier_mm <- tiered_l[[tf]]$Mouse_specific
  
  if (is.null(tier_hg) || is.null(tier_mm)) {
    return(NA)
  }
  
  n_specific_hg <- nrow(tier_hg)
  n_not_ortho_hg <- length(setdiff(tier_hg$Symbol, pc_ortho$Symbol_hg))
  n_specific_mm <- nrow(tier_mm)
  n_not_ortho_mm <- length(setdiff(tier_mm$Symbol, pc_ortho$Symbol_mm))
  
  data.frame(
    Symbol = tf,
    N_specific_hg = n_specific_hg,
    N_ortho_hg = (n_specific_hg - n_not_ortho_hg),
    N_not_ortho_hg = n_not_ortho_hg,
    # Prop_not_ortho_hg = (n_not_ortho_hg / n_specific_hg),
    N_specific_mm = n_specific_mm,
    N_ortho_mm = (n_specific_mm - n_not_ortho_mm),
    N_not_ortho_mm = n_not_ortho_mm
    # Prop_not_ortho_mm = (n_not_ortho_mm / n_specific_mm)
  )
})


not_ortho <- do.call(rbind, not_ortho[!is.na(not_ortho)])



# Reporting
# plot_hist(n_not_ortho, stat_col = "N_not_ortho_hg")
# plot_hist(n_not_ortho, stat_col = "N_not_ortho_mm")
# summary(n_not_ortho$N_not_ortho_hg)
# summary(n_not_ortho$N_not_ortho_mm)
# setdiff(tiered_l$STAT1$Human_specific$Symbol, pc_ortho$Symbol_hg)
# setdiff(tiered_l$STAT1$Mouse_specific$Symbol, pc_ortho$Symbol_mm)



# Candidate evolutionary divergent interactions
# ------------------------------------------------------------------------------


# Select for genes that are top k in one species and bottom half in other species
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


# AP1 members typically have greatest similarity and size of reproducible sets.
# Look for genes ranked favourably across all
# ------------------------------------------------------------------------------


ap1_genes <- c("FOS", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND")

ap1_hg <- lapply(ap1_genes, function(x) {
  rank_hg[[x]] %>% 
    arrange(match(Symbol, pc_hg$Symbol)) %>%
    select(Symbol, Rank_bind, Rank_RSR, RP) %>% 
    rename_with(~paste0(., "_", x), -c("Symbol")) 
})


ap1_hg <- plyr::join_all(ap1_hg, by = "Symbol")
ap1_hg <- ap1_hg[order(rowSums(select_if(ap1_hg, is.integer))), ]



# Saving out objects
# ------------------------------------------------------------------------------


tiered_path <- "/space/scratch/amorin/R_objects/tiered_evidence_list.RDS"


# TODO: yet another reshape of the tiered list. must reconsider
# Get a list grouped by tier, with each element containing a single table of
# the TR-gene pairs at that tier


tiered_output <- lapply(tier_names, function(tier) {
  
  l <- lapply(names(tiered_l), function(x) {
    tf_by_tier <- tiered_l[[x]][[tier]]
    if (length(tf_by_tier) == 0 || nrow(tf_by_tier) == 0) return(NA)
    data.frame(TR = x, tf_by_tier)
  })
  
  l <- l[!is.na(l)]
  do.call(rbind, l)
})
names(tiered_output) <- tier_names


saveRDS(tiered_output, tiered_path)
  








# Plots
# ------------------------------------------------------------------------------


# Stacked barchart of tiered evidence counts


# Using only ortho

plot_df1a <- n_interactions %>% 
  pivot_longer(cols = c("Stringent", "Elevated", "Mixed"),
               names_to = "Scheme",
               values_to = "Count") %>% 
  filter(!is.na(Count)) %>% 
  mutate(Scheme = factor(Scheme, levels = c("Mixed", "Elevated", "Stringent")))
  

plot_df1a$Symbol <- factor(plot_df1a$Symbol, levels = arrange(n_interactions, Stringent)$Symbol)


p1a <- ggplot(
  plot_df1a, 
  aes(x = reorder(Symbol, Count, FUN = sum), y = Count, fill = Scheme, colour = Scheme)) +
  # aes(x = Symbol, y = Count, fill = Scheme, colour = Scheme)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Count of reproducible interactions") +
  xlab("Transcription regulator") +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))



# Including species specific

plot_df1b <- n_interactions %>% 
  pivot_longer(cols = c("Stringent", "Elevated", "Human_specific", "Mouse_specific", "Mixed"),
               names_to = "Scheme",
               values_to = "Count") %>% 
  filter(!is.na(Count)) %>%
  mutate(Scheme = factor(Scheme, levels = c("Mixed", "Mouse_specific", "Human_specific","Elevated", "Stringent")))
  

plot_df1b$Symbol <- factor(plot_df1b$Symbol, levels = arrange(n_interactions, Stringent)$Symbol)


p1b <- ggplot(
  plot_df1b, 
  aes(x = reorder(Symbol, Count, FUN = sum), y = Count, fill = Scheme, colour = Scheme)) +
  # aes(x = Symbol, y = Count, fill = Scheme, colour = Scheme)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Count of reproducible interactions") +
  xlab("Transcription regulator") +
  scale_fill_manual(values = c("#1b9e77", "goldenrod", "royalblue", "#d95f02", "#7570b3")) +
  scale_colour_manual(values = c("#1b9e77", "goldenrod", "royalblue", "#d95f02", "#7570b3")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))



# Plotting interactions for each species grouped by ortho status


plot_df2a <- n_interactions %>% 
  mutate(
    Ortho = !is.na(Mouse),
    # Orthologous = rowSums(n_interactions[, c("Stringent", "Elevated", "Mixed")], na.rm = TRUE),
    Orthologous = rowSums(n_interactions[, c("Stringent", "Elevated")], na.rm = TRUE),
    Orthologous = ifelse(Orthologous == 0, NA, Orthologous),
    Human_specific = ifelse(Ortho, Human_specific, Human)
    ) %>% 
  pivot_longer(cols = c("Human_specific", "Orthologous"),
               names_to = "Scheme",
               values_to = "Count") %>% 
  filter(!is.na(Count))


plot_df2b <- n_interactions %>% 
  mutate(
    Ortho = !is.na(Human),
    # Orthologous = rowSums(n_interactions[, c("Stringent", "Elevated", "Mixed")], na.rm = TRUE),
    Orthologous = rowSums(n_interactions[, c("Stringent", "Elevated")], na.rm = TRUE),
    Orthologous = ifelse(Orthologous == 0, NA, Orthologous),
    Mouse_specific = ifelse(Ortho, Mouse_specific, Mouse)
    ) %>% 
  pivot_longer(cols = c("Mouse_specific", "Orthologous"),
               names_to = "Scheme",
               values_to = "Count") %>% 
  filter(!is.na(Count))


p2a <- ggplot(
  plot_df2a, 
  aes(x = reorder(Symbol, Count, FUN = median), y = Count, fill = Scheme, colour = Scheme)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~Ortho, scales = "free") +
  ylab("Count of reproducible interactions") +
  xlab("Transcription regulator") +
  scale_fill_manual(values = c("royalblue", "darkgrey"), labels = c("Human specific", "Orthologous")) +
  scale_colour_manual(values = c("royalblue", "darkgrey")) +
  guides(colour = "none") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        # axis.title = element_text(size = 20),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        legend.position = c(0.75, 0.85),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        plot.margin = margin(c(10, 20, 10, 10)))


p2b <- ggplot(
  plot_df2b, 
  aes(x = reorder(Symbol, Count, FUN = median), y = Count, fill = Scheme, colour = Scheme)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~Ortho, scales = "free") +
  ylab("Count of reproducible interactions") +
  xlab("Transcription regulator") +
  scale_fill_manual(values = c("goldenrod", "darkgrey"), labels = c("Mouse specific", "Orthologous")) +
  scale_colour_manual(values = c("goldenrod", "darkgrey")) +
  guides(colour = "none") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        # axis.title = element_text(size = 20),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        legend.position = c(0.75, 0.85),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        plot.margin = margin(c(10, 20, 10, 10)))


# Combining species plots with common label
# https://stackoverflow.com/questions/33114380/centered-x-axis-label-for-muliplot-using-cowplot-package
p2 <- plot_grid(p2a, p2b, ncol = 1) 


y_grob <- textGrob("Count of reproducible interactions", 
                   gp = gpar(fontsize = 20), rot = 90)

x_grob <- textGrob("Transcription regulator", 
                   gp = gpar(fontsize = 20))

p2 <- grid.arrange(arrangeGrob(p2, left = y_grob, bottom = x_grob))


# Combine ortho with species specific
p3a <- plot_grid(p1a, p2, rel_widths = c(1, 0.5))
p3b <- plot_grid(p1b, p2, rel_widths = c(1, 0.5))


ggsave(p3a, height = 8, width = 20, device = "png", dpi = 300, bg = "white",
       filename = file.path(plot_dir, "ortho_reproducible_interaction_counts.png"))





# Heatmap of tiered evidence for a given TF. For species-specific interactions,
# join the equivalent ortho ranking

plot_tf <- "ASCL1"

rank_l <- tiered_l[[plot_tf]]

rank_l$Human_specific <- rank_l$Human_specific %>%
  dplyr::rename(Rank_RSR_hg = Rank_RSR, Rank_bind_hg = Rank_bind) %>%
  mutate(Symbol_hg = Symbol) %>% 
  left_join(pc_ortho, by = "Symbol_hg") %>% 
  left_join(rank_mm[[str_to_title(plot_tf)]][, c("Symbol", "Rank_bind", "Rank_RSR")],
  by = c("Symbol_mm" = "Symbol")) %>%
  dplyr::rename(Rank_RSR_mm = Rank_RSR, Rank_bind_mm = Rank_bind)
  

rank_l$Mouse_specific <- rank_l$Mouse_specific %>%
  dplyr::rename(Rank_RSR_mm = Rank_RSR, Rank_bind_mm = Rank_bind) %>%
  mutate(Symbol_mm = Symbol) %>% 
  left_join(pc_ortho, by = "Symbol_mm") %>% 
  left_join(rank_hg[[plot_tf]][, c("Symbol", "Rank_bind", "Rank_RSR")],
  by = c("Symbol_hg" = "Symbol")) %>%
  dplyr::rename(Rank_RSR_hg = Rank_RSR, Rank_bind_hg = Rank_bind) %>% 
  mutate(Symbol_hg = ifelse(is.na(Symbol_hg), Symbol_mm, Symbol_hg))


rank_l <- rank_l[c("Stringent", "Elevated", "Human_specific", "Mouse_specific", "Mixed")]


plot_df3 <- lapply(rank_l, `[`, c("Rank_RSR_hg", "Rank_RSR_mm", "Rank_bind_hg", "Rank_bind_mm")) %>%
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


rownames(plot_df3) <- unlist(lapply(rank_l, pluck, "Symbol_hg"))


# Padding between species and genes in evidence tiers 

gap_genes <- rep(
  head(cumsum(unlist(lapply(rank_l, nrow))), -1),
  each = 4)

# gap_evidence <- rep(1:4, each = 4)
gap_evidence <- c(rep(1, 4), rep(2, 8), rep(3, 4), rep(4, 2))


# With data types/species labels
pheatmap(t(plot_df3),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rev(c('white', "#a6bddb",'#045a8d')),
         border_color = "black",
         na_col = "grey",
         gaps_col = gap_genes,
         gaps_row = gap_evidence,
         cellheight = 10,
         legend = FALSE,
         angle_col = 90,
         width = 14,
         filename = file.path(plot_dir, paste0(plot_tf, "_tiered_evidence_targets_heatmap_label.png"))
)


# Without data types/species labels (for matching curated heatmap in illustrator)
pheatmap(t(plot_df3),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rev(c('white', "#a6bddb",'#045a8d')),
         border_color = "black",
         na_col = "grey",
         gaps_col = gap_genes,
         gaps_row = gap_evidence,
         show_rownames = FALSE,
         cellheight = 10,
         legend = FALSE,
         angle_col = 90,
         width = 14,
         filename = file.path(plot_dir, paste0(plot_tf, "_tiered_evidence_targets_heatmap_nolabel.png"))
)


# Binary heatmap of curation status

labels_curated <- get_curated_labels(tf = plot_tf, 
                                     curated_df = curated, 
                                     ortho_df = pc_ortho,
                                     pc_df = pc_hg, 
                                     species = "Human", 
                                     remove_self = TRUE)


curated_vec <- setNames(as.integer(rownames(plot_df3) %in% labels_curated), rownames(plot_df3))


# With gene symbols
pheatmap(t(curated_vec),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "black"),
         border_color = "black",
         gaps_col = gap_genes,
         cellheight = 10,
         legend = FALSE,
         angle_col = 90,
         width = 14,
         filename = file.path(plot_dir, paste0(plot_tf, "_curation_heatmap_label.png"))
)


# Without gene symbols (for matching tiered evidence heatmap in illustrator)
pheatmap(t(curated_vec),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "black"),
         border_color = "black",
         gaps_col = gap_genes,
         show_colnames = FALSE,
         cellheight = 10,
         legend = FALSE,
         angle_col = 90,
         width = 14,
         filename = file.path(plot_dir, paste0(plot_tf, "_curation_heatmap_nolabel.png"))
)



# Barchart of count of candidate evo divergent interactions

p4 <- 
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


# Non-orthologous genes gained in species-specific comparison


p_df5a <- not_ortho %>% 
  dplyr::select(Symbol, N_ortho_hg, N_not_ortho_hg) %>% 
  pivot_longer(cols = c("N_ortho_hg", "N_not_ortho_hg"),
               names_to = "Has_ortholog",
               values_to = "Count") %>% 
  mutate(Has_ortholog = ifelse(Has_ortholog == "N_ortho_hg", TRUE, FALSE))


p_df5b <- not_ortho %>% 
  dplyr::select(Symbol, N_ortho_mm, N_not_ortho_mm) %>% 
  pivot_longer(cols = c("N_ortho_mm", "N_not_ortho_mm"),
               names_to = "Has_ortholog",
               values_to = "Count") %>% 
  mutate(Has_ortholog = ifelse(Has_ortholog == "N_ortho_mm", TRUE, FALSE))


p5a <- 
  ggplot(p_df5a, aes(x = reorder(Symbol, Count, fun = sum), y = Count, fill = Has_ortholog, colour = Has_ortholog)) +
  geom_bar(position = "stack", stat = "identity") +
  ggtitle("Human") +
  scale_fill_manual(values = c("royalblue", "darkgrey")) +
  scale_colour_manual(values = c("royalblue", "darkgrey")) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = c(0.7, 0.8),
        plot.margin = margin(c(10, 20, 10, 10)))
  

p5b <- 
  ggplot(p_df5b, aes(x = reorder(Symbol, Count, fun = sum), y = Count, fill = Has_ortholog, colour = Has_ortholog)) +
  geom_bar(position = "stack", stat = "identity") +
  ggtitle("Mouse") +
  scale_fill_manual(values = c("goldenrod", "darkgrey")) +
  scale_colour_manual(values = c("goldenrod", "darkgrey")) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = c(0.7, 0.8),
        plot.margin = margin(c(10, 20, 10, 10)))


p5 <- plot_grid(p5a, p5b, ncol = 1)


ggsave(p5, height = 15, width = 18, device = "png", dpi = 300, bg = "white",
       filename = file.path(plot_dir, "species_specific_genes_gained.png"))



# Missing data type in elevated set


head(missing_prop)

p6a <- missing_prop %>% 
  filter(N >= 5) %>% 
  pivot_longer(cols = c("Prop_human", "Prop_mouse"),
               names_to = "Species",
               values_to = "Proportion") %>% 
  mutate(Species = ifelse(Species == "Prop_human", "Human", "Mouse")) %>% 
  ggplot(aes(., x = Species, y = Proportion)) +
  geom_boxplot(width = 0.3, fill = "slategrey") +
  ylab("Proportion Missing") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_blank())


p6b <- missing_prop %>% 
  filter(N >= 5) %>% 
  pivot_longer(cols = c("Prop_coexpr_hg", "Prop_bind_hg", "Prop_coexpr_mm", "Prop_bind_mm"),
               names_to = "Rank",
               values_to = "Proportion") %>% 
  mutate(Rank = case_when(
    Rank == "Prop_bind_hg" ~ "Binding human",
    Rank == "Prop_coexpr_hg" ~ "Coexpr. human",
    Rank == "Prop_bind_mm" ~ "Binding mouse",
    Rank == "Prop_coexpr_mm" ~ "Coexpr. mouse")
  ) %>% 
  ggplot(aes(., x = Rank, y = Proportion)) +
  geom_boxplot(width = 0.3, fill = "slategrey") +
  ylab("Proportion Missing") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_blank())


ggsave(p6a, height = 6, width = 6, device = "png", dpi = 300, bg = "white",
       filename = file.path(plot_dir, "species_elevated_proportion_missing.png"))


ggsave(p6b, height = 6, width = 12, device = "png", dpi = 300, bg = "white",
       filename = file.path(plot_dir, "rank_elevated_proportion_missing.png"))
