## Exploring the tiered integrated rankings, and focusing on ASCL1's targets
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(cowplot)
library(ggrepel)
library(grid)
library(gridExtra)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 500

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Integrated rankings
rank_hg <- readRDS(rank_int_hg_path)
rank_mm <- readRDS(rank_int_mm_path)
rank_ortho <- readRDS(rank_int_ortho_path)

# Tiered list (grouped by TFs)
tier_l <- readRDS(tiered_evidence_path)

# Tiered list (grouped by tiers)
tier_flat <- readRDS(tiered_evidence_flat_path)


tier_names <- names(tier_l[[1]])


# Inspecting the collections
# ------------------------------------------------------------------------------


# Tally each TF's interactions at each tier, setting missing to NA

tally_interactions <- function(tier_l) {
  
  n_intr <- lapply(tier_l, function(tf) {
    unlist(lapply(tf, function(x) ifelse(is.null(x), NA, nrow(x))))
  }) 
  
  n_intr <- do.call(rbind, n_intr) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Symbol")
  
  return(n_intr)
}


n_intr <- tally_interactions(tier_l)



# Summary of total interactions at each tier
# Paper focus on elevated and stringent: summ_tier$Elevated

summ_tier <- lapply(tier_names, function(x) {
  list(
    N_pairs = sum(n_intr[, x], na.rm = TRUE),
    N_TFs = sum(!is.na(n_intr[, x]) & n_intr[, x] > 0),
    Summ = summary(n_intr[n_intr[, x] > 0, x]))
})
names(summ_tier) <- tier_names



# The Species-specific tier can include interactions for TRs that had ChIP-seq
# in both species. Isolate counts for TRs with only ChIP-seq in one species

n_intr_bind_in_one_species <- n_intr %>% 
  filter(is.na(Stringent)) %>%   # removes TRs eligible for ortho (ChIP-seq in both species)
  summarise(
    N_human_TRs = sum(!is.na(Human)),
    Max_human = slice_max(., Human, n = 1)$Symbol,
    N_mouse_TRs = sum(!is.na(Mouse)),
    Max_mouse = slice_max(., Mouse, n = 1)$Symbol
  )


# Genes can reoccur across different TFs reproducible interactions. Tally these
# reoccurences, and group which TFs are affiliated with each unique gene
# ------------------------------------------------------------------------------


tally_symbols <- lapply(tier_flat, function(x) {
  sort(table(x[, intersect(colnames(x), c("Symbol_hg", "Symbol"))]))
})


# Get the count of genes affiliated with more than 1 TF at each tier
n_multi <- lapply(tally_symbols, function(x) sum(x > 1))


# Get the corresponding TFs for the genes belonging in the most TF sets
max_genes <- lapply(tally_symbols, function(x) names(x[x == max(x)]))

max_tfs <- lapply(tier_names, function(x) {
  symbol <- intersect(colnames(tier_flat[[x]]), c("Symbol", "Symbol_hg"))
  unique(filter(tier_flat[[x]], !!sym(symbol) %in% max_genes[[x]])$TR)
})
names(max_tfs) <- tier_names



# Which data type is typically absent from the elevated collection?
# ------------------------------------------------------------------------------



elevated_missing <- function(tier_l, k) {
  
  missing_counts <- lapply(tier_l, function(x) {
    
    if (is.null(x$Elevated)) {
      return(NA)
    }
    
    df <- data.frame(
      N = nrow(x$Elevated),
      N_coexpr_hg = sum(x$Elevated$Rank_aggr_coexpr_hg > k),
      N_coexpr_mm = sum(x$Elevated$Rank_aggr_coexpr_mm > k),
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
  
  missing_prop[, 2:ncol(missing_prop)] <- 
    t(apply(missing_counts, 1, function(x) x[2:ncol(missing_counts)] / x["N"]))
  
  missing_prop <- round(missing_prop, 3)
  
  return(list(Count = missing_counts, Prop = missing_prop))
}



missing <- elevated_missing(tier_l, k)


# Paper focuses on MAFB, which had the largest elevated collection of TFs
# with no stringent interactions. Find human binding is usually absent

mafb <- missing$Prop["MAFB", ]

# The higher the number, the more time it was the one ranking under k=500
summ_counts <- summary(missing$Count)
summ_prop <- summary(missing$Prop)



# How many interactions are in the species-specific sets because they lack an
# ortholog?
# ------------------------------------------------------------------------------


tally_not_ortholog <- function(tier_l, pc_ortho) {
  
  not_ortho <- lapply(names(tier_l), function(tf) {
    
    tier_hg <- tier_l[[tf]]$Human_specific
    tier_mm <- tier_l[[tf]]$Mouse_specific
    
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
      Prop_not_ortho_hg = (n_not_ortho_hg / n_specific_hg),
      N_specific_mm = n_specific_mm,
      N_ortho_mm = (n_specific_mm - n_not_ortho_mm),
      N_not_ortho_mm = n_not_ortho_mm,
      Prop_not_ortho_mm = (n_not_ortho_mm / n_specific_mm)
    )
  })
  
  
  not_ortho <- do.call(rbind, not_ortho[!is.na(not_ortho)])
  return(not_ortho)
}



not_ortho <- tally_not_ortholog(tier_l, pc_ortho)


# Note that immune TFs seem to have most interactions that are well supported
# but lack an ortholog

top_not_ortho_hg <- not_ortho %>% 
  filter(N_not_ortho_hg >= 3) %>%   # prevent small n
  slice_max(Prop_not_ortho_hg, n = 10)


top_not_ortho_mm <- not_ortho %>% 
  filter(N_not_ortho_mm >= 3) %>% 
  slice_max(Prop_not_ortho_mm, n = 10)



# Candidate evolutionary divergent interactions for each TF in ortho set
# ------------------------------------------------------------------------------


# Select genes that are top k in one species and ranked bottom half the in other

filter_evo_candidates <- function(rank_ortho, k) {
  
  mid_ortho <- nrow(rank_ortho[[1]]) / 2
  
  evo_div <- lapply(rank_ortho, function(rank_df) {
    
    # Decrease k if the current k includes tied values
    k_bind_hg <- check_k(sort(rank_df$Bind_score_hg, decreasing = TRUE), k = k)
    k_bind_mm <- check_k(sort(rank_df$Bind_score_mm, decreasing = TRUE), k = k)
    k_coexpr_hg <- check_k(sort(rank_df$Avg_aggr_coexpr_hg, decreasing = TRUE), k = k)
    k_coexpr_mm <- check_k(sort(rank_df$Avg_aggr_coexpr_mm, decreasing = TRUE), k = k)
    
    human <- 
      filter(rank_df,
            (Rank_aggr_coexpr_hg <= k_coexpr_hg & Rank_bind_hg <= k_bind_hg) &
            (Rank_aggr_coexpr_mm >= mid_ortho & Rank_bind_mm >= mid_ortho))
    
    mouse <- 
      filter(rank_df,
            (Rank_aggr_coexpr_mm <= k_coexpr_mm & Rank_bind_mm <= k_bind_mm) &
            (Rank_aggr_coexpr_hg >= mid_ortho & Rank_bind_hg >= mid_ortho))
    
    list(Human = human, Mouse = mouse)
    
  })
  
  return(evo_div)
}



evo_div <- filter_evo_candidates(rank_ortho, k)


# Count candidates
species_count <- function(species_l) {
  data.frame(Human = nrow(species_l$Human), Mouse = nrow(species_l$Mouse))
}

n_evo_div <- do.call(rbind, lapply(evo_div, species_count)) %>% 
  rownames_to_column(var = "Symbol") %>% 
  arrange(desc(Human))



# AP1 members typically have greatest similarity and size of reproducible sets.
# Look for genes ranked highly across all
# ------------------------------------------------------------------------------


ap1_genes <- c("FOS", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND")


collect_gene_ranks <- function(genes, rank_l, pc_df) {
  
  gene_ranks <- lapply(genes, function(x) {
    rank_l[[x]] %>% 
      arrange(match(Symbol, pc_df$Symbol)) %>%
      select(Symbol, Rank_bind, Rank_aggr_coexpr, Rank_integrated) %>% 
      rename_with(~paste0(., "_", x), -c("Symbol")) 
  })
  
  gene_ranks <- plyr::join_all(gene_ranks, by = "Symbol")
  gene_ranks <- gene_ranks[order(rowSums(select_if(gene_ranks, is.integer))), ]
  return(gene_ranks)
}


ap1_hg <- collect_gene_ranks(ap1_genes, rank_hg, pc_hg)



# Plots
# ------------------------------------------------------------------------------


# Stacked barchart of tiered evidence counts. Each bar is a count of unique
# interactions that were gained in the new tier (eg, gained in elevated that
# were not in stringent)


# Using only ortho (used in paper)

pdf1a <- n_intr %>% 
  pivot_longer(cols = c("Stringent", "Elevated", "Mixed"),
               names_to = "Scheme",
               values_to = "Count") %>% 
  filter(!is.na(Count)) %>% 
  mutate(Scheme = factor(Scheme, levels = c("Mixed", "Elevated", "Stringent")))
  

pdf1a$Symbol <- factor(pdf1a$Symbol, levels = arrange(n_intr, Stringent)$Symbol)


p1a <- ggplot(
  pdf1a, 
  aes(x = reorder(Symbol, Count, FUN = sum), 
      y = Count, fill = Scheme, 
      colour = Scheme)) +
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



# Including species specific (too cluttered)

plot_df1b <- n_intr %>% 
  pivot_longer(cols = c(
    "Stringent", "Elevated", "Human_specific", "Mouse_specific", "Mixed"),
    names_to = "Scheme",
    values_to = "Count") %>% 
  filter(!is.na(Count)) %>%
  mutate(Scheme = factor(
    Scheme,
    levels = c(
      "Mixed",
      "Mouse_specific",
      "Human_specific",
      "Elevated",
      "Stringent"
    )
  ))
  

plot_df1b$Symbol <- factor(plot_df1b$Symbol, 
                           levels = arrange(n_intr, Stringent)$Symbol)


p1b <- ggplot(
  plot_df1b, 
  aes(x = reorder(Symbol, Count, FUN = sum), y = Count, 
      fill = Scheme, colour = Scheme)) +
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



# Plotting stacked bar of interactions for each species grouped by ortho status:
# Split TRs by whether they have binding data in both species or not.
# If not, they were not eligible for consideration in the ortho tiers (ie,
# stringent/elevated). All resulting interactions are species specific. For 
# those with data in both, coloured bar indicates unique species-specific
# interactions gained relative to the stringent/elevated class

hist_species <- function(n_intr, species, colour) {
  
  ortho_counts <- rowSums(n_intr[, c("Stringent", "Elevated")], na.rm = TRUE)
  ortho_counts <- ifelse(ortho_counts == 0, NA, ortho_counts)
  
  # If Human check if Mouse counts exist (implies ChIP-seq); same with reverse
  in_ortho <- !is.na(n_intr[[setdiff(c("Human", "Mouse"), species)]])
  
  # uniquely gained relative to stringent and elevated
  species_specific <- n_intr[[paste0(species, "_specific")]]
  species_all <- n_intr[[species]]
  species_counts <- ifelse(in_ortho, species_specific, species_all)
  
  pdf <- n_intr %>% 
    mutate(In_ortho = in_ortho,
           Orthologous = ortho_counts,
           Species_specific = species_counts) %>% 
    pivot_longer(cols = c("Species_specific", "Orthologous"),
                 names_to = "Scheme",
                 values_to = "Count") %>% 
    filter(!is.na(Count)) %>% 
    mutate(Scheme = factor(Scheme, 
                           levels = c("Species_specific", "Orthologous")))

  
  ggplot(pdf,
         aes(x = reorder(Symbol, Count, FUN = median), 
             y = Count,
             fill = Scheme,
             colour = Scheme)) +
    geom_bar(position = "stack", stat = "identity") +
    facet_wrap(~In_ortho, scales = "free") +
    ylab("Count of reproducible interactions") +
    xlab("Transcription regulator") +
    scale_fill_manual(values = c(colour, "darkgrey"), 
                      labels = c(paste0(species, " specific"), "Orthologous")) +
    scale_colour_manual(values = c(colour, "darkgrey")) +
    guides(colour = "none") +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 20),
          legend.position = c(0.75, 0.85),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 15),
          plot.margin = margin(c(10, 20, 10, 10)))
}



p2a <- hist_species(n_intr, "Human", "royalblue")
p2b <- hist_species(n_intr, "Mouse", "goldenrod")


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



# Barchart of count of candidate evo divergent interactions

p4 <- 
  n_evo_div %>% 
  pivot_longer(cols = c("Human", "Mouse"),
               names_to = "Species",
               values_to = "Count") %>% 
  ggplot(
  aes(x = reorder(Symbol, Count, FUN = median), y = Count, fill = Species)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  ylab("Count of candidate interactions") +
  xlab("Transcription factor") +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))


# Genes with no orthologs that were gained in species-specific comparison


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

p6a <- missing$Prop %>% 
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


p6b <- missing$Prop %>% 
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
