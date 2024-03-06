## Combine the coexpression and binding rankings to nominate TF-gene interactions
## with reproducible evidence across species and methods
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 500
force_resave <- TRUE

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Saved list of the aggregate coexpression ranks
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)
rank_tf_ortho <- readRDS(rank_tf_mm_path)

# Average bind scores
bind_summary <- readRDS(bind_summary_path)


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
      mutate(
        Rank_bind = rank(-Bind_score, ties.method = "min"),
        Rank_integrated = rank(Rank_aggr_coexpr * Rank_bind,
                               ties.method = "min",
                               na.last = "keep")
      ) %>% 
      relocate(c(Rank_integrated, Rank_aggr_coexpr, Rank_bind),
               .after = Symbol) %>% 
      arrange(Rank_integrated)
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
  
  # Iteratively subset to ortho genes, join human and mouse, and re-rank
  ortho_l <- mclapply(1:nrow(tf_ortho), function(x) {
    
    
    df_hg <- left_join(rank_hg[[tf_ortho$Symbol_hg[x]]],
                       pc_ortho[, c("Symbol_hg", "ID")],
                       by = c("Symbol" = "Symbol_hg")) %>%
      dplyr::select(-Rank_integrated) %>%
      filter(!is.na(ID))
    
    
    df_mm <- left_join(rank_mm[[tf_ortho$Symbol_mm[x]]],
                       pc_ortho[, c("Symbol_mm", "ID")],
                       by = c("Symbol" = "Symbol_mm")) %>%
      dplyr::select(-Rank_integrated) %>%
      filter(!is.na(ID))
    
    
    df_ortho <- left_join(df_hg, df_mm,
                          by = "ID",
                          suffix = c("_hg", "_mm")) %>%
      # filter(!is.na(Avg_aggr_coexpr_hg) & !is.na(Avg_aggr_coexpr_mm)) %>%
      dplyr::select(-ID) %>% 
      mutate(
        Rank_aggr_coexpr_hg = rank(-Avg_aggr_coexpr_hg, 
                                   ties.method = "min",
                                   na.last = "keep"),
        Rank_aggr_coexpr_mm = rank(-Avg_aggr_coexpr_mm, 
                                   ties.method = "min",
                                   na.last = "keep"),
        Rank_bind_hg = rank(-Bind_score_hg, ties.method = "min"),
        Rank_bind_mm = rank(-Bind_score_mm, ties.method = "min"),
        Rank_integrated = rank(
          log(Rank_aggr_coexpr_hg + Rank_aggr_coexpr_mm + Rank_bind_hg + Rank_bind_mm),
          ties.method = "min",
          na.last = "keep")
      ) %>% 
        relocate(c(Symbol_mm, 
                   Rank_integrated, 
                   Rank_aggr_coexpr_hg, 
                   Rank_aggr_coexpr_mm, 
                   Rank_bind_hg,
                   Rank_bind_mm),
                 .after = Symbol_hg) %>% 
      arrange(Rank_integrated)
  
  }, mc.cores = ncores)
  
  names(ortho_l) <- tf_ortho$Symbol_hg
  return(ortho_l)
}



# Taking the list of ortho rankings, subset each rank df into tiers of
# evidence: stringent (top K in all 4 rankings), elevated (top K in 3/4), 
# species-specific (topk K in coexpr and binding within species only), and
# mixed-species, which has top K in both data types but one in each species

subset_tiered_evidence <- function(rank_ortho, 
                                   rank_hg,
                                   rank_mm,
                                   pc_ortho,
                                   k, 
                                   ncores = 1) {
  # Ortho symbol map
  tf_ortho <- filter(pc_ortho, 
                Symbol_hg %in% names(rank_hg) |
                Symbol_mm %in% names(rank_mm))
  
  
  tiered_l <- mclapply(1:nrow(tf_ortho), function(x) {
    
    # NULL if no data across species
    specific_hg <- specific_mm <- stringent <- elevated <- mixed <- NULL
    
    # Species and ortho rankings for the given TF
    df_ortho <- rank_ortho[[tf_ortho$Symbol_hg[x]]]
    df_hg <- rank_hg[[tf_ortho$Symbol_hg[x]]]
    df_mm <- rank_mm[[tf_ortho$Symbol_mm[x]]]
    
    # Convert ortho ranks to logical status for cut-off. Filter species to genes
    # meeting cut-off in both data types. May need to decrease k if the current 
    # k includes tied values
    
    if (!is.null(df_hg)) {
      
      k_coexpr_hg <- check_k(sort(df_hg$Avg_aggr_coexpr, decreasing = TRUE), k = k)
      k_bind_hg <- check_k(sort(df_hg$Bind_score, decreasing = TRUE), k = k)
      df_hg <- filter(df_hg, Rank_aggr_coexpr <= k_coexpr_hg & Rank_bind <= k_bind_hg)
      
    }
    
    if (!is.null(df_mm)) {
      
      k_coexpr_mm <- check_k(sort(df_mm$Avg_aggr_coexpr, decreasing = TRUE), k = k)
      k_bind_mm <- check_k(sort(df_mm$Bind_score, decreasing = TRUE), k = k)
      df_mm <- filter(df_mm, Rank_aggr_coexpr <= k_coexpr_mm & Rank_bind <= k_bind_mm)
      
    }
    
    if (!is.null(df_ortho)) {
      
      k_coexpr_hg <- check_k(sort(df_ortho$Avg_aggr_coexpr_hg, decreasing = TRUE), k = k)
      k_bind_hg <- check_k(sort(df_ortho$Bind_score_hg, decreasing = TRUE), k = k)
      k_coexpr_mm <- check_k(sort(df_ortho$Avg_aggr_coexpr_mm, decreasing = TRUE), k = k)
      k_bind_mm <- check_k(sort(df_ortho$Bind_score_mm, decreasing = TRUE), k = k)
      
      df_ortho <- mutate(df_ortho,
                         Cutoff_aggr_coexpr_hg = Rank_aggr_coexpr_hg <= k_coexpr_hg,
                         Cutoff_aggr_coexpr_mm = Rank_aggr_coexpr_mm <= k_coexpr_mm,
                         Cutoff_bind_hg = Rank_bind_hg <= k_bind_hg,
                         Cutoff_bind_mm = Rank_bind_mm <= k_bind_mm)
      
      # Stringent requires top k in all four of human and mouse coexpr and binding 
      stringent <- df_ortho %>%  
        filter(
          Cutoff_aggr_coexpr_hg & 
          Cutoff_aggr_coexpr_mm & 
          Cutoff_bind_hg & 
          Cutoff_bind_mm
        ) %>%
        dplyr::select(-contains("Cutoff")) %>%
        arrange(Rank_integrated)
      
      # Elevated requires evidence in 3/4 rankings
      elevated <- df_ortho %>%
        mutate(n = rowSums(.[grep("Cutoff", names(.))])) %>%
        filter(n == 3) %>%
        dplyr::select(-c(contains("Cutoff"), n)) %>%
        arrange(Rank_integrated)
      
      # Species specific
      specific_hg <- filter(df_hg, Symbol %!in% c(stringent$Symbol_hg, elevated$Symbol_hg))
      specific_mm <- filter(df_mm, Symbol %!in% c(stringent$Symbol_mm, elevated$Symbol_mm))
      
      # Mixed species requires evidence in one data type for each species
      mixed <- df_ortho %>%
        filter(
          (Cutoff_aggr_coexpr_hg & Cutoff_bind_mm) | 
          (Cutoff_aggr_coexpr_mm & Cutoff_bind_hg)
        ) %>% 
        filter(
          Symbol_hg %!in% c(stringent$Symbol_hg, elevated$Symbol_hg, specific_hg$Symbol) &
          Symbol_mm %!in% specific_mm$Symbol
        ) %>% 
        select(-contains("Cutoff")) %>% 
        arrange(Rank_integrated)
      
    } 
    
    list(Stringent = stringent, 
         Elevated = elevated,
         Human = df_hg,
         Human_specific = specific_hg,
         Mouse = df_mm,
         Mouse_specific = specific_mm,
         Mixed = mixed)
    
  }, mc.cores = ncores)
  
  names(tiered_l) <- tf_ortho$Symbol_hg
  return(tiered_l)
}



# Collapse each tier of evidence into a single dataframe of TR-gene pairs

flatten_tiered_evidence <- function(tiered_evidence) {
  
  tier_names <- names(tiered_evidence[[1]]) 
  
  tiered_output <- lapply(tier_names, function(tier) {
    
    l <- lapply(names(tiered_evidence), function(x) {
      
      tf_by_tier <- tiered_evidence[[x]][[tier]]
      if (length(tf_by_tier) == 0 || nrow(tf_by_tier) == 0) return(NA)
      
      data.frame(TR = x, tf_by_tier)
      
    })
    
    l <- l[!is.na(l)]
    do.call(rbind, l)
  })
  
  names(tiered_output) <- tier_names
  return(tiered_output)
}



# Generate and save objects
# ------------------------------------------------------------------------------


# Human TFs
save_function_results(
  path = rank_int_hg_path,
  fun = join_coexpr_bind,
  args = list(
    coexpr_l = rank_tf_hg,
    bind_mat = bind_summary$Human_TF
  ),
  force_resave = force_resave
)



# Mouse TFs
save_function_results(
  path = rank_int_mm_path,
  fun = join_coexpr_bind,
  args = list(
    coexpr_l = rank_tf_mm,
    bind_mat = bind_summary$Mouse_TF
  ),
  force_resave = force_resave
)



# Ortho
rank_hg <- readRDS(rank_int_hg_path)
rank_mm <- readRDS(rank_int_mm_path)

save_function_results(
  path = rank_int_ortho_path,
  fun = join_and_rank_ortho,
  args = list(
    rank_hg = rank_hg,
    rank_mm = rank_mm,
    pc_ortho = pc_ortho,
    ncores = ncore
  ),
  force_resave = force_resave
)



# Tiered evidence
rank_ortho <- readRDS(rank_int_ortho_path)

save_function_results(
  path = tiered_evidence_path,
  fun = subset_tiered_evidence,
  args = list(
    rank_ortho = rank_ortho,
    rank_hg = rank_hg,
    rank_mm = rank_mm,
    pc_ortho = pc_ortho,
    k = k,
    ncores = ncore
  ),
  force_resave = force_resave
)



# TR-gene pairs for each tier
tiered_evidence <- readRDS(tiered_evidence_path)

save_function_results(
  path = tiered_evidence_flat_path,
  fun = flatten_tiered_evidence,
  args = list(tiered_evidence = tiered_evidence),
  force_resave = force_resave
)
