## Perform GO enrichment of biological processes on all aggregate coexpression
## profiles using ermineR's precision recall implementation
## https://github.com/PavlidisLab/ermineR
## -----------------------------------------------------------------------------

library(ermineR)
library(rJava)
library(tidyverse)
library(parallel)
source("R/utils/functions.R")
source("R/00_config.R")

set.seed(5)


# Download GO terms and human/mouse gene to GO map

if (!file.exists(go_path)) {
  goAtDate(go_path, go_date, overwrite = FALSE)
}



if (!file.exists(anno_hg_path)) {
  invisible(
    gemma.R::get_platform_annotations(platform = "Generic_human_ncbiIds",
                                      file = anno_hg_path,
                                      unzip = TRUE)
  )
}


if (!file.exists(anno_mm_path)) {
  invisible(
    gemma.R::get_platform_annotations(platform = "Generic_mouse_ncbiIds",
                                      file = anno_mm_path,
                                      unzip = TRUE)
  )
}



# Assumes dat_l as named list of TF coexpression rankings
# Assumes summary df has Rank_aggr_coexpr and Symbol in column names
# Iterates through each aggregate ranking in dat_l, performing GO enrichment
# of biological process terms using ermineR::precRecall
# ------------------------------------------------------------------------------


erminer_enrich_list <- function(dat_l, 
                                score_col,
                                go_path,
                                anno_path,
                                ncores = 1) {
  
  anno <- read.delim(anno_path)
  anno$ElementName <- anno$GeneSymbols  # ncbi IDs to symbols
  
  
  tfs <- names(dat_l)
  
  erminer_l <- mclapply(tfs, function(tf) {
    
    df <- dat_l[[tf]]
    rownames(df) <- df$Symbol
    
    tryCatch(
      {
        res <- precRecall(annotation = anno,
                          scores = df,
                          scoreColumn = score_col,
                          geneSetDescription = go_path,
                          aspects = "B",  # biological process
                          iterations = 10000,
                          bigIsBetter = FALSE)  # integer rank smaller better
        
        res$results
        
      },
      error = function(e) {
        NA
      }
    )
    
  }, mc.cores = ncores)
  names(erminer_l) <- tfs
  
  return(erminer_l)
}




# Human

if (!file.exists(erminer_coexpr_hg_path)) {
  
  message(paste("Human", Sys.time()))
  
  coexpr_hg <- readRDS(rank_tf_hg_path)
  
  save_function_results(
    path = erminer_coexpr_hg_path,
    fun = erminer_enrich_list,
    args = list(
      dat_l = coexpr_hg,
      score_col = "Rank_aggr_coexpr",
      go_path = go_path,
      anno_path = anno_hg_path,
      ncores = ncore
    )
  )
}



# Mouse

if (!file.exists(erminer_coexpr_mm_path)) {
  
  message(paste("Mouse", Sys.time()))
  
  coexpr_mm <- readRDS(rank_tf_mm_path)
  
  save_function_results(
    path = erminer_coexpr_mm_path,
    fun = erminer_enrich_list,
    args = list(
      dat_l = coexpr_mm,
      score_col = "Rank_aggr_coexpr",
      go_path = go_path,
      anno_path = anno_mm_path,
      ncores = ncore
    )
  )
}


# Note that for ortho (ranking placing equal weight on mouse and human aggregate)
# uses human symbols/annotations. Inspection using mouse anno gave similar results

if (!file.exists(erminer_coexpr_ortho_path)) {
  
  message(paste("Ortho", Sys.time()))
  
  coexpr_ortho <- readRDS(rank_tf_ortho_path)
  
  # Coerce to expected column names
  coexpr_ortho <- lapply(coexpr_ortho, function(x) {
    dplyr::rename(x, Symbol = Symbol_hg, Rank_aggr_coexpr = Rank_aggr_coexpr_ortho)
  })
  
  save_function_results(
    path = erminer_coexpr_ortho_path,
    fun = erminer_enrich_list,
    args = list(
      dat_l = coexpr_ortho,
      score_col = "Rank_aggr_coexpr",
      go_path = go_path,
      anno_path = anno_hg_path,
      ncores = ncore
    )
  )
}
