## TODO
## -----------------------------------------------------------------------------

# https://github.com/PavlidisLab/ermineR
# devtools::install_github('PavlidisLab/ermineR')
# install.packages("rJava")
library(ermineR)
library(rJava)
library(tidyverse)
library(parallel)
source("R/00_config.R")

erminer_coexpr_hg_path <- "/space/scratch/amorin/R_objects/coexpr_erminer_pr_hg.RDS"
erminer_coexpr_mm_path <- "/space/scratch/amorin/R_objects/coexpr_erminer_pr_mm.RDS"
erminer_coexpr_ortho_path <- "/space/scratch/amorin/R_objects/coexpr_erminer_pr_ortho.RDS"


go_date <- "2024-09-08"
go_path <- paste0("/space/scratch/amorin/R_objects/go_terms_", go_date, ".xml")

if (!file.exists(go_path)) {
  goAtDate(go_path, go_date, overwrite = FALSE)
}


anno_hg_path <- "/space/scratch/amorin/R_objects/gemma_generic_human_anno"
anno_mm_path <- "/space/scratch/amorin/R_objects/gemma_generic_mouse_anno"


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
                          aspects = "B",
                          iterations = 10000,
                          bigIsBetter = FALSE)
        
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



# Note that ortho (ranking placing equal weight on mouse and human aggregate)
# uses human symbols

if (!file.exists(erminer_coexpr_hg_path)) {
  
  message(paste("Human", Sys.time()))
  
  coexpr_hg <- readRDS(rank_tf_hg_path)
  
  erminer_coexpr_hg <- erminer_enrich_list(dat_l = coexpr_hg, 
                                           score_col = "Rank_aggr_coexpr",
                                           go_path = go_path,
                                           anno_path = anno_hg_path,
                                           ncores = ncore)
  
  
  # Note I did this interactively to inspect failures (which are random)
  which_na <- which(is.na(erminer_coexpr_hg))
  
  
  erminer_coexpr_hg[which_na] <- erminer_enrich_list(dat_l = coexpr_hg[which_na], 
                                                     score_col = "Rank_aggr_coexpr",
                                                     go_path = go_path,
                                                     anno_path = anno_hg_path,
                                                     ncores = ncore)
  
  
  saveRDS(erminer_coexpr_hg, erminer_coexpr_hg_path)
}




if (!file.exists(erminer_coexpr_mm_path)) {
  
  message(paste("Mouse", Sys.time()))
  
  coexpr_mm <- readRDS(rank_tf_mm_path)
  
  erminer_coexpr_mm <- erminer_enrich_list(dat_l = coexpr_mm, 
                                           score_col = "Rank_aggr_coexpr",
                                           go_path = go_path,
                                           anno_path = anno_mm_path,
                                           ncores = ncore)
  
  
  # Note I did this interactively to inspect failures (which are random)
  which_na <- which(is.na(erminer_coexpr_mm))
  
  
  erminer_coexpr_mm[which_na] <- erminer_enrich_list(dat_l = coexpr_mm[which_na], 
                                                     score_col = "Rank_aggr_coexpr",
                                                     go_path = go_path,
                                                     anno_path = anno_mm_path,
                                                     ncores = ncore)
  
  
  saveRDS(erminer_coexpr_mm, erminer_coexpr_mm_path)
}




if (!file.exists(erminer_coexpr_ortho_path)) {
  
  message(paste("Ortho", Sys.time()))
  
  coexpr_ortho <- readRDS(rank_tf_ortho_path)
  
  # Coerce to expected column names
  coexpr_ortho <- lapply(coexpr_ortho, function(x) {
    dplyr::rename(x, Symbol = Symbol_hg, Rank_aggr_coexpr = Rank_aggr_coexpr_ortho)
  })
  
  erminer_coexpr_ortho <- erminer_enrich_list(dat_l = coexpr_ortho, 
                                              score_col = "Rank_aggr_coexpr",
                                              go_path = go_path,
                                              anno_path = anno_hg_path,
                                              ncores = ncore)
  
  # Note I did this interactively to inspect failures (which are random)
  which_na <- which(is.na(erminer_coexpr_ortho))
  
  
  erminer_coexpr_ortho[which_na] <- erminer_enrich_list(dat_l = coexpr_ortho[which_na], 
                                                        score_col = "Rank_aggr_coexpr",
                                                        go_path = go_path,
                                                        anno_path = anno_hg_path,
                                                        ncores = ncore)
  
  saveRDS(erminer_coexpr_ortho, erminer_coexpr_ortho_path)
}
