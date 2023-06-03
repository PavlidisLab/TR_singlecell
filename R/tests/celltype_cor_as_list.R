## Generate cell type correlation into a list for inspection of intermediates of
## aggregation.




# Return a df of the max cor for each cell type for the given gene, removing NAs
# and the gene itself to prevent cor=1

max_cor_df <- function(cmat_list, tf) {
  
  cor_max <- lapply(names(cmat_list), function(x) {
    
    cor_mat <- cmat_list[[x]]
    
    vec <- cor_mat[tf, setdiff(colnames(cor_mat), tf)]
    
    if (all(is.na(vec))) {
      return(NA)
    }
    
    data.frame(
      Cell_type = x,
      Value = max(vec, na.rm = TRUE),
      Symbol = names(vec)[which.max(vec)])
  })
  
  cor_max <- cor_max[!is.na(cor_max)]
  
  cor_max <- data.frame(do.call(rbind, cor_max)) %>% 
    arrange(desc(Value))
  
  return(cor_max)
}