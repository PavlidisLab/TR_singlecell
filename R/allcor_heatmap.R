library(pheatmap)
library(tidyverse)

cor_l <- readRDS("~/Plots/TR_singlecell/All_corplots/Mouse_Ascl1_Dll3.RDS")



# Why does this not work?
# cor_df <- rbind(do.call, lapply(names(cor_l), function(x) {
#   data.frame(ID = x,
#              Cell_type = cor_l[[x]],
#              Cor = cor_l[[x]])
# }))


cor_df <- lapply(names(cor_l), function(x) {
  data.frame(ID = x,
             Cell_type = names(cor_l[[x]]),
             Cor = cor_l[[x]])
}) %>% 
  do.call(rbind, .) %>% 
  arrange(Cor)


pheatmap(t(cor_df[, "Cor", drop = FALSE]),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("#0571b0", "white", "#ca0020"))(100),
         breaks = seq(-1, 1, length.out =  100),
         border_color = NA,
         cellheight = 10)
