source("R/utils/vector_comparison_functions.R")

test_gene <- "RPL3"
test_mat <- agg_hg$GSE180928
test_vec <- agg_hg$GSE216019[, test_gene]




test_gene_similarity_matrix <- 
  gene_similarity_matrix(agg_l = agg_hg,
                         gene = test_gene,
                         msr = "Topk",
                         ncores = ncore)

