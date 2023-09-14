# Sample only the experiments in which a given TF was measured

source("R/05b_save_sampled_null_topk.R")


tf <- "Ascl1"

tf_exp <- names(which(msr_mm[tf, ] == 1))


set.seed(5)


topk_mm <- lapply(1:n_samps, function(x) {
  
  message(paste("Sample #", x, Sys.time()))
  
  sample_topk_intersect(agg_l = agg_tf_mm[tf_exp],
                        genes = tfs_mm$Symbol,
                        msr_mat = msr_mm,
                        k = k)
})


saveRDS(topk_mm, "~/scratch/R_objects/03-09-2023_topk_size_matched_null.RDS")