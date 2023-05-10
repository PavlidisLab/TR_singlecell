## https://tanaylab.github.io/metacell/articles/a-basic_pbmc8k.html
## TODO: Understand scdb loading structure


library("metacell")
library("tidyverse")

options(timeout = 1000)  # default 60 was insufficient for download time


# Creates a dir for a project database 
test_dir <- "/space/scratch/amorin/metacell_test/"
if (!dir.exists(test_dir)) dir.create(test_dir)
scdb_init(test_dir, force_reinit = TRUE)


# Default behavior creates a figure dir where subsequent plot calls autosave
if (!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")


# Loading a test data set (10X 8k PBMCs) as 'tgScMat' object
mcell_import_scmat_10x("test", base_dir = "http://www.wisdom.weizmann.ac.il/~atanay/metac_data/pbmc_8k/")
mat <- scdb_mat("test")
print(dim(mat@mat)) # 28826  8276


# Histogram of UMI counts per cell for use in filtering poor cells
mcell_plot_umis_per_cell("test")


# The following is a demonstration of removing genes that may prove problematic
# for subsequent clustering

nms <- c(rownames(mat@mat), rownames(mat@ignore_gmat))

ig_genes <- c(grep("^IGJ", nms, value = TRUE), 
             grep("^IGH", nms, value = TRUE),
             grep("^IGK", nms, value = TRUE), 
             grep("^IGL", nms, value = TRUE))

bad_genes <-  unique(c(
  grep("^MT-", nms, value = TRUE),
  grep("^MTMR", nms, value = TRUE),
  grep("^MTND", nms, value = TRUE),
  "NEAT1",
  "TMSB4X",
  "TMSB10",
  ig_genes
))


# TODO: This is editing an object and replacing the matrix slot. How to access
# main object?

mcell_mat_ignore_genes(new_mat_id = "test",
                       mat_id = "test", 
                       ig_genes = bad_genes) 


# Filtering matrix in object to remove low UMI cells
# TODO: auto threshold rule for min UMI? also check histo for current 

mcell_mat_ignore_small_cells(new_mat_id = "test", 
                             mat_id = "test", 
                             min_umis = 800)

# TODO: what the hell is this? what is force and why is default FALSE

mcell_add_gene_stat(gstat_id = "test", 
                    mat_id = "test", 
                    force = TRUE)

# TODO: T_vm exact definition
# Selecting highly variable genes for downstream clustering

mcell_gset_filter_varmean(gset_id = "test_feats", 
                          gstat_id = "test", 
                          T_vm = 0.08, 
                          force_new = TRUE)

# TODO: T_tot, T_top3 exact definition

mcell_gset_filter_cov(gset_id = "test_feats", 
                      gstat_id = "test", 
                      T_tot = 100, 
                      T_top3 = 2)


# TODO: What is sz correlation?

mcell_plot_gstats(gstat_id = "test", gset_id = "test_feats")

# Balanced K-nn using highly variable genes ("test_feats")
# TODO: How to access variable genes? And firm up downsample

mcell_add_cgraph_from_mat_bknn(mat_id = "test", 
                               gset_id = "test_feats", 
                               graph_id = "test_graph",
                               K = 100,
                               dsamp = TRUE)

# Generate metacells

mcell_coclust_from_graph_resamp(
  coc_id = "test_coc500", 
  graph_id = "test_graph",
  min_mc_size = 20, 
  p_resamp = 0.75, 
  n_resamp = 500)


mcell_mc_from_coclust_balanced(
  coc_id = "test_coc500", 
  mat_id = "test",
  mc_id = "test_mc", 
  K = 30, 
  min_mc_size = 30, 
  alpha = 2)


mcell_plot_outlier_heatmap(mc_id="test_mc", mat_id = "test", T_lfc=3)


mcell_mc_split_filt(new_mc_id="test_mc_f", 
                    mc_id="test_mc", 
                    mat_id="test",
                    T_lfc=3, 
                    plot_mats=F)



mcell_gset_from_mc_markers(gset_id="test_markers", mc_id="test_mc_f")


marks_colors = read.table(system.file("extdata", "pbmc_mc_colorize.txt", package="metacell"), sep="\t", h=T, stringsAsFactors=F)
mc_colorize("test_mc_f", marker_colors=marks_colors)
mc = scdb_mc("test_mc_f")

table(mc@colors)

mcell_mc_plot_marks(mc_id="test_mc_f", gset_id="test_markers", mat_id="test")

lfp = log2(mc@mc_fp)
tail(sort(lfp["RUNX1",]))

#
gene_df <- data.frame(mc@e_gc) %>% rownames_to_column(var = "Symbol")
fc_df <- data.frame(log2(mc@mc_fp)) %>% rownames_to_column(var = "Symbol")

#

mcell_mc2d_force_knn(mc2d_id="test_2dproj",mc_id="test_mc_f", graph_id="test_graph")
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="test_2dproj")


mc_hc = mcell_mc_hclust_confu(mc_id="test_mc_f", graph_id="test_graph")


mc_sup = mcell_mc_hierarchy(mc_id="test_mc_f", mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="test_mc_f", 
                        graph_id="test_graph", 
                        mc_order=mc_hc$order, 
                        sup_mc = mc_sup, 
                        width=2800, 
                        heigh=2000, 
                        min_nmc=2)
