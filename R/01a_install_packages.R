## Collecting packages used. Note that most would have been installed individually.
## -----------------------------------------------------------------------------



packages <- c(
  "testthat",
  "parallel",
  "tidyverse",
  "BiocManager",
  "remotes",
  "Seurat",
  "WGCNA",
  "umap",
  "pheatmap",
  "RcppTOML",
  "igraph",
  "RANN",
  "WeightedCluster",
  "corpcor",
  "weights",
  "Hmisc",
  "Matrix",
  "patchwork",
  "plyr",
  "irlba",
  "DescTools",
  "googlesheets4",
  "ROCR",
  "data.table",
  "ggrepel",
  "proxyC",
  "qlcMatrix"
)


installed_packages <- packages %in% rownames(installed.packages())


if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}



BiocManager::install("preprocessCore")
BiocManager::install("GenomicRanges")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("ComplexHeatmap")
