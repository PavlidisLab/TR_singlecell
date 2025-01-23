## Collecting packages used. Note that most would have been installed individually.
## -----------------------------------------------------------------------------


options(repos = "http://cran.rstudio.com/")


packages <- c(
  "testthat",
  "parallel",
  "tidyverse",
  "BiocManager",
  "remotes",
  "Seurat",
  "umap",
  "pheatmap",
  "RcppTOML",
  "igraph",
  "RANN",
  "WeightedCluster",
  "corpcor",
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
  "qlcMatrix",
  "microbenchmark",
  "credentials",
  "rJava",
  "egg"
)


installed_packages <- packages %in% rownames(installed.packages())


if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}



BiocManager::install("WGCNA")
BiocManager::install("preprocessCore")
BiocManager::install("GenomicRanges")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("ComplexHeatmap")
BiocManager::install("biomaRt")

devtools::install_github('PavlidisLab/ermineR')