## Collecting packages used. Note that most would have been installed individually.
## -----------------------------------------------------------------------------

if (!require("BiocManager")) install.packages('BiocManager')

if (!requireNamespace("remotes")) install.packages("remotes")


# BiocManager::install("tanaylab/metacell")
 
# remotes::install_github("GfellerLab/SuperCell")

install.packages("Seurat")
install.packages("WGCNA")
install.packages("umap")
install.packages("RcppTOML")
install.packages("igraph")
install.packages("RANN")
install.packages("WeightedCluster")
install.packages("corpcor")
install.packages("weights")
install.packages("Hmisc")
install.packages("Matrix")
install.packages("patchwork")
install.packages("plyr")
install.packages("irlba")
install.packages("DescTools")
install.packages("lobstr")
install.packages("googlesheets4")
install.packages("ROCR")
install.packages("data.table")
install.packages("ggrepel")
