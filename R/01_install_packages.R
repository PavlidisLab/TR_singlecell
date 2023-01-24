if (!require("BiocManager")) install.packages('BiocManager')
BiocManager::install("tanaylab/metacell")


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


if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("GfellerLab/SuperCell")

