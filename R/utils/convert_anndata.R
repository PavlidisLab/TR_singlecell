# TODO: can't be done yet because hdf5 not on servers

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html