library(sceasy)

h5adPath <- commandArgs(trailingOnly = TRUE)[1]

sceasy::convertFormat(h5adPath, from="anndata", to="seurat", outFile = gsub(".h5ad", ".rds", h5adPath), main_layer = "data")
