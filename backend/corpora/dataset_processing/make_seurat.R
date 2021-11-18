library(sceasy)

require(devtools)
install_version("SeuratObject", version = "4.0.2", repos = "http://cran.rstudio.com/")

h5adPath <- commandArgs(trailingOnly = TRUE)[1]

sceasy::convertFormat(h5adPath, from="anndata", to="seurat", outFile = gsub(".h5ad", ".rds", h5adPath), main_layer = "data")
