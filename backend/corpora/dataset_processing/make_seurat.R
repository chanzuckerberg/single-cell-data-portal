library(sceasy)

require(devtools)

h5adPath <- commandArgs(trailingOnly = TRUE)[1]

sceasy::convertFormat(h5adPath,
                      from="anndata",
                      to="seurat",
                      outFile = gsub(".h5ad", ".rds", h5adPath),
                      main_layer = "data",
                      target_uns_keys = c("schema_version", "title", "batch_condition", "default_embedding", "X_approximate_distribution"))
