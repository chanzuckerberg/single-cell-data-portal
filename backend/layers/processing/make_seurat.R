install.packages('devtools')
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rhdf5")
devtools::install_github("cellgeni/schard")

library(sceasy)

require(devtools)

h5adPath <- '/Users/trentsmith/workspace/single-cell-curation/cellxgene_schema_cli/tests/fixtures/h5ads/example_valid.h5ad'

# load h5ad as Seurat
converted = schard::h5ad2seurat('/Users/trentsmith/workspace/single-cell-curation/cellxgene_schema_cli/tests/fixtures/h5ads/example_valid.h5ad')
saveRDS(converted, '/Users/trentsmith/workspace/single-cell-curation/cellxgene_schema_cli/tests/fixtures/h5ads/example_valid.rds')


target_uns_keys <- c("schema_version",
                     "title",
                     "batch_condition",
                     "default_embedding",
                     "X_approximate_distribution",
                     "citation",
                     "schema_reference"
                    )

sceasy::convertFormat(h5adPath,
                      from="anndata",
                      to="seurat",
                      outFile = gsub(".h5ad", ".rds", h5adPath),
                      main_layer = "data",
                      target_uns_keys = target_uns_keys)

