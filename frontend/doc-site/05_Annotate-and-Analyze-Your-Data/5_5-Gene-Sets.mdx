# Gene Sets

Cellxgene Desktop enables cellxgene users to pre-load gene sets for exploration, and to create or modify gene sets while exploring data. This is particularly useful for recording cell type markers when annotating cells THIS SHOULD LINK TO THE DESKTOP ANNOTATION CELLS PAGE. This section explains how to format gene sets data for preloading of gene sets and where to find gene sets information after making changes in the explorer.

To read more about how to use gene sets in the cellxgene interface, please see Exploring Gene Expression. LINK TO THIS GITHUB PAGE IS NOT WORKING FOR ME

## Gene Sets Requirements

A gene set must have a unique name. A given gene may be included in multiple gene sets. A gene set may optionally include a description for the entire gene set, and descriptions for each gene included the gene set.


## `cellxgene` CLI Command

Users MAY use `--gene-sets-file` `name_of_file.csv` to designate a CSV with preexisting gene sets you would like to view in cellxgene. The CSV MUST follow the format below.

Users MAY use `--user-generated-data-dir` `file_name` to designate a file where both gene sets and annotations CSVs will be saved. File names for gene sets and annotations will contain a pseudo-session ID by default. `--user-generated-data-dir` is incompatible with `--gene-sets-file` and `--annotations-file` and will error if used together. See `--help` for more details.

## `cellxgene` Gene Set CSV Data Format

The cellxgene gene set data format is a Tidy CSV (comma-separated values) THIS GITHUB LINK IS NOT WORKING FOR ME file using ASCII encoding. Multiple gene sets MAY be included in the file.

The first row MUST contain column headers in the following order:
`gene_set_name`
`gene_set_description`
`gene_symbol`
`gene_description`

Each row represents a gene in a gene set and must be unique.

Example:

EMBED TABLE


Users MAY include additional columns. Once gene sets are edited from cellxgene, the additional columns will no longer be stored in the file designated using the command `--gene-sets-file`.

## `gene_set_name`

The` gene_set_name` column MUST contain a value and MUST NOT contain the following ASCII characters or sequences:

control characters (decimal 0-31)
DEL (decimal 127)
leading spaces (decimal 32) in a field - " This is an example"
trailing spaces (decimal 32) in a field - "This is an example "
multiple spaces (decimal 32) "internal" to a field - "This is an example"

Note: If gene_symbol(s) for a gene_set_name exist on noncontiguous rows, they will be added to the existing gene set. For example, CD163 below is added to the club.cell gene set:

EMBED TABLE


## `gene_set_description`
Populating `gene_set_description` is optional. The first instance where `gene_set_description` is populated for a specific `gene_set_name` will be surfaced when a user hovers over `gene_set_name` in cellxgene. All other instances are ignored in subsequent rows for the same `gene_set_name`.

## `gene_symbol`
A given `gene_symbol` may only be added once to a gene set and exist as a VAR in the underlying anndata file.

## `gene_description`
Populating `gene_description` is optional and will be surfaced when a user hovers on `gene_symbol` in cellxgene.
