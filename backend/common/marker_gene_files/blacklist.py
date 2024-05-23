import os

# This file contains blacklisted marker genes that match the following criteria:
# Ensembl 104 gene IDs that are annotated with the GO terms: ['GO:0005840', 'GO:0005739']
# These GO terms correspond to the cellular components "ribosome" and "mitochondrion".
file_dir = os.path.dirname(os.path.realpath(__file__))
MARKER_GENE_BLACKLIST_FILENAME = os.path.join(file_dir, "marker_gene_blacklist.txt")

marker_gene_blacklist = []
with open(MARKER_GENE_BLACKLIST_FILENAME, "r") as f:
    marker_gene_blacklist = f.read().split(",")
