import os
from dataclasses import dataclass

import pandas as pd

file_dir = os.path.dirname(os.path.realpath(__file__))
ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME = os.path.join(file_dir, "ensembl_gene_ids_to_descriptions.tsv.gz")


@dataclass
class GeneMetadata:
    gene_id_to_name: dict[str, str]
    gene_id_to_symbol: dict[str, str]


# This file contains blacklisted marker genes that match the following criteria:
# Ensembl 104 gene IDs that are annotated with the GO terms: ['GO:0005840', 'GO:0005739']
# These GO terms correspond to the cellular components "ribosome" and "mitochondrion".
file_dir = os.path.dirname(os.path.realpath(__file__))
MARKER_GENE_BLACKLIST_FILENAME = os.path.join(file_dir, "marker_gene_blacklist.txt")

marker_gene_blacklist = []
with open(MARKER_GENE_BLACKLIST_FILENAME, "r") as f:
    marker_gene_blacklist = f.read().split(",")


def get_gene_id_to_name_and_symbol() -> GeneMetadata:
    """
    This function reads a file containing gene metadata and returns a GeneMetadata object.
    The GeneMetadata object contains two dictionaries: gene_id_to_name and gene_id_to_symbol.
    gene_id_to_name is a dictionary where the keys are Ensembl GeneIDs and the values are gene descriptions.
    gene_id_to_symbol is a dictionary where the keys are Ensembl GeneIDs and the values are gene symbols.

    Returns
    -------
    GeneMetadata
        An object containing two dictionaries: gene_id_to_name and gene_id_to_symbol.
    """

    gene_metadata = pd.read_csv(ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME, sep="\t")
    gene_id_to_name = gene_metadata.set_index("Ensembl GeneIDs")["Description"].to_dict()
    gene_id_to_symbol = gene_metadata.set_index("Ensembl GeneIDs")["Symbols"].to_dict()
    return GeneMetadata(gene_id_to_name=gene_id_to_name, gene_id_to_symbol=gene_id_to_symbol)
