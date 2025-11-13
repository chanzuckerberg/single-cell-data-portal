import csv
import gzip
import pathlib
from typing import IO, Optional

from backend.common.census_cube.utils import ontology_parser

genes_files = [
    "genes_ercc.csv.gz",
    "genes_sars_cov_2.csv.gz",
    "genes_homo_sapiens.csv.gz",
    "genes_mus_musculus.csv.gz",
]


def ontology_term_label(ontology_term_id: str) -> Optional[str]:
    """
    Returns the label for an ontology term, given its id. This excludes gene ontology term, which are handled
    separately by gene_gene_term_label(). Return None if ontology term id is invalid.
    """
    try:
        return ontology_parser.get_term_label(ontology_term_id)
    # If the ontology term id is invalid or not found in the ontology schema,
    # return the ontology term id itself. This is useful for cases like publication
    # citation and newly added cell types not yet in the ontology schema.
    except (ValueError, KeyError):
        return ontology_term_id


def is_ontology_term_deprecated(ontology_term_id: str) -> bool:
    """
    Returns whether an ontology term is deprecated. Returns False for unknown terms.
    """
    try:
        return ontology_parser.is_term_deprecated(ontology_term_id)
    except (ValueError, KeyError):
        # Assume unknown terms are not deprecated
        return False


def ontology_term_description(ontology_term_id: str) -> Optional[str]:
    """
    Returns the description for an ontology term, given its id. Returns None if the term is not found.
    """
    try:
        return ontology_parser.get_term_description(ontology_term_id)
    except (ValueError, KeyError):
        # Return None for unknown terms (no description available)
        return None


def ontology_term_synonyms(ontology_term_id: str) -> list[str]:
    """
    Returns the synonyms for an ontology term, given its id. Returns empty list if the term is not found.
    """
    try:
        return ontology_parser.get_term_synonyms(ontology_term_id)
    except (ValueError, KeyError):
        # Return empty list for unknown terms (no synonyms available)
        return []


def gene_term_label(gene_ontology_term_id: str) -> Optional[str]:
    """
    Returns the label for a gene ontology term, given its id. Return None if ontology term id is invalid.
    """
    global gene_term_id_labels
    return gene_term_id_labels.get(gene_ontology_term_id, gene_ontology_term_id)


def __load_genes() -> None:
    global gene_term_id_labels

    gene_term_id_labels = {}
    for genes_file in genes_files:
        for row in csv.DictReader(
            __open_ontology_resource(genes_file).read().decode("utf-8").split("\n"),
            fieldnames=["term_id", "label", "version"],
        ):
            gene_term_id_labels[row["term_id"]] = row["label"]


def __open_ontology_resource(file) -> IO:
    curr_path = pathlib.Path(__file__).parent.absolute()
    root_path = curr_path.parent.parent.parent
    file_path = root_path.joinpath("common", "ontology_files", file)
    return gzip.open(file_path)


__load_genes()
