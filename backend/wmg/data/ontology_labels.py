import csv
import gzip
import json
import pathlib
from typing import IO, Optional

# TODO: Place this module into a common ontology util package with development_stage_ontology_mapping.py and
#  extract_ontology_terms_from_owl.py. https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues
#  /chanzuckerberg/single-cell-data-portal/2133

all_ontologies_json_file = "all_ontology.json.gz"

genes_files = [
    "genes_ercc.csv.gz",
    "genes_sars_cov_2.csv.gz",
    "genes_homo_sapiens.csv.gz",
    "genes_mus_musculus.csv.gz",
]

SUFFIXES_TO_STRIP = ["organoid", "cell culture"]


def ontology_term_label(ontology_term_id: str) -> Optional[str]:
    """
    Returns the label for an ontology term, given its id. This excludes gene ontology term, which are handled
    separately by gene_gene_term_label(). Return None if ontology term id is invalid.
    """
    global ontology_term_id_labels

    # this catches the organoid tissue edge case (e.g. UBERON:0000995 (organoid)) or the cell culture edge case
    # (e.g. CL:0000082 (cell culture))
    suffix_to_strip = None
    for suffix in SUFFIXES_TO_STRIP:
        if suffix in ontology_term_id:
            suffix_to_strip = suffix
            break

    if suffix_to_strip:
        ontology_term_id = ontology_term_id.split(f"({suffix_to_strip})")[0].strip()
    ontology_term_id_label = ontology_term_id_labels.get(ontology_term_id)
    if suffix_to_strip:
        ontology_term_id_label = ontology_term_id_label + f" ({suffix_to_strip})"

    return ontology_term_id_label


def gene_term_label(gene_ontology_term_id: str) -> Optional[str]:
    """
    Returns the label for a gene ontology term, given its id. Return None if ontology term id is invalid.
    """
    global gene_term_id_labels
    return gene_term_id_labels.get(gene_ontology_term_id)


def __load_ontologies() -> None:
    global ontology_term_id_labels

    ontology_term_id_labels = {}
    all_ontologies = json.loads(__open_ontology_resource(all_ontologies_json_file).read().decode("utf-8"))
    for _, ontology_dict in all_ontologies.items():
        for term_id, term in ontology_dict.items():
            ontology_term_id_labels[term_id] = term["label"]


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
    root_path = curr_path.parent.parent
    file_path = root_path.joinpath("common", "ontology_files", file)
    return gzip.open(file_path)


__load_ontologies()
__load_genes()
