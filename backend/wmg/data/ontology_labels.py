import csv
import gzip
import json
import logging
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

logger = logging.getLogger("wmg-ontology-labels")


def ontology_term_label(ontology_term_id: str) -> Optional[str]:
    """
    Returns the label for an ontology term, given its id. This excludes gene ontology term
    and self reported ethnicity term, which are handled separately by gene_term_label()
    and ethnicity_term_label() respectively. Return None if ontology term id is invalid.
    """

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


def ethnicity_term_label(self_reported_ethnicity_ontology_term_id: str) -> str:
    """
    Return the label for self reported ethnicity ontology term id.

    This function is compatible with schema-4 and schema-3 format of
    `self_reported_ethnicity_ontology_term_id`.

    NOTE: It is assumed that self_reported_ethnicity_ontology_term_id has been validated against
    a schema (ex: schema-3 or schema-4) before this function is called. Therefore,
    no error checking on the format of the input is done.

    Parameters
    ----------
    self_reported_ethnicity_ontology_term_id: An ontology_term_id string that can
                                              encode multiple ethnicities with comma delimiter
                                              in Schema-4.
                                              NOTE: Schema-3 does not support comma delimited
                                              ontology term IDs to encode multiple ethnicities
                                              in a single string.

    Returns
    -------
    ethnicity_term_label: A comma delimited ethnicity label corresponding to the input ethnicity_term_id
    """
    # In schema-4, self_reported_ethnicity_ontology_term_id can be a comma
    # separated string to denote multiple ethnicities.
    individual_term_ids = self_reported_ethnicity_ontology_term_id.split(",")
    logger.info(
        f"PRATHAP! ethnicity_term_label() func - self_reported_ethnicity_ontology_term_id: {self_reported_ethnicity_ontology_term_id}"
    )
    logger.info(f"PRATHAP! ethnicity_term_label() func - individual_term_ids: {individual_term_ids}")
    term_labels = [ontology_term_id_labels.get(term_id) for term_id in individual_term_ids]
    logger.info(f"PRATHAP! ethnicity_term_label() func - term_labels: {term_labels}")
    ethnicity_term_label = ",".join(term_labels)
    return ethnicity_term_label


def gene_term_label(gene_ontology_term_id: str) -> Optional[str]:
    """
    Returns the label for a gene ontology term, given its id. Return None if ontology term id is invalid.
    """
    return gene_term_id_labels.get(gene_ontology_term_id)


def __load_ontologies() -> None:
    global ontology_term_id_labels

    all_ontologies = json.loads(__open_ontology_resource(all_ontologies_json_file).read().decode("utf-8"))
    for _, ontology_dict in all_ontologies.items():
        for term_id, term in ontology_dict.items():
            ontology_term_id_labels[term_id] = term["label"]


def __load_genes() -> None:
    global gene_term_id_labels

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


ontology_term_id_labels = {}
gene_term_id_labels = {}

__load_ontologies()
__load_genes()
