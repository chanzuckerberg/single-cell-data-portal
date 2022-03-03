import csv
import gzip
import json
from collections import defaultdict
from typing import IO, Dict

# TODO: Place this module into a common ontology util package with ontology_mapping.py and
#  extract_ontology_terms_from_owl.py

all_ontologies_json_file = "all_ontology.json.gz"

genes_files = [
    "genes_ercc.csv.gz",
    "genes_sars_cov_2.csv.gz",
    "genes_homo_sapiens.csv.gz",
    "genes_mus_musculus.csv.gz",
]
ontology_term_id_labels: Dict[str, str] = None
gene_term_id_labels: Dict[str, str] = None


def ontology_term_label(ontology_term_id: str) -> str:
    if ontology_term_id_labels is None:
        __load_ontologies()

    return ontology_term_id_labels.get(ontology_term_id)


def gene_term_label(gene_ontology_term_id: str) -> str:
    if gene_term_id_labels is None:
        __load_genes()

    return gene_term_id_labels.get(gene_ontology_term_id)


def __load_ontologies() -> None:
    global ontology_term_id_labels

    ontology_term_id_labels = {}
    all_ontologies = json.loads(__open_ontology_resource(all_ontologies_json_file).read().decode('utf-8'))
    for ontology_name, ontology_dict in all_ontologies.items():
        for term_id, term in ontology_dict.items():
            ontology_term_id_labels[term_id] = term['label']


def __load_genes() -> None:
    global gene_term_id_labels

    gene_term_id_labels = {}
    for genes_file in genes_files:
        for row in csv.DictReader(__open_ontology_resource(genes_file).read().decode('utf-8').split('\n'),
                                  fieldnames=["term_id", "label", "version"]):
            gene_term_id_labels[row['term_id']] = row['label']


def __open_ontology_resource(file) -> IO:
    from importlib import resources

    path = resources.files('cellxgene_schema').joinpath(f'ontology_files/{file}')
    return gzip.open(path)
