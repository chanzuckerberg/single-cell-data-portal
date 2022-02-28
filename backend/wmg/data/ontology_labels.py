import csv
import gzip
import json
from typing import IO

# TODO: Place this module into a common ontology util package with ontology_mapping.py and
#  extract_ontology_terms_from_owl.py

all_ontologies_json_file = "all_ontology.json.gz"

genes_files = [
    "genes_ercc.csv.gz",
    "genes_sars_cov_2.csv.gz",
    "genes_homo_sapiens.csv.gz",
    "genes_mus_musculus.csv.gz",
]


def load_ontologies() -> dict:
    all_ontologies = json.loads(open_ontology_resource(all_ontologies_json_file).read().decode('utf-8'))
    ontology_term_id_labels = {}
    for _, ontology in all_ontologies.items():
        for term_id, term in ontology.items():
            ontology_term_id_labels[term_id] = term['label']
    return ontology_term_id_labels


def load_genes() -> dict:
    gene_term_id_labels = {}
    for genes_file in genes_files:
        for row in csv.DictReader(open_ontology_resource(genes_file).read().decode('utf-8').split('\n'),
                                  fieldnames=["term_id", "label", "version"]):
            gene_term_id_labels[row['term_id']] = row['label']
    return gene_term_id_labels


def open_ontology_resource(file) -> IO:
    from importlib import resources

    path = resources.files('cellxgene_schema').joinpath(f'ontology_files/{file}')
    return gzip.open(path)
