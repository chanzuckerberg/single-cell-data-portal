import tiledb

from backend.wmg.data.ontology_labels import get_ontology_term_label, get_gene_term_label
from typing import Dict, List, Callable
import json

from backend.wmg.data.snapshot import (
    EXPRESSION_SUMMARY_CUBE_NAME,
    PRIMARY_FILTER_DIMENSIONS_FILENAME,
)


def list_primary_filter_dimension_term_ids(cube, primary_dim_name: str):
    return cube.query(attrs=[], dims=[primary_dim_name]).df[:].groupby([primary_dim_name]).first().index.tolist()


def list_grouped_primary_filter_dimensions_term_ids(
    cube, primary_dim_name: str, group_by_dim: str
) -> Dict[str, List[str]]:
    return (
        cube.query(attrs=[], dims=[primary_dim_name, group_by_dim])
        .df[:]
        .drop_duplicates()
        .groupby(group_by_dim)
        .agg(list)
        .to_dict()[primary_dim_name]
    )


def build_id_label_mapping(ontology_term_ids: List[str], label_lookup_function: Callable) -> List[dict]:
    ontology_term_id_label_map = []
    for ontology_term_id in ontology_term_ids:
        readable_label = label_lookup_function(ontology_term_id)
        ontology_term_id_label_map.append({ontology_term_id: readable_label})
    return ontology_term_id_label_map


def generate_gene_terms_by_organism(cube) -> dict[str, List[str]]:
    """
        Gene terms are grouped by organism, and represented as a nested lists in dict, keyed by organism
    {
      "organism_0_id": [{"gene_0_id": "gene_0_label"}, {"gene_1_id": "gene_1_label"}, {"gene_2_id": "gene_2_label"}],
      "organism_1_id": [{"gene_3_id": "gene_3_label"}, {"gene_4_id": "gene_4_label"}, {"gene_5_id": "gene_5_label"}],
      "organism_2_id": [{"gene_7_id": "gene_7_label"}, {"gene_8_id": "gene_8_label"}, {"gene_9_id": "gene_9_label"}]
      }
    """
    organism_gene_terms = []
    organism_gene_ids = list_grouped_primary_filter_dimensions_term_ids(
        cube, "gene_ontology_term_id", group_by_dim="organism_ontology_term_id"
    )

    for organism_term_id, gene_term_ids in organism_gene_ids.items():
        organism_gene_terms.append({organism_term_id: build_id_label_mapping(gene_term_ids, get_gene_term_label)})
    return organism_gene_terms


def generate_tissue_terms_by_organism(cube) -> dict[str, List[str]]:
    """
        Tissue terms are grouped by organism, and represented as a nested lists in dict, keyed by organism
    {
      "organism_0_id": [{"tissue_0_id": "tissue_0_label"}, {"tissue_1_id": "tissue_1_label"}],
      "organism_1_id": [{"tissue_3_id": "tissue_3_label"}, {"tissue_4_id": "tissue_4_label"}],
      "organism_2_id": [{"tissue_7_id": "tissue_7_label"}, {"tissue_8_id": "tissue_8_label"}]
      }
    """
    organism_tissue_terms = []
    organism_tissue_ids = list_grouped_primary_filter_dimensions_term_ids(
        cube, "tissue_ontology_term_id", group_by_dim="organism_ontology_term_id"
    )
    for organism_term_id, tissue_term_ids in organism_tissue_ids.items():
        organism_tissue_terms.append(
            {organism_term_id: build_id_label_mapping(tissue_term_ids, get_ontology_term_label)}
        )
    return organism_tissue_terms


def generate_primary_filter_dimensions(snapshot_path: str, corpus_name: str, snapshot_id: int):
    with tiledb.open(f"{snapshot_path}/{corpus_name}/{EXPRESSION_SUMMARY_CUBE_NAME}") as cube:
        organism_term_ids = list_primary_filter_dimension_term_ids(cube, "organism_ontology_term_id")
        organism_terms = build_id_label_mapping(organism_term_ids, get_ontology_term_label)

        organism_gene_terms = generate_gene_terms_by_organism(cube)
        organism_tissue_terms = generate_tissue_terms_by_organism(cube)

        result = dict(
            snapshot_id=snapshot_id,
            organism_terms=organism_terms,
            tissue_terms=organism_tissue_terms,
            gene_terms=organism_gene_terms,
        )

        with open(f"{snapshot_path}/{PRIMARY_FILTER_DIMENSIONS_FILENAME}", "w") as f:
            json.dump(result, f)
