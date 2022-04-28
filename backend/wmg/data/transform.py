import tiledb

from backend.wmg.data.ontology_labels import ontology_term_label, gene_term_label
from typing import Dict, List, Iterable
import json

from backend.wmg.data.snapshot import (
    EXPRESSION_SUMMARY_CUBE_NAME,
    PRIMARY_FILTER_DIMENSIONS_FILENAME,
)


def generate_primary_filter_dimensions(snapshot_path: str, corpus_name: str, snapshot_id: int):

    # TODO: remove them from WmgQuery (next 4 following functions)
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

    def build_gene_id_label_mapping(gene_ontology_term_ids: List[str]) -> List[dict]:
        return [
            {gene_ontology_term_id: gene_term_label(gene_ontology_term_id)}
            for gene_ontology_term_id in gene_ontology_term_ids
        ]

    def build_ontology_term_id_label_mapping(ontology_term_ids: Iterable[str]) -> List[dict]:
        return [{ontology_term_id: ontology_term_label(ontology_term_id)} for ontology_term_id in ontology_term_ids]

    with tiledb.open(f"{snapshot_path}/{corpus_name}/{EXPRESSION_SUMMARY_CUBE_NAME}") as cube:

        # gene terms are grouped by organism, and represented as a nested lists in dict, keyed by organism
        organism_gene_ids: dict[str, List[str]] = list_grouped_primary_filter_dimensions_term_ids(
            cube, "gene_ontology_term_id", group_by_dim="organism_ontology_term_id"
        )
        organism_gene_terms = {
            organism_term_id: build_gene_id_label_mapping(gene_term_ids)
            for organism_term_id, gene_term_ids in organism_gene_ids.items()
        }

        # tissue terms are grouped by organism, and represented as a nested lists in dict, keyed by organism
        organism_tissue_ids: dict[str, List[str]] = list_grouped_primary_filter_dimensions_term_ids(
            cube, "tissue_ontology_term_id", group_by_dim="organism_ontology_term_id"
        )
        organism_tissue_terms = {
            organism_term_id: build_ontology_term_id_label_mapping(tissue_term_ids)
            for organism_term_id, tissue_term_ids in organism_tissue_ids.items()
        }

        result = dict(
            snapshot_id=snapshot_id,
            organism_terms=build_ontology_term_id_label_mapping(
                list_primary_filter_dimension_term_ids(cube, "organism_ontology_term_id")
            ),
            tissue_terms=organism_tissue_terms,
            gene_terms=organism_gene_terms,
        )

        with open(f"{snapshot_path}/{PRIMARY_FILTER_DIMENSIONS_FILENAME}", "w") as f:
            json.dump(result, f)
