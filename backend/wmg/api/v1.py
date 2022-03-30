from collections import defaultdict
from typing import Dict, List, Any

import connexion
from flask import jsonify
from pandas import DataFrame

from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.wmg.data.cube import load_cubes, WmgCubes
from backend.wmg.data.query import (
    build_gene_id_label_mapping,
    build_ontology_term_id_label_mapping,
    WmgQuery,
    WmgQueryCriteria,
    build_dot_plot_matrix,
)
from backend.wmg.data.schema import cube_non_indexed_dims


# TODO: add cache directives: no-cache (i.e. revalidate); impl etag
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2132
def primary_filter_dimensions():
    cubes: WmgCubes = load_cubes()
    qry = WmgQuery(cubes)

    # gene terms are grouped by organism, and represented as a nested lists in dict, keyed by organism
    organism_gene_ids: dict[str, List[str]] = qry.list_grouped_primary_filter_dimensions_term_ids(
        "gene_ontology_term_id", group_by_dim="organism_ontology_term_id"
    )
    organism_gene_terms = {
        organism_term_id: build_gene_id_label_mapping(gene_term_ids)
        for organism_term_id, gene_term_ids in organism_gene_ids.items()
    }

    result = dict(
        snapshot_id=cubes.snapshot_identifier,
        organism_terms=build_ontology_term_id_label_mapping(
            qry.list_primary_filter_dimension_term_ids("organism_ontology_term_id")
        ),
        tissue_terms=build_ontology_term_id_label_mapping(
            qry.list_primary_filter_dimension_term_ids("tissue_ontology_term_id")
        ),
        gene_terms=organism_gene_terms,
    )
    return jsonify(result)


def query():
    request = connexion.request.json
    criteria = WmgQueryCriteria(**request["filter"])

    cubes: WmgCubes = load_cubes()
    query = WmgQuery(cubes)
    query_result = query.expression_summary(criteria)
    cell_counts = query.cell_counts(criteria)
    dot_plot_matrix_df = build_dot_plot_matrix(query_result, cell_counts)
    all_filter_dims_values = extract_filter_dims_values(query_result)

    include_filter_dims = request.get("include_filter_dims", False)
    response_filter_dims_values = build_filter_dims_values(all_filter_dims_values) if include_filter_dims else {}

    cell_type_term_ids = all_filter_dims_values["cell_type_ontology_term_id"]

    return jsonify(
        dict(
            snapshot_id=cubes.snapshot_identifier,
            expression_summary=build_expression_summary(dot_plot_matrix_df),
            term_id_labels=dict(
                genes=build_gene_id_label_mapping(criteria.gene_ontology_term_ids),
                cell_types=build_ontology_term_id_label_mapping(cell_type_term_ids),
            ),
            filter_dims=response_filter_dims_values,
        )
    )


# TODO: Read this from generated data artifact instead of DB.
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2086. This code is without a unit test, but we are intending to replace it.
def fetch_datasets_metadata(dataset_ids: List[str]) -> List[Dict]:
    # We return a DTO because the db entity can't access its attributes after the db session ends,
    # and we want to keep session management out of the calling method

    with db_session_manager() as session:

        def get_dataset(dataset_id_):
            dataset = Dataset.get(session, dataset_id_)
            if dataset is None:
                # Handle possible missing dataset due to db state evolving past wmg snapshot
                return dict(id=dataset_id_, label="", collection_id="", collection_label="")
            return dict(
                id=dataset.id,
                label=dataset.name,
                collection_id=dataset.collection.id,
                collection_label=dataset.collection.name,
            )

        return [get_dataset(dataset_id) for dataset_id in dataset_ids]


def build_filter_dims_values(all_filter_dims_values):
    response_filter_dims_values = dict(
        datasets=fetch_datasets_metadata(all_filter_dims_values["dataset_id"]),
        disease_terms=build_ontology_term_id_label_mapping(all_filter_dims_values["disease_ontology_term_id"]),
        sex_terms=build_ontology_term_id_label_mapping(all_filter_dims_values["sex_ontology_term_id"]),
        development_stage_terms=build_ontology_term_id_label_mapping(
            all_filter_dims_values["development_stage_ontology_term_id"]
        ),
        ethnicity_terms=build_ontology_term_id_label_mapping(all_filter_dims_values["ethnicity_ontology_term_id"]),
        # excluded per product requirements, but keeping in, commented-out, to reduce future head-scratching
        # assay_ontology_terms=build_ontology_term_id_label_mapping(all_filter_dims_values["assay_ontology_term_id"]),
    )
    return response_filter_dims_values


def build_expression_summary(query_result: DataFrame) -> dict:
    # Create nested dicts with gene_ontology_term_id, tissue_ontology_term_id keys, respectively
    structured_result: Dict[str, Dict[str, List[Dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for row in query_result.itertuples(index=False):
        structured_result[row.gene_ontology_term_id][row.tissue_ontology_term_id].append(
            dict(
                id=row.cell_type_ontology_term_id,
                n=row.nnz,
                me=row.sum / row.nnz,
                pc=row.nnz / row.n_cells_cell_type,
                tpc=row.nnz / row.n_cells_tissue,
            )
        )
    return structured_result


def extract_filter_dims_values(query_result: DataFrame) -> dict:
    """
    Return unique values for each dimension in the specified query result
    """
    dims: set = set(query_result.columns) & set(cube_non_indexed_dims)
    dim_uniq_values: dict = {dim: query_result.groupby(dim).groups.keys() for dim in dims}
    return dim_uniq_values
