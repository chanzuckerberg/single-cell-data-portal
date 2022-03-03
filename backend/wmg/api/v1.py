import logging
from collections import defaultdict
from typing import Dict, List
from uuid import uuid4

import connexion
import tiledb
from flask import jsonify
from pandas import DataFrame

from backend.wmg.data.config import fast_config, create_ctx
from backend.wmg.data.query import (
    build_gene_id_label_mapping,
    build_ontology_term_id_label_mapping,
    WmgQuery,
    WmgQueryCriteria, build_dot_plot_matrix,
)
# TODO: Replace with real snapshot uuid
from backend.wmg.data.schema import cube_non_indexed_dims

DUMMY_SNAPSHOT_UUID = uuid4().hex


# TODO: add cache directives
def primary_filter_dimensions():
    organism_terms = [dict(oid1="olbl1"), dict(oid2="olbl2")]
    tissue_terms = [dict(ttid1="ttlbl1"), dict(ttid2="ttlbl2")]
    result = dict(snapshot_id=DUMMY_SNAPSHOT_UUID, organism_terms=organism_terms, tissue_terms=tissue_terms)
    return jsonify(result)


cube = None


def find_cube_latest_snapshot():
    global cube

    if cube is None:
        # TODO: Remove tiledb dependency from this module
        # TODO: Okay to keep open indefinitely? Is it faster than re-opening each request?
        cube_uri = build_cube_uri()
        logging.info(f"Opening WMG cube at {cube_uri}")
        cube = tiledb.open(cube_uri, ctx=create_ctx(fast_config()))

    return cube


def build_cube_uri():
    # TODO: Retrieve from app config
    cube_base_uri = "s3://wmg-dev"
    # TODO: Retrieve from s3://wmg-<env>/latest_snapshot_uuid
    cube_latest_snapshot = "dummy-snapshot"
    cube_uri = f"{cube_base_uri}/{cube_latest_snapshot}/cube/"
    return cube_uri


def query():
    request = connexion.request.json

    criteria = WmgQueryCriteria(**request["filter"])
    query_result = WmgQuery(find_cube_latest_snapshot()).execute(criteria)
    dot_plot_matrix_df = build_dot_plot_matrix(query_result)
    all_filter_dims_values = extract_filter_dims_values(query_result)
    response_filter_dims_values = build_filter_dims_values(all_filter_dims_values)

    cell_type_term_ids = all_filter_dims_values['cell_type_ontology_term_id']

    return jsonify(
        dict(
            snapshot_id=DUMMY_SNAPSHOT_UUID,
            expression_summary=build_expression_summary(dot_plot_matrix_df),
            term_id_labels=dict(
                genes=build_gene_id_label_mapping(criteria.gene_ontology_term_ids),
                cell_types=build_ontology_term_id_label_mapping(cell_type_term_ids),
            ),
            filter_dims=response_filter_dims_values
        )
    )


def build_datasets(dataset_ids: List[str]) -> List[Dict]:
    return [dict(id=dataset_id,
                 # TODO: retrieve from db or API
                 label='',
                 # TODO: retrieve from db or API
                 collection_label='',
                 # TODO: retrieve from db or API
                 collection_url=''
                 ) for dataset_id in dataset_ids]


def build_filter_dims_values(all_filter_dims_values):
    response_filter_dims_values = dict(
            datasets=build_datasets(all_filter_dims_values["dataset_id"]),
            disease_terms=build_ontology_term_id_label_mapping(all_filter_dims_values["disease_ontology_term_id"]),
            sex_terms=build_ontology_term_id_label_mapping(all_filter_dims_values["sex_ontology_term_id"]),
            development_stage_terms=build_ontology_term_id_label_mapping(
                    all_filter_dims_values["development_stage_ontology_term_id"]),
            ethnicity_terms=build_ontology_term_id_label_mapping(all_filter_dims_values["ethnicity_ontology_term_id"]),
            # excluded per product requirements, but keeping in, commented-out, to reduce future head-scratching
            # assay_ontology_terms=build_ontology_term_id_label_mapping(all_filter_dims_values["assay_ontology_term_id"]),
    )
    return response_filter_dims_values


def build_expression_summary(query_result: DataFrame) -> dict:
    # Create nested dicts with gene_ontology_term_id, tissue_ontology_term_id keys, respectively
    structured_result = defaultdict(lambda: defaultdict(list))
    for group_by_key, cell_type_stats in query_result.to_dict("index").items():
        gene_ontology_term_id, tissue_ontology_term_id, cell_type_ontology_term_id = [s for s in group_by_key]
        structured_result[gene_ontology_term_id][tissue_ontology_term_id].append(
            dict(
                id=cell_type_ontology_term_id,
                n=cell_type_stats["nnz"],
                me=cell_type_stats["sum"] / cell_type_stats["n_cells"],
                # TODO
                pc=0.0,
                # TODO
                tpc=0.0,
            )
        )
    return structured_result


def extract_filter_dims_values(query_result: DataFrame) -> dict:
    """
    Return unique values for each dimension in the specified query result
    """
    return {col: query_result.groupby(col).groups.keys()
            for col in set(query_result.columns) & set(cube_non_indexed_dims)}
