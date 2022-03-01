import json
from collections import defaultdict
from tempfile import mkdtemp
from uuid import uuid4

import connexion
import tiledb
from flask import jsonify
from pandas import DataFrame

from backend.wmg.data.config import fast_config, create_ctx
from backend.wmg.data.query import (
    build_gene_id_label_mapping,
    build_cell_type_id_label_mapping,
    WmgQuery,
    WmgQueryCriteria
)

# TODO: Replace with real snapshot uuid
from unit.backend.wmg.fixtures.cube import create_cube

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
        # TODO: Replace with S3 location
        cube_tmp_dir = mkdtemp()
        create_cube(cube_dir=cube_tmp_dir, dim_size=3)

        # TODO: Remove tiledb dependency from this module
        # TODO: Okay to keep open indefinitely? Is it faster than re-opening each request?
        cube = tiledb.open(cube_tmp_dir, ctx=create_ctx(fast_config()))

    return cube


def query():
    request = connexion.request.json

    criteria = WmgQueryCriteria(**request["filter"])
    query_result = WmgQuery(find_cube_latest_snapshot()).execute(criteria)

    cell_type_term_ids = {key[2] for key in query_result.index}

    return jsonify(
        dict(
            snapshot_id=DUMMY_SNAPSHOT_UUID,
            expression_summary=build_expression_summary_json(query_result),
            term_id_labels=dict(
                genes=build_gene_id_label_mapping(criteria.gene_ontology_term_ids),
                cell_types=build_cell_type_id_label_mapping(cell_type_term_ids),
            ),
        )
    )


def build_expression_summary_json(query_result: DataFrame) -> str:
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
