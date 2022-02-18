from typing import List, Dict
from uuid import uuid4

import connexion
from flask import jsonify

from backend.wmg.data.query import build_genes, build_gene_id_label_mapping, build_cell_type_id_label_mapping, \
    cell_type_term_ids

DUMMY_SNAPSHOT_UUID = uuid4().hex


def primary_filter_dimensions():
    organism_terms = [dict(oid1="olbl1"), dict(oid2="olbl2")]
    tissue_terms = [dict(ttid1="ttlbl1"), dict(ttid2="ttlbl2")]
    result = dict(snapshot_id=DUMMY_SNAPSHOT_UUID, organism_terms=organism_terms, tissue_terms=tissue_terms)
    return jsonify(result)


def query():
    request = connexion.request.json

    gene_term_ids = request["filter"]["gene_term_ids"]
    tissue_term_ids = request["filter"]["tissue_term_ids"]

    return jsonify(
        dict(
            snapshot_id=DUMMY_SNAPSHOT_UUID,
            expression_summary=build_genes(gene_term_ids, tissue_term_ids),
            term_id_labels=dict(
                genes=build_gene_id_label_mapping(gene_term_ids),
                cell_types=build_cell_type_id_label_mapping(cell_type_term_ids),
            ),
        )
    )
