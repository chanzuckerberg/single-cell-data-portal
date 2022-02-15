from uuid import UUID, uuid4

import connexion
from flask import jsonify

_logger = None


def latest_snapshot():
    return jsonify(snapshot_id=uuid4().hex)


def primary_filter_dimensions(snapshot_id: UUID):
    organism_terms = [dict(oid1='olbl1'),
                      dict(oid2='olbl2')]
    tissue_type_terms = [dict(ttid1='ttlbl1'),
                         dict(ttid2='ttlbl2')]
    result = dict(organism_terms=organism_terms,
                  tissue_type_terms=tissue_type_terms)
    return jsonify(result)


def query(snapshot_id: UUID):
    request = connexion.request.json
    print(snapshot_id)
    print(request['filter'])
    print(request['response_option'])
    return jsonify(dict(expression_summary={}, term_id_labels=dict(genes=[], cell_types=[])))

