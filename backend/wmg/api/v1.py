from uuid import UUID, uuid4

from flask import jsonify

_logger = None


def latest_snapshot():
    return jsonify(snapshot_id=uuid4().hex)


def primary_filter_dimensions(snapshot_id: UUID):
    organism_terms = [dict(id='oid1', label='olbl1'),
                      dict(id='oid2', label='olbl2')]
    tissue_type_terms = [dict(id='ttid1', label='ttlbl1'),
                         dict(id='ttid2', label='ttlbl2')]
    result = dict(organism_terms=organism_terms,
                  tissue_type_terms=tissue_type_terms)
    return jsonify(result)
