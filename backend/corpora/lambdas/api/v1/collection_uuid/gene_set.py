from flask import make_response, g, jsonify
from sqlalchemy.exc import IntegrityError

from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Collection
from backend.corpora.common.entities.geneset import Geneset
from backend.corpora.common.utils.exceptions import ForbiddenHTTPException, InvalidParametersHTTPException


def post(collection_uuid: str, body: dict, user: str):
    db_session = g.db_session
    collection = Collection.if_owner(db_session, collection_uuid, CollectionVisibility.PRIVATE, user)
    if not collection:
        raise ForbiddenHTTPException
    gene_sets = body["gene_sets"]
    for gene_set in gene_sets:
        try:
            Geneset.create(
                db_session, name=gene_set["gene_set_name"], description=gene_set["gene_set_description"],
                gene_symbols=gene_set["genes"], collection=collection
            )
        except IntegrityError:
            db_session.rollback()
            raise InvalidParametersHTTPException("Duplicate geneset name")

    result = Geneset.retrieve_all_genesets_for_a_collection(db_session, collection_uuid)

    return make_response(jsonify(result), 200)
