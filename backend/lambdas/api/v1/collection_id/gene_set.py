from flask import make_response, g, jsonify
from sqlalchemy.exc import IntegrityError

from backend.common.corpora_orm import CollectionVisibility
from backend.common.entities import Collection
from backend.common.entities.geneset import Geneset
from backend.api_server.db import dbconnect
from backend.common.utils.http_exceptions import ForbiddenHTTPException, InvalidParametersHTTPException
from backend.lambdas.api.v1.authorization import owner_or_allowed


@dbconnect
def post(collection_id: str, body: dict, token_info: dict):
    db_session = g.db_session
    collection = Collection.get_collection(
        db_session,
        collection_id,
        CollectionVisibility.PRIVATE.name,
        owner=owner_or_allowed(token_info),
    )
    if not collection:
        raise ForbiddenHTTPException
    gene_sets = body["gene_sets"]
    for gene_set in gene_sets:
        try:
            Geneset.create(
                db_session,
                name=gene_set["gene_set_name"],
                description=gene_set["gene_set_description"],
                genes=gene_set["genes"],
                collection_id=collection.id,
            )
        except IntegrityError:
            db_session.rollback()
            raise InvalidParametersHTTPException(detail="Duplicate geneset name")

    result = Geneset.retrieve_all_genesets_for_a_collection(db_session, collection_id)

    return make_response(jsonify(result), 200)
