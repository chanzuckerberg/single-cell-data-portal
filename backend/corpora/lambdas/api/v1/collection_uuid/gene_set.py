from flask import make_response, g, jsonify
from sqlalchemy.exc import IntegrityError

from .....common.corpora_orm import CollectionVisibility
from .....common.entities import Collection
from .....common.entities.geneset import Geneset
from .....api_server.db import dbconnect
from .....common.utils.exceptions import ForbiddenHTTPException, InvalidParametersHTTPException
from backend.corpora.lambdas.api.v1.collection import _owner_or_allowed


@dbconnect
def post(collection_uuid: str, body: dict, user: str):
    db_session = g.db_session
    collection = Collection.get_collection(db_session, collection_uuid, CollectionVisibility.PRIVATE.name, owner=_owner_or_allowed(user))
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
                collection_visibility=collection.visibility.name,
            )
        except IntegrityError:
            db_session.rollback()
            raise InvalidParametersHTTPException("Duplicate geneset name")

    result = Geneset.retrieve_all_genesets_for_a_collection(db_session, collection_uuid)

    return make_response(jsonify(result), 200)
