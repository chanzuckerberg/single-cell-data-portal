from flask import make_response, g

from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Collection
from backend.corpora.common.entities.geneset import Geneset
from backend.corpora.common.utils.exceptions import ForbiddenHTTPException


def post(collection_uuid: str, body: dict, user: str):
    db_session = g.db_sesson
    collection = Collection.if_owner(db_session, collection_uuid, CollectionVisibility.PRIVATE, user)
    if not collection:
        raise ForbiddenHTTPException
    gene_sets = body["gene_sets"]
    for gene_set in gene_sets:
        Geneset.create(db_session, name=gene_set[0], description=gene_set[1], gene_symbols=gene_set[2],
                       collection=collection)

    gene_sets = collection.genesets

