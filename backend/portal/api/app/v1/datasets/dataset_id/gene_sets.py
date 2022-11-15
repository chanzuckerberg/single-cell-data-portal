from flask import g, make_response, jsonify

from backend.api_server.db import dbconnect
from backend.common.corpora_orm import CollectionVisibility
from backend.common.entities import Dataset, Collection
from backend.common.entities.geneset import GenesetDatasetLink
from backend.common.utils.exceptions import CorporaException
from backend.common.utils.http_exceptions import ForbiddenHTTPException, NotFoundHTTPException
from backend.portal.api.app.v1.authorization import owner_or_allowed


def validate_genesets_in_collection_and_linked_to_dataset(dataset, collection, update_list):
    dataset_geneset_ids = [x.id for x in dataset.genesets]
    collection_geneset_ids = [x.id for x in collection.genesets]

    add_list_in_collection = all(item in collection_geneset_ids for item in update_list["add"])
    remove_list_in_collection = all(item in collection_geneset_ids for item in update_list["remove"])
    if not (add_list_in_collection and remove_list_in_collection):
        raise NotFoundHTTPException()
    remove_list_in_dataset = all(item in dataset_geneset_ids for item in update_list["remove"])
    if not remove_list_in_dataset:
        raise NotFoundHTTPException()
    add_list_in_dataset = any(item in dataset_geneset_ids for item in update_list["add"])
    if add_list_in_dataset:
        raise NotFoundHTTPException()


@dbconnect
def post(dataset_id: str, body: object, token_info: dict):
    db_session = g.db_session
    dataset = Dataset.get(db_session, dataset_id)
    if not dataset:
        raise ForbiddenHTTPException()
    collection = Collection.get_collection(
        db_session, dataset.collection.id, CollectionVisibility.PRIVATE.name, owner=owner_or_allowed(token_info)
    )
    if not collection:
        raise ForbiddenHTTPException()
    validate_genesets_in_collection_and_linked_to_dataset(dataset, collection, body)
    try:
        GenesetDatasetLink.update_links_for_a_dataset(db_session, dataset_id, add=body["add"], remove=body["remove"])
    except CorporaException:
        raise NotFoundHTTPException()
    gene_sets = [
        x.to_dict(
            remove_attr=[
                "collection",
                "collection_id",
                "created_at",
                "updated_at",
                "genes",
            ]
        )
        for x in dataset.genesets
    ]
    return make_response(jsonify(gene_sets), 202)
