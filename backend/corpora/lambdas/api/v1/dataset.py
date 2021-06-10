from flask import make_response, jsonify, g

from ....common.corpora_orm import CollectionVisibility
from ....common.entities import Dataset, Collection
from ....common.entities.geneset import GenesetDatasetLink
from ....api_server.db import dbconnect
from ....common.utils.exceptions import (
    NotFoundHTTPException,
    ServerErrorHTTPException,
    ForbiddenHTTPException,
    CorporaException,
)


@dbconnect
def post_dataset_asset(dataset_uuid: str, asset_uuid: str):
    db_session = g.db_session
    # retrieve the dataset
    dataset = Dataset.get(db_session, dataset_uuid)
    if not dataset:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}' not found.")

    # retrieve the artifact
    asset = dataset.get_asset(asset_uuid)
    if not asset:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}/asset/{asset_uuid}' not found.")

    # Retrieve S3 metadata
    file_size = asset.get_file_size()
    if not file_size:
        raise ServerErrorHTTPException()

    # Generate pre-signed URL
    presigned_url = asset.generate_file_url()
    if not presigned_url:
        raise ServerErrorHTTPException()

    return make_response(
        jsonify(
            dataset_id=dataset_uuid,
            file_name=asset.filename,
            file_size=file_size,
            presigned_url=presigned_url,
        ),
        200,
    )


@dbconnect
def get_status(dataset_uuid: str, user: str):
    db_session = g.db_session
    dataset = Dataset.get(db_session, dataset_uuid)
    if not dataset:
        raise ForbiddenHTTPException()
    if not Collection.get_collection(db_session, dataset.collection.id, dataset.collection.visibility, owner=user):
        raise ForbiddenHTTPException()
    status = dataset.processing_status.to_dict(remove_none=True)
    for remove in ["dataset", "created_at", "updated_at"]:
        status.pop(remove)
    return make_response(jsonify(status), 200)


@dbconnect
def delete_dataset(dataset_uuid: str, user: str):
    """
    Deletes an existing dataset or cancels an in progress upload.
    """
    db_session = g.db_session
    dataset = Dataset.get(db_session, dataset_uuid, include_tombstones=True)
    if not dataset:
        raise ForbiddenHTTPException()
    if not Collection.get_collection(db_session, dataset.collection.id, dataset.collection.visibility, owner=user):
        raise ForbiddenHTTPException()
    if dataset.collection_visibility == CollectionVisibility.PUBLIC:
        return make_response(jsonify("Can not delete a public dataset"), 405)
    if dataset.tombstone is False:
        if dataset.published:
            dataset.update(tombstone=True)
        else:
            dataset.asset_deletion()
            dataset.delete()
    return "", 202


@dbconnect
def post_dataset_gene_sets(dataset_uuid: str, body: object, user: str):
    db_session = g.db_session
    dataset = Dataset.get(db_session, dataset_uuid)
    if not dataset:
        raise ForbiddenHTTPException()
    collection = Collection.get_collection(
        db_session, dataset.collection.id, CollectionVisibility.PRIVATE.name, owner=user
    )
    if not collection:
        raise ForbiddenHTTPException()
    validate_genesets_in_collection_and_linked_to_dataset(dataset, collection, body)
    try:
        GenesetDatasetLink.update_links_for_a_dataset(db_session, dataset_uuid, add=body["add"], remove=body["remove"])
    except CorporaException:
        raise NotFoundHTTPException()
    gene_sets = [
        x.to_dict(
            remove_attr=[
                "collection",
                "collection_visibility",
                "collection_id",
                "created_at",
                "updated_at",
                "genes",
            ]
        )
        for x in dataset.genesets
    ]
    return make_response(jsonify(gene_sets), 202)


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
