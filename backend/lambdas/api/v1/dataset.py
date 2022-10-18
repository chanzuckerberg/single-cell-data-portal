from flask import make_response, jsonify, g

from backend.lambdas.api.v1.common import delete_dataset_common, get_collection_else_forbidden
from backend.common.corpora_orm import CollectionVisibility, DatasetArtifactFileType
from backend.common.entities import Dataset, Collection
from backend.common.entities.geneset import GenesetDatasetLink
from backend.api_server.db import dbconnect
from backend.common.utils.http_exceptions import (
    NotFoundHTTPException,
    ServerErrorHTTPException,
    ForbiddenHTTPException,
)
from backend.common.utils.exceptions import CorporaException
from backend.lambdas.api.v1.authorization import owner_or_allowed


@dbconnect
def post_dataset_asset(dataset_id: str, asset_id: str):
    db_session = g.db_session
    # retrieve the dataset
    dataset = Dataset.get(db_session, dataset_id)
    if not dataset:
        raise NotFoundHTTPException(detail=f"'dataset/{dataset_id}' not found.")

    # retrieve the artifact
    asset = dataset.get_asset(asset_id)
    if not asset:
        raise NotFoundHTTPException(detail=f"'dataset/{dataset_id}/asset/{asset_id}' not found.")

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
            dataset_id=dataset_id,
            file_name=asset.filename,
            file_size=file_size,
            presigned_url=presigned_url,
        ),
        200,
    )


@dbconnect
def get_dataset_assets(dataset_id: str):
    db_session = g.db_session
    # retrieve the dataset
    dataset = Dataset.get(db_session, dataset_id)
    assets = dataset.get_assets()
    return make_response(jsonify(assets=assets))


@dbconnect
def get_status(dataset_id: str, token_info: dict):
    db_session = g.db_session
    dataset = Dataset.get(db_session, dataset_id)
    if not dataset:
        raise ForbiddenHTTPException()
    get_collection_else_forbidden(
        db_session,
        dataset.collection.id,
        owner=owner_or_allowed(token_info),
    )
    status = dataset.processing_status.to_dict(remove_none=True, remove_attr=["dataset", "created_at", "updated_at"])
    return make_response(jsonify(status), 200)


@dbconnect
def get_datasets_index():
    db_session = g.db_session
    datasets = Collection.list_public_datasets_for_index(db_session)
    return make_response(jsonify(datasets), 200)


@dbconnect
def delete_dataset(dataset_id: str, token_info: dict):
    """
    Deletes an existing dataset or cancels an in progress upload.
    """
    db_session = g.db_session
    dataset = Dataset.get(db_session, dataset_id, include_tombstones=True)
    delete_dataset_common(db_session, dataset, token_info)
    return "", 202


@dbconnect
def get_dataset_identifiers(url: str):
    db_session = g.db_session
    dataset = Dataset.get_by_explorer_url(db_session, url)
    if not dataset:
        raise NotFoundHTTPException()
    artifact = dataset.get_most_recent_artifact(filetype=DatasetArtifactFileType.CXG)
    s3_uri = artifact.s3_uri if artifact else None

    dataset_identifiers = {
        "s3_uri": s3_uri,
        "dataset_id": dataset.id,
        "collection_id": dataset.collection_id,
        "collection_visibility": dataset.collection.visibility,
        "tombstoned": dataset.tombstone,
    }
    return make_response(jsonify(dataset_identifiers), 200)


@dbconnect
def post_dataset_gene_sets(dataset_id: str, body: object, token_info: dict):
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
