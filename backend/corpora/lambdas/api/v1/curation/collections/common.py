import typing
from sqlalchemy.orm import Session

from ...authorization import is_user_owner_or_allowed
from ......common.corpora_config import CorporaConfig
from ......common.corpora_orm import (
    CollectionVisibility,
    DbCollectionLink,
    DbCollection,
    DbDataset,
    DbDatasetProcessingStatus,
    DbDatasetArtifact,
    DatasetArtifactFileType,
    Base,
    ProcessingStatus,
)


DATASET_ONTOLOGY_ELEMENTS = (
    "sex",
    "ethnicity",
    "development_stage",
    "cell_type",
    "tissue",
    "assay",
    "disease",
    "organism",
)

DATASET_ONTOLOGY_ELEMENTS_PREVIEW = (
    "tissue",
    "assay",
    "disease",
    "organism",
)


def reshape_for_curation_api_and_is_allowed(collection, token_info, id_provided=False, preview=False):
    """
    Reshape Collection data for the Curation API response. Remove tombstoned Datasets.
    :param collection: the Collection being returned in the API response
    :param token_info: user access token
    :param id_provided: bool - whether or not the collection uuid was provided by the user, for access purposes
    :param preview: boool - whether the dataset is in preview form or not.
    :return: whether or not the Collection should be included in the response per ownership/access rules
    """

    owner = collection.pop("owner")  # Don't actually want to return 'owner' in response
    if is_user_owner_or_allowed(token_info, owner):
        collection["access_type"] = "WRITE"
    elif not id_provided and collection["visibility"] == CollectionVisibility.PRIVATE:
        # User neither provided the uuid for access nor are they authorized by their access token
        return False
    elif token_info:
        # Access token was provided but user is not authorized
        collection["access_type"] = "READ"
    collection["collection_url"] = f"{CorporaConfig().collections_base_url}/collections/{collection['id']}"

    if datasets := collection.get("datasets"):
        collection["datasets"] = reshape_datasets_for_curation_api(datasets, preview)
    return True


def reshape_datasets_for_curation_api(datasets: typing.List[dict], preview=False) -> typing.List[dict]:
    active_datasets = []
    for dataset in datasets:
        if dataset.get("tombstone"):
            continue  # Remove tombstoned Datasets
        active_datasets.append(reshape_dataset_for_curation_api(dataset, preview))
    return active_datasets


def reshape_dataset_for_curation_api(dataset: dict, preview=False) -> dict:
    if artifacts := dataset.get("artifacts"):
        dataset["dataset_assets"] = []
        for asset in artifacts:
            if asset["filetype"] in (DatasetArtifactFileType.H5AD, DatasetArtifactFileType.RDS):
                dataset["dataset_assets"].append(asset)
        del dataset["artifacts"]
    if dataset.get("processing_status"):
        dataset["processing_status"] = dataset["processing_status"]["processing_status"]
    dataset_ontology_elements = DATASET_ONTOLOGY_ELEMENTS_PREVIEW if preview else DATASET_ONTOLOGY_ELEMENTS
    for ontology_element in dataset_ontology_elements:
        if dataset_ontology_element := dataset.get(ontology_element):
            if not isinstance(dataset_ontology_element, list):
                # Package in array
                dataset[ontology_element] = [dataset_ontology_element]
        else:
            dataset[ontology_element] = []
    return dataset


def list_collections_curation(
    session: Session, collection_columns: typing.Dict[Base, typing.List[str]], visibility: str = None
) -> typing.List[dict]:
    """
    Get a subset of columns, in dict form, for all Collections with the specified visibility. If visibility is None,
    return *all* Collections that are *not* tombstoned.
    :param session: the SQLAlchemy session
    :param collection_columns: the list of columns to be returned (see usage by TransformingBase::to_dict_keep)
    :param visibility: the CollectionVisibility string name
    @return: a list of dict representations of Collections
    """
    filters = [DbCollection.tombstone == False]  # noqa
    if visibility == CollectionVisibility.PUBLIC.name:
        filters.append(DbCollection.visibility == CollectionVisibility.PUBLIC)
    elif visibility == CollectionVisibility.PRIVATE.name:
        filters.append(DbCollection.visibility == CollectionVisibility.PRIVATE)

    resp_collections = []
    for collection in session.query(DbCollection).filter(*filters).all():
        resp_collection = collection.to_dict_keep(collection_columns)
        resp_collection["processing_status"] = add_collection_level_processing_status(collection)
        resp_collections.append(resp_collection)
    return resp_collections


def add_collection_level_processing_status(collection: DbCollection):
    # Add a Collection-level processing status
    status = None
    has_statuses = False
    for dataset in collection.datasets:
        processing_status = dataset.processing_status
        if processing_status:
            has_statuses = True
            if processing_status.processing_status == ProcessingStatus.PENDING:
                status = ProcessingStatus.PENDING
            elif processing_status.processing_status == ProcessingStatus.FAILURE:
                status = ProcessingStatus.FAILURE
                break
    if has_statuses and not status:  # At least one dataset processing status exists, and all were SUCCESS
        status = ProcessingStatus.SUCCESS
    return status


class EntityColumns:

    collections_cols = [
        "id",
        "name",
        "visibility",
        "tombstone",
        "contact_name",
        "contact_email",
        "curator_name",
        "revised_at",
        "created_at",
        "published_at",
        "description",
        "publisher_metadata",
        "revision_of",
        "tombstone",
        "owner",  # Needed for determining view permissions
        "links",
        "datasets",
    ]

    link_cols = [
        "link_name",
        "link_url",
        "link_type",
    ]

    dataset_preview_cols = [
        "id",
        "curator_tag",
        "tissue",
        "assay",
        "disease",
        "organism",
        "tombstone",
    ]

    dataset_cols = [
        *dataset_preview_cols,
        "name",
        "revision",
        "revised_at",
        "is_primary_data",
        "x_normalization",
        "artifacts",
        "sex",
        "ethnicity",
        "development_stage",
        "explorer_url",
        "cell_type",
        "cell_count",
        "x_approximate_distribution",
        # "batch_condition",  # TODO: https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data-portal/1461  # noqa: E501
        "mean_genes_per_cell",
        "schema_version",
        "processing_status",
    ]

    dataset_asset_cols = [
        "filetype",
        "filename",
    ]

    dataset_processing_status_cols = [
        "processing_status",
    ]

    columns_for_collections = {
        DbCollectionLink: link_cols,
        DbCollection: collections_cols,
        DbDataset: dataset_preview_cols,
        DbDatasetProcessingStatus: dataset_processing_status_cols,
    }

    columns_for_collection_id = {
        DbCollectionLink: link_cols,
        DbCollection: collections_cols,
        DbDataset: dataset_cols,
        DbDatasetArtifact: dataset_asset_cols,
        DbDatasetProcessingStatus: dataset_processing_status_cols,
    }

    columns_for_dataset = {
        DbDataset: dataset_cols,
        DbDatasetArtifact: dataset_asset_cols,
        DbDatasetProcessingStatus: dataset_processing_status_cols,
    }
