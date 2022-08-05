import typing

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


def reshape_for_curation_api_and_is_allowed(collection, token_info, id_provided=False):
    """
    Reshape Collection data for the Curation API response. Remove tombstoned Datasets.
    :param collection: the Collection being returned in the API response
    :param token_info: user access token
    :param id_provided: bool - whether or not the collection uuid was provided by the user, for access purposes
    :return: whether or not the Collection should be included in the response per ownership/access rules
    """
    owner = collection["owner"]
    if is_user_owner_or_allowed(token_info, owner):
        collection["access_type"] = "WRITE"
    elif not id_provided and collection["visibility"] == CollectionVisibility.PRIVATE:
        # User neither provided the uuid for access nor are they authorized by their access token
        return False
    elif token_info:
        # Access token was provided but user is not authorized
        collection["access_type"] = "READ"

    del collection["owner"]  # Don't actually want to return 'owner' in response
    collection["collection_url"] = f"{CorporaConfig().collections_base_url}/collections/{collection['id']}"

    if datasets := collection.get("datasets"):
        collection["datasets"] = reshape_datasets_for_curation_api(datasets)
    return True


def reshape_datasets_for_curation_api(datasets: typing.List[dict]) -> typing.List[dict]:
    active_datasets = []
    for dataset in datasets:
        if dataset.get("tombstone"):
            continue  # Remove tombstoned Datasets
        active_datasets.append(reshape_dataset_for_curation_api(dataset))
    return active_datasets


def reshape_dataset_for_curation_api(dataset: dict) -> dict:
    if artifacts := dataset.get("artifacts"):
        dataset["dataset_assets"] = []
        for asset in artifacts:
            if asset["filetype"] in (DatasetArtifactFileType.H5AD, DatasetArtifactFileType.RDS):
                dataset["dataset_assets"].append(asset)
        del dataset["artifacts"]
    if dataset.get("processing_status"):
        dataset["processing_status"] = dataset["processing_status"]["processing_status"]
    for ontology_element in DATASET_ONTOLOGY_ELEMENTS:
        if dataset_ontology_element := dataset.get(ontology_element):
            if not isinstance(dataset_ontology_element, list):
                # Package in array
                dataset[ontology_element] = [dataset_ontology_element]
        else:
            dataset[ontology_element] = []
    return dataset


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
