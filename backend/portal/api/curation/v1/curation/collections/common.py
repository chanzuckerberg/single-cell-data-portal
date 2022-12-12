from dataclasses import asdict
from typing import List, Optional, Tuple
from urllib.parse import urlparse

from backend.common.corpora_config import CorporaConfig
from backend.common.corpora_orm import (
    DatasetArtifactFileType,
    DbCollection,
    DbCollectionLink,
    DbDataset,
    DbDatasetArtifact,
    DbDatasetProcessingStatus,
    ProcessingStatus,
    ValidationStatus,
)
from backend.common.utils.http_exceptions import ForbiddenHTTPException
from backend.layers.api.router import get_business_logic
from backend.layers.auth.user_info import UserInfo
from backend.layers.common.entities import (
    CollectionId,
    CollectionVersion,
    CollectionVersionId,
    DatasetProcessingStatus,
    DatasetVersion,
    Link,
)

allowed_dataset_asset_types = ("H5AD", "RDS")

DATASET_ONTOLOGY_ELEMENTS = (
    "sex",
    "self_reported_ethnicity",
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


def extract_doi_from_links(links: List[Link]) -> Tuple[Optional[str], List[dict]]:
    """
    Pull out the DOI from the 'links' list and return it along with the altered links array
    :param links: a list of links.
    :return: the DOI link and the remaining links
    """
    doi, dict_links = None, []
    for link in links:
        if link["type"] == "DOI":
            doi = urlparse(link["uri"]).path.strip("/")
        else:
            dict_links.append(dict(link_name=link.name, link_url=link.uri, link_type=link.type))
    return doi, links


def reshape_for_curation_api(collection_version: CollectionVersion, user_info: UserInfo, preview: bool = False) -> dict:
    """
    Reshape Collection data for the Curation API response. Remove tombstoned Datasets.
    :param collection: the Collection being returned in the API response
    :param user_info:
    :param preview: boool - whether the dataset is in preview form or not.
    :return: the response.
    """
    # get collectoin attributes based on published status
    if collection_version.published_at is None:
        # Unpublished
        collection_id = collection_version.version_id
        revision_of = collection_version.collection_id
        revising_in = None
    else:
        # Published
        collection_id = collection_version.collection_id
        revision_of = None
        revising_in = (
            None
            if not user_info.is_user_owner_or_allowed(collection_version.owner)
            else get_business_logic().get_unplublished_collection_version_from_canonical(collection_id)
        )

    # get collection dataset attributes
    response_datasets = []
    collection_level_processing_status = None  # Noe if no datasets
    for dataset_version_id in collection_version.datasets:
        dataset_version = get_business_logic().get_dataset_version(dataset_version_id)
        ds = asdict(dataset_version.metadata)

        # get dataset asset attributes
        assets = []
        for artifact in dataset_version.artifacts:
            if artifact.type in allowed_dataset_asset_types:
                assets.append(dict(filetype=artifact.type, filename=artifact.uri.split("/")[-1]))
        ds["dataset_assets"] = assets
        ds["processing_status_detail"] = dataset_version.status.validation_message
        ds["revised_at"] = dataset_version.canonical_dataset.revised_at
        response_datasets.append(ds)

        # get the Collection-level processing status
        dataset_processing_status = dataset_version.status.processing_status
        if dataset_processing_status:
            if dataset_processing_status in (DatasetProcessingStatus.PENDING, DatasetProcessingStatus.INITIALIZED):
                collection_level_processing_status = DatasetProcessingStatus.PENDING
            elif dataset_processing_status == DatasetProcessingStatus.FAILURE:
                collection_level_processing_status = dataset_processing_status

    # build response
    doi, links = extract_doi_from_links(collection_version.metadata.links)
    response = dict(
        collection_url=f"{CorporaConfig().collections_base_url}/collections/{collection_id.id}",
        contact_email=collection_version.metadata.contact_email,
        contact_name=collection_version.metadata.contact_name,
        created_at=collection_version.created_at,
        curator_name=collection_version.owner,
        datasets=response_datasets,
        description=collection_version.metadata.description,
        dio=doi,
        id=collection_id.id,
        links=links,
        name=collection_version.metadata.name,
        processing_status=collection_level_processing_status,
        published_at=collection_version.canonical_collection.originally_published_at,
        publisher_metadata=collection_version.publisher_metadata,
        revised_at=get_business_logic()
        .get_published_collection_version(collection_version.canonical_collection.id)
        .published_at,
        revising_in=revising_in,
        revision_of=revision_of,
        visibility=get_visibility(collection_version),
    )
    return response


def reshape_datasets_for_curation_api(datasets: List[dict], preview=False) -> List[dict]:
    active_datasets = []
    for dataset in datasets:
        if dataset.get("tombstone"):
            continue  # Remove tombstoned Datasets
        active_datasets.append(reshape_dataset_for_curation_api(dataset, preview))
    return active_datasets


def reshape_dataset_for_curation_api(dataset: dict, preview=False) -> dict:
    if artifacts := dataset.pop("artifacts", []):
        dataset["dataset_assets"] = []
        for asset in artifacts:
            if asset["filetype"] in (DatasetArtifactFileType.H5AD, DatasetArtifactFileType.RDS):
                dataset["dataset_assets"].append(asset)
    if processing_status := dataset.pop("processing_status", None):
        if processing_status["processing_status"] == ProcessingStatus.FAILURE:
            if processing_status["validation_status"] == ValidationStatus.INVALID:
                dataset["processing_status_detail"] = processing_status["validation_message"]
                dataset["processing_status"] = "VALIDATION_FAILURE"
            else:
                dataset["processing_status"] = "PIPELINE_FAILURE"
        else:
            dataset["processing_status"] = processing_status["processing_status"]
    dataset_ontology_elements = DATASET_ONTOLOGY_ELEMENTS_PREVIEW if preview else DATASET_ONTOLOGY_ELEMENTS
    for ontology_element in dataset_ontology_elements:
        if dataset_ontology_element := dataset.get(ontology_element):
            if not isinstance(dataset_ontology_element, list):
                # Package in array
                dataset[ontology_element] = [dataset_ontology_element]
        else:
            dataset[ontology_element] = []

    if not preview:  # Add these fields only to full (and not preview) Dataset metadata response
        dataset["revision_of"] = dataset.pop("original_id", None)
        dataset["title"] = dataset.pop("name", None)
        if value := dataset.pop("is_primary_data", None):
            dataset["is_primary_data"] = is_primary_data_mapping.get(value, [])

    return dataset


is_primary_data_mapping = {
    "PRIMARY": [True],
    "SECONDARY": [False],
    "BOTH": [True, False],
}


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
        "revising_in",
    ]

    link_cols = [
        "link_name",
        "link_url",
        "link_type",
    ]

    dataset_preview_cols = [
        "id",
        "tissue",
        "assay",
        "disease",
        "organism",
        "tombstone",
        "suspension_type",
    ]

    dataset_cols = [
        *dataset_preview_cols,
        "original_id",
        "name",
        "revision",
        "revised_at",
        "is_primary_data",
        "artifacts",
        "sex",
        "self_reported_ethnicity",
        "development_stage",
        "explorer_url",
        "cell_type",
        "cell_count",
        "x_approximate_distribution",
        "batch_condition",
        "mean_genes_per_cell",
        "schema_version",
        "processing_status",
        "donor_id",
    ]

    dataset_asset_cols = [
        "filetype",
        "filename",
    ]

    dataset_processing_status_cols = ["processing_status", "validation_message", "validation_status"]

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


def get_visibility(collection_version: CollectionVersion) -> str:
    return "PUBLIC" if collection_version.published_at else "PRIVATE"


def get_collection_level_processing_status(datasets: List[DatasetVersion]) -> str:
    # get the Collection-level processing status
    if not datasets:  # Return None if no datasets.
        return None
    return_status = DatasetProcessingStatus.SUCCESS
    for dataset in datasets:
        status = dataset.status.processing_status
        if status:
            if status in (DatasetProcessingStatus.PENDING, DatasetProcessingStatus.INITIALIZED):
                return_status = DatasetProcessingStatus.PENDING
            elif status == DatasetProcessingStatus.FAILURE:
                return status
    return return_status


def get_infered_collection_version_else_forbidden(collection_id: str) -> Optional[CollectionVersion]:
    """
    Infer the collection version from either a CollectionId or a CollectionVersionId and return the CollectionVersion.
    :param collection_id: identifies the collection version
    :return: The CollectionVersion if it exists.
    """
    version = get_business_logic().get_published_collection_version(CollectionId(collection_id))
    if version is None:
        version = get_business_logic().get_collection_version(CollectionVersionId(collection_id))
    if version is None or version.canonical_collection.tombstoned is True:
        raise ForbiddenHTTPException()
    return version


def is_owner_or_allowed_else_forbidden(collection_version, user_info):
    if not user_info.is_user_owner_or_allowed(collection_version.owner):
        raise ForbiddenHTTPException()
