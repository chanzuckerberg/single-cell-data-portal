from dataclasses import asdict
from typing import List, Optional, Tuple, Union
from urllib.parse import urlparse

from backend.common.corpora_config import CorporaConfig
from backend.common.corpora_orm import (
    DbCollection,
    DbCollectionLink,
    DbDataset,
    DbDatasetArtifact,
    DbDatasetProcessingStatus,
)
from backend.common.utils.http_exceptions import ForbiddenHTTPException
from backend.layers.api.explorer_url import generate as generate_explorer_url
from backend.layers.api.providers import get_business_logic
from backend.layers.auth.user_info import UserInfo
from backend.layers.common.entities import (
    CollectionId,
    CollectionVersion,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    DatasetArtifactType,
    DatasetId,
    DatasetProcessingStatus,
    DatasetValidationStatus,
    DatasetVersion,
    DatasetVersionId,
    Link,
    OntologyTermId,
)

allowed_dataset_asset_types = (DatasetArtifactType.H5AD, DatasetArtifactType.RDS)


def get_collections_base_url():
    return CorporaConfig().collections_base_url


def extract_doi_from_links(links: List[Link]) -> Tuple[Optional[str], List[dict]]:
    """
    Pull out the DOI from the 'links' list and return it along with the altered links array
    :param links: a list of links.
    :return: the DOI link and the remaining links
    """
    doi, dict_links = None, []
    for link in links:
        if link.type == "DOI":
            doi = urlparse(link.uri).path.strip("/")
        else:
            dict_links.append(dict(link_name=link.name, link_url=link.uri, link_type=link.type))
    return doi, dict_links


def reshape_for_curation_api(
    collection_version: Union[CollectionVersion, CollectionVersionWithDatasets],
    user_info: UserInfo,
    preview: bool = False,
) -> dict:
    """
    Reshape Collection data for the Curation API response. Remove tombstoned Datasets.
    :param collection: the Collection being returned in the API response
    :param user_info:
    :param preview: boool - whether the dataset is in preview form or not.
    :return: the response.
    """
    business_logic = get_business_logic()
    is_published = collection_version.published_at is not None
    # get collection attributes based on published status
    if is_published:
        # Published
        collection_id = collection_version.collection_id
        collection_url = f"{get_collections_base_url()}/collections/{collection_id.id}"
        revision_of = None
        if not user_info.is_user_owner_or_allowed(collection_version.owner):
            _revising_in = None
        else:
            _revising_in = business_logic.get_unpublished_collection_version_from_canonical(
                collection_version.collection_id
            )
        revising_in = _revising_in.version_id.id if _revising_in else None
        use_canonical_id = True
    else:
        # Unpublished - need to determine if it's a revision or first time collection
        # For that, we look at whether the canonical collection is published
        is_revision = collection_version.canonical_collection.originally_published_at is not None
        if is_revision:
            # If it's a revision, both collection_id and collection_url need to point to the version_id,
            # and datasets should expose the private url (based on version_id)
            collection_id = collection_version.version_id
            collection_url = f"{get_collections_base_url()}/collections/{collection_id.id}"
            revision_of = collection_version.collection_id.id
            use_canonical_id = False
        else:
            # If it's an unpublished, unrevised collection, then collection_url will point to the permalink
            # (aka the link to the canonical_id) and the collection_id will point to version_id.
            # Also, revision_of should be None, and the datasets should expose the canonical url
            collection_id = collection_version.collection_id
            collection_url = f"{get_collections_base_url()}/collections/{collection_version.collection_id}"
            revision_of = None
            use_canonical_id = True
        revising_in = None

    # get collection dataset attributes
    response_datasets = reshape_datasets_for_curation_api(
        collection_version.datasets, is_published, use_canonical_id, preview
    )

    # build response
    doi, links = extract_doi_from_links(collection_version.metadata.links)
    published_version = business_logic.get_published_collection_version(collection_version.canonical_collection.id)
    if published_version is not None:
        revised_at = published_version.published_at
    else:
        revised_at = None
    response = dict(
        collection_url=collection_url,
        consortia=collection_version.metadata.consortia,
        contact_email=collection_version.metadata.contact_email,
        contact_name=collection_version.metadata.contact_name,
        created_at=collection_version.created_at,
        curator_name=collection_version.curator_name,
        datasets=response_datasets,
        description=collection_version.metadata.description,
        doi=doi,
        id=collection_id.id,
        links=links,
        name=collection_version.metadata.name,
        processing_status=get_collection_level_processing_status(collection_version.datasets),
        published_at=collection_version.canonical_collection.originally_published_at,
        publisher_metadata=collection_version.publisher_metadata,
        revised_at=revised_at,
        revising_in=revising_in,
        revision_of=revision_of,
        visibility=get_visibility(collection_version),
    )
    return response


def reshape_datasets_for_curation_api(
    datasets: List[Union[DatasetVersionId, DatasetVersion]],
    is_published: bool,
    use_canonical_id: bool,
    preview: bool = False,
) -> List[dict]:
    active_datasets = []
    for dv in datasets:
        dataset_version = get_business_logic().get_dataset_version(dv) if isinstance(dv, DatasetVersionId) else dv
        active_datasets.append(
            reshape_dataset_for_curation_api(dataset_version, is_published, use_canonical_id, preview)
        )
    return active_datasets


def reshape_dataset_for_curation_api(
    dataset_version: DatasetVersion, is_published: bool, use_canonical_id: bool, preview=False
) -> dict:
    ds = dict()

    # Determine what columns to include from the dataset
    if preview:
        columns = EntityColumns.dataset_metadata_preview_cols
    else:
        columns = EntityColumns.dataset_metadata_cols

    # Get dataset metadata fields.
    # Metadata can be None if the dataset isn't still fully processed, so we account for that
    if dataset_version.metadata is not None:
        for column in columns:
            col = getattr(dataset_version.metadata, column)
            if isinstance(col, OntologyTermId):
                col = [asdict(col)]
            elif isinstance(col, list) and len(col) != 0 and isinstance(col[0], OntologyTermId):
                col = [asdict(i) for i in col]
            ds[column] = col

    # Get none preview specific dataset fields
    if not preview:
        # get dataset asset attributes
        assets = []
        for artifact in dataset_version.artifacts:
            if artifact.type in allowed_dataset_asset_types:
                assets.append(dict(filetype=artifact.type.upper(), filename=artifact.uri.split("/")[-1]))

        ds["dataset_assets"] = assets
        ds["processing_status_detail"] = dataset_version.status.validation_message
        ds["revised_at"] = dataset_version.canonical_dataset.revised_at
        ds["revision_of"] = (
            None
            if dataset_version.canonical_dataset.dataset_version_id == dataset_version.version_id
            else dataset_version.dataset_id.id
        )
        ds["revision"] = 0  # TODO this should be the number of times this dataset has been revised and published
        ds["title"] = ds.pop("name", None)
        ds["explorer_url"] = generate_explorer_url(dataset_version, use_canonical_id)
        ds["tombstone"] = False  # TODO this will always be false. Remove in the future
        if dataset_version.metadata is not None:
            ds["is_primary_data"] = is_primary_data_mapping.get(ds.pop("is_primary_data"), [])
            if ds["x_approximate_distribution"]:
                ds["x_approximate_distribution"] = ds["x_approximate_distribution"].upper()
        if status := dataset_version.status:
            if status.processing_status == DatasetProcessingStatus.FAILURE:
                if status.validation_status == DatasetValidationStatus.INVALID:
                    ds["processing_status_detail"] = status.validation_message
                    ds["processing_status"] = "VALIDATION_FAILURE"
                else:
                    ds["processing_status"] = "PIPELINE_FAILURE"
            else:
                ds["processing_status"] = status.processing_status
    if use_canonical_id:
        ds["id"] = dataset_version.dataset_id.id
    else:
        ds["id"] = dataset_version.version_id.id
    return ds


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
        "owner",  # Needed for determining view permissions
        "links",
        "datasets",
        "revising_in",
        "consortia",
    ]

    link_cols = [
        "link_name",
        "link_url",
        "link_type",
    ]

    dataset_metadata_preview_cols = [
        "tissue",
        "assay",
        "disease",
        "organism",
        "suspension_type",
    ]

    dataset_metadata_cols = [
        *dataset_metadata_preview_cols,
        "name",
        "is_primary_data",
        "sex",
        "self_reported_ethnicity",
        "development_stage",
        "cell_type",
        "cell_count",
        "x_approximate_distribution",
        "batch_condition",
        "mean_genes_per_cell",
        "schema_version",
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
        DbDataset: dataset_metadata_preview_cols,
        DbDatasetProcessingStatus: dataset_processing_status_cols,
    }

    columns_for_collection_id = {
        DbCollectionLink: link_cols,
        DbCollection: collections_cols,
        DbDataset: dataset_metadata_cols,
        DbDatasetArtifact: dataset_asset_cols,
        DbDatasetProcessingStatus: dataset_processing_status_cols,
    }

    columns_for_dataset = {
        DbDataset: dataset_metadata_cols,
        DbDatasetArtifact: dataset_asset_cols,
        DbDatasetProcessingStatus: dataset_processing_status_cols,
    }


def get_visibility(collection_version: CollectionVersion) -> str:
    return "PUBLIC" if collection_version.published_at else "PRIVATE"


def get_collection_level_processing_status(datasets: List[DatasetVersion]) -> str:
    if not datasets:  # Return None if no datasets.
        return None

    return_status = DatasetProcessingStatus.SUCCESS
    for dv in datasets:
        version = get_business_logic().get_dataset_version(dv) if isinstance(dv, DatasetVersionId) else dv
        status = version.status.processing_status
        if status:
            if status in (DatasetProcessingStatus.PENDING, DatasetProcessingStatus.INITIALIZED):
                return_status = DatasetProcessingStatus.PENDING
            elif status == DatasetProcessingStatus.FAILURE:
                return_status = status
    return return_status


def get_infered_collection_version_else_forbidden(collection_id: str) -> CollectionVersionWithDatasets:
    """
    Infer the collection version from either a CollectionId or a CollectionVersionId and return the CollectionVersion.
    :param collection_id: identifies the collection version
    :return: The CollectionVersion if it exists.
    """
    version = get_business_logic().get_published_collection_version(CollectionId(collection_id))
    if version is None:
        version = get_business_logic().get_collection_version(CollectionVersionId(collection_id))
    if version is None:
        version = get_business_logic().get_unpublished_collection_version_from_canonical(CollectionId(collection_id))
    if version is None or version.canonical_collection.tombstoned is True:
        raise ForbiddenHTTPException()
    return version


def get_infered_dataset_version(dataset_id: str) -> Optional[DatasetVersion]:
    """
    Infer the dataset version from either a DatasetId or a DatasetVersionId and return the DatasetVersion.
    :param dataset_id: identifies the dataset version
    :return: The DatasetVersion if it exists.
    """
    version = get_business_logic().get_dataset_version(DatasetVersionId(dataset_id))
    if version is None:
        version = get_business_logic().get_dataset_version_from_canonical(DatasetId(dataset_id))
    return version


def is_owner_or_allowed_else_forbidden(collection_version, user_info):
    if not user_info.is_user_owner_or_allowed(collection_version.owner):
        raise ForbiddenHTTPException()
