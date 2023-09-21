from dataclasses import asdict
from typing import List, Optional, Tuple, Union
from urllib.parse import urlparse
from uuid import UUID

from backend.common.corpora_config import CorporaConfig
from backend.common.utils.http_exceptions import ForbiddenHTTPException, GoneHTTPException, NotFoundHTTPException
from backend.layers.auth.user_info import UserInfo
from backend.layers.common.entities import (
    CollectionId,
    CollectionVersion,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    CollectionVersionWithPublishedDatasets,
    DatasetArtifact,
    DatasetArtifactType,
    DatasetId,
    DatasetProcessingStatus,
    DatasetValidationStatus,
    DatasetVersion,
    DatasetVersionId,
    Link,
    OntologyTermId,
    PublishedDatasetVersion,
)
from backend.portal.api.explorer_url import generate as generate_explorer_url
from backend.portal.api.providers import get_business_logic

allowed_dataset_asset_types = (DatasetArtifactType.H5AD, DatasetArtifactType.RDS)


def get_collections_base_url():
    return CorporaConfig().collections_base_url


def extract_dataset_assets(dataset_version: DatasetVersion):
    asset_list = list()
    for asset in dataset_version.artifacts:
        if asset.type not in allowed_dataset_asset_types:
            continue
        filesize = get_business_logic().s3_provider.get_file_size(asset.uri)
        if filesize is None:
            filesize = -1
        url = get_business_logic().generate_permanent_url(dataset_version.version_id, asset.type)
        result = {
            "filesize": filesize,
            "filetype": asset.type.upper(),
            "url": url,
        }
        asset_list.append(result)
    return _with_duplicates_removed(asset_list)  # TODO: de-dupe on DatasetArtifact insertion via unique constraint


def _with_duplicates_removed(asset_list: List[dict]) -> List[dict]:
    asset_types = set()
    assets = []
    for asset in asset_list:
        if asset["filetype"] in asset_types:
            continue
        asset_types.add(asset["filetype"])
        assets.append(asset)
    return assets


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
    collection_version: Union[CollectionVersion, CollectionVersionWithDatasets, CollectionVersionWithPublishedDatasets],
    user_info: UserInfo = None,
    reshape_for_version_endpoint: bool = False,
    preview: bool = False,
) -> dict:
    """
    Reshape Collection data for the Curation API response. Remove tombstoned Datasets.
    :param collection_version: the Collection Version being returned in the API response
    :param user_info:
    :param reshape_for_version_endpoint CollectionVersion is being returned in a version endpoint
    :param preview: bool - whether the dataset is in preview form or not.
    :return: the response.
    """
    business_logic = get_business_logic()
    is_published = collection_version.published_at is not None
    # get collection attributes based on endpoint type and published status
    collection_id = collection_version.collection_id.id
    if not reshape_for_version_endpoint:
        published_at = collection_version.canonical_collection.originally_published_at
        if is_published:
            # Published
            collection_url = f"{get_collections_base_url()}/collections/{collection_version.collection_id.id}"
            revision_of = None
            if not user_info or not user_info.is_user_owner_or_allowed(collection_version.owner):
                _revising_in = None
            else:
                _revising_in = business_logic.get_unpublished_collection_version_from_canonical(
                    collection_version.collection_id
                )
            revising_in = _revising_in.version_id.id if _revising_in else None
            use_canonical_url = True
        else:
            # Unpublished - need to determine if it's a revision or first time collection
            # For that, we look at whether the canonical collection is published
            is_revision = collection_version.canonical_collection.originally_published_at is not None
            if is_revision:
                # If it's a revision, both collection_id and collection_url need to point to the version_id,
                # and datasets should expose the private url (based on version_id)
                collection_id = collection_version.version_id.id
                collection_url = f"{get_collections_base_url()}/collections/{collection_version.version_id.id}"
                revision_of = collection_version.collection_id.id
                use_canonical_url = False
            else:
                # If it's an unpublished, unrevised collection, then collection_url will point to the permalink
                # (aka the link to the canonical_id) and the collection_id will point to version_id.
                # Also, revision_of should be None, and the datasets should expose the canonical url
                collection_url = f"{get_collections_base_url()}/collections/{collection_version.collection_id.id}"
                revision_of = None
                use_canonical_url = True
            revising_in = None
    else:
        collection_url = f"{get_collections_base_url()}/collections/{collection_version.version_id.id}"
        use_canonical_url = False
        published_at = collection_version.published_at
        revision_of = collection_version.collection_id.id
        revising_in = None

    # get collection dataset attributes
    response_datasets = sorted(
        reshape_datasets_for_curation_api(
            collection_version.datasets,
            use_canonical_url,
            preview,
            as_version=reshape_for_version_endpoint,
            is_published=is_published,
        ),
        key=lambda d: d["dataset_id"],  # For stable ordering
        reverse=True,  # To stay consistent with Datasets index endpoint sorting
    )

    # build response
    doi, links = extract_doi_from_links(collection_version.metadata.links)
    response = dict(
        collection_id=collection_id,
        collection_url=collection_url,
        collection_version_id=collection_version.version_id.id,
        consortia=collection_version.metadata.consortia,
        contact_email=collection_version.metadata.contact_email,
        contact_name=collection_version.metadata.contact_name,
        created_at=collection_version.created_at,
        curator_name=collection_version.curator_name,
        description=collection_version.metadata.description,
        doi=doi,
        links=links,
        name=collection_version.metadata.name,
        published_at=published_at,
        publisher_metadata=collection_version.publisher_metadata,
        visibility=get_visibility(collection_version),
    )
    if reshape_for_version_endpoint:
        response["dataset_versions"] = response_datasets
    else:
        response["datasets"] = response_datasets
        response.update(
            revised_at=collection_version.canonical_collection.revised_at,
            revising_in=revising_in,
            revision_of=revision_of,
        )
    if not is_published:
        response.update(processing_status=get_collection_level_processing_status(collection_version.datasets))
    return response


def reshape_datasets_for_curation_api(
    datasets: List[Union[DatasetVersionId, DatasetVersion]],
    use_canonical_url: bool,
    preview: bool = False,
    as_version: bool = False,
    is_published: bool = False,
) -> List[dict]:
    active_datasets = []
    for dv in datasets:
        dataset_version = get_business_logic().get_dataset_version(dv) if isinstance(dv, DatasetVersionId) else dv
        reshaped_dataset = reshape_dataset_for_curation_api(
            dataset_version, use_canonical_url, preview, as_canonical=not as_version, is_published=is_published
        )
        active_datasets.append(reshaped_dataset)
    return active_datasets


def reshape_dataset_for_curation_api(
    dataset_version: DatasetVersion,
    use_canonical_url: bool,
    preview=False,
    as_canonical=True,
    is_published=False,
) -> dict:
    ds = dict()

    # Determine what columns to include from the dataset
    columns = EntityColumns.dataset_metadata_preview_cols if preview else EntityColumns.dataset_metadata_cols
    # Get dataset metadata fields.
    # Metadata can be None if the dataset isn't still fully processed, so we account for that
    if dataset_version.metadata is not None:
        for column in columns:
            col = getattr(dataset_version.metadata, column, None)
            if isinstance(col, OntologyTermId):
                col = [asdict(col)]
            elif isinstance(col, list) and len(col) != 0 and isinstance(col[0], OntologyTermId):
                col = [asdict(i) for i in col]
            ds[column] = col

    ds["dataset_id"] = dataset_version.dataset_id.id
    ds["dataset_version_id"] = dataset_version.version_id.id
    # Get none preview specific dataset fields
    if not preview:
        ds["assets"] = extract_dataset_assets(dataset_version)
        ds["title"] = ds.pop("name", None)
        ds["tombstone"] = False  # TODO this will always be false. Remove in the future
        if dataset_version.metadata is not None:
            ds["is_primary_data"] = is_primary_data_mapping.get(ds.pop("is_primary_data"), [])
            if ds["x_approximate_distribution"]:
                ds["x_approximate_distribution"] = ds["x_approximate_distribution"].upper()
        if not is_published and (status := dataset_version.status):
            if status.processing_status == DatasetProcessingStatus.FAILURE:
                if status.validation_status == DatasetValidationStatus.INVALID:
                    ds["processing_status_detail"] = status.validation_message
                    ds["processing_status"] = "VALIDATION_FAILURE"
                else:
                    ds["processing_status"] = "PIPELINE_FAILURE"
            else:
                ds["processing_status"] = status.processing_status
        if isinstance(dataset_version, PublishedDatasetVersion):
            ds["collection_id"] = dataset_version.collection_id.id
            ds["collection_version_id"] = dataset_version.collection_version_id.id
            ds["published_at"] = dataset_version.published_at
        if as_canonical:
            ds["published_at"] = dataset_version.canonical_dataset.published_at
            ds["revised_at"] = dataset_version.canonical_dataset.revised_at
            ds["explorer_url"] = generate_explorer_url(dataset_version, use_canonical_url)
    return ds


is_primary_data_mapping = {
    "PRIMARY": [True],
    "SECONDARY": [False],
    "BOTH": [True, False],
}


class EntityColumns:
    collections_cols = [
        "collection_id",
        "collection_version_id",
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
        "primary_cell_count",
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


def get_visibility(collection_version: CollectionVersion) -> str:
    return "PUBLIC" if collection_version.published_at else "PRIVATE"


def validate_uuid_else_forbidden(_id: str):
    try:
        UUID(_id)
    except ValueError as e:
        raise ForbiddenHTTPException() from e


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


def get_inferred_collection_version(_id: str) -> CollectionVersionWithDatasets:
    """
    Infer the collection version from either a CollectionId or a CollectionVersionId and return the CollectionVersion,
    if currently mapped published version or open unpublished version.
    :param _id: identifies the collection version, whether by canonical id or version id
    :return: The CollectionVersion if it exists.
    """
    validate_uuid_else_forbidden(_id)
    business_logic = get_business_logic()

    # gets currently mapped collection version, or unpublished version if never published
    version = business_logic.get_collection_version_from_canonical(CollectionId(_id))
    if version is None:
        version = business_logic.get_collection_version(CollectionVersionId(_id), get_tombstoned=True)
        if version is None:
            raise NotFoundHTTPException()
        # Only allow fetch by Collection Version ID if unpublished revision of published collection
        if version.published_at is not None or version.canonical_collection.originally_published_at is None:
            raise ForbiddenHTTPException()

    if version.canonical_collection.tombstoned is True:
        raise GoneHTTPException()
    return version


def get_dataset_version_from_canonical_id(dataset_id: str) -> Optional[DatasetVersion]:
    """
    Get the dataset version from a DatasetId and return the DatasetVersion.
    :param dataset_id: identifies the dataset version
    :return: The DatasetVersion if it exists.
    """
    validate_uuid_else_forbidden(dataset_id)
    return get_business_logic().get_dataset_version_from_canonical(DatasetId(dataset_id))


def is_owner_or_allowed_else_forbidden(collection_version, user_info):
    if not user_info.is_user_owner_or_allowed(collection_version.owner):
        raise ForbiddenHTTPException()
