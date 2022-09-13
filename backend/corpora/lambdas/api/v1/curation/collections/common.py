import typing

from sqlalchemy.orm import Session
from urllib.parse import urlparse

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
    ProcessingStatus,
    Base,
    IsPrimaryData,
    ProjectLinkType,
    ValidationStatus,
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


def extract_doi_from_links(collection: dict):
    """
    Pull out the DOI from the 'links' array and return it along with the altered links array
    :param collection: the Collection
    :return: None
    """
    doi, resp_links = None, []
    for link in collection.get("links", []):
        if link["link_type"] == ProjectLinkType.DOI:
            doi = urlparse(link["link_url"]).path.strip("/")
        else:
            resp_links.append(link)

    return doi, resp_links


def reshape_for_curation_api(
    db_session: Session, collection: DbCollection, token_info: dict, preview: bool = False
) -> dict:
    """
    Reshape Collection data for the Curation API response. Remove tombstoned Datasets.
    :param db_session: the db Session
    :param collection: the Collection being returned in the API response
    :param token_info: user access token
    :param preview: boool - whether the dataset is in preview form or not.
    :return: whether or not the Collection should be included in the response per ownership/access rules
    """
    entity_columns = EntityColumns.columns_for_collections if preview else EntityColumns.columns_for_collection_id
    resp_collection = collection.to_dict_keep(entity_columns)
    resp_collection["processing_status"] = add_collection_level_processing_status(collection)

    owner = resp_collection.pop("owner")  # Don't actually want to return 'owner' in response
    resp_collection["revising_in"] = get_revising_in(db_session, collection, token_info, owner)
    resp_collection["collection_url"] = f"{CorporaConfig().collections_base_url}/collections/{collection.id}"
    if datasets := resp_collection.get("datasets"):
        resp_collection["datasets"] = reshape_datasets_for_curation_api(datasets, preview)

    resp_collection["doi"], resp_collection["links"] = extract_doi_from_links(resp_collection)

    return resp_collection


def get_revising_in(
    db_session: Session, collection: DbCollection, token_info: dict, owner: str
) -> typing.Optional[str]:
    """
    If the Collection is public AND the user is authorized use a database call to populate 'revising_in' attribute.
    None -> 1) revision does not exist, 2) the Collection is private, or 3) the user is not authorized
    "<revision_id>" -> user is authorized and revision exists
    :param db_session: the db Session
    :param collection: the Collection
    :param token_info: the user's access token info
    :param owner: the owner of the Collection
    :return: None
    """
    if collection.visibility == CollectionVisibility.PUBLIC and is_user_owner_or_allowed(token_info, owner):
        if result := db_session.query(DbCollection.id).filter(DbCollection.revision_of == collection.id).one_or_none():
            return result.id
    return None


def reshape_datasets_for_curation_api(datasets: typing.List[dict], preview=False) -> typing.List[dict]:
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
    IsPrimaryData.PRIMARY: [True],
    IsPrimaryData.SECONDARY: [False],
    IsPrimaryData.BOTH: [True, False],
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
        "curator_tag",
        "tissue",
        "assay",
        "disease",
        "organism",
        "tombstone",
    ]

    dataset_cols = [
        *dataset_preview_cols,
        "original_id",
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
        "batch_condition",
        "mean_genes_per_cell",
        "schema_version",
        "processing_status",
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


def add_collection_level_processing_status(collection: DbCollection) -> str:
    # Add a Collection-level processing status
    if not collection.datasets:  # Return None if the collection has no datasets.
        return None
    return_status = ProcessingStatus.SUCCESS
    for dataset in collection.datasets:
        if processing_status := dataset.processing_status:
            status = processing_status.processing_status
            if status == ProcessingStatus.PENDING:
                return_status = status
            elif status == ProcessingStatus.FAILURE:
                return status
    return return_status
