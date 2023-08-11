from datetime import datetime
from typing import List, Optional, Tuple, Union

from backend.layers.business.exceptions import DatasetVersionNotFoundException
from backend.layers.common.entities import (
    CollectionVersion,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    DatasetVersion,
    PublishedDatasetVersion,
)


def get_published_at_and_collection_version_id_else_not_found(
    dataset_version: DatasetVersion,
    collection_versions: Union[List[CollectionVersionWithDatasets], List[CollectionVersion]],
) -> Tuple[Optional[datetime], CollectionVersionId]:
    """
    Determines the published_at date for a DatasetVersion and the id of the CollectionVersion under which it was initially published.
    """
    for collection_version in sorted(
        #  Alphabetically sorts 'False' before 'True'; iterates over published versions first
        collection_versions,
        key=lambda cv: (cv.published_at is None, cv.published_at),
    ):
        if collection_version.published_at is not None:
            if isinstance(collection_version, CollectionVersionWithDatasets) and dataset_version.version_id.id in {
                dv.version_id.id for dv in collection_version.datasets
            }:
                return collection_version.published_at, collection_version.version_id
            elif isinstance(collection_version, CollectionVersion) and dataset_version.version_id.id in {
                dv.id for dv in collection_version.datasets
            }:
                return collection_version.published_at, collection_version.version_id
        else:
            break
    raise DatasetVersionNotFoundException("No such published Dataset version")


def set_revised_at_field(dataset_versions: List[DatasetVersion], collection_versions: List[CollectionVersion]) -> None:
    """
    Sets the `revised_at` field on the CanonicalDataset object for each DatasetVersion object
    """
    for dataset_version in dataset_versions:
        try:
            version_published_at, collection_version_id = get_published_at_and_collection_version_id_else_not_found(
                dataset_version, collection_versions
            )
            if version_published_at > dataset_version.canonical_dataset.published_at:  # type: ignore
                # Dataset has been revised
                dataset_version.canonical_dataset.revised_at = version_published_at
        except DatasetVersionNotFoundException:
            # Dataset has never been published
            pass


def get_dataset_versions_with_published_at_and_collection_version_id(
    dataset_versions: List[DatasetVersion],
    all_collection_versions: Union[List[CollectionVersionWithDatasets], List[CollectionVersion]],
) -> List[PublishedDatasetVersion]:
    """
    Takes a list of DatasetVersion entities and determines the initial data of publish as well as the id of
    the CollectionVersion for that initial publication, then returns a list of the corresponding
    PublishedDatasetVersion which contains these two extra values.
    """
    published_datasets_for_collection: List[PublishedDatasetVersion] = []
    for dataset_version in dataset_versions:
        published_at, collection_version_id = get_published_at_and_collection_version_id_else_not_found(
            dataset_version, all_collection_versions
        )
        published_datasets_for_collection.append(
            PublishedDatasetVersion(
                collection_version_id=collection_version_id,
                published_at=published_at,  # type: ignore
                **vars(dataset_version),
            )
        )
    return published_datasets_for_collection
