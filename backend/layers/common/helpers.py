from datetime import datetime
from typing import List, Optional, Tuple, Union

from backend.layers.business.exceptions import DatasetVersionNotFoundException
from backend.layers.common.entities import (
    CollectionVersion,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    DatasetVersion,
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
