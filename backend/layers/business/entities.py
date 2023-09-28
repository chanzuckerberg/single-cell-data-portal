from dataclasses import dataclass
from typing import List, Optional

from backend.layers.common.entities import DatasetArtifactType, Link


@dataclass
class CollectionQueryFilter:
    is_published: Optional[bool] = None
    owner: Optional[str] = None
    curator_name: Optional[str] = None
    # TODO: add list of fields to be returned (if needed)


@dataclass
class DatasetArtifactDownloadData:
    file_size: int
    url: str


@dataclass
class CollectionMetadataUpdate:
    """
    This class can be used to issue an update to the collection metadata.
    Since we support partial updates, i.e. missing fields will be ignored,
    all the fields are marked an optional
    """

    name: Optional[str]
    description: Optional[str]
    contact_name: Optional[str]
    contact_email: Optional[str]
    links: Optional[List[Link]]
    consortia: Optional[List[str]]
