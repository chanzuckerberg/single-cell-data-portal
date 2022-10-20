from dataclasses import dataclass
from datetime import datetime
from typing import List


@dataclass
class DatasetStatus:
    status: str  # TODO: use an enum


@dataclass
class DatasetArtifact:
    type: str
    uri: str


@dataclass
class DatasetMetadata:
    organism: str
    tissue: str
    assay: str
    disease: str
    sex: str
    self_reported_ethnicity: str
    development_stage: str
    cell_type: str
    cell_count: int


@dataclass
class Dataset:
    id: str
    status: DatasetStatus
    metadata: DatasetMetadata
    artifacts: List[DatasetArtifact]


@dataclass
class DatasetVersion:
    version_id: str
    dataset: Dataset


@dataclass
class CollectionMetadata:
    name: str
    description: str
    owner: str
    contact_name: str
    contact_email: str
    links: List[dict]  # TODO: use a dataclass


@dataclass
class Collection:
    id: str
    metadata: CollectionMetadata
    publisher_metadata: dict  # TODO: use a dataclass
    datasets: List[DatasetVersion]


@dataclass
class CollectionVersion:
    version_id: str
    collection: Collection
    published_at: datetime
