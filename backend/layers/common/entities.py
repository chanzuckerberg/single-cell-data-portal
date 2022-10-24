from dataclasses import dataclass
from datetime import datetime
from typing import List, Optional


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
class DatasetVersion:
    dataset_id: str
    version_id: str
    processing_status: Optional[DatasetStatus]
    metadata: DatasetMetadata
    artifacts: List[DatasetArtifact]


@dataclass
class CollectionMetadata:
    name: str
    description: str
    owner: str
    contact_name: str
    contact_email: str
    links: List[dict]  # TODO: use a dataclass


@dataclass
class CollectionVersion:
    collection_id: str
    version_id: str
    owner: str
    metadata: CollectionMetadata
    publisher_metadata: Optional[dict]  # TODO: use a dataclass
    datasets: List[DatasetVersion] 
    published_at: Optional[datetime]
