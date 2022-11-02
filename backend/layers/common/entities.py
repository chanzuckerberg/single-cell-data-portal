from dataclasses import dataclass
from datetime import datetime
from typing import List, Optional
from enum import Enum


# TODO: copy and paste the docs for these

class DatasetStatusGeneric:
    pass

class DatasetUploadStatus(DatasetStatusGeneric, Enum):
    NA = "N/A"
    WAITING = "Waiting"
    UPLOADING = "Uploading"
    UPLOADED = "Uploaded"
    FAILED = "Failed"
    CANCEL_PENDING = "Cancel pending"
    CANCELED = "Canceled"

class DatasetValidationStatus(DatasetStatusGeneric, Enum):
    NA = "N/A"
    VALIDATING = "Validating"
    VALID = "Valid"
    INVALID = "Invalid"

class DatasetConversionStatus(DatasetStatusGeneric, Enum):
    NA = "N/A"
    CONVERTING = "Converting"
    CONVERTED = "Converted"
    UPLOADING = "Uploading"
    UPLOADED = "Uploaded"
    FAILED = "Failed"
    SKIPPED = "Skipped"

class DatasetProcessingStatus(DatasetStatusGeneric, Enum):
    INITIALIZED = "INITIALIZED"
    PENDING = "PENDING"
    SUCCESS = "SUCCESS"
    FAILURE = "FAILURE"

@dataclass
class DatasetStatus:
    upload_status: DatasetUploadStatus
    validation_status: DatasetValidationStatus
    cxg_status: DatasetConversionStatus
    rds_status: DatasetConversionStatus
    h5ad_status: DatasetConversionStatus
    processing_status: DatasetProcessingStatus

@dataclass
class CollectionId:
    id: str
@dataclass
class CollectionVersionId:
    id: str

@dataclass
class DatasetId:
    id: str

@dataclass
class DatasetVersionId:
    id: str

@dataclass
class DatasetArtifact:
    id: str
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
    dataset_id: DatasetId
    version_id: DatasetVersionId
    status: DatasetStatus
    metadata: DatasetMetadata
    artifacts: List[DatasetArtifact]


@dataclass
class Link:
    name: str
    type: str
    uri: str
    
@dataclass
class CollectionMetadata:
    name: str
    description: str
    contact_name: str
    contact_email: str
    links: List[Link]


@dataclass
class CollectionVersion:
    collection_id: CollectionId
    version_id: CollectionVersionId
    owner: str
    metadata: CollectionMetadata
    publisher_metadata: Optional[dict]  # TODO: use a dataclass
    datasets: List[DatasetVersion] 
    published_at: Optional[datetime]

class CollectionLinkType(Enum):
    DOI = "doi"
    RAW_DATA = "raw_data"
    PROTOCOL = "protocol"
    LAB_WEBSITE = "lab_website"
    OTHER = "other"
    DATA_SOURCE = "data_source"
