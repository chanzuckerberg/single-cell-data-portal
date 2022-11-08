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
    upload_status: Optional[DatasetUploadStatus]
    validation_status: Optional[DatasetValidationStatus]
    cxg_status: Optional[DatasetConversionStatus]
    rds_status: Optional[DatasetConversionStatus]
    h5ad_status: Optional[DatasetConversionStatus]
    processing_status: Optional[DatasetProcessingStatus]

    @staticmethod 
    def empty():
        return DatasetStatus(None, None, None, None, None, None)

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
class OntologyTermId:
    label: str
    ontology_term_id: str

@dataclass
class DatasetMetadata:
    organism: List[OntologyTermId]
    tissue: List[OntologyTermId]
    assay: List[OntologyTermId]
    disease: List[OntologyTermId]
    sex: List[OntologyTermId]
    self_reported_ethnicity: List[OntologyTermId]
    development_stage: List[OntologyTermId]
    cell_type: List[OntologyTermId]
    cell_count: int
    schema_version: str
    mean_genes_per_cell: float
    batch_condition: List[str]
    suspension_type: List[str]
    donor_id: List[str]
    is_primary_data: str

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
