from dataclasses import dataclass
from datetime import datetime
from typing import List, Optional
from enum import Enum

import json
from dataclasses_json import dataclass_json


# TODO: copy and paste the docs for these


class DatasetStatusKey(Enum):
    UPLOAD = "upload"
    VALIDATION = "validation"
    CXG = "cxg"
    RDS = "rds"
    H5AD = "h5ad"
    PROCESSING = "processing"


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


@dataclass_json
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

    def to_json(self):
        return json.dumps(self, default=lambda obj: obj.__dict__)


@dataclass(eq=True, frozen=True)
class CollectionId:
    id: str

    def __repr__(self) -> str:
        return self.id


@dataclass
class CollectionVersionId:
    id: str

    def __repr__(self) -> str:
        return self.id


@dataclass
class DatasetId:
    id: str

    def __repr__(self) -> str:
        return self.id


@dataclass
class DatasetVersionId:
    id: str

    def __repr__(self) -> str:
        return self.id


@dataclass
class DatasetArtifactId:
    id: str

    def __repr__(self) -> str:
        return self.id


@dataclass
class DatasetArtifact:
    id: DatasetArtifactId
    type: str
    uri: str


@dataclass
class OntologyTermId:
    label: str
    ontology_term_id: str


@dataclass_json
@dataclass
class DatasetMetadata:
    name: str
    schema_version: str
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
    x_approximate_distribution: str

    def to_json(self):
        return json.dumps(self, default=lambda obj: obj.__dict__)


@dataclass
class CanonicalDataset:
    dataset_id: DatasetId
    dataset_version_id: DatasetVersionId
    published_at: Optional[datetime]


@dataclass
class DatasetVersion:
    dataset_id: DatasetId
    version_id: DatasetVersionId
    collection_id: CollectionId  # Pointer to the canonical collection id this dataset belongs to
    status: DatasetStatus
    metadata: DatasetMetadata
    artifacts: List[DatasetArtifact]
    created_at: datetime
    canonical_dataset: CanonicalDataset


@dataclass
class Link:
    name: Optional[str]
    type: str
    uri: str


@dataclass_json
@dataclass
class CollectionMetadata:
    name: str
    description: str
    contact_name: str
    contact_email: str
    links: List[Link]

@dataclass
class CanonicalCollection:
    id: CollectionId
    version_id: CollectionVersionId # Needs to be optional, or not exist
    originally_published_at: Optional[datetime]
    tombstoned: bool


@dataclass
class CollectionVersionBase:
    collection_id: CollectionId
    version_id: CollectionVersionId
    owner: str
    metadata: CollectionMetadata
    publisher_metadata: Optional[dict]  # TODO: use a dataclass
    published_at: Optional[datetime]
    created_at: datetime
    canonical_collection: CanonicalCollection

@dataclass
class CollectionVersion(CollectionVersionBase):
    datasets: List[DatasetVersionId]

@dataclass
class CollectionVersionWithDatasets(CollectionVersionBase):
    datasets: List[DatasetVersion]


class CollectionLinkType(Enum):
    DOI = "doi"
    RAW_DATA = "raw_data"
    PROTOCOL = "protocol"
    LAB_WEBSITE = "lab_website"
    OTHER = "other"
    DATA_SOURCE = "data_source"


class DatasetArtifactType(Enum):
    """
    Enumerates DatasetArtifact file types.

    H5AD - An AnnData object describing an expression matrix, post-processing by cellxgene pipeline.
        Uses the .h5ad extension.
    RAW_H5AD - An AnnData object describing an expression matrix, as directly uploaded by users.
        Uses the .h5ad extension.
    RDS - A Seurat file object describing an expression matrix. Uses the .rds extension.
    CXG - A TileDb object describing a cellxgene object. Uses .cxg extension.
    """

    RAW_H5AD = "raw_h5ad"
    H5AD = "h5ad"
    RDS = "rds"
    CXG = "cxg"
