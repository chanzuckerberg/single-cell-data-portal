import uuid
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import List, Optional
from urllib.parse import urlparse

from dataclasses_json import dataclass_json

# TODO: copy and paste the docs for these


class DatasetStatusKey(str, Enum):
    UPLOAD = "upload"
    VALIDATION = "validation"
    CXG = "cxg"
    RDS = "rds"
    H5AD = "h5ad"
    PROCESSING = "processing"


class DatasetStatusGeneric:
    pass


class DatasetProcessingStatus(DatasetStatusGeneric, Enum):
    """
    Enumerates the status of processing a dataset.

    INITIALIZED = Dataset id created, and awaiting upload.
    PENDING = Processing has not started
    SUCCESS = Processing succeeded
    FAILURE = Processing failed
    """

    INITIALIZED = "INITIALIZED"
    PENDING = "PENDING"
    SUCCESS = "SUCCESS"
    FAILURE = "FAILURE"


class DatasetUploadStatus(DatasetStatusGeneric, Enum):
    NA = "NA"
    WAITING = "WAITING"
    UPLOADING = "UPLOADING"
    UPLOADED = "UPLOADED"
    FAILED = "FAILED"
    CANCEL_PENDING = "CANCEL PENDING"
    CANCELED = "CANCELED"


class DatasetValidationStatus(DatasetStatusGeneric, Enum):
    NA = "NA"
    VALIDATING = "VALIDATING"
    VALID = "VALID"
    INVALID = "INVALID"


class DatasetConversionStatus(DatasetStatusGeneric, Enum):
    NA = "NA"
    CONVERTING = "CONVERTING"
    CONVERTED = "CONVERTED"
    UPLOADING = "UPLOADING"
    UPLOADED = "UPLOADED"
    FAILED = "FAILED"
    SKIPPED = "SKIPPED"


class CollectionLinkType(str, Enum):
    DOI = "doi"
    RAW_DATA = "raw_data"
    PROTOCOL = "protocol"
    LAB_WEBSITE = "lab_website"
    OTHER = "other"
    DATA_SOURCE = "data_source"


class DatasetArtifactType(str, Enum):
    RAW_H5AD = "raw_h5ad"
    H5AD = "h5ad"
    RDS = "rds"
    CXG = "cxg"


class CollectionVisibility(Enum):
    """
    Describes a DbCollection's visibility.
    At most, one LIVE and one EDIT entry of a Collection may exist at a time.

    PUBLIC - a published and publicly viewable Collection.
    PRIVATE - an open Submission, i.e an unpublished and non-public Collection.
    """

    PUBLIC = "Public"
    PRIVATE = "Private"


@dataclass_json
@dataclass
class DatasetStatus:
    upload_status: Optional[DatasetUploadStatus]
    validation_status: Optional[DatasetValidationStatus]
    cxg_status: Optional[DatasetConversionStatus]
    rds_status: Optional[DatasetConversionStatus]
    h5ad_status: Optional[DatasetConversionStatus]
    processing_status: Optional[DatasetProcessingStatus]
    validation_message: Optional[str] = None

    @staticmethod
    def empty():
        return DatasetStatus(None, None, None, None, None, None)


@dataclass
class EntityId:
    id: str

    def __init__(self, entity_id: str = None):
        self.id = str(entity_id) if entity_id is not None else str(uuid.uuid4())

    def __repr__(self) -> str:
        return self.id


class CollectionId(EntityId):
    pass


class CollectionVersionId(EntityId):
    pass


class DatasetId(EntityId):
    pass


class DatasetVersionId(EntityId):
    pass


class DatasetArtifactId(EntityId):
    pass


@dataclass
class DatasetArtifact:
    id: DatasetArtifactId
    type: DatasetArtifactType
    uri: str

    def get_file_name(self):
        return urlparse(self.uri).path.split("/")[-1]


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
    mean_genes_per_cell: float
    batch_condition: List[str]
    suspension_type: List[str]
    donor_id: List[str]
    is_primary_data: str
    x_approximate_distribution: Optional[str]
    primary_cell_count: Optional[int] = None
    feature_count: Optional[int]
    feature_biotype: Optional[List[str]]
    feature_reference: Optional[List[str]]
    default_embedding: Optional[str]
    embeddings: Optional[List[str]]


@dataclass
class CanonicalDataset:
    dataset_id: DatasetId
    dataset_version_id: Optional[DatasetVersionId]
    tombstoned: bool
    published_at: Optional[datetime] = None
    revised_at: Optional[datetime] = None  # The last time this Dataset Version was Published


@dataclass
class DatasetVersion:
    dataset_id: DatasetId
    version_id: DatasetVersionId
    collection_id: CollectionId  # Pointer to the canonical collection id this dataset belongs to
    status: DatasetStatus
    metadata: Optional[DatasetMetadata]
    artifacts: List[DatasetArtifact]
    created_at: datetime
    canonical_dataset: CanonicalDataset


@dataclass
class PublishedDatasetVersion(DatasetVersion):
    collection_version_id: CollectionVersionId  # Pointer to collection version it was originally published under
    published_at: datetime
    revised_at: datetime = None


@dataclass
class Link:
    name: Optional[str]
    type: str
    uri: str

    def strip_fields(self):
        if self.name:
            self.name = self.name.strip()
        self.type = self.type.strip()
        self.uri = self.uri.strip()


@dataclass_json
@dataclass
class CollectionMetadata:
    name: str
    description: str
    contact_name: str
    contact_email: str
    links: List[Link]
    consortia: List[str] = field(default_factory=list)


@dataclass
class CanonicalCollection:
    id: CollectionId
    version_id: Optional[CollectionVersionId]  # Needs to be optional, or not exist
    originally_published_at: Optional[datetime]
    revised_at: Optional[datetime]
    tombstoned: bool


@dataclass
class CollectionVersionBase:
    collection_id: CollectionId
    version_id: CollectionVersionId
    owner: str
    curator_name: str
    metadata: CollectionMetadata
    publisher_metadata: Optional[dict]  # TODO: use a dataclass
    published_at: Optional[datetime]
    created_at: datetime
    schema_version: str
    canonical_collection: CanonicalCollection

    def is_published(self) -> bool:
        """
        This collection version has been published.
        TODO: After old API code is removed consider moving closer to API layer
        """
        return self.published_at is not None and self.canonical_collection.originally_published_at is not None

    def is_unpublished_version(self) -> bool:
        """
        The collection has been published, and this is a unpublished version of the collection.
        TODO: After old API code is removed consider moving closer to API layer
        """
        return self.published_at is None and self.canonical_collection.originally_published_at is not None

    def is_initial_unpublished_version(self) -> bool:
        """
        The collection is unpublished, this version is unpublished, and no previous versions have been
        published.
        TODO: After old API code is removed consider moving closer to API layer
        """
        return self.published_at is None and self.canonical_collection.originally_published_at is None


@dataclass
class CollectionVersion(CollectionVersionBase):
    datasets: List[DatasetVersionId]


@dataclass
class CollectionVersionWithDatasets(CollectionVersionBase):
    datasets: List[DatasetVersion]


@dataclass
class CollectionVersionWithPublishedDatasets(CollectionVersionBase):
    datasets: List[PublishedDatasetVersion]
