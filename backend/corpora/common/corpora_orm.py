import enum

from datetime import datetime
from sqlalchemy import (
    Boolean,
    Column,
    DateTime,
    Enum,
    Float,
    ForeignKey,
    ForeignKeyConstraint,
    Integer,
    String,
    UniqueConstraint,
    types,
)
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.ext.declarative import declarative_base, DeclarativeMeta
from sqlalchemy.orm import relationship
from typing import Optional, List
from uuid import uuid4

from .utils.exceptions import CorporaException


def generate_uuid():
    return str(uuid4())


class StrippedString(types.TypeDecorator):
    """
    Returns a string with spaces stripped.
    """

    impl = types.String

    def process_result_value(self, value, dialect):
        """
        Strip the trailing spaces on resulting values.
        If value is false, we return it as-is; it might be none
        for nullable columns.
        """
        return value.strip() if value else value

    def copy(self):
        """
        Make a copy of this type.
        """
        return StrippedString(self.impl.length)


class TransformingBase(object):
    """
    Add functionality to transform a Base object, and recursively transform its linked entities.
    """

    def __iter__(self):
        return iter(self.to_dict().items())

    def to_dict(
        self,
        backref: List["Base"] = None,
        remove_none: bool = False,
        remove_attr: Optional[List[str]] = None,
        remove_relationships: bool = False,
    ) -> dict:
        """
        Converts the columns and relationships of a SQLAlchemy Base object into a python dictionary.

        :param backref: used to avoid recursively looping between two tables.
        :param remove_none: If true, removes keys that are none from the result.
        :param remove_attr: Attributes not to convert.
        :param remove_relationships: Ignore relationships.
        :return: a dictionary representation of the database object.
        """
        result = dict()
        remove_attr = remove_attr if remove_attr else []

        # Populate result with columns.
        for attr, column in self.__mapper__.c.items():
            if column.key in remove_attr:
                continue
            if remove_none:
                if getattr(self, attr) is not None:
                    result[column.key] = getattr(self, attr)
                else:
                    continue
            else:
                result[column.key] = getattr(self, attr)

        # Populate result with relationships.
        if not remove_relationships:
            if backref:
                backref.append(self.__table__)
            else:
                backref = [self.__table__]

            for attr, relation in self.__mapper__.relationships.items():
                if attr in remove_attr:
                    continue
                # Avoid recursive loop between multiple tables.
                if relation.target in backref:
                    continue
                value = getattr(self, attr)
                if value is None:
                    if not remove_none:
                        result[relation.key] = None
                elif isinstance(value.__class__, DeclarativeMeta):
                    result[relation.key] = value.to_dict(backref=backref, remove_none=remove_none)
                elif isinstance(value, list):
                    result[relation.key] = [i.to_dict(backref=backref, remove_none=remove_none) for i in value]
                else:
                    raise CorporaException(f"Unable to convert to dictionary. Unexpected type: {type(value)}.")
            backref.pop()
        return result

    id = Column(String, primary_key=True, default=generate_uuid)


class AuditMixin(object):
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow)


class TimestampMixin(object):
    published_at = Column(DateTime, nullable=True)
    revised_at = Column(DateTime, nullable=True)


Base = declarative_base(cls=TransformingBase)


class CollectionVisibility(enum.Enum):
    """
    Describes a DbCollection's visibility.
    At most, one LIVE and one EDIT entry of a Collection may exist at a time.

    PUBLIC - a published and publicly viewable Collection.
    PRIVATE - an open Submission, i.e an unpublished and non-public Collection.
    """

    PUBLIC = "Public"
    PRIVATE = "Private"


class XApproximateDistribution(enum.Enum):
    """
    Describes a DbDataset's x_approximate_distribution.

    COUNT - for data whose distributions are best approximated by counting distributions
            like Poisson, Binomial, or Negative Binomial.
    NORMAL - for data whose distributions are best approximated by the Gaussian distribution.
    """

    COUNT = "count"
    NORMAL = "normal"


class IsPrimaryData(enum.Enum):
    """
    Describes a DbDataset's is_primary_data.

    PRIMARY - when all observation values are True.
    SECONDARY - when all observation values are False.
    BOTH - when observation values are either True or False.
    """

    PRIMARY = "primary"
    SECONDARY = "secondary"
    BOTH = "both"


class ProjectLinkType(enum.Enum):
    """
    Enumerates DbCollection external web link types.

    PROTOCOL - A link to a sequencing protocol.
    RAW_DATA - A link to a raw data repository.
    OTHER - Other.
    """

    DOI = "doi"
    RAW_DATA = "raw_data"
    PROTOCOL = "protocol"
    LAB_WEBSITE = "lab_website"
    OTHER = "other"
    DATA_SOURCE = "data_source"


# provide a consistent name
CollectionLinkType = ProjectLinkType


class DatasetArtifactFileType(enum.Enum):
    """
    Enumerates DatasetArtifact file types.

    H5AD - An AnnData object describing an expression matrix. Uses the .h5ad extension.
    RDS - A Seurat file object describing an expression matrix. Uses the .rds extension.
    LOOM - Removed. No longer supported. (#1427)
    CXG - A TileDb object describing a cellxgene object. Uses .cxg extension.
    """

    H5AD = "h5ad"
    RDS = "rds"
    CXG = "cxg"


class DatasetArtifactType(enum.Enum):
    """
    Enumerates DatasetArtifact types.

    ORIGINAL - A data artifact that adheres to the minimal metadata schema requirements.
    REMIX - A data artifact that adheres to the Corpora metadata schema requirements.
    """

    ORIGINAL = "Original"
    REMIX = "Remix"


class DbCollection(Base, AuditMixin, TimestampMixin):
    """
    A Corpora collection represents an in progress or live submission of a lab experiment.
    DbCollections are associated with one or more single-cell datasets and links to external repositories.
    """

    # the tablename is "project" instead of "collection" to avoid migrating the database
    __tablename__ = "project"

    visibility = Column(Enum(CollectionVisibility), primary_key=True, nullable=False)
    owner = Column(StrippedString, nullable=False)
    name = Column(StrippedString)
    description = Column(StrippedString)
    obfuscated_uuid = Column(String, default="")
    contact_name = Column(StrippedString, default="")
    contact_email = Column(StrippedString, default="")
    data_submission_policy_version = Column(StrippedString, nullable=True)
    tombstone = Column(Boolean, default=False, nullable=False)

    # Relationships
    links = relationship("DbProjectLink", back_populates="collection", cascade="all, delete-orphan")
    datasets = relationship("DbDataset", back_populates="collection", cascade="all, delete-orphan")
    genesets = relationship("DbGeneset", back_populates="collection", cascade="all, delete-orphan")


class DbProjectLink(Base, AuditMixin):
    """
    Represents an external web link for DbCollections such as protocols and supplementary data repositories.
    """

    # the tablename is "project_link" instead of "collection_link" to avoid migrating the database
    __tablename__ = "project_link"

    collection_id = Column(String, nullable=False)
    collection_visibility = Column(Enum(CollectionVisibility), nullable=False)
    link_name = Column(StrippedString)
    link_url = Column(StrippedString)
    link_type = Column(Enum(CollectionLinkType))

    # Relationships
    collection = relationship("DbCollection", uselist=False, back_populates="links")

    # Composite FK
    __table_args__ = (
        ForeignKeyConstraint([collection_id, collection_visibility], [DbCollection.id, DbCollection.visibility]),
        {},
    )


# provide a consistent name
DbCollectionLink = DbProjectLink


class DbDataset(Base, AuditMixin, TimestampMixin):
    """
    Models a single experiment uploaded and processed by Corpora.
    Describes experiment metadata such as specimen and assay data.
    Related data files are represented by DbDataArtifacts.
    """

    __tablename__ = "dataset"

    revision = Column(Integer)
    name = Column(String)
    organism = Column(JSONB)
    tissue = Column(JSONB)
    assay = Column(JSONB)
    disease = Column(JSONB)
    sex = Column(JSONB)
    ethnicity = Column(JSONB)
    development_stage = Column(JSONB)
    cell_type = Column(JSONB)
    cell_count = Column(Integer)
    is_valid = Column(Boolean, default=False)
    is_primary_data = Column(Enum(IsPrimaryData))
    collection_id = Column(String, nullable=False)
    collection_visibility = Column(Enum(CollectionVisibility), nullable=False)
    tombstone = Column(Boolean, default=False, nullable=False)
    original_id = Column(String)
    published = Column(Boolean, default=False)
    explorer_url = Column(String, index=True)
    x_normalization = Column(String)
    x_approximate_distribution = Column(Enum(XApproximateDistribution))
    mean_genes_per_cell = Column(Float, default=0.0)
    schema_version = Column(String)

    # Relationships
    collection = relationship("DbCollection", uselist=False, back_populates="datasets")
    artifacts = relationship("DbDatasetArtifact", back_populates="dataset", cascade="all, delete-orphan")
    processing_status = relationship(
        "DbDatasetProcessingStatus", back_populates="dataset", cascade="all, delete-orphan", uselist=False
    )
    genesets = relationship("DbGeneset", secondary="geneset_dataset_link", back_populates="datasets")

    # Composite FK
    __table_args__ = (
        ForeignKeyConstraint([collection_id, collection_visibility], [DbCollection.id, DbCollection.visibility]),
        {},
    )

    def to_dict(self, *args, **kwargs):
        kwargs["remove_attr"] = kwargs.get("remove_attr", []) + ["genesets"]
        result = super(Base, self).to_dict(*args, **kwargs)
        result["linked_genesets"] = [gs.id for gs in self.genesets]
        return result


class DbDatasetArtifact(Base, AuditMixin):
    """
    Represents a user uploaded or Corpora generated file linked to a DbDataset.
    All matrices and cellxgene objects are examples of a DbDatasetArtifact.
    """

    __tablename__ = "dataset_artifact"
    __mapper_args__ = {"confirm_deleted_rows": False}

    dataset_id = Column(ForeignKey("dataset.id"), nullable=False)
    filename = Column(String)
    filetype = Column(Enum(DatasetArtifactFileType))
    type = Column(Enum(DatasetArtifactType))
    user_submitted = Column(Boolean)
    s3_uri = Column(String)

    # Relationships
    dataset = relationship("DbDataset", uselist=False, back_populates="artifacts")


class UploadStatus(enum.Enum):
    """
    Enumerates the status of an upload

    NA - No associated upload with the dataset
    WAITING - The upload is enqueued, and waiting for the upload container
    UPLOADING - The file is actively being uploaded
    UPLOADED - The upload was completed successfully
    FAILED - The upload has failed
    CANCEL_PENDING - The upload is in the process of being canceled
    CANCELED - The upload has been canceled
    """

    NA = "N/A"
    WAITING = "Waiting"
    UPLOADING = "Uploading"
    UPLOADED = "Uploaded"
    FAILED = "Failed"
    CANCEL_PENDING = "Cancel pending"
    CANCELED = "Canceled"


class ValidationStatus(enum.Enum):
    """
    Enumerates the status of validation of an uploaded dataset file

    NA - No associated validation with the dataset
    VALIDATING - The validation script is running
    VALID - The uploaded file successfully passed validation
    INVALID - The uploaded file failed validation
    """

    NA = "N/A"
    VALIDATING = "Validating"
    VALID = "Valid"
    INVALID = "Invalid"


class ConversionStatus(enum.Enum):
    """
    Enumerates the status of conversion of a valid uploaded file into another file format

    NA - No associated conversion with the dataset, perhaps because the uploaded dataset file
         was already in this format.
    CONVERTING = The conversion script is running
    CONVERTED - Conversion completed
    UPLOADING - The file is being uploaded to the S3 artifact bucket
    UPLOADED - The file was successfully uploaded to the S3 artifact bucket and the dataset artifact was updated
    FAILED - Conversion failed
    """

    NA = "N/A"
    CONVERTING = "Converting"
    CONVERTED = "Converted"
    UPLOADING = "Uploading"
    UPLOADED = "Uploaded"
    FAILED = "Failed"


class ProcessingStatus(enum.Enum):
    """
    Enumerates the status of processing a dataset.

    PENDING = Processing has not started
    SUCCESS - Processing succeeded
    FAILURE - Processing failed
    """

    PENDING = "PENDING"
    SUCCESS = "SUCCESS"
    FAILURE = "FAILURE"


class DbDatasetProcessingStatus(Base, AuditMixin):
    """
    Represents progress and status of user-initiated upload, validation, and conversion.
    """

    __tablename__ = "dataset_processing_status"

    dataset_id = Column(ForeignKey("dataset.id"), nullable=False)
    upload_status = Column(Enum(UploadStatus))
    upload_progress = Column(Float)
    upload_message = Column(String)
    validation_status = Column(Enum(ValidationStatus))
    validation_message = Column(String)
    rds_status = Column(Enum(ConversionStatus))
    cxg_status = Column(Enum(ConversionStatus))
    h5ad_status = Column(Enum(ConversionStatus))
    processing_status = Column(Enum(ProcessingStatus))

    # Relationships
    dataset = relationship("DbDataset", back_populates="processing_status")


class DbGeneset(Base, AuditMixin):
    """
    Represents a geneset linking a list of genes to a collection and specific datasets within that collection
    """

    __tablename__ = "geneset"

    name = Column(String, nullable=False)
    description = Column(String)
    genes = Column(JSONB)
    collection_id = Column(String, nullable=False)
    collection_visibility = Column(Enum(CollectionVisibility), nullable=False)
    collection = relationship("DbCollection", uselist=False, back_populates="genesets")
    datasets = relationship("DbDataset", secondary="geneset_dataset_link", back_populates="genesets")

    __table_args__ = (
        ForeignKeyConstraint([collection_id, collection_visibility], [DbCollection.id, DbCollection.visibility]),
        UniqueConstraint("name", "collection_id", "collection_visibility", name="_geneset_name__collection_uc"),
    )

    def to_dict(self, *args, **kwargs):
        kwargs["remove_attr"] = kwargs.get("remove_attr", []) + ["datasets"]
        result = super(Base, self).to_dict(*args, **kwargs)
        result["linked_datasets"] = [ds.id for ds in self.datasets]
        return result


class DbGenesetDatasetLink(Base, AuditMixin):
    """
    Represents a link between a geneset and a dataset supporting a many to many relationship
    """

    __tablename__ = "geneset_dataset_link"

    geneset_id = Column(String, ForeignKey("geneset.id"), index=True)
    dataset_id = Column(String, ForeignKey("dataset.id"), index=True)
