import enum
import os
import sys
from datetime import datetime

from sqlalchemy import (
    Boolean,
    Column,
    create_engine,
    DateTime,
    Enum,
    ForeignKey,
    ForeignKeyConstraint,
    Integer,
    String,
)
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.ext.declarative import declarative_base, DeclarativeMeta
from sqlalchemy.orm import relationship, sessionmaker

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from .corpora_config import CorporaDbConfig
from .utils.exceptions import CorporaException


class TransformingBase(object):
    """
    Add functionality to transform a Base object, and recursively transform its linked entities.
    """

    def __iter__(self):
        return iter(self.to_dict().items())

    def to_dict(self, backref: "Base" = None) -> dict:
        """
        Converts the columns and relationships of a SQLAlchemy Base object into a python dictionary.

        :param backref: used to avoid recursively looping between two tables.
        :return: a dictionary representation of the database object.
        """

        # Populate result with columns.
        result = {column.key: getattr(self, attr) for attr, column in self.__mapper__.c.items()}

        # Populate result with relationships.
        for attr, relation in self.__mapper__.relationships.items():
            # Avoid recursive loop between two tables.
            if backref == relation.target:
                continue
            value = getattr(self, attr)
            if value is None:
                result[relation.key] = None
            elif isinstance(value.__class__, DeclarativeMeta):
                result[relation.key] = value.to_dict(backref=self.__table__)
            elif isinstance(value, list):
                result[relation.key] = [i.to_dict(backref=self.__table__) for i in value]
            else:
                raise CorporaException(f"Unable to convert to dictionary. Unexpected type: {type(value)}.")
        return result


Base = declarative_base(cls=TransformingBase)


class DBSessionMaker:
    def __init__(self):
        self.engine = create_engine(CorporaDbConfig().database_uri, connect_args={"connect_timeout": 5})
        self.session_maker = sessionmaker(bind=self.engine)

    def session(self, **kwargs):
        return self.session_maker(**kwargs)


class CollectionVisibility(enum.Enum):
    """
    Describes a DbCollection's visibility.
    At most, one LIVE and one EDIT entry of a Collection may exist at a time.

    PUBLIC - a published and publicly viewable Collection.
    PRIVATE - an open Submission, i.e an unpublished and non-public Collection.
    """

    PUBLIC = "Public"
    PRIVATE = "Private"


class ProcessingState(enum.Enum):
    """
    Enumerates DbCollection states in the data processing pipeline from upload to deployment.

    NA - Not in the data processing pipeline which can represent pre or post completion of the pipeline.
    IN_VALIDATION - Following submission, validate datasets for required metadata and absence of PII.
    IN_ARTIFACT_CREATION - Following validation, create all Original + Remix matrix formats and cellxgene objects.
    IN_DEPLOYMENT - The final stage in the pipeline: deploying artifacts to Data Portal and cellxgene applications.
    """

    NA = "N/A"
    IN_VALIDATION = "In validation"
    IN_ARTIFACT_CREATION = "In artifact creation"
    IN_DEPLOYMENT = "In deployment"


class ValidationState(enum.Enum):
    """
    Enumerates DbCollection validation states.

    NOT_VALIDATED - Validation not performed yet.
    VALID - Collection is valid.
    INVALID - Collection is invalid.
    """

    NOT_VALIDATED = "Not validated"
    VALID = "Valid"
    INVALID = "Invalid"


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


# provide a consistent name
CollectionLinkType = ProjectLinkType


class DatasetArtifactFileType(enum.Enum):
    """
    Enumerates DatasetArtifact file types.

    H5AD - An AnnData object describing an expression matrix. Uses the .h5ad extension.
    RDS - A Seurat file object describing an expression matrix. Uses the .rds extension.
    LOOM - A AnnData object describing an expression matrix. Uses the .loom extension.
    CXG - A TileDb object describing a cellxgene object. Uses .cxg extension.
    """

    H5AD = "h5ad"
    RDS = "rds"
    LOOM = "loom"
    CXG = "cxg"


class DatasetArtifactType(enum.Enum):
    """
    Enumerates DatasetArtifact types.

    ORIGINAL - A data artifact that adheres to the minimal metadata schema requirements.
    REMIX - A data artifact that adheres to the Corpora metadata schema requirements.
    """

    ORIGINAL = "Original"
    REMIX = "Remix"


class DbCollection(Base):
    """
    A Corpora collection represents an in progress or live submission of a lab experiment.
    DbCollections are associated with one or more single-cell datasets and links to external repositories.
    """

    # the tablename is "project" instead of "collection" to avoid migrating the database
    __tablename__ = "project"

    id = Column(String, primary_key=True)
    visibility = Column(
        Enum(CollectionVisibility), primary_key=True, nullable=False
    )  # Enum(CollectionVisibility). Enum type unsupported for composite FKs.
    owner = Column(String, nullable=False)
    name = Column(String)
    description = Column(String)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow)
    obfuscated_uuid = Column(String, default="")
    contact_name = Column(String, default="")
    contact_email = Column(String, default="")
    data_submission_policy_version = Column(String, nullable=False)

    # Relationships
    links = relationship("DbProjectLink", back_populates="collection", cascade="all, delete-orphan")
    datasets = relationship("DbDataset", back_populates="collection", cascade="all, delete-orphan")


class DbProjectLink(Base):
    """
    Represents an external web link for DbCollections such as protocols and supplementary data repositories.
    """

    # the tablename is "project_link" instead of "collection_link" to avoid migrating the database
    __tablename__ = "project_link"

    id = Column(String, primary_key=True)
    collection_id = Column(String, nullable=False)
    collection_visibility = Column(Enum(CollectionVisibility), nullable=False)
    link_name = Column(String)
    link_url = Column(String)
    link_type = Column(Enum(CollectionLinkType))
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow)

    # Relationships
    collection = relationship("DbCollection", uselist=False, back_populates="links")

    # Composite FK
    __table_args__ = (
        ForeignKeyConstraint([collection_id, collection_visibility], [DbCollection.id, DbCollection.visibility]),
        {},
    )


# provide a consistent name
DbCollectionLink = DbProjectLink


class DbDataset(Base):
    """
    Models a single experiment uploaded and processed by Corpora.
    Describes experiment metadata such as specimen and assay data.
    Related data files are represented by DbDataArtifacts.
    """

    __tablename__ = "dataset"

    id = Column(String, primary_key=True)
    revision = Column(Integer)
    name = Column(String)
    organism = Column(JSONB)
    tissue = Column(JSONB)
    assay = Column(JSONB)
    disease = Column(JSONB)
    sex = Column(JSONB)
    ethnicity = Column(JSONB)
    development_stage = Column(JSONB)
    cell_count = Column(Integer)
    is_valid = Column(Boolean, default=False)
    collection_id = Column(String, nullable=False)
    collection_visibility = Column(Enum(CollectionVisibility), nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow)

    # Relationships
    collection = relationship("DbCollection", uselist=False, back_populates="datasets")
    artifacts = relationship("DbDatasetArtifact", back_populates="dataset", cascade="all, delete-orphan")
    deployment_directories = relationship(
        "DbDeploymentDirectory", back_populates="dataset", cascade="all, delete-orphan"
    )

    # Composite FK
    __table_args__ = (
        ForeignKeyConstraint([collection_id, collection_visibility], [DbCollection.id, DbCollection.visibility]),
        {},
    )


class DbDatasetArtifact(Base):
    """
    Represents a user uploaded or Corpora generated file linked to a DbDataset.
    All matrices and cellxgene objects are examples of a DbDatasetArtifact.
    """

    __tablename__ = "dataset_artifact"

    id = Column(String, primary_key=True)
    dataset_id = Column(ForeignKey("dataset.id"), nullable=False)
    filename = Column(String)
    filetype = Column(Enum(DatasetArtifactFileType))
    type = Column(Enum(DatasetArtifactType))
    user_submitted = Column(Boolean)
    s3_uri = Column(String)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow)

    # Relationships
    dataset = relationship("DbDataset", uselist=False, back_populates="artifacts")


class DbDeploymentDirectory(Base):
    """
    Represents the deployment of a dataset to a Corpora application.
    This entity only supports cellxgene deployments.
    """

    __tablename__ = "deployment_directory"

    id = Column(String, primary_key=True)
    dataset_id = Column(ForeignKey("dataset.id"), nullable=False)
    url = Column(String)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow)

    # Relationships
    dataset = relationship("DbDataset", uselist=False, back_populates="deployment_directories")
