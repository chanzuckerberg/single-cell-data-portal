import enum
import os
import sys

from sqlalchemy import (
    Boolean,
    Column,
    create_engine,
    DateTime,
    Enum,
    ForeignKey,
    Integer,
    String,
    text,
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from common.corpora_config import CorporaDbConfig


Base = declarative_base()
deployment_stage = os.environ["DEPLOYMENT_STAGE"]

# TODO: check syntax (SQL -> Postgres)
# The in-memory SQLite database used in unit tests does not support now() SQL syntax
# DEFAULT_DATETIME = "2000-01-01 00:00:00" if deployment_stage == "test" else text("now()")
DEFAULT_DATETIME = text("now()")


class DBSessionMaker:
    def __init__(self):
        connection = "sqlite:///:memory:" if deployment_stage == "test" else CorporaDbConfig().database_uri
        self.engine = create_engine(connection)
        self.session_maker = sessionmaker(bind=self.engine)

    def session(self, **kwargs):
        return self.session_maker(**kwargs)


class ProjectStatus(enum.Enum):
    """
    Describes a DbProject's status.
    At most, one LIVE and one EDITING entry of a Project may exist at a time.

    LIVE - a published and publicly viewable Project.
    EDITING - an open Submission, i.e an unpublished and non-public Project.
    """
    LIVE = "Live"
    EDITING = "Editing"


# TODO: Rethink.
#  Should this be per dataset?
#  Should it be expanded to represent individual datasets?
#  Or reported state represents slowest Dataset in the Project?
class ProcessingState(enum.Enum):
    """
    Enumerates DbProject states in the data processing pipeline from upload to deployment.

    NA - Not in the data processing pipeline which can represent pre or post completion of the pipeline.
    IN_VALIDATION - Following submission, validate datasets for required metadata and absence of PII.
    IN_ARTIFACT_CREATION - Following validation, create all Original+Remix matrix formats + cellxgene objects.
    IN_DEPLOYMENT - The final stage in the pipeline: deploying artifacts to Data Portal and cellxgene applications.
    """
    NA = "N/A"
    IN_VALIDATION = "In validation"
    IN_ARTIFACT_CREATION = "In artifact creation"
    IN_DEPLOYMENT = "In deployment"


class ValidationState(enum.Enum):
    """
    Enumerates DbProject validation states.

    NOT_VALIDATED - Validation not performed yet.
    VALID - Project is valid.
    INVALID - Project is invalid.
    """
    NOT_VALIDATED = "Not validated"
    VALID = "Valid"
    INVALID = "Invalid"


class ProjectLinkType(enum.Enum):
    """
    Enumerates DbProject external web link types.

    PROTOCOL - A link to a sequencing protocol.
    RAW_DATA - A link to a raw data repository.
    OTHER - Other.
    """
    PROTOCOL = "Protocol"
    RAW_DATA = "Raw data"
    OTHER = "Other"


class DatasetArtifactType(enum.Enum):
    """
    Enumerates DatasetArtifact types.

    ORIGINAL - A data artifact that adheres to the minimal metadata schema requirements.
    REMIX - A data artifact that adheres to the Corpora metadata schema requirements.
    """
    ORIGINAL = "Original"
    REMIX = "Remix"


class DbUser(Base):
    """
    A registered Corpora user.
    Maintains user details such as contact information and access control settings.
    """
    __tablename__ = "user"

    id = Column(String, primary_key=True)
    name = Column(String)
    email = Column(String)
    created_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)


class DbProject(Base):
    """
    A Corpora project represents an in progress or live submission of a lab experiment.
    DbProjects are associated with one or more single-cell datasets and links to external repositories.
    """
    __tablename__ = "project"

    id = Column(String, primary_key=True)
    status = Column(Enum(ProjectStatus), primary_key=True)
    submitter = Column(ForeignKey("user.id"), nullable=False)
    name = Column(String)
    description = Column(String)
    s3_bucket = Column(String)
    tc_uri = Column(String)
    needs_attestation = Column(Boolean)
    processing_state = Column(Enum(ProcessingState))
    validation_state = Column(Enum(ValidationState))
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)

    user = relationship("DbUser")


class DbProjectDataset(Base):
    """
    Associates a DbProject with a DbDataset.
    A DbProject may link to several DbDatasets.
    A DbDataset must belong to one DbProject.
    """
    __tablename__ = "project_dataset"

    id = Column(String, primary_key=True)
    project_id = Column(ForeignKey("project.id"), nullable=False)
    project_status = Column(ForeignKey("project.status"), nullable=False)
    dataset_id = Column(ForeignKey("dataset.id"), nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)

    project = relationship("DbProject")
    dataset = relationship("DbDataset")


class DbProjectLink(Base):
    """
    Represents an external web link for DbProjects such as protocols and supplementary data repositories.
    """
    __tablename__ = "project_link"

    id = Column(String, primary_key=True)
    project_id = Column(ForeignKey("project.id"), nullable=False)
    project_status = Column(ForeignKey("project.status"), nullable=False)
    link_url = Column(String)
    link_type = Column(Enum(ProjectLinkType))
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)

    project = relationship("DbProject")


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
    organism = Column(String)
    organism_ontology = Column(String)
    tissue = Column(String)
    tissue_ontology = Column(String)
    assay = Column(String)
    assay_ontology = Column(String)
    disease = Column(String)
    disease_ontology = Column(String)
    sex = Column(String)
    ethnicity = Column(String)
    ethnicity_ontology = Column(String)
    source_data_location = Column(String)
    preprint_doi = Column(String)
    publication_doi = Column(String)
    created_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)


class DbDatasetArtifact(Base):
    """
    Represents a DbUser uploaded or Corpora generated file linked to a DbDataset.
    All matrices and cellxgene objects are examples of a DbDatasetArtifact.
    """
    __tablename__ = "dataset_artifact"

    id = Column(String, primary_key=True)
    dataset_id = Column(ForeignKey("dataset.id"), nullable=False)
    filename = Column(String)
    filetype = Column(String)
    type = Column(Enum(DatasetArtifactType))
    s3_uri = Column(String)
    created_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)

    dataset = relationship("DbDataset")


class DbDeploymentDirectory(Base):
    """
    Represents the deployment of a dataset to a Corpora application.
    This entity only supports cellxgene deployments.
    """
    __tablename__ = "deployment_directory"

    id = Column(String, primary_key=True)
    dataset_id = Column(ForeignKey("dataset.id"), nullable=False)
    environment = Column(String)
    url = Column(String)
    created_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)

    dataset = relationship("DbDataset")


class DbDatasetContributor(Base):
    """
    Associates a DbDataset with a DbContributor.
    DbDatasets may have many DbContributors.
    DbContributors may have many DbDatasets.
    """
    __tablename__ = "dataset_contributor"

    id = Column(String, primary_key=True)
    contributor_id = Column(ForeignKey("contributor.id"), nullable=False)
    dataset_id = Column(ForeignKey("dataset.id"), nullable=False)
    created_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)


class DbContributor(Base):
    """
    A data contributor. Typically a researcher associated with an institution.
    """
    __tablename__ = "contributor"

    id = Column(String, primary_key=True)
    name = Column(String)
    institution = Column(String)
    email = Column(String)
    created_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)
