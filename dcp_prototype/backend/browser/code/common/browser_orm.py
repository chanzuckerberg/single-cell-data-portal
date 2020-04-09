import os
import sys

from sqlalchemy import (
    create_engine,
    Boolean,
    Column,
    DateTime,
    ForeignKey,
    Integer,
    String,
    text,
    BigInteger,
)
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.ext.declarative import declarative_base

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from browser.code.config.db_config import BrowserDbConfig


Base = declarative_base()
deployment_stage = os.environ["DEPLOYMENT_STAGE"]

# The in-memory SQLite database used in
# unit tests does not support now() SQL syntax
DEFAULT_DATETIME = "2000-01-01 00:00:00" if deployment_stage == "test" else text("now()")


class DBSessionMaker:
    def __init__(self):
        connection = "sqlite:///:memory:" if deployment_stage == "test" else BrowserDbConfig().database_uri
        self.engine = create_engine(connection)
        self.session_maker = sessionmaker(bind=self.engine)

    def session(self, **kwargs):
        return self.session_maker(**kwargs)


class Project(Base):
    __tablename__ = "project"

    id = Column(String(64), primary_key=True)
    title = Column(String(255))
    label = Column(String(255))
    description = Column(String(3000))
    biosample_categories = Column(String(32))
    development_stages = Column(String(150))
    diseases = Column(String(1000))
    cell_isolation_methods = Column(String(100))
    cell_types = Column(String(150))
    cell_count = Column(Integer)
    paired_end = Column(String(16))
    nucleic_acid_sources = Column(String(64))
    input_nucleic_acid_molecules = Column(String(100))
    publication_title = Column(String(200))
    publication_doi = Column(String(32))
    cxg_enabled = Column(Boolean)
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)


class File(Base):
    __tablename__ = "file"

    id = Column(String(64), primary_key=True)
    project_id = Column(ForeignKey("project.id"), nullable=False)
    filename = Column(String(80))
    file_format = Column(String(8))
    file_size = Column(BigInteger)
    file_type = Column(String(20))
    s3_uri = Column(String(200))
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)


class LibraryConstructionMethod(Base):
    __tablename__ = "library_construction_method"

    id = Column(Integer, primary_key=True)
    name = Column(String(64))
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)


class Organ(Base):
    __tablename__ = "organ"

    id = Column(Integer, primary_key=True)
    name = Column(String(32))
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)


class Species(Base):
    __tablename__ = "species"

    id = Column(Integer, primary_key=True)
    name = Column(String(16))
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)


class DataRepository(Base):
    __tablename__ = "data_repository"

    id = Column(Integer, primary_key=True)
    name = Column(String(32))
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)


class Contributor(Base):
    __tablename__ = "contributor"

    id = Column(Integer, primary_key=True)
    first_name = Column(String(16))
    middle_name = Column(String(16))
    last_name = Column(String(32))
    institution = Column(String(100))
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)


class ExternalAccession(Base):
    __tablename__ = "external_accession"

    id = Column(Integer, primary_key=True)
    project_id = Column(ForeignKey("project.id"), nullable=False)
    data_repository_id = Column(ForeignKey("data_repository.id"), nullable=False)
    accession = Column(String(32))
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)

    project = relationship("Project")
    data_repository = relationship("DataRepository")


class LibraryConstructionMethodJoinProject(Base):
    __tablename__ = "library_construction_method_join_project"

    id = Column(Integer, primary_key=True)
    library_construction_method_id = Column(ForeignKey("library_construction_method.id"), nullable=False)
    project_id = Column(ForeignKey("project.id"), nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)

    library_construction_method = relationship("LibraryConstructionMethod")
    project = relationship("Project")


class OrganJoinProject(Base):
    __tablename__ = "organ_join_project"

    id = Column(Integer, primary_key=True)
    organ_id = Column(ForeignKey("organ.id"), nullable=False)
    project_id = Column(ForeignKey("project.id"), nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)

    organ = relationship("Organ")
    project = relationship("Project")


class SpeciesJoinProject(Base):
    __tablename__ = "species_join_project"

    id = Column(Integer, primary_key=True)
    species_id = Column(ForeignKey("species.id"), nullable=False)
    project_id = Column(ForeignKey("project.id"), nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)

    species = relationship("Species")
    project = relationship("Project")


class ContributorJoinProject(Base):
    __tablename__ = "contributor_join_project"

    id = Column(Integer, primary_key=True)
    contributor_id = Column(ForeignKey("contributor.id"), nullable=False)
    project_id = Column(ForeignKey("project.id"), nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime(True), nullable=False, server_default=DEFAULT_DATETIME)

    contributor = relationship("Contributor")
    project = relationship("Project")
