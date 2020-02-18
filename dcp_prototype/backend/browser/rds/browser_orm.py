import os
import sys

from sqlalchemy import create_engine, Table, Column, ForeignKey, Integer, String, DateTime
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.ext.declarative import declarative_base

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from browser.rds.db_config import BrowserDbConfig


stage = os.environ['DEPLOYMENT_STAGE']
db_name = f"browser_{stage}"

Base = declarative_base()
engine = create_engine(BrowserDbConfig().database_uri)
conn = engine.connect()
engine.execute(f"USE {db_name}")


class DBSessionMaker:
    def __init__(self):
        Base.metadata.bind = engine
        self.session_maker = sessionmaker()
        self.session_maker.bind = engine

    def session(self, **kwargs):
        return self.session_maker(**kwargs)


class Project(Base):
    __tablename__ = "project"

    id = Column(String(64), primary_key=True)
    title = Column(String(255))
    label = Column(String(100))  # NA
    description = Column(String(3000))  # NA
    category = Column(String(32))
    developmental_stage = Column(String(32))
    disease_ontology = Column(String(16))  # multiple
    sample_type = Column(String(32))  # NA
    organ_part = Column(String(32))  # NA
    analysis_protocol = Column(String(64))  # NA
    cell_count = Column(Integer)  # NA
    donor_count = Column(Integer)  # NA
    publication_title = Column(String(128))
    publication_doi = Column(String(32))
    contact_name = Column(String(32))
    contact_institution = Column(String(64))
    contact_email = Column(String(32))


class File(Base):
    __tablename__ = "file"

    id = Column(String(64), primary_key=True)
    project_id = Column(ForeignKey("project.id"), nullable=False)
    filename = Column(String(80))
    file_format = Column(String(8))
    file_size = Column(Integer)  # NA
    flowcell_id = Column(String(64))  # NA
    lane_index = Column(String(8))  # NA
    read_index = Column(String(8))
    s3_uri = Column(String(200))
    created_at = Column(DateTime(True), nullable=False)
    updated_at = Column(DateTime(True), nullable=False)
    species = Column(String(16))
    library_construction_method_ontology = Column(String(16))
    tissue_ontology = Column(String(16))
    paired_end_sequencing = Column(String(16))
    instrument_manufacturer_model_ontology = Column(String(32))
    file_output_type = Column(String(32))
    software = Column(String(32))
    algorithm = Column(String(32))
    genome_reference = Column(String(32))
    genomic_annotation = Column(String(32))
    genomic_annotation_biotypes = Column(String(32))
    quantification_software = Column(String(32))


class LibraryPrepProtocol(Base):
    __tablename__ = "library_prep_protocol"

    id = Column(Integer, primary_key=True)
    construction_method_ontology = Column(String(16))
    end_bias = Column(String(16))
    nucleic_acid_source = Column(String(16))


class Tissue(Base):
    __tablename__ = "tissue"

    id = Column(Integer, primary_key=True)
    tissue_ontology = Column(String(16))


class Species(Base):
    __tablename__ = "species"

    id = Column(Integer, primary_key=True)
    species_ontology = Column(String(16))


class DataRepository(Base):
    __tablename__ = "data_repository"

    id = Column(Integer, primary_key=True)
    name = Column(String(32))


class Contributor(Base):
    __tablename__ = "contributor"

    id = Column(Integer, primary_key=True)
    name = Column(String(32))
    institution = Column(String(32))
    email = Column(String(32))


class ExternalAccession(Base):
    __tablename__ = "external_accession"

    id = Column(Integer, primary_key=True)
    project_id = Column(ForeignKey("project.id"), nullable=False)
    data_repository_id = Column(ForeignKey("data_repository.id"), nullable=False)
    accession = Column(String(32))

    project = relationship("Project")
    data_repository = relationship("DataRepository")


class LibraryPrepProtocolJoinProject(Base):
    __tablename__ = "library_prep_protocol_join_project"

    id = Column(Integer, primary_key=True)
    library_prep_protocol_id = Column(ForeignKey("library_prep_protocol.id"), nullable=False)
    project_id = Column(ForeignKey("project.id"), nullable=False)

    library_prep_protocol = relationship("LibraryPrepProtocol")
    project = relationship("Project")


class TissueJoinProject(Base):
    __tablename__ = "tissue_join_project"

    id = Column(Integer, primary_key=True)
    tissue_id = Column(ForeignKey("tissue.id"), nullable=False)
    project_id = Column(ForeignKey("project.id"), nullable=False)

    tissue = relationship("Tissue")
    project = relationship("Project")


class SpeciesJoinProject(Base):
    __tablename__ = "species_join_project"

    id = Column(Integer, primary_key=True)
    species_id = Column(ForeignKey("species.id"), nullable=False)
    project_id = Column(ForeignKey("project.id"), nullable=False)

    species = relationship("Species")
    project = relationship("Project")


class ContributorJoinProject(Base):
    __tablename__ = "contributor_join_project"

    id = Column(Integer, primary_key=True)
    contributor_id = Column(ForeignKey("contributor.id"), nullable=False)
    project_id = Column(ForeignKey("project.id"), nullable=False)

    contributor = relationship("Contributor")
    project = relationship("Project")


# uncomment to drop and recreate all tables
# engine.execute("DROP TABLE tissue_join_project, species_join_project, library_prep_protocol_join_project, contributor_join_project, external_accession, file, library_prep_protocol, project, species, tissue, contributor, data_repository;")
# Base.metadata.create_all(engine)
# tables = engine.execute("SHOW TABLES;")
# tables = [t[0] for t in tables]
# print(tables)
