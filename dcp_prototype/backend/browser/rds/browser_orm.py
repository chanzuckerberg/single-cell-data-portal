import os
import sys

from sqlalchemy import create_engine, Table, Column, ForeignKey, Integer, String, DateTime
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.ext.declarative import declarative_base

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from browser.config.db_config import BrowserDbConfig


Base = declarative_base()


class DBSessionMaker:
    def __init__(self):
        engine = create_engine(BrowserDbConfig().database_uri)
        Base.metadata.bind = engine
        self.session_maker = sessionmaker()
        self.session_maker.bind = engine

    def session(self, **kwargs):
        return self.session_maker(**kwargs)


class Project(Base):
    __tablename__ = "project"

    id = Column(String(64), primary_key=True)
    title = Column(String(255))
    label = Column(String(100))
    description = Column(String(3000))
    category = Column(String(32))
    developmental_stage = Column(String(32))
    disease_ontology = Column(String(16))
    sample_type = Column(String(32))
    organ_part = Column(String(32))
    analysis_protocol = Column(String(64))
    cell_count = Column(Integer)
    donor_count = Column(Integer)
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
    file_size = Column(Integer)
    flowcell_id = Column(String(64))
    lane_index = Column(String(8))
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
