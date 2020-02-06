import os
import sys
from datetime import datetime

from sqlalchemy import (
    create_engine,
    Column,
    String,
    DateTime,
    ForeignKey
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa
from config.ledger_config import LedgerDbConfig

Base = declarative_base()


class DBSessionMaker:
    def __init__(self):
        engine = create_engine(LedgerDbConfig().database_uri)
        Base.metadata.bind = engine
        self.session_maker = sessionmaker()
        self.session_maker.bind = engine

    def session(self, **kwargs):
        return self.session_maker(**kwargs)


class Library(Base):
    __tablename__ = "library"
    id = Column(String(), primary_key=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(
        DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow
    )
    library_prep_protocol_id = Column(
        String(), ForeignKey("library_prep_protocol.id"), nullable=False
    )
    library_prep_protocol = relationship(
        "LibraryPrepProtocol", back_populates="libraries"
    )
    sequencing_protocol_id = Column(
        String(), ForeignKey("sequencing_protocol.id"), nullable=False
    )
    sequencing_protocol = relationship("SequencingProtocol", back_populates="libraries")
    project_id = Column(String(), ForeignKey("project.id"), nullable=False)
    project = relationship("Project", back_populates="libraries")


class LibraryPrepProtocol(Base):
    __tablename__ = "library_prep_protocol"
    id = Column(String(), primary_key=True)
    input_nucleic_acid_molecule = Column(String(), nullable=True)
    library_construction_method_ontology = Column(String(), nullable=True)
    nucleic_acid_source = Column(String(), nullable=True)
    end_bias = Column(String(), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(
        DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow
    )

    biosample_prep_id = Column(
        String(), ForeignKey("biosample_prep.id"), nullable=False
    )
    biosample_prep = relationship(
        "BiosamplePrep", back_populates="library_prep_protocols"
    )

    libraries = relationship(
        "Library",
        order_by=Library.id,
        back_populates="library_prep_protocol",
        cascade="all, delete, delete-orphan",
    )


class BiosamplePrep(Base):
    __tablename__ = "biosample_prep"
    id = Column(String(), primary_key=True)
    category = Column(String(), nullable=True)
    organ_ontology = Column(String(), nullable=True)
    developmental_stage = Column(String(), nullable=True)
    disease_ontology_label = Column(String(), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(
        DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow
    )

    library_prep_protocols = relationship(
        "LibraryPrepProtocol",
        order_by=LibraryPrepProtocol.id,
        back_populates="biosample_prep",
    )


class Contributor(Base):
    __tablename__ = "contributor"
    id = Column(String(), primary_key=True)
    name = Column(String(), nullable=True)
    email = Column(String(), nullable=True)
    phone_number = Column(String(), nullable=True)
    corresponding_contributor = Column(String(), nullable=True)
    lab = Column(String(), nullable=True)
    street_address = Column(String(), nullable=True)
    country = Column(String(), nullable=True)
    contributor_role_ontology = Column(String(), nullable=True)
    orcid_id = Column(String(), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(
        DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow
    )

    project_id = Column(String(), ForeignKey("project.id"), nullable=False)
    project = relationship("Project", back_populates="contributors")


class Project(Base):
    __tablename__ = "project"
    id = Column(String(), primary_key=True)
    project_short_name = Column(String(), nullable=True)
    publication_title = Column(String(), nullable=True)
    publication_doi = Column(String(), nullable=True)
    external_accessions = Column(String(), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(
        DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow
    )

    contributors = relationship(
        "Contributor",
        order_by=Contributor.id,
        back_populates="project",
        cascade="all, delete, delete-orphan",
    )
    libraries = relationship(
        "Library",
        order_by=Library.id,
        back_populates="project",
        cascade="all, delete, delete-orphan",
    )


class SequenceFile(Base):
    __tablename__ = "sequence_file"
    id = Column(String(), primary_key=True)
    filename = Column(String(), nullable=True)
    file_format = Column(String(), nullable=True)
    flowcell_id = Column(String(), nullable=True)
    lane_index = Column(String(), nullable=True)
    read_index = Column(String(), nullable=True)
    s3_uri = Column(String(), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(
        DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow
    )

    sequencing_protocol_id = Column(
        String(), ForeignKey("sequencing_protocol.id"), nullable=False
    )
    sequencing_protocol = relationship(
        "SequencingProtocol", back_populates="sequence_files"
    )


class SequencingProtocol(Base):
    __tablename__ = "sequencing_protocol"
    id = Column(String(), primary_key=True)
    paired_end_sequencing = Column(String(), nullable=True)
    instrument_manufacturer_model = Column(String(), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(
        DateTime, default=datetime.utcnow, nullable=False, onupdate=datetime.utcnow
    )

    libraries = relationship(
        "Library",
        order_by=Library.id,
        back_populates="sequencing_protocol",
        cascade="all, delete, delete-orphan",
    )
    sequence_files = relationship(
        "SequenceFile",
        order_by=SequenceFile.id,
        back_populates="sequencing_protocol",
        cascade="all, delete, delete-orphan",
    )
