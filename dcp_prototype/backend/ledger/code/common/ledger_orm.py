# coding: utf-8
import os
import sys

from sqlalchemy import ARRAY, BigInteger, Column, create_engine, DateTime, Enum, ForeignKey, String, text
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.ext.declarative import declarative_base

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


class AlignmentProtocol(Base):
    __tablename__ = 'alignment_protocol'

    id = Column(String, primary_key=True)
    software = Column(String)
    algorithm = Column(String)
    genome_reference = Column(String)
    genomic_annotation = Column(String)
    genomic_annotation_biotypes = Column(ARRAY(String()))
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class BiosamplePrep(Base):
    __tablename__ = 'biosample_prep'

    id = Column(String, primary_key=True)
    donor_development_stage_at_collection = Column(String)
    donor_species = Column(String)
    organ = Column(String)
    system = Column(String)
    category = Column(String)
    name = Column(String)
    diseases = Column(ARRAY(String()))
    biosample_diseases = Column(ARRAY(String()))
    selected_cell_markers = Column(ARRAY(String()))
    selected_cell_types = Column(ARRAY(String()))
    cell_isolated_method = Column(String)
    protocols_used = Column(ARRAY(String()))
    biosample_summary = Column(String)
    collection_method = Column(String)
    differentiation_method = Column(String)
    dissociation_method = Column(String)
    donor_accession = Column(String)
    donor_ethnicity = Column(String)
    donor_sex = Column(String)
    donor_strain = Column(String)
    enrichment_method = Column(String)
    induction_method = Column(String)
    percent_cell_viability = Column(String)
    reprogramming_factors = Column(ARRAY(String()))
    suspension_type = Column(String)
    target_pathway = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class Contributor(Base):
    __tablename__ = 'contributor'

    id = Column(String, primary_key=True)
    name = Column(String)
    institution = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class File(Base):
    __tablename__ = 'file'

    id = Column(String, primary_key=True)
    type = Column(Enum('EXPRESSION', 'ANALYSIS', 'SEQUENCE', name='file_type_enum'), nullable=False)
    filename = Column(String, nullable=False)
    file_format = Column(String, nullable=False)
    file_size = Column(BigInteger, nullable=False)
    s3_uri = Column(String, nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class LibraryPrepProtocol(Base):
    __tablename__ = 'library_prep_protocol'

    id = Column(String, primary_key=True)
    input_nucleic_acid_molecule_ontology = Column(String, nullable=False)
    library_construction_method_ontology = Column(String, nullable=False)
    nucleic_acid_source = Column(String, nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class Process(Base):
    __tablename__ = 'process'

    id = Column(String, primary_key=True)
    input_id = Column(String, nullable=False)
    output_id = Column(String, nullable=False)
    protocol_id = Column(String, nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class Project(Base):
    __tablename__ = 'project'

    id = Column(String, primary_key=True)
    project_title = Column(String)
    project_description = Column(String)
    array_express_accessions = Column(ARRAY(String()))
    biostudies_accessions = Column(ARRAY(String()))
    geo_series_accessions = Column(ARRAY(String()))
    insdc_project_accessions = Column(ARRAY(String()))
    publication_id = Column(String, nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class Publication(Base):
    __tablename__ = 'publication'

    id = Column(String, primary_key=True)
    doi = Column(String)
    pmid = Column(String)
    title = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class QuantificationProtocol(Base):
    __tablename__ = 'quantification_protocol'

    id = Column(String, primary_key=True)
    quantification_software = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class SequencingProtocol(Base):
    __tablename__ = 'sequencing_protocol'

    id = Column(String, primary_key=True)
    paired_end = Column(String)
    instrument_manufacturer_model = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class Library(Base):
    __tablename__ = 'library'

    id = Column(String, primary_key=True)
    cell_count = Column(String, nullable=False)
    assay_category = Column(String, nullable=False)
    project_id = Column(ForeignKey('project.id'), nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))

    project = relationship('Project')


class ProjectContributorJoin(Base):
    __tablename__ = 'project_contributor_join'

    id = Column(String, primary_key=True)
    contributor_id = Column(ForeignKey('contributor.id'), nullable=False)
    project_id = Column(ForeignKey('project.id'), nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))

    contributor = relationship('Contributor')
    project = relationship('Project')