# coding: utf-8
import os
import sys

from sqlalchemy import create_engine, Column, DateTime, ForeignKey, String, text
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
    __tablename__ = "alignment_protocol"

    id = Column(String, primary_key=True)
    software = Column(String)
    algorithm = Column(String)
    genome_reference = Column(String)
    genomic_annotation = Column(String)
    genomic_annotation_biotypes = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class AnalysisFile(Base):
    __tablename__ = "analysis_file"

    id = Column(String, primary_key=True)
    filename = Column(String)
    file_format = Column(String)
    file_output_type = Column(String)
    s3_uri = Column(String, nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class BiosamplePrep(Base):
    __tablename__ = "biosample_prep"

    id = Column(String, primary_key=True)
    category = Column(String)
    organ_ontology = Column(String)
    developmental_stage = Column(String)
    disease_ontology = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class Contributor(Base):
    __tablename__ = "contributor"

    id = Column(String, primary_key=True)
    name = Column(String)
    email = Column(String)
    phone_number = Column(String)
    corresponding_contributor = Column(String)
    lab = Column(String)
    street_address = Column(String)
    country = Column(String)
    contributor_role_ontology = Column(String)
    orcid_id = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class ExpressionFile(Base):
    __tablename__ = "expression_file"

    id = Column(String, primary_key=True)
    filename = Column(String)
    file_format = Column(String)
    file_size = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class LibraryPrepProtocol(Base):
    __tablename__ = "library_prep_protocol"

    id = Column(String, primary_key=True)
    input_nucleic_acid_molecule = Column(String)
    library_construction_method_ontology = Column(String)
    nucleic_acid_source = Column(String)
    end_bias = Column(String)
    barcoded_read = Column(String)
    barcoded_offset = Column(String)
    barcoded_length = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class Project(Base):
    __tablename__ = "project"

    id = Column(String, primary_key=True)
    project_title = Column(String)
    publication_title = Column(String)
    publication_doi = Column(String)
    external_accessions = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class QuantificationProtocol(Base):
    __tablename__ = "quantification_protocol"

    id = Column(String, primary_key=True)
    quantification_software = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class SequenceFile(Base):
    __tablename__ = "sequence_file"

    id = Column(String, primary_key=True)
    filename = Column(String)
    file_format = Column(String)
    file_size = Column(String)
    flowcell_id = Column(String)
    lane_index = Column(String)
    read_index = Column(String)
    s3_uri = Column(String, nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class SequencingProtocol(Base):
    __tablename__ = "sequencing_protocol"

    id = Column(String, primary_key=True)
    paired_end_sequencing = Column(String)
    instrument_manufacturer_model = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))


class AnalysisFileQuantificationProtocolExpressionFileProcessJ(Base):
    __tablename__ = "analysis_file_quantification_protocol_expression_file_process_j"

    id = Column(String, primary_key=True)
    analysis_file_id = Column(ForeignKey("analysis_file.id"), nullable=False)
    quantification_protocol_id = Column(
        ForeignKey("quantification_protocol.id"), nullable=False
    )
    expression_file_id = Column(ForeignKey("expression_file.id"), nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))

    analysis_file = relationship("AnalysisFile")
    expression_file = relationship("ExpressionFile")
    quantification_protocol = relationship("QuantificationProtocol")


class BamQcMetric(Base):
    __tablename__ = "bam_qc_metric"

    id = Column(String, primary_key=True)
    analysis_file_id = Column(ForeignKey("analysis_file.id"), nullable=False)
    percent_mapped_to_genome = Column(String)
    percent_reads_mapped_to_intergenic = Column(String)
    percent_reads_mapped_uniquely = Column(String)
    percent_reads_mapped_multiple = Column(String)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))

    analysis_file = relationship("AnalysisFile")


class Library(Base):
    __tablename__ = "library"

    id = Column(String, primary_key=True)
    project_id = Column(ForeignKey("project.id"), nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))

    project = relationship("Project")


class ProjectContributorJoin(Base):
    __tablename__ = "project_contributor_join"

    id = Column(String, primary_key=True)
    contributor_id = Column(ForeignKey("contributor.id"), nullable=False)
    project_id = Column(ForeignKey("project.id"), nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))

    contributor = relationship("Contributor")
    project = relationship("Project")


class SequenceFileAlignmentProtocolAnalysisFileProcessJoin(Base):
    __tablename__ = "sequence_file_alignment_protocol_analysis_file_process_join"

    id = Column(String, primary_key=True)
    sequence_file_id = Column(ForeignKey("sequence_file.id"), nullable=False)
    analysis_file_id = Column(ForeignKey("analysis_file.id"), nullable=False)
    alignment_protocol_id = Column(ForeignKey("alignment_protocol.id"), nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))

    alignment_protocol = relationship("AlignmentProtocol")
    analysis_file = relationship("AnalysisFile")
    sequence_file = relationship("SequenceFile")


class BiosamplePrepLibraryLibraryPrepProtocolProcessJoin(Base):
    __tablename__ = "biosample_prep_library_library_prep_protocol_process_join"

    id = Column(String, primary_key=True)
    biosample_prep_id = Column(ForeignKey("biosample_prep.id"), nullable=False)
    library_prep_protocol_id = Column(
        ForeignKey("library_prep_protocol.id"), nullable=False
    )
    library_id = Column(ForeignKey("library.id"), nullable=False)
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))

    biosample_prep = relationship("BiosamplePrep")
    library = relationship("Library")
    library_prep_protocol = relationship("LibraryPrepProtocol")


class LibrarySequenceFileSequencingProtocolProcessJoin(Base):
    __tablename__ = "library_sequence_file_sequencing_protocol_process_join"

    id = Column(String, primary_key=True)
    library_id = Column(ForeignKey("library.id"), nullable=False)
    sequence_file_id = Column(ForeignKey("sequence_file.id"), nullable=False)
    sequencing_protocol_id = Column(
        ForeignKey("sequencing_protocol.id"), nullable=False
    )
    created_at = Column(DateTime(True), nullable=False, server_default=text("now()"))
    updated_at = Column(DateTime(True), nullable=False, server_default=text("now()"))

    library = relationship("Library")
    sequence_file = relationship("SequenceFile")
    sequencing_protocol = relationship("SequencingProtocol")
