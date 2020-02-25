"""initial_table_setup

Revision ID: 407ab1d9843a
Revises:
Create Date: 2020-02-04 14:23:22.436636

"""
from alembic import op
from sqlalchemy import ARRAY, Column, String
from sqlalchemy import text
from sqlalchemy.types import DateTime
from sqlalchemy.dialects.postgresql import BIGINT, ENUM


# revision identifiers, used by Alembic.
revision = "407ab1d9843a"
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    # Biosample Prep
    op.create_table(
        "biosample_prep",
        Column("id", String, nullable=False, primary_key=True),
        Column("donor_development_stage_at_collection", String, nullable=True),
        Column("donor_species", String, nullable=True),
        Column("organ", String, nullable=True),
        Column("system", String, nullable=True),
        Column("category", String, nullable=True),
        Column("name", String, nullable=True),
        Column("diseases", ARRAY(String), nullable=True),
        Column("donor_diseases", ARRAY(String), nullable=True),
        Column("selected_cell_markers", ARRAY(String), nullable=True),
        Column("selected_cell_types", ARRAY(String), nullable=True),
        Column("cell_isolation_method", String, nullable=True),
        Column("protocols_used", ARRAY(String), nullable=True),
        Column("biosample_summary", String, nullable=True),
        Column("collection_method", String, nullable=True),
        Column("differentiation_method", String, nullable=True),
        Column("dissociation_method", String, nullable=True),
        Column("donor_accession", String, nullable=True),
        Column("donor_ethnicity", String, nullable=True),
        Column("donor_sex", String, nullable=True),
        Column("donor_strain", String, nullable=True),
        Column("enrichment_method", String, nullable=True),
        Column("induction_method", String, nullable=True),
        Column("percent_cell_viability", String, nullable=True),
        Column("reprogramming_factors", ARRAY(String), nullable=True),
        Column("suspension_type", String, nullable=True),
        Column("target_pathway", String, nullable=True),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Contributor
    op.create_table(
        "contributor",
        Column("id", String, nullable=False, primary_key=True),
        Column("name", String, nullable=True),
        Column("institution", String, nullable=True),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Publication
    op.create_table(
        "publication",
        Column("id", String, nullable=False, primary_key=True),
        Column("doi", String, nullable=True),
        Column("pmid", String, nullable=True),
        Column("title", String, nullable=True),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Library
    op.create_table(
        "library",
        Column("id", String, nullable=False, primary_key=True),
        Column("cell_count", String, nullable=True),
        Column("assay_category", String, nullable=True),
        Column("project_id", String, nullable=False),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Library Prep Protocol
    op.create_table(
        "library_prep_protocol",
        Column("id", String, nullable=False, primary_key=True),
        Column("input_nucleic_acid_molecule_ontology", String, nullable=False),
        Column("library_construction_method_ontology", String, nullable=False),
        Column("nucleic_acid_source", String, nullable=False),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Project
    op.create_table(
        "project",
        Column("id", String, nullable=False, primary_key=True),
        Column("title", String, nullable=True),
        Column("description", String, nullable=True),
        Column("array_express_accessions", ARRAY(String), nullable=True),
        Column("biostudies_accessions", ARRAY(String), nullable=True),
        Column("geo_series_accessions", ARRAY(String), nullable=True),
        Column("insdc_project_accessions", ARRAY(String), nullable=True),
        Column("publication_id", String, nullable=False),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # File
    file_type_enum = ENUM("EXPRESSION", "ANALYSIS", "SEQUENCE", name="file_type_enum", create_type=False)
    file_type_enum.create(op.get_bind(), checkfirst=True)
    op.create_table(
        "file",
        Column("id", String, nullable=False, primary_key=True),
        Column("type", file_type_enum, nullable=False),
        Column("filename", String, nullable=False),
        Column("file_format", String, nullable=False),
        Column("file_size", BIGINT, nullable=False),
        Column("s3_uri", String, nullable=False),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Sequencing Protocol
    op.create_table(
        "sequencing_protocol",
        Column("id", String, nullable=False, primary_key=True),
        Column("paired_end", String, nullable=True),
        Column("instrument_manufacturer_model", String, nullable=True),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Project x Contributor Join
    op.create_table(
        "project_contributor_join",
        Column("contributor_id", String, nullable=False),
        Column("project_id", String, nullable=False),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Alignment Protocol
    op.create_table(
        "alignment_protocol",
        Column("id", String, nullable=False, primary_key=True),
        Column("software", String, nullable=True),
        Column("algorithm", String, nullable=True),
        Column("genome_reference", String, nullable=True),
        Column("genomic_annotation", String, nullable=True),
        Column("genomic_annotation_biotypes", ARRAY(String), nullable=True),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Quantification Protocol
    op.create_table(
        "quantification_protocol",
        Column("id", String, nullable=False, primary_key=True),
        Column("quantification_software", String, nullable=True),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Biosample Prep x Library x Library Prep Protocol
    op.create_table(
        "biosample_prep_library_library_prep_protocol_process_join",
        Column("biosample_prep_id", String, nullable=False),
        Column("library_id", String, nullable=False),
        Column("library_prep_protocol_id", String, nullable=False),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Library x Sequence File x Sequencing Protocol
    op.create_table(
        "library_sequence_file_sequencing_protocol_process_join",
        Column("library_id", String, nullable=False),
        Column("sequence_file_id", String, nullable=False),
        Column("sequencing_protocol_id", String, nullable=False),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Sequence File x Analysis File x Alignment Protocol
    op.create_table(
        "sequence_file_analysis_file_alignment_protocol_process_join",
        Column("analysis_file_id", String, nullable=False),
        Column("sequence_file_id", String, nullable=False),
        Column("alignment_protocol_id", String, nullable=False),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )

    # Analysis File x Expression File x Quantification Protocol
    op.create_table(
        "analysis_file_expression_file_quantification_protocol_process_join",
        Column("analysis_file_id", String, nullable=False),
        Column("expression_file_id", String, nullable=False),
        Column("quantification_protocol_id", String, nullable=False),
        Column("created_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
        Column("updated_at", DateTime(timezone=True), nullable=False, server_default=text("now()"),),
    )


def downgrade():
    pass
