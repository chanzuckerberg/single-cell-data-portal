"""initial_table_setup

Revision ID: 407ab1d9843a
Revises:
Create Date: 2020-02-04 14:23:22.436636

"""
from alembic import op
from sqlalchemy import Column, String
from sqlalchemy import text
from sqlalchemy.types import DateTime


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
        Column("category", String, nullable=True),
        Column("organ_ontology", String, nullable=True),
        Column("developmental_stage", String, nullable=True),
        Column("disease_ontology", String, nullable=True),
        Column(
            "created_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
        Column(
            "updated_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
    )

    # Contributor
    op.create_table(
        "contributor",
        Column("id", String, nullable=False, primary_key=True),
        Column("name", String, nullable=True),
        Column("email", String, nullable=True),
        Column("phone_number", String, nullable=True),
        Column("corresponding_contributor", String, nullable=True),
        Column("lab", String, nullable=True),
        Column("street_address", String, nullable=True),
        Column("country", String, nullable=True),
        Column("contributor_role_ontology", String, nullable=True),
        Column("orcid_id", String, nullable=True),
        Column("project_id", String, nullable=False),
        Column(
            "created_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
        Column(
            "updated_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
    )

    # Library
    op.create_table(
        "library",
        Column("id", String, nullable=False, primary_key=True),
        Column("library_prep_protocol_id", String, nullable=False),
        Column("sequencing_protocol_id", String, nullable=False),
        Column("project_id", String, nullable=False),
        Column(
            "created_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
        Column(
            "updated_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
    )

    # Library Prep Protocol
    op.create_table(
        "library_prep_protocol",
        Column("id", String, nullable=False, primary_key=True),
        Column("input_nucleic_acid_molecule", String, nullable=True),
        Column("library_construction_method_ontology", String, nullable=True),
        Column("nucleic_acid_source", String, nullable=True),
        Column("end_bias", String, nullable=True),
        Column("biosample_prep_id", String, nullable=False),
        Column(
            "created_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
        Column(
            "updated_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
    )

    # Project
    op.create_table(
        "project",
        Column("id", String, nullable=False, primary_key=True),
        Column("project_short_name", String, nullable=True),
        Column("publication_title", String, nullable=True),
        Column("publication_doi", String, nullable=True),
        Column("external_accessions", String, nullable=True),
        Column(
            "created_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
        Column(
            "updated_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
    )

    # Sequence File
    op.create_table(
        "sequence_file",
        Column("id", String, nullable=False, primary_key=True),
        Column("filename", String, nullable=True),
        Column("file_format", String, nullable=True),
        Column("flowcell_id", String, nullable=True),
        Column("lane_index", String, nullable=True),
        Column("read_index", String, nullable=True),
        Column("s3_uri", String, nullable=False),
        Column("sequencing_protocol_id", String, nullable=False),
        Column(
            "created_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
        Column(
            "updated_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
    )

    # Sequencing Protocol
    op.create_table(
        "sequencing_protocol",
        Column("id", String, nullable=False, primary_key=True),
        Column("paired_end_sequencing", String, nullable=True),
        Column("instrument_manufacturer_model", String, nullable=True),
        Column(
            "created_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
        Column(
            "updated_at",
            DateTime(timezone=True),
            nullable=False,
            server_default=text("now()"),
        ),
    )


def downgrade():
    pass
