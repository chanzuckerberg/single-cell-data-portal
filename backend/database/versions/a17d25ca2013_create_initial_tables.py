"""create_initial_tables

Revision ID: a17d25ca2013
Revises:
Create Date: 2020-07-17 21:57:32.177139

"""

from alembic import op
from sqlalchemy import Column, String, text, Boolean, Integer, func, TIMESTAMP
from sqlalchemy.dialects.postgresql import ENUM
from sqlalchemy.types import DateTime

# revision identifiers, used by Alembic.
revision = "a17d25ca2013"
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    # User table
    op.create_table(
        "user",
        Column("id", String, nullable=False, primary_key=True),
        Column("name", String, nullable=True),
        Column("email", String, nullable=True),
        Column("created_at", TIMESTAMP, nullable=False, server_default=func.now()),
        Column("updated_at", TIMESTAMP, nullable=False, server_default=func.now()),
    )

    # Project table
    processing_state_enum = ENUM(
        "NA", "IN_VALIDATION", "IN_ARTIFACT_CREATION", "IN_DEPLOYMENT", name="processingstate", create_type=False
    )
    processing_state_enum.create(op.get_bind(), checkfirst=True)
    validation_state_enum = ENUM(
        "NA", "IN_VALIDATION", "IN_ARTIFACT_CREATION", "IN_DEPLOYMENT", name="validationstate", create_type=False
    )
    validation_state_enum.create(op.get_bind(), checkfirst=True)

    op.create_table(
        "project",
        Column("id", String, nullable=False, primary_key=True),
        Column("status", String, nullable=False, primary_key=True),
        Column("owner", String, nullable=False),
        Column("name", String, nullable=True),
        Column("description", String, nullable=True),
        Column("s3_bucket", String, nullable=True),
        Column("tc_uri", String, nullable=True),
        Column("needs_attestation", Boolean, nullable=False),
        Column("processing_state", processing_state_enum, nullable=False),
        Column("validation_state", validation_state_enum, nullable=False),
        Column("created_at", TIMESTAMP, nullable=False, server_default=func.now()),
        Column("updated_at", TIMESTAMP, nullable=False, server_default=func.now()),
    )

    op.create_foreign_key(
        "fk_project_owner", "project", "user", ["owner"], ["id"],
    )

    # ProjectLink table
    project_link_enum = ENUM("PROTOCOL", "RAW_DATA", "OTHER", name="projectlink", create_type=False)
    project_link_enum.create(op.get_bind(), checkfirst=True)

    op.create_table(
        "project_link",
        Column("id", String, nullable=False, primary_key=True),
        Column("project_id", String, nullable=False),
        Column("project_status", String, nullable=False),
        Column("link_url", String, nullable=True),
        Column("link_type", project_link_enum, nullable=True),
        Column("created_at", TIMESTAMP, nullable=False, server_default=func.now()),
        Column("updated_at", TIMESTAMP, nullable=False, server_default=func.now()),
    )

    op.create_foreign_key(
        "fk_project", "project_link", "project", ["project_id", "project_status"], ["id", "status"],
    )

    # Dataset table
    op.create_table(
        "dataset",
        Column("id", String, nullable=False, primary_key=True),
        Column("revision", Integer, nullable=False),
        Column("name", String, nullable=True),
        Column("organism", String, nullable=True),
        Column("organism_ontology", String, nullable=True),
        Column("tissue", String, nullable=True),
        Column("tissue_ontology", String, nullable=True),
        Column("assay", String, nullable=True),
        Column("assay_ontology", String, nullable=True),
        Column("disease", String, nullable=True),
        Column("disease_ontology", String, nullable=True),
        Column("sex", String, nullable=True),
        Column("ethnicity", String, nullable=True),
        Column("ethnicity_ontology", String, nullable=True),
        Column("source_data_location", String, nullable=True),
        Column("preprint_doi", String, nullable=True),
        Column("publication_doi", String, nullable=True),
        Column("created_at", TIMESTAMP, nullable=False, server_default=func.now()),
        Column("updated_at", TIMESTAMP, nullable=False, server_default=func.now()),
    )

    # DatasetArtifact table
    dataset_artifact_file_type_enum = ENUM(
        "H5AD", "RDS", "LOOM", "CXG", name="datasetartifactfiletype", create_type=False
    )
    dataset_artifact_file_type_enum.create(op.get_bind(), checkfirst=True)
    dataset_artifact_type_enum = ENUM("ORIGINAL", "REMIX", name="datasetartifacttype", create_type=False)
    dataset_artifact_type_enum.create(op.get_bind(), checkfirst=True)

    op.create_table(
        "dataset_artifact",
        Column("id", String, nullable=False, primary_key=True),
        Column("dataset_id", String, nullable=False),
        Column("filename", String, nullable=True),
        Column("filetype", dataset_artifact_file_type_enum, nullable=True),
        Column("type", dataset_artifact_type_enum, nullable=True),
        Column("user_submitted", Boolean, nullable=True),
        Column("s3_uri", String, nullable=True),
        Column("created_at", TIMESTAMP, nullable=False, server_default=func.now()),
        Column("updated_at", TIMESTAMP, nullable=False, server_default=func.now()),
    )

    op.create_foreign_key(
        "fk_dataset_id", "dataset_artifact", "dataset", ["dataset_id"], ["id"],
    )

    # DatasetDeployment table
    op.create_table(
        "dataset_deployment",
        Column("id", String, nullable=False, primary_key=True),
        Column("dataset_id", String, nullable=False),
        Column("environment", String, nullable=True),
        Column("url", String, nullable=True),
        Column("created_at", TIMESTAMP, nullable=False, server_default=func.now()),
        Column("updated_at", TIMESTAMP, nullable=False, server_default=func.now()),
    )

    op.create_foreign_key(
        "fk_dataset_id", "dataset_deployment", "dataset", ["dataset_id"], ["id"],
    )

    # ProjectDataset table
    op.create_table(
        "project_dataset",
        Column("id", String, nullable=False, primary_key=True),
        Column("project_id", String, nullable=False),
        Column("project_status", String, nullable=False),
        Column("dataset_id", String, nullable=False),
        Column("created_at", TIMESTAMP, nullable=False, server_default=func.now()),
        Column("updated_at", TIMESTAMP, nullable=False, server_default=func.now()),
    )

    op.create_foreign_key(
        "fk_project", "project_dataset", "project", ["project_id", "project_status"], ["id", "status"],
    )
    op.create_foreign_key(
        "fk_dataset_id", "project_dataset", "dataset", ["dataset_id"], ["id"],
    )

    # Contributor table
    op.create_table(
        "contributor",
        Column("id", String, nullable=False, primary_key=True),
        Column("name", String, nullable=True),
        Column("institution", String, nullable=True),
        Column("email", String, nullable=True),
        Column("created_at", TIMESTAMP, nullable=False, server_default=func.now()),
        Column("updated_at", TIMESTAMP, nullable=False, server_default=func.now()),
    )

    # DatasetContributor table
    op.create_table(
        "dataset_contributor",
        Column("id", String, nullable=False, primary_key=True),
        Column("contributor_id", String, nullable=False),
        Column("dataset_id", String, nullable=False),
        Column("created_at", TIMESTAMP, nullable=False, server_default=func.now()),
        Column("updated_at", TIMESTAMP, nullable=False, server_default=func.now()),
    )

    op.create_foreign_key(
        "fk_contributor_id", "dataset_contributor", "contributor", ["contributor_id"], ["id"],
    )
    op.create_foreign_key(
        "fk_dataset_id", "dataset_contributor", "dataset", ["dataset_id"], ["id"],
    )


def downgrade():
    op.drop_table("dataset_contributor")
    op.drop_table("contributor")
    op.drop_table("project_dataset")
    op.drop_table("dataset_deployment")
    op.drop_table("dataset_artifact")
    op.drop_table("dataset")
    op.drop_table("project_link")
    op.drop_table("project")
    op.drop_table("user")
