"""updating_project_dataset_relationship

Revision ID: 899656d37baa
Revises: a17d25ca2013
Create Date: 2020-07-13 21:11:12.683976

"""
from alembic import op
from sqlalchemy import Column, String, TIMESTAMP, func, INTEGER

# revision identifiers, used by Alembic.
revision = "899656d37baa"
down_revision = "a17d25ca2013"
branch_labels = None
depends_on = None


def upgrade():
    # For more information about the changes executed here, please visit the design document and the approved
    # modifications to the database UML:
    # https://docs.google.com/document/d/1d8tv2Ub5b3E7Il85adOAUcG8P05N6UBZJ3XbhJSRrFs/edit#heading=h.xbr2m3xe2mwj

    # Add two new columns to the Dataset table to reflect new metadata in Corpora Schema.
    op.add_column("dataset", Column("development_stage", String()))
    op.add_column("dataset", Column("development_stage_ontology", String()))

    # Add a new column to the ProjectLink table to reflect a new piece of metadata about project links
    op.add_column("project_link", Column("link_name", String()))

    # Add a foreign key to the Dataset table to reflect the one and only project that the dataset is in.
    op.add_column("dataset", Column("project_id", String(), nullable=True))
    op.add_column("dataset", Column("project_status", String(), nullable=True))
    op.create_foreign_key(
        "fk_project",
        "dataset",
        "project",
        ["project_id", "project_status"],
        ["id", "status"],
    )

    # Drop the ProjectDataset table
    op.drop_table("project_dataset")


def downgrade():
    # Remove columns from Dataset table
    op.drop_column("dataset", "development_stage")
    op.drop_column("dataset", "development_stage_ontology")

    # Remove column from ProjectLink table
    op.drop_column("project_link", "link_name")

    # Create ProjectDataset join table
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
        "fk_project",
        "project_dataset",
        "project",
        ["project_id", "project_status"],
        ["id", "status"],
    )
    op.create_foreign_key(
        "fk_dataset_id",
        "project_dataset",
        "dataset",
        ["dataset_id"],
        ["id"],
    )
    op.create_table(
        "project_dataset",
        Column("id", INTEGER, primary_key=True),
        Column("project_id", String(), nullable=False),
        Column("project_status", String(), nullable=False),
        Column("dataset_id", String(), nullable=False),
        Column("created_at", TIMESTAMP, server_default=func.now()),
        Column("updated_at", TIMESTAMP, server_default=func.now()),
    )
    op.create_foreign_key("fk_project_id", "project", "project_dataset", ["id"], ["project_id"])
    op.create_foreign_key("fk_project_status", "project", "project_dataset", ["status"], ["project_status"])
    op.create_foreign_key("fk_dataset_id", "dataset", "project_dataset", ["id"], ["dataset_id"])

    # Remove foreign key in Dataset table that references the project ID
    op.drop_column("dataset", "project_id")
    op.drop_column("dataset", "project_status")
