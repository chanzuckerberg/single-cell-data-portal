"""update_for_collections

Revision ID: f6c4f66e8db6
Revises: 9023d0bf61ab
Create Date: 2020-10-27 16:14:31.084060

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = "f6c4f66e8db6"
down_revision = "9023d0bf61ab"
branch_labels = None
depends_on = None

# Enum 'type' for PostgreSQL
enum_projectlink = "projectlink"
# Set temporary enum 'type' for PostgreSQL
tmp_enum_projectlink = "tmp_" + enum_projectlink

# Options for Enum
old_options = ("RAW_DATA", "PROTOCOL", "SUMMARY", "OTHER")
new_options = sorted(["DOI", "RAW_DATA", "PROTOCOL", "LAB_WEBSITE", "OTHER"])

# Create enum fields
old_type = sa.Enum(*old_options, name=enum_projectlink)
new_type = sa.Enum(*new_options, name=enum_projectlink)


def add_nullable_false(table: str, column: str, type, default):
    op.add_column(table, sa.Column(column, type, nullable=True))
    op.execute(f"UPDATE {table} SET {column} = {default}")
    op.alter_column(table, column, nullable=False)


def add_enum(name: str, values: list):
    values.sort()
    _enum = postgresql.ENUM(*values, name=name.lower())
    _enum.create(op.get_bind())


def upgrade():
    # For more information about the changes executed here, please visit the design document and the approved
    # modifications to the database UML:
    # https://github.com/chanzuckerberg/single-cell/blob/main/rfcs/0001-data-portal-architecture/text.md

    def _add_visibility(table):
        prefix1 = "" if table == "project" else "collection_"
        prefix2 = "" if table == "project" else "project_"
        op.add_column(
            table,
            sa.Column(f"{prefix1}visibility", sa.Enum("PUBLIC", "PRIVATE", name="collectionvisibility"), nullable=True),
        )
        op.execute(f"UPDATE {table} SET {prefix1}visibility = 'PUBLIC' WHERE {prefix2}status = 'LIVE'")
        op.execute(f"UPDATE {table} SET {prefix1}visibility = 'PRIVATE' WHERE {prefix2}status = 'EDIT'")
        op.alter_column(table, f"{prefix1}visibility", nullable=False)
        if table != "project":
            op.alter_column(table, f"{prefix2}id", nullable=False, new_column_name=f"{prefix1}id")
        op.drop_column(table, f"{prefix2}status")

    # Create Visibility Enum
    add_enum("collectionvisibility", ["PUBLIC", "PRIVATE"])

    # Remove (id, status) Primary and Foreign Key Constraint
    op.execute("ALTER TABLE project DROP CONSTRAINT project_pkey CASCADE")

    # Update Project Table
    _add_visibility("project")
    op.add_column("project", sa.Column("contact_email", sa.String(), nullable=True))
    op.add_column("project", sa.Column("contact_name", sa.String(), nullable=True))
    add_nullable_false("project", "data_submission_policy_version", sa.String(), default="0")
    op.add_column("project", sa.Column("obfuscated_uuid", sa.String(), nullable=True))
    op.create_primary_key("project_pkey", "project", ["id", "visibility"])
    op.drop_constraint("project_owner_fkey", "project", type_="foreignkey")
    op.drop_column("project", "processing_state")
    op.drop_column("project", "s3_bucket")
    op.drop_column("project", "needs_attestation")
    op.drop_column("project", "tc_uri")
    op.drop_column("project", "validation_state")

    # Update Project Link Table
    _add_visibility("project_link")
    op.create_foreign_key(
        None, "project_link", "project", ["collection_id", "collection_visibility"], ["id", "visibility"]
    )
    # Rename current enum type to tmp_
    op.execute("DELETE FROM project_link WHERE project_link.link_type = 'SUMMARY'")
    op.execute("ALTER TYPE " + enum_projectlink + " RENAME TO " + tmp_enum_projectlink)
    # Create new enum type in db
    new_type.create(op.get_bind())
    # Update column to use new enum type
    op.execute(
        "ALTER TABLE project_link ALTER COLUMN link_type TYPE "
        + enum_projectlink
        + " USING link_type::text::"
        + enum_projectlink
    )
    # Drop old enum type
    op.execute("DROP TYPE " + tmp_enum_projectlink)

    # Update Dataset Table
    _add_visibility("dataset")
    op.add_column("dataset", sa.Column("cell_count", sa.Integer(), nullable=True))
    op.add_column("dataset", sa.Column("is_valid", sa.Boolean(), nullable=True))
    op.create_foreign_key(None, "dataset", "project", ["collection_id", "collection_visibility"], ["id", "visibility"])
    op.drop_column("dataset", "publication_doi")
    op.drop_column("dataset", "source_data_location")
    op.drop_column("dataset", "preprint_doi")

    # Update Deployment_Directory Table
    op.drop_column("deployment_directory", "environment")

    # Remove Tables
    op.drop_table("dataset_contributor")
    op.drop_table("contributor")
    op.drop_table("user")

    # Drop Unused Enums
    sa.Enum(name="projectstatus").drop(op.get_bind(), checkfirst=False)
    sa.Enum(name="validationstate").drop(op.get_bind(), checkfirst=False)
    sa.Enum(name="processingstate").drop(op.get_bind(), checkfirst=False)
    # ### end Alembic commands ###


def downgrade():
    def _add_status(table):
        prefix1 = "" if table == "project" else "project_"
        prefix2 = "" if table == "project" else "collection_"
        op.add_column(
            table, sa.Column(f"{prefix1}status", sa.Enum("PUBLIC", "PRIVATE", name="projectstatus"), nullable=True)
        )
        op.execute(f"UPDATE {table} SET {prefix1}status = 'LIVE' WHERE {prefix2}visibility = 'PUBLIC'")
        op.execute(f"UPDATE {table} SET {prefix1}status = 'EDIT' WHERE {prefix2}visibility = 'PRIVATE'")
        op.alter_column(table, f"{prefix1}status", nullable=False)
        if table != "project":
            op.alter_column(table, f"{prefix2}id", nullable=False, new_column_name=f"{prefix1}id")
        op.drop_column(table, f"{prefix2}visibility")

    # Create Status Enum
    add_enum("processingState", ["NA", "IN_VALIDATION", "IN_ARTIFACT_CREATION", "IN_DEPLOYMENT"])
    add_enum("projectstatus", ["LIVE", "EDIT"])
    add_enum("validationstate", ["NOT_VALIDATED", "VALID", "INVALID"])

    op.create_table(
        "user",
        sa.Column("id", sa.VARCHAR(), autoincrement=False, nullable=False),
        sa.Column("name", sa.VARCHAR(), autoincrement=False, nullable=True),
        sa.Column("email", sa.VARCHAR(), autoincrement=False, nullable=True),
        sa.Column(
            "created_at", postgresql.TIMESTAMP(), server_default=sa.text("now()"), autoincrement=False, nullable=False
        ),
        sa.Column(
            "updated_at", postgresql.TIMESTAMP(), server_default=sa.text("now()"), autoincrement=False, nullable=False
        ),
        sa.PrimaryKeyConstraint("id", name="user_pkey"),
    )
    op.execute("""INSERT INTO public.user (id) SELECT DISTINCT owner FROM project""")

    # Create Contributors Table
    op.create_table(
        "contributor",
        sa.Column("id", sa.VARCHAR(), autoincrement=False, nullable=False),
        sa.Column("name", sa.VARCHAR(), autoincrement=False, nullable=True),
        sa.Column("institution", sa.VARCHAR(), autoincrement=False, nullable=True),
        sa.Column("email", sa.VARCHAR(), autoincrement=False, nullable=True),
        sa.Column(
            "created_at", postgresql.TIMESTAMP(), server_default=sa.text("now()"), autoincrement=False, nullable=False
        ),
        sa.Column(
            "updated_at", postgresql.TIMESTAMP(), server_default=sa.text("now()"), autoincrement=False, nullable=False
        ),
        sa.PrimaryKeyConstraint("id", name="contributor_pkey"),
        postgresql_ignore_search_path=False,
    )

    # Create Contributor Dataset Table
    op.create_table(
        "dataset_contributor",
        sa.Column("id", sa.VARCHAR(), autoincrement=False, nullable=False),
        sa.Column("contributor_id", sa.VARCHAR(), autoincrement=False, nullable=False),
        sa.Column("dataset_id", sa.VARCHAR(), autoincrement=False, nullable=False),
        sa.Column(
            "created_at", postgresql.TIMESTAMP(), server_default=sa.text("now()"), autoincrement=False, nullable=False
        ),
        sa.Column(
            "updated_at", postgresql.TIMESTAMP(), server_default=sa.text("now()"), autoincrement=False, nullable=False
        ),
        sa.ForeignKeyConstraint(["contributor_id"], ["contributor.id"], name="dataset_contributor_contributor_id_fkey"),
        sa.ForeignKeyConstraint(["dataset_id"], ["dataset.id"], name="dataset_contributor_dataset_id_fkey"),
        sa.PrimaryKeyConstraint("id", name="dataset_contributor_pkey"),
    )

    # Remove (id, status) Primary and Foreign Key Constraint
    op.execute("ALTER TABLE project DROP CONSTRAINT project_pkey CASCADE")

    # Update Project Table
    _add_status("project")
    op.add_column(
        "project",
        sa.Column(
            "validation_state",
            postgresql.ENUM("NOT_VALIDATED", "VALID", "INVALID", name="validationstate"),
            autoincrement=False,
            nullable=True,
        ),
    )

    op.add_column(
        "project",
        sa.Column(
            "processing_state",
            postgresql.ENUM("NA", "IN_VALIDATION", "IN_ARTIFACT_CREATION", "IN_DEPLOYMENT", name="processingstate"),
            autoincrement=False,
            nullable=True,
        ),
    )
    op.add_column("project", sa.Column("tc_uri", sa.VARCHAR(), autoincrement=False, nullable=True))
    op.add_column("project", sa.Column("needs_attestation", sa.BOOLEAN(), autoincrement=False, nullable=True))
    op.add_column("project", sa.Column("s3_bucket", sa.VARCHAR(), autoincrement=False, nullable=True))
    op.create_primary_key("project_pkey", "project", ["id", "status"])
    op.create_foreign_key("project_owner_fkey", "project", "user", ["owner"], ["id"])
    op.drop_column("project", "obfuscated_uuid")
    op.drop_column("project", "data_submission_policy_version")
    op.drop_column("project", "contact_name")
    op.drop_column("project", "contact_email")

    # Update Project Link Table
    _add_status("project_link")
    op.create_foreign_key(None, "project_link", "project", ["project_id", "project_status"], ["id", "status"])
    op.execute("DELETE FROM project_link WHERE project_link.link_type = 'DOI'")
    op.execute("DELETE FROM project_link WHERE project_link.link_type = 'LAB_WEBSITE'")
    op.execute("ALTER TYPE " + enum_projectlink + " RENAME TO " + tmp_enum_projectlink)
    # Create new enum type in db
    old_type.create(op.get_bind())
    # Update column to use new enum type
    op.execute(
        "ALTER TABLE project_link ALTER COLUMN link_type TYPE "
        + enum_projectlink
        + " USING link_type::text::"
        + enum_projectlink
    )
    # Drop old enum type
    op.execute("DROP TYPE " + tmp_enum_projectlink)

    # Update Deployment Directory Table
    op.add_column("deployment_directory", sa.Column("environment", sa.VARCHAR(), autoincrement=False, nullable=True))

    # Update Dataset Table
    _add_status("dataset")
    op.add_column("dataset", sa.Column("preprint_doi", sa.VARCHAR(), autoincrement=False, nullable=True))
    op.add_column("dataset", sa.Column("source_data_location", sa.VARCHAR(), autoincrement=False, nullable=True))
    op.add_column("dataset", sa.Column("publication_doi", sa.VARCHAR(), autoincrement=False, nullable=True))
    op.create_foreign_key("fk_project", "dataset", "project", ["project_id", "project_status"], ["id", "status"])
    op.drop_column("dataset", "is_valid")
    op.drop_column("dataset", "cell_count")

    # Drop collection visibility Enum
    sa.Enum(name="collectionvisibility").drop(op.get_bind(), checkfirst=False)
    # ### end Alembic commands ###
