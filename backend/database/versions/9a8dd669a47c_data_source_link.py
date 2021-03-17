"""data_source_link

Revision ID: 9a8dd669a47c
Revises: bf83d170bb7f
Create Date: 2021-03-08 16:34:52.801864

"""
import sqlalchemy as sa

# revision identifiers, used by Alembic.
from alembic import op

revision = "9a8dd669a47c"
down_revision = "bf83d170bb7f"
branch_labels = None
depends_on = None

# Enum 'type' for PostgreSQL
enum_projectlink = "projectlink"
# Set temporary enum 'type' for PostgreSQL
tmp_enum_projectlink = "tmp_" + enum_projectlink

# Options for Enum
old_options = sorted(["DOI", "RAW_DATA", "PROTOCOL", "LAB_WEBSITE", "OTHER"])
new_options = sorted(["DOI", "RAW_DATA", "PROTOCOL", "LAB_WEBSITE", "OTHER", "DATA_SOURCE"])

# Create enum fields
old_type = sa.Enum(*old_options, name=enum_projectlink)
new_type = sa.Enum(*new_options, name=enum_projectlink)


def upgrade():
    # Rename current enum type to tmp_
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


def downgrade():
    op.execute("DELETE FROM project_link WHERE project_link.link_type = 'DATA_SOURCE'")
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
