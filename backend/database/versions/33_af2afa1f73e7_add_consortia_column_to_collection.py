"""add_consortia_column_to_collection

Revision ID: 33_af2afa1f73e7
Revises: 32_c27083d1a76d
Create Date: 2022-11-28 23:20:52.397594

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql


# revision identifiers, used by Alembic.
revision = "33_af2afa1f73e7"
down_revision = "32_c27083d1a76d"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("project", sa.Column("consortia", postgresql.ARRAY(sa.String()), nullable=True))


def downgrade():
    op.drop_column("project", "consortia")
