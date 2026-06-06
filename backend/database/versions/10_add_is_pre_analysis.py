"""add is_pre_analysis flag to CollectionVersion

Revision ID: 10_add_is_pre_analysis
Revises: 09_add_artifact_filesize
Create Date: 2026-04-07

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "10_add_is_pre_analysis"
down_revision = "09_add_artifact_filesize"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column(
        "CollectionVersion",
        sa.Column("is_pre_analysis", sa.BOOLEAN(), nullable=True),
        schema="persistence_schema",
    )


def downgrade():
    op.drop_column("CollectionVersion", "is_pre_analysis", schema="persistence_schema")
