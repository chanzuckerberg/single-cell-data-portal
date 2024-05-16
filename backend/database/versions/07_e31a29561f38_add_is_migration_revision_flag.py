"""add is_auto_version flag

Revision ID: 07_e31a29561f38
Revises: 06_8bced1b1470b
Create Date: 2024-04-23 10:32:09.778855

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "07_e31a29561f38"
down_revision = "06_8bced1b1470b"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column(
        "CollectionVersion",
        sa.Column("is_auto_version", sa.BOOLEAN(), nullable=True),
        schema="persistence_schema",
    )


def downgrade():
    op.drop_column("CollectionVersion", "is_auto_version", schema="persistence_schema")
