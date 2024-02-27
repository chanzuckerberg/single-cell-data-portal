"""add data submission policy version

Revision ID: 06_8bced1b1470b
Revises: 05_67c4015bce73
Create Date: 2024-02-14 15:42:20.549707

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "06_8bced1b1470b"
down_revision = "05_67c4015bce73"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column(
        "CollectionVersion",
        sa.Column("data_submission_policy_version", sa.String(), nullable=True),
        schema="persistence_schema",
    )


def downgrade():
    op.drop_column("CollectionVersion", "data_submission_policy_version", schema="persistence_schema")
