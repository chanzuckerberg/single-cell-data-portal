"""add custom dataset order

Revision ID: 05_67c4015bce73
Revises: 04_2f30f3bcc9aa
Create Date: 2024-01-25 11:42:46.126701

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "05_67c4015bce73"
down_revision = "04_2f30f3bcc9aa"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column(
        "CollectionVersion",
        sa.Column("has_custom_dataset_order", sa.BOOLEAN(), nullable=True),
        schema="persistence_schema",
    )


def downgrade():
    op.drop_column("CollectionVersion", "has_custom_dataset_order", schema="persistence_schema")
