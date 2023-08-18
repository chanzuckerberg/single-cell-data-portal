"""add_collection_version_schema_field

Revision ID: 04_2f30f3bcc9aa
Revises: 03_98cf7564214f
Create Date: 2023-07-24 10:09:10.690501

"""
import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "04_2f30f3bcc9aa"
down_revision = "03_98cf7564214f"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column(
        "CollectionVersion", sa.Column("schema_version", sa.String(), nullable=True), schema="persistence_schema"
    )


def downgrade():
    op.drop_column("CollectionVersion", "schema_version", schema="persistence_schema")
