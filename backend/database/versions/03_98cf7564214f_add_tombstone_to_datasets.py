"""add_tombstone_to_datasets

Revision ID: 03_98cf7564214f
Revises: 02_0a5021e09eff
Create Date: 2023-06-02 18:25:16.188088

"""
import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "03_98cf7564214f"
down_revision = "02_0a5021e09eff"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("Dataset", sa.Column("tombstone", sa.BOOLEAN(), nullable=True), schema="persistence_schema")


def downgrade():
    op.drop_column("Dataset", "tombstone", schema="persistence_schema")
