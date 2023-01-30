"""adding revised_at for collections

Revision ID: 35_0a5021e09eff
Revises: 34_2be441104b48
Create Date: 2023-01-30 11:31:09.335634

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = "35_0a5021e09eff"
down_revision = "34_2be441104b48"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("Collection", sa.Column("revised_at", sa.DateTime, nullable=True), schema="persistence_schema")


def downgrade():
    op.drop_column("Collection", "revised_at", schema="persistence_schema")
