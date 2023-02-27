"""adding revised_at for collections

Revision ID: 02_0a5021e09eff
Revises: 01_2be441104b48
Create Date: 2023-01-30 11:31:09.335634

"""
import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "02_0a5021e09eff"
down_revision = "01_2be441104b48"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("Collection", sa.Column("revised_at", sa.DateTime, nullable=True), schema="persistence_schema")


def downgrade():
    op.drop_column("Collection", "revised_at", schema="persistence_schema")
