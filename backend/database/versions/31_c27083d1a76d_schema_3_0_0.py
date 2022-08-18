"""schema 3.0.0

Revision ID: c27083d1a76d
Revises: 30_26c54abcaac9
Create Date: 2022-08-18 11:15:51.141347

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql


# revision identifiers, used by Alembic.
revision = '31_c27083d1a76d'
down_revision = '30_26c54abcaac9'
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("dataset", sa.Column("donor_id", postgresql.ARRAY(sa.String()), nullable=True))


def downgrade():
    op.drop_column("dataset", "donor_id")