"""schema 3.0.0

Revision ID: 32_c27083d1a76d
Revises: 31_253d1d67ea4a
Create Date: 2022-08-18 11:15:51.141347

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql


# revision identifiers, used by Alembic.
revision = "32_c27083d1a76d"
down_revision = "31_253d1d67ea4a"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("dataset", sa.Column("suspension_type", postgresql.ARRAY(sa.String()), nullable=True))
    op.add_column("dataset", sa.Column("donor_id", postgresql.ARRAY(sa.String()), nullable=True))
    op.alter_column("dataset", "ethnicity", new_column_name="self_reported_ethnicity")
    op.drop_column("dataset", "x_normalization")


def downgrade():
    op.drop_column("dataset", "donor_id")
    op.drop_column("dataset", "suspension_type")
    op.alter_column("dataset", "self_reported_ethnicity", new_column_name="ethnicity")
    op.add_column("dataset", sa.Column("x_normalization", sa.String, nullable=True))
