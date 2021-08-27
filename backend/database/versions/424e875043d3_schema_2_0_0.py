"""schema-2.0.0

Revision ID: 424e875043d3
Revises: b65dad1f2a7e
Create Date: 2021-08-12 10:22:45.891595

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB, ENUM


# revision identifiers, used by Alembic.
revision = "424e875043d3"
down_revision = "b65dad1f2a7e"
branch_labels = None
depends_on = None


def create_enum(name: str, values: list):
    values.sort()
    _enum = ENUM(*values, name=name.lower(), create_type=False)
    _enum.create(op.get_bind(), checkfirst=True)
    return _enum


def upgrade():
    x_approximate_distribution_enum = create_enum("x_approximate_distribution", ["COUNT", "NORMAL"])
    is_primary_data_enum = create_enum("is_primary_data", ["PRIMARY", "SECONDARY", "BOTH"])

    op.add_column("dataset", sa.Column("cell_type", JSONB, nullable=True))
    op.add_column("dataset", sa.Column("is_primary_data", is_primary_data_enum, nullable=True))
    op.add_column("dataset", sa.Column("x_normalization", sa.String, nullable=True))
    op.add_column("dataset", sa.Column("x_approximate_distribution", x_approximate_distribution_enum, nullable=True))
    op.add_column("dataset", sa.Column("schema_version", sa.String, nullable=True))
    op.add_column("dataset", sa.Column("mean_genes_per_cell", sa.Float, server_default="0"))


def downgrade():
    op.drop_column("dataset", "mean_genes_per_cell")
    op.drop_column("dataset", "schema_version")
    op.drop_column("dataset", "x_approximate_distribution")
    op.drop_column("dataset", "x_normalization")
    op.drop_column("dataset", "is_primary_data")
    op.drop_column("dataset", "cell_type")
