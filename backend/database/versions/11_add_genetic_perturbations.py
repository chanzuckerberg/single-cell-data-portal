"""add genetic_perturbations column to DatasetVersion

Revision ID: 11_add_genetic_perturbations
Revises: 10_add_is_pre_analysis
Create Date: 2026-04-09

Stores the full uns['genetic_perturbations'] dictionary as a separate nullable JSON column
to avoid bloating the existing dataset_metadata column (~6.5 MB serialized per perturbation dataset).
"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "11_add_genetic_perturbations"
down_revision = "10_add_is_pre_analysis"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column(
        "DatasetVersion",
        sa.Column("genetic_perturbations", sa.JSON(), nullable=True),
        schema="persistence_schema",
    )


def downgrade():
    op.drop_column("DatasetVersion", "genetic_perturbations", schema="persistence_schema")
