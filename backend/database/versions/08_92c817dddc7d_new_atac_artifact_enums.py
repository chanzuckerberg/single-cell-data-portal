"""new-atac-artifact-enums

Changing DatasetArtifactTable.type from an enum to a string to allow for new ATAC artifact types.

Revision ID: 08_92c817dddc7d
Revises: 07_e31a29561f38
Create Date: 2025-03-06 12:31:12.249307

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "08_92c817dddc7d"
down_revision = "07_e31a29561f38"
branch_labels = None
depends_on = None


def upgrade():
    op.alter_column(
        "DatasetArtifact",
        "type",
        type_=sa.String(),
        existing_type=sa.Enum("CXG", "RDS", "H5AD", "RAW_H5AD", name="datasetartifacttype"),
        schema="persistence_schema",
    )


def downgrade():
    op.alter_column(
        "DatasetArtifact",
        "type",
        type_=sa.Enum("CXG", "RDS", "H5AD", "RAW_H5AD", name="datasetartifacttype"),
        existing_type=sa.String(),
        schema="persistence_schema",
    )
