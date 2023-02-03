"""Remove DatasetArtifactType

Revision ID: 28_d7976d87cc2f
Revises: 27_4d70e8c321e3
Create Date: 2022-04-03 15:01:26.748802

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = "28_d7976d87cc2f"
down_revision = "27_4d70e8c321e3"
branch_labels = None
depends_on = None


def upgrade():
    op.drop_column("dataset_artifact", "type")
    sa.Enum(name="datasetartifacttype").drop(op.get_bind())


def downgrade():
    dataset_artifact_type_enum = postgresql.ENUM("ORIGINAL", "REMIX", name="datasetartifacttype")
    dataset_artifact_type_enum.create(op.get_bind())
    op.add_column("dataset_artifact", sa.Column("type", dataset_artifact_type_enum, autoincrement=False, nullable=True))
