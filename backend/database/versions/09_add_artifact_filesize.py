"""Add filesize column to DatasetArtifact

Stores artifact file size in the database so it can be served directly without
a live S3 HEAD request per artifact. Nullable to allow zero-downtime rollout
and backfill of existing rows.

Revision ID: 09_add_artifact_filesize
Revises: 08_92c817dddc7d
Create Date: 2026-04-10

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "09_add_artifact_filesize"
down_revision = "08_92c817dddc7d"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column(
        "DatasetArtifact",
        sa.Column("filesize", sa.BigInteger(), nullable=True),
        schema="persistence_schema",
    )


def downgrade():
    op.drop_column("DatasetArtifact", "filesize", schema="persistence_schema")
