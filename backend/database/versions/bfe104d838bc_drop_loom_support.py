"""drop_loom_support

Revision ID: 20_bfe104d838bc
Revises: 19_23bade351f23
Create Date: 2021-10-08 12:32:38.874124

"""
from alembic import op

# revision identifiers, used by Alembic.
revision = "bfe104d838bc"
down_revision = "23bade351f23"
branch_labels = None
depends_on = None


def upgrade():
    op.drop_column("dataset_processing_status", "conversion_loom_status")


def downgrade():
    op.execute(
        """
        ALTER TABLE dataset_processing_status
        ADD COLUMN conversion_loom_status conversionstatus
    """
    )
