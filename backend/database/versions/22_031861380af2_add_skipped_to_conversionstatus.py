"""add_SKIPPED_to_conversionstatus

Revision ID: 22_031861380af2
Revises: 21_c830a9e1f874
Create Date: 2021-11-22 15:43:00.821321

"""
from alembic import op

# revision identifiers, used by Alembic.
revision = "22_031861380af2"
down_revision = "21_c830a9e1f874"
branch_labels = None
depends_on = None


def upgrade():
    with op.get_context().autocommit_block():
        op.execute("ALTER TYPE conversionstatus ADD VALUE 'SKIPPED'")


def downgrade():
    op.execute("UPDATE dataset_processing_status SET h5ad_status = 'NA' WHERE h5ad_status = 'SKIPPED';")
    op.execute("UPDATE dataset_processing_status SET rds_status = 'NA' WHERE rds_status = 'SKIPPED';")
    op.execute("UPDATE dataset_processing_status SET cxg_status = 'NA' WHERE cxg_status = 'SKIPPED';")

    op.execute("ALTER TYPE conversionstatus RENAME TO conversionstatus_old;")
    op.execute(
        "CREATE TYPE conversionstatus AS ENUM('NA', 'CONVERTING', 'CONVERTED', 'FAILED', 'UPLOADING', 'UPLOADED');"
    )
    op.execute(
        "ALTER TABLE dataset_processing_status ALTER COLUMN cxg_status TYPE conversionstatus USING "
        "cxg_status::text::conversionstatus"
    )
    op.execute(
        "ALTER TABLE dataset_processing_status ALTER COLUMN h5ad_status TYPE conversionstatus USING "
        "h5ad_status::text::conversionstatus"
    )
    op.execute(
        "ALTER TABLE dataset_processing_status ALTER COLUMN rds_status TYPE conversionstatus USING "
        "rds_status::text::conversionstatus"
    )

    op.execute("DROP TYPE conversionstatus_old")
