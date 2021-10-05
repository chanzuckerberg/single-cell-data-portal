"""update_conversion_status

Revision ID: c830a9e1f874
Revises: 424e875043d3
Create Date: 2021-10-04 22:56:28.324384

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = 'c830a9e1f874'
down_revision = '424e875043d3'
branch_labels = None
depends_on = None


def upgrade():
    with op.get_context().autocommit_block():
        op.alter_column('dataset_processing_status', 'conversion_anndata_status', new_column_name='anndata_status')
        op.alter_column('dataset_processing_status', 'conversion_rds_status', new_column_name='rds_status')
        op.alter_column('dataset_processing_status', 'conversion_cxg_status', new_column_name='cxg_status')
        op.alter_column('dataset_processing_status', 'conversion_loom_status', new_column_name='loom_status')

        op.execute("ALTER TYPE conversionstatus ADD VALUE 'UPLOADING'")
        op.execute("ALTER TYPE conversionstatus ADD VALUE 'UPLOADED'")


def downgrade():
    op.execute("UPDATE dataset_processing_status SET anndata_status = 'NA' WHERE anndata_status in ('UPLOADING', 'UPLOADED');")
    op.execute("UPDATE dataset_processing_status SET rds_status = 'NA' WHERE rds_status in ('UPLOADING', 'UPLOADED');")
    op.execute("UPDATE dataset_processing_status SET cxg_status = 'NA' WHERE cxg_status in ('UPLOADING', 'UPLOADED');")
    op.execute("UPDATE dataset_processing_status SET loom_status = 'NA' WHERE loom_status in ('UPLOADING', 'UPLOADED');")


    op.execute("ALTER TYPE conversionstatus RENAME TO old_conversionstatus;")
    op.execute("CREATE TYPE conversionstatus AS ENUM('NA', 'CONVERTING', 'CONVERTED', 'FAILED');")
    op.execute(
        "ALTER TABLE dataset_processing_status ALTER COLUMN cxg_status TYPE conversionstatus USING "
        "cxg_status::text::conversionstatus"
    )
    op.execute(
        "ALTER TABLE dataset_processing_status ALTER COLUMN anndata_status TYPE conversionstatus USING "
        "anndata_status::text::conversionstatus"
    )
    op.execute(
        "ALTER TABLE dataset_processing_status ALTER COLUMN rds_status TYPE conversionstatus USING "
        "rds_status::text::conversionstatus"
    )
    op.execute(
        "ALTER TABLE dataset_processing_status ALTER COLUMN loom_status TYPE conversionstatus USING "
        "cxg_status::text::conversionstatus"
    )

    # op.execute("DROP TYPE old_conversionstatus") TODO uncomment out before merging

    op.alter_column('dataset_processing_status', 'anndata_status', new_column_name='conversion_anndata_status')
    op.alter_column('dataset_processing_status', 'rds_status', new_column_name='conversion_rds_status')
    op.alter_column('dataset_processing_status', 'cxg_status', new_column_name='conversion_cxg_status')
    op.alter_column('dataset_processing_status', 'loom_status', new_column_name='conversion_loom_status')

