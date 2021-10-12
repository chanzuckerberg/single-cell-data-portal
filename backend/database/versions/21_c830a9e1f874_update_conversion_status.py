"""update_conversion_status

Revision ID: 21_c830a9e1f874
Revises: 20_bfe104d838bc
Create Date: 2021-10-04 22:56:28.324384

"""
from alembic import op

# revision identifiers, used by Alembic.
revision = "21_c830a9e1f874"
down_revision = "bfe104d838bc"
branch_labels = None
depends_on = None


def upgrade():
    with op.get_context().autocommit_block():
        op.alter_column("dataset_processing_status", "conversion_anndata_status", new_column_name="h5ad_status")
        op.alter_column("dataset_processing_status", "conversion_rds_status", new_column_name="rds_status")
        op.alter_column("dataset_processing_status", "conversion_cxg_status", new_column_name="cxg_status")

        op.execute("ALTER TYPE conversionstatus ADD VALUE 'UPLOADING'")
        op.execute("ALTER TYPE conversionstatus ADD VALUE 'UPLOADED'")


def downgrade():
    op.execute("UPDATE dataset_processing_status SET h5ad_status = 'NA' WHERE h5ad_status in ('UPLOADING', 'UPLOADED');") # noqa E501
    op.execute("UPDATE dataset_processing_status SET rds_status = 'NA' WHERE rds_status in ('UPLOADING', 'UPLOADED');")
    op.execute("UPDATE dataset_processing_status SET cxg_status = 'NA' WHERE cxg_status in ('UPLOADING', 'UPLOADED');")

    op.execute("ALTER TYPE conversionstatus RENAME TO conversionstatus_old;")
    op.execute("CREATE TYPE conversionstatus AS ENUM('NA', 'CONVERTING', 'CONVERTED', 'FAILED');")
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

    op.alter_column("dataset_processing_status", "h5ad_status", new_column_name="conversion_anndata_status")
    op.alter_column("dataset_processing_status", "rds_status", new_column_name="conversion_rds_status")
    op.alter_column("dataset_processing_status", "cxg_status", new_column_name="conversion_cxg_status")
