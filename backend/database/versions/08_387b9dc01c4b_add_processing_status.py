"""add_processing_status

Revision ID: 387b9dc01c4b
Revises: 07_f6c4f66e8db6
Create Date: 2020-11-12 10:24:40.435293

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ENUM


# revision identifiers, used by Alembic.
revision = "08_387b9dc01c4b"
down_revision = "07_f6c4f66e8db6"
branch_labels = None
depends_on = None


def create_enum(name: str, values: list):
    values.sort()
    _enum = ENUM(*values, name=name.lower(), create_type=False)
    _enum.create(op.get_bind(), checkfirst=True)
    return _enum


def upgrade():

    upload_status_enum = create_enum(
        "uploadstatus", ["NA", "WAITING", "UPLOADING", "UPLOADED", "FAILED", "CANCEL_PENDING", "CANCELED"]
    )
    validation_status_enum = create_enum("validationstatus", ["NA", "VALIDATING", "VALID", "INVALID"])
    conversion_status_enum = create_enum("conversionstatus", ["NA", "CONVERTING", "CONVERTED", "FAILED"])

    op.create_table(
        "dataset_processing_status",
        sa.Column("id", sa.String(), nullable=False),
        sa.Column("dataset_id", sa.String(), nullable=False),
        sa.Column("upload_status", upload_status_enum, nullable=True),
        sa.Column("upload_progress", sa.Float(), nullable=True),
        sa.Column("upload_message", sa.String(), nullable=True),
        sa.Column("validation_status", validation_status_enum, nullable=True),
        sa.Column("validation_message", sa.String(), nullable=True),
        sa.Column("conversion_loom_status", conversion_status_enum, nullable=True),
        sa.Column("conversion_rds_status", conversion_status_enum, nullable=True),
        sa.Column("conversion_cxg_status", conversion_status_enum, nullable=True),
        sa.Column("conversion_anndata_status", conversion_status_enum, nullable=True),
        sa.ForeignKeyConstraint(["dataset_id"], ["dataset.id"], name="dataset_processing_status_dataset_id_fkey"),
        sa.PrimaryKeyConstraint("id", name="dataset_processing_status_pkey"),
    )


def downgrade():
    op.drop_table("dataset_processing_status")
    sa.Enum(name="uploadstatus").drop(op.get_bind(), checkfirst=False)
    sa.Enum(name="validationstatus").drop(op.get_bind(), checkfirst=False)
    sa.Enum(name="conversionstatus").drop(op.get_bind(), checkfirst=False)
