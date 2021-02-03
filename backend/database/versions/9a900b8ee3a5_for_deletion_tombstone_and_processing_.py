"""For_deletion—tombstone_and_processing_status

Revision ID: 9a900b8ee3a5
Revises: 7794b1ea430f
Create Date: 2021-01-25 15:44:47.054499

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ENUM


# revision identifiers, used by Alembic.
revision = "9a900b8ee3a5"
down_revision = "7794b1ea430f"
branch_labels = None
depends_on = None


def create_enum(name: str, values: list):
    values.sort()
    _enum = ENUM(*values, name=name.lower(), create_type=False)
    _enum.create(op.get_bind(), checkfirst=True)
    return _enum


def upgrade():
    # Create processing_status enum
    processing_status_enum = create_enum("processingstatus", ["PENDING", "SUCCESS", "FAILURE"])

    # add processing_status_enum to processing_status column in dataset_process_status table
    op.add_column("dataset_processing_status", sa.Column("processing_status", processing_status_enum))

    # Add tombstone column to dataset with default false
    op.add_column("dataset", sa.Column("tombstone", sa.Boolean, default=False))

    # Add tombstone column to collection with default false
    op.add_column("project", sa.Column("tombstone", sa.Boolean, default=False))


def downgrade():
    # Remove processing_status from dataset_processing_status table
    op.drop_column("dataset_processing_status", "processing_status")

    # Remove tombstone column from dataset
    op.drop_column("dataset", "tombstone")

    # remove tombstone column from collection
    op.drop_column("project", "tombstone")

    # Remove processing_status_enum
    sa.Enum(name="processingstatus").drop(op.get_bind(), checkfirst=False)
