"""remove_curator_tag and add ProcessingStatus Initialize

Revision ID: 31_253d1d67ea4a
Revises: 30_26c54abcaac9
Create Date: 2022-09-07 13:39:39.401203

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = "31_253d1d67ea4a"
down_revision = "30_26c54abcaac9"
branch_labels = None
depends_on = None

# Enum 'type' for PostgreSQL
name = "processingstatus"
# Set temporary enum 'type' for PostgreSQL
temp_name = name + "_temp"


def upgrade():

    # Add INITIALIZE as a processing status
    op.execute(f"ALTER TYPE {name} RENAME TO {temp_name};")
    op.execute(f"CREATE TYPE {name} AS ENUM('PENDING', 'SUCCESS', 'FAILURE', 'INITIALIZED');")
    op.execute(
        f"ALTER TABLE dataset_processing_status ALTER COLUMN processing_status TYPE {name} USING "
        f"processing_status::text::{name}"
    )
    op.execute(f"DROP TYPE {temp_name};")

    # Remove curator_tag
    op.drop_constraint("_dataset__collection_id_curator_tag", "dataset", type_="unique")
    op.drop_column("dataset", "curator_tag")


def downgrade():
    # Remove INITIALIZE as a processing status
    op.execute(f"ALTER TYPE {name} RENAME TO {temp_name};")
    op.execute(f"CREATE TYPE {name} AS ENUM('PENDING', 'SUCCESS', 'FAILURE');")
    op.execute(
        f"ALTER TABLE dataset_processing_status ALTER COLUMN processing_status TYPE {name} USING "
        f"processing_status::text::{name}"
    )
    op.execute(f"DROP TYPE {temp_name};")

    # Add curator_tag
    op.add_column("dataset", sa.Column("curator_tag", sa.VARCHAR(), autoincrement=False, nullable=True))
    op.create_unique_constraint("_dataset__collection_id_curator_tag", "dataset", ["collection_id", "curator_tag"])
