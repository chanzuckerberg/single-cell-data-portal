"""remove_loom_enum

Revision ID: 22_449eea8f5aeb
Revises: 21_c830a9e1f874
Create Date: 2021-10-13 15:53:34.183438

"""
from alembic import op

# revision identifiers, used by Alembic.
revision = '22_449eea8f5aeb'
down_revision = '21_c830a9e1f874'
branch_labels = None
depends_on = None


def upgrade():
    op.execute("ALTER TYPE datasetartifactfiletype RENAME TO datasetartifactfiletype_old;")
    op.execute("CREATE TYPE datasetartifactfiletype AS ENUM('H5AD', 'RDS', 'CXG');")
    op.execute(
        "ALTER TABLE dataset_artifact ALTER COLUMN filetype TYPE datasetartifactfiletype USING "
        "filetype::text::datasetartifactfiletype"
    )
    op.execute("DROP TYPE datasetartifactfiletype_old;")


def downgrade():
    #  db/test_migration  fails if the enums are out of order so it is necessary to recreate the enum
    # instead of just adding the value via:
    # op.execute("ALTER TYPE datasetartifactfiletype ADD VALUE 'LOOM';")

    op.execute("ALTER TYPE datasetartifactfiletype RENAME TO datasetartifactfiletype_old;")
    op.execute("CREATE TYPE datasetartifactfiletype AS ENUM('H5AD', 'RDS', 'LOOM', 'CXG');")
    op.execute(
        "ALTER TABLE dataset_artifact ALTER COLUMN filetype TYPE datasetartifactfiletype USING "
        "filetype::text::datasetartifactfiletype"
    )
    op.execute("DROP TYPE datasetartifactfiletype_old;")
