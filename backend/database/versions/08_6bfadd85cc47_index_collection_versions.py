"""index_collection_versions

Revision ID: 6bfadd85cc47
Revises: 07_e31a29561f38
Create Date: 2024-06-17 09:32:59.945987

"""

from alembic import op

# revision identifiers, used by Alembic.
revision = "08_6bfadd85cc47"
down_revision = "07_e31a29561f38"
branch_labels = None
depends_on = None


def upgrade():
    op.create_index("CollectionVersionIndex", "CollectionVersion", ["id"])


def downgrade():
    op.drop_index("CollectionVersionIndex", "CollectionVersion")
