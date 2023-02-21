"""deprecate_obfuscated_uuid

Revision ID: 26_913a1204bf96
Revises: 25_43e20ba9f9f0
Create Date: 2022-03-07 14:04:09.806504

"""
import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "26_913a1204bf96"
down_revision = "25_43e20ba9f9f0"
branch_labels = None
depends_on = None


def upgrade():
    op.drop_column("project", "obfuscated_uuid")


def downgrade():
    op.add_column("project", sa.Column("obfuscated_uuid", sa.VARCHAR(), autoincrement=False, nullable=True))
