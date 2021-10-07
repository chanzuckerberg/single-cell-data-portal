"""update_data_submission_policy_version

Revision ID: 23bade351f23
Revises: a8cd0dc08805
Create Date: 2021-09-24 16:41:33.297342

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = "23bade351f23"
down_revision = "a8cd0dc08805"
branch_labels = None
depends_on = None


def upgrade():
    op.alter_column("project", "data_submission_policy_version", existing_type=sa.VARCHAR(), nullable=True)


def downgrade():
    op.execute("UPDATE project SET data_submission_policy_version='' WHERE data_submission_policy_version=null")

    op.alter_column("project", "data_submission_policy_version", existing_type=sa.VARCHAR(), nullable=False)
