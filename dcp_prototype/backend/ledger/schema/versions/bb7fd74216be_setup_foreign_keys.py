"""setup_foreign_keys

Revision ID: bb7fd74216be
Revises: 407ab1d9843a
Create Date: 2020-02-04 14:51:37.477811

"""
from alembic import op


# revision identifiers, used by Alembic.
revision = "bb7fd74216be"
down_revision = "407ab1d9843a"
branch_labels = None
depends_on = None


def upgrade():
    op.create_foreign_key("library_project", "library", "project", ["project_id"], ["id"])

    op.create_foreign_key(
        "project_contributor_join_project", "project_contributor_join", "project", ["project_id"], ["id"],
    )

    op.create_foreign_key(
        "project_contributor_join_contributor", "project_contributor_join", "contributor", ["contributor_id"], ["id"],
    )


def downgrade():
    pass
