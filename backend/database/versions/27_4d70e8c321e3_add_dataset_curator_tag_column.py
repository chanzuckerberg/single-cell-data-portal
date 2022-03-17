"""Add dataset.curator_tag column

Revision ID: 27_4d70e8c321e3
Revises: 26_913a1204bf96
Create Date: 2022-03-14 16:07:56.963055

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = "27_4d70e8c321e3"
down_revision = "26_913a1204bf96"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("dataset", sa.Column("curator_tag", sa.String(), nullable=True))
    op.create_unique_constraint("_dataset__collection_id_curator_tag", "dataset", ["collection_id", "curator_tag"])


def downgrade():
    op.drop_constraint("_dataset__collection_id_curator_tag", "dataset", type_="unique")
    op.drop_column("dataset", "curator_tag")
