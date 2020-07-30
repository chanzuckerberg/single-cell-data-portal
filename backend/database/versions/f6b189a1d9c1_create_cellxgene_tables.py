"""create_cellxgene_tables

Revision ID: f6b189a1d9c1
Revises: 899656d37baa
Create Date: 2020-07-30 12:03:18.275767

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
from sqlalchemy import Column, String, TIMESTAMP, func

revision = 'f6b189a1d9c1'
down_revision = '899656d37baa'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "cxguser",
        Column("id", String, nullable=False, primary_key=True),
        Column("created_at", TIMESTAMP, nullable=False, server_default=func.now()),
        Column("updated_at", TIMESTAMP, nullable=False, server_default=func.now()),
    )
    op.create_table(
        "cxgdataset",
        Column("id", String, nullable=False, primary_key=True),
        Column("name", String),
        Column("created_at", TIMESTAMP, nullable=False, server_default=func.now()),
        Column("updated_at", TIMESTAMP, nullable=False, server_default=func.now()),
    )
    op.create_table(
        "annotation",
        Column("id", String, nullable=False, primary_key=True),
        Column("user_id", String, nullable=False),
        Column("dataset_id", String, nullable=False),
        Column("tiledb_uri", String, nullable=False),
        Column("created_at", TIMESTAMP, nullable=False, server_default=func.now()),
    )
    op.create_foreign_key(
        "fk_annotation_user", "annotation", "cxguser", ["user_id"], ["id"],
    )
    op.create_foreign_key(
        "fk_annotation_dataset", "annotation", "cxgdataset", ["dataset_id"], ["id"],
    )



def downgrade():
    op.drop_table("cxguser")
    op.drop_table("cxgdataset")
    op.drop_table("annotation")

