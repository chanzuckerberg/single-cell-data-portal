"""rename_gene_symbols_field_to_genes

Revision ID: 7685265932f6
Revises: 9a8dd669a47c
Create Date: 2021-04-12 17:37:38.641573

"""
from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = '7685265932f6'
down_revision = '9a8dd669a47c'
branch_labels = None
depends_on = None


def upgrade():
    op.alter_column("geneset", "gene_symbols", new_column_name="genes")


def downgrade():
    op.alter_column("geneset", "genes", new_column_name="gene_symbols")
