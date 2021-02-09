"""tombstone_false

Revision ID: 5a363594dd06
Revises: 9a900b8ee3a5
Create Date: 2021-02-08 16:45:25.402786

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '5a363594dd06'
down_revision = '9a900b8ee3a5'
branch_labels = None
depends_on = None


def upgrade():
    def tombstone_false(table: str, column: str, default):
        op.execute(f"UPDATE {table} SET {column} = {default}")
        op.alter_column(table, column, nullable=False)

    tombstone_false('dataset', 'tombstone', False)
    tombstone_false('project', 'tombstone', False)


def downgrade():
    op.alter_column('dataset', 'tombstone', nullable=True)
    op.alter_column('project', 'tombstone', nullable=True)
