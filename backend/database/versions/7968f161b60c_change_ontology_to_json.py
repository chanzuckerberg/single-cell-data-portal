"""change_ontology_to_json

Revision ID: 7968f161b60c
Revises: 899656d37baa
Create Date: 2020-08-20 08:58:21.531154

"""
from alembic import op
from sqlalchemy import Column, String
from sqlalchemy.dialects.postgresql import JSONB


# revision identifiers, used by Alembic.
revision = "7968f161b60c"
down_revision = "899656d37baa"
branch_labels = None
depends_on = None


def upgrade():
    op.alter_column(table_name="dataset", column_name="organism", type_=JSONB())
    op.drop_column(table_name="dataset", column_name="organism_ontology")

    op.alter_column(table_name="dataset", column_name="tissue", type_=JSONB())
    op.drop_column(table_name="dataset", column_name="tissue_ontology")

    op.alter_column(table_name="dataset", column_name="assay", type_=JSONB())
    op.drop_column(table_name="dataset", column_name="assay_ontology")

    op.alter_column(table_name="dataset", column_name="disease", type_=JSONB())
    op.drop_column(table_name="dataset", column_name="disease_ontology")

    op.alter_column(table_name="dataset", column_name="ethnicity", type_=JSONB())
    op.drop_column(table_name="dataset", column_name="ethnicity_ontology")

    op.alter_column(table_name="dataset", column_name="development_stage", type_=JSONB())
    op.drop_column(table_name="dataset", column_name="development_stage_ontology")

    op.alter_column(table_name="dataset", column_name="sex", type_=JSONB())


def downgrade():
    op.alter_column(table_name="dataset", column_name="organism", type_=String())
    op.add_column("dataset", Column("organism_ontology", String()))

    op.alter_column(table_name="dataset", column_name="tissue", type_=String())
    op.add_column("dataset", Column("tissue_ontology", String()))

    op.alter_column(table_name="dataset", column_name="assay", type_=String())
    op.add_column("dataset", Column("assay_ontology", String()))

    op.alter_column(table_name="dataset", column_name="disease", type_=String())
    op.add_column("dataset", Column("disease_ontology", String()))

    op.alter_column(table_name="dataset", column_name="ethnicity", type_=String())
    op.add_column("dataset", Column("ethnicity_ontology", String()))

    op.alter_column(table_name="dataset", column_name="development_stage", type_=String())
    op.add_column("dataset", Column("development_stage_ontology", String()))

    op.alter_column(table_name="dataset", column_name="sex", type_=String())
