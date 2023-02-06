"""change_ontology_to_json

Revision ID: 02_7968f161b60c
Revises: 01_899656d37baa
Create Date: 2020-08-20 08:58:21.531154

"""
from alembic import op
from sqlalchemy import Column, String
from sqlalchemy.dialects.postgresql import JSONB

# revision identifiers, used by Alembic.
revision = "02_7968f161b60c"
down_revision = "01_899656d37baa"
branch_labels = None
depends_on = None


def upgrade():
    op.alter_column(
        table_name="dataset",
        column_name="organism",
        type_=JSONB(),
        postgresql_using="json_build_object('ontology_term_id', organism_ontology, 'label', organism)",
    )
    op.drop_column(table_name="dataset", column_name="organism_ontology")

    op.alter_column(
        table_name="dataset",
        column_name="tissue",
        type_=JSONB(),
        postgresql_using="json_build_array(json_build_object('ontology_term_id', tissue_ontology, 'label', tissue))",
    )
    op.drop_column(table_name="dataset", column_name="tissue_ontology")

    op.alter_column(
        table_name="dataset",
        column_name="assay",
        type_=JSONB(),
        postgresql_using="json_build_array(json_build_object('ontology_term_id', assay_ontology, 'label', assay))",
    )
    op.drop_column(table_name="dataset", column_name="assay_ontology")

    op.alter_column(
        table_name="dataset",
        column_name="disease",
        type_=JSONB(),
        postgresql_using="json_build_array(json_build_object('ontology_term_id', disease_ontology, 'label', disease))",
    )
    op.drop_column(table_name="dataset", column_name="disease_ontology")

    op.alter_column(
        table_name="dataset",
        column_name="ethnicity",
        type_=JSONB(),
        postgresql_using=(
            "json_build_array(json_build_object(" "'ontology_term_id', ethnicity_ontology, 'label', ethnicity))"
        ),
    )
    op.drop_column(table_name="dataset", column_name="ethnicity_ontology")

    op.alter_column(
        table_name="dataset",
        column_name="development_stage",
        type_=JSONB(),
        postgresql_using=(
            "json_build_array(json_build_object("
            "'ontology_term_id', development_stage_ontology, 'label', development_stage))"
        ),
    )
    op.drop_column(table_name="dataset", column_name="development_stage_ontology")

    op.alter_column(table_name="dataset", column_name="sex", type_=JSONB(), postgresql_using="json_build_array(sex)")


def downgrade():
    """
    This will cleanly downgrade the _schema_, but it's not possible to losslessly downgrade
    the _data_. So rather than discard elements in the json arrays, it just stringifies the
    whole thing. After running this, you'll have a bunch of values that need to be parsed
    and decisions about which elements to get rid off will need to be made.
    """
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
