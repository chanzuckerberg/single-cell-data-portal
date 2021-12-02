"""add_summary_link

Revision ID: 03_fa0950f6df19
Revises: 02_7968f161b60c
Create Date: 2020-08-24 11:04:30.640496

"""
from alembic import op
from sqlalchemy.dialects.postgresql import ENUM


# revision identifiers, used by Alembic.
revision = "03_fa0950f6df19"
down_revision = "02_7968f161b60c"
branch_labels = None
depends_on = None


def upgrade():
    project_link_enum = ENUM(name="projectlink")
    project_link_enum.drop(op.get_bind(), checkfirst=True)
    project_link_enum = ENUM("PROTOCOL", "SUMMARY", "RAW_DATA", "OTHER", name="projectlink", create_type=False)
    project_link_enum.create(op.get_bind(), checkfirst=True)

    op.alter_column(
        table_name="project_link",
        column_name="link_type",
        type_=project_link_enum,
        nullable=True,
        postgresql_using=(
            "CASE WHEN link_type = 'PROTOCOL' THEN 'PROTOCOL'::projectlink "
            "     WHEN link_type = 'RAW_DATA' THEN 'RAW_DATA'::projectlink "
            "     WHEN link_type = 'OTHER' THEN 'OTHER'::projectlink "
            "     END"
        ),
    )


def downgrade():
    project_link_enum = ENUM(name="projectlink")
    project_link_enum.drop(op.get_bind(), checkfirst=True)
    project_link_enum = ENUM("PROTOCOL", "RAW_DATA", "OTHER", name="projectlink", create_type=False)
    project_link_enum.create(op.get_bind(), checkfirst=True)

    op.alter_column(
        table_name="project_link",
        column_name="link_type",
        type_=project_link_enum,
        nullable=True,
        postgresql_using=(
            "CASE WHEN link_type = 'PROTOCOL' THEN 'PROTOCOL'::projectlink "
            "     WHEN link_type = 'RAW_DATA' THEN 'RAW_DATA'::projectlink "
            "     WHEN link_type = 'OTHER' THEN 'OTHER'::projectlink "
            "     WHEN link_type = 'SUMMARY' THEN 'OTHER'::projectlink "
            "     END"
        ),
    )
