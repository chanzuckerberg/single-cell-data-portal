"""make_collection_id_pkey

Revision ID: 28_e08b6b67f076
Revises: 27_4d70e8c321e3
Create Date: 2022-03-04 21:54:51.795657

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = "29_e08b6b67f076"
down_revision = "28_d7976d87cc2f"
branch_labels = None
depends_on = None


def upgrade():
    # Drop foreign keys that reference composite key [project.id, project.visibility]
    op.drop_constraint("dataset_collection_id_fkey", "dataset", type_="foreignkey")
    op.drop_column("dataset", "collection_visibility")
    op.drop_constraint("project_link_collection_id_fkey", "project_link", type_="foreignkey")
    op.drop_column("project_link", "collection_visibility")

    # Drop foreign key that references composite key [project.id, project.visibility]
    op.drop_constraint("geneset_collection_id_fkey", "geneset", type_="foreignkey")
    # Drop unique key that uses collection_visibility
    op.drop_constraint("_geneset_name__collection_uc", "geneset", type_="unique")
    op.drop_column("geneset", "collection_visibility")
    # Recreate unique key without collection_visibility
    op.create_unique_constraint(op.f("_geneset_name__collection_uc"), "geneset", ["name", "collection_id"])

    # Add project.revision_of
    op.add_column("project", sa.Column("revision_of", sa.String(), nullable=True, unique=True))

    # Drop project composite primary key; add project.id as primary key
    op.drop_constraint("project_pkey", "project")
    op.create_primary_key(op.f("pk_project"), "project", ["id"])

    # Add foreign keys for project.id
    op.create_foreign_key(op.f("fk_dataset_collection_id_project"), "dataset", "project", ["collection_id"], ["id"])
    op.create_foreign_key(op.f("fk_geneset_collection_id_project"), "geneset", "project", ["collection_id"], ["id"])
    op.create_foreign_key(op.f("fk_project_revision_of_project"), "project", "project", ["revision_of"], ["id"])
    op.create_foreign_key(
        op.f("fk_project_link_collection_id_project"), "project_link", "project", ["collection_id"], ["id"]
    )

    # Check constraint: published Collections must not be a revision_of another Collection
    op.create_check_constraint(
        op.f("ck_project_private_revision"), "project", "revision_of IS NULL OR visibility = 'PRIVATE'"
    )


def downgrade():
    def remake_collection_id_and_collection_visibility_composite_foreign_key(table: str):
        op.add_column(
            table,
            sa.Column(
                "collection_visibility",
                postgresql.ENUM("PRIVATE", "PUBLIC", name="collectionvisibility"),
                autoincrement=False,
            ),
        )
        op.execute(
            f"""
            UPDATE {table}
                SET collection_visibility = (
                    SELECT visibility
                    FROM project
                    WHERE project.id = {table}.collection_id
                )
            """
        )
        op.alter_column(table, "collection_visibility", nullable=False)
        op.create_foreign_key(
            f"{table}_collection_id_fkey",
            table,
            "project",
            ["collection_id", "collection_visibility"],
            ["id", "visibility"],
        )

    # Drop foreign keys for project.id
    op.drop_constraint(op.f("fk_project_link_collection_id_project"), "project_link", type_="foreignkey")
    op.drop_constraint(op.f("fk_project_revision_of_project"), "project", type_="foreignkey")
    op.drop_constraint(op.f("fk_geneset_collection_id_project"), "geneset", type_="foreignkey")
    op.drop_constraint(op.f("fk_dataset_collection_id_project"), "dataset", type_="foreignkey")

    # Drop project.id as primary key
    op.drop_constraint(op.f("pk_project"), "project")

    # Reassign id for revisions
    op.execute(
        """
        UPDATE project
            SET id = revision_of, visibility = 'PRIVATE'
        WHERE revision_of IS NOT NULL
        """
    )

    # Remove project.revision_of
    op.drop_constraint(op.f("ck_project_private_revision"), "project", type_="check")
    op.drop_column("project", "revision_of")

    # Add project composite primary key
    op.create_primary_key("project_pkey", "project", ["id", "visibility"])

    # Create foreign keys that reference composite key [project.id, project.visibility]
    remake_collection_id_and_collection_visibility_composite_foreign_key("project_link")
    remake_collection_id_and_collection_visibility_composite_foreign_key("geneset")
    remake_collection_id_and_collection_visibility_composite_foreign_key("dataset")

    # Reset unique constraint
    op.drop_constraint("_geneset_name__collection_uc", "geneset", type_="unique")
    op.create_unique_constraint(
        "_geneset_name__collection_uc", "geneset", ["name", "collection_id", "collection_visibility"]
    )
