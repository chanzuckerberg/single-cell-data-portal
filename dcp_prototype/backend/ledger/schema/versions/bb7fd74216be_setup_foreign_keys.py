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
    op.create_foreign_key(
        "library_project", "library", "project", ["project_id"], ["id"]
    )

    op.create_foreign_key(
        "library_process_to_biosample_prep",
        "biosample_prep_library_library_prep_protocol_process_join",
        "biosample_prep",
        ["biosample_prep_id"],
        ["id"],
    )

    op.create_foreign_key(
        "library_process_to_library_prep_protocol",
        "biosample_prep_library_library_prep_protocol_process_join",
        "library_prep_protocol",
        ["library_prep_protocol_id"],
        ["id"],
    )

    op.create_foreign_key(
        "library_process_to_library_prep_protocol_library",
        "biosample_prep_library_library_prep_protocol_process_join",
        "library",
        ["library_id"],
        ["id"],
    )

    op.create_foreign_key(
        "project_contributor_join_project",
        "project_contributor_join",
        "project",
        ["project_id"],
        ["id"],
    )

    op.create_foreign_key(
        "project_contributor_join_contributor",
        "project_contributor_join",
        "contributor",
        ["contributor_id"],
        ["id"],
    )

    op.create_foreign_key(
        "sequencing_process_to_library",
        "library_sequence_file_sequencing_protocol_process_join",
        "library",
        ["library_id"],
        ["id"],
    )

    op.create_foreign_key(
        "sequencing_process_to_sequence_file",
        "library_sequence_file_sequencing_protocol_process_join",
        "sequence_file",
        ["sequence_file_id"],
        ["id"],
    )

    op.create_foreign_key(
        "sequencing_process_to_sequencing_protocol",
        "library_sequence_file_sequencing_protocol_process_join",
        "sequencing_protocol",
        ["sequencing_protocol_id"],
        ["id"],
    )

    op.create_foreign_key(
        "bam_qc_metric_analysis_file",
        "bam_qc_metric",
        "analysis_file",
        ["analysis_file_id"],
        ["id"],
    )

    op.create_foreign_key(
        "alignment_process_to_analysis_file",
        "sequence_file_alignment_protocol_analysis_file_process_join",
        "analysis_file",
        ["analysis_file_id"],
        ["id"],
    )

    op.create_foreign_key(
        "alignment_process_to_sequence_file",
        "sequence_file_alignment_protocol_analysis_file_process_join",
        "sequence_file",
        ["sequence_file_id"],
        ["id"],
    )

    op.create_foreign_key(
        "alignment_process_to_alignment_protocol",
        "sequence_file_alignment_protocol_analysis_file_process_join",
        "alignment_protocol",
        ["alignment_protocol_id"],
        ["id"],
    )

    op.create_foreign_key(
        "quantification_process_to_analysis_file",
        "analysis_file_quantification_protocol_expression_file_process_join",
        "analysis_file",
        ["analysis_file_id"],
        ["id"],
    )

    op.create_foreign_key(
        "quantification_process_to_quantification_protocol",
        "analysis_file_quantification_protocol_expression_file_process_join",
        "quantification_protocol",
        ["quantification_protocol_id"],
        ["id"],
    )

    op.create_foreign_key(
        "quantification_process_to_expression_file",
        "analysis_file_quantification_protocol_expression_file_process_join",
        "expression_file",
        ["expression_file_id"],
        ["id"],
    )


def downgrade():
    pass
