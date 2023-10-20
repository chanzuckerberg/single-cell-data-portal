import logging
import os
import pathlib
import sys
import time
from typing import Optional

import tiledb

from backend.common.utils.result_notification import (
    format_failed_batch_issue_slack_alert,
    gen_wmg_pipeline_failure_message,
    gen_wmg_pipeline_success_message,
    notify_slack,
)
from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME
from backend.wmg.pipeline.cell_type_ordering import create_cell_type_ordering
from backend.wmg.pipeline.constants import (
    CELL_TYPE_ORDERING_CREATED_FLAG,
    DATASET_METADATA_CREATED_FLAG,
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG,
    FILTER_RELATIONSHIPS_CREATED_FLAG,
    MARKER_GENES_CUBE_CREATED_FLAG,
    PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG,
    WMG_DATA_SCHEMA_VERSION,
    WMG_PIPELINE_TEST_RUN_KEY,
)
from backend.wmg.pipeline.dataset_metadata import create_dataset_metadata
from backend.wmg.pipeline.errors import PipelineStepMissing
from backend.wmg.pipeline.expression_summary_and_cell_counts import create_expression_summary_and_cell_counts_cubes
from backend.wmg.pipeline.expression_summary_default import create_expression_summary_default_cube
from backend.wmg.pipeline.filter_relationships import create_filter_relationships_graph
from backend.wmg.pipeline.load_cube import upload_artifacts_to_s3
from backend.wmg.pipeline.marker_genes import create_marker_genes_cube
from backend.wmg.pipeline.primary_filter_dimensions import create_primary_filter_dimensions
from backend.wmg.pipeline.utils import load_pipeline_state, log_func_runtime
from backend.wmg.pipeline.validation.validation import Validation

logger = logging.getLogger(__name__)

PIPELINE_STEPS = [
    {
        "flag": EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
        "step": create_expression_summary_and_cell_counts_cubes,
    },
    {"flag": EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG, "step": create_expression_summary_default_cube},
    {"flag": FILTER_RELATIONSHIPS_CREATED_FLAG, "step": create_filter_relationships_graph},
    {"flag": PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG, "step": create_primary_filter_dimensions},
    {"flag": MARKER_GENES_CUBE_CREATED_FLAG, "step": create_marker_genes_cube},
    {"flag": DATASET_METADATA_CREATED_FLAG, "step": create_dataset_metadata},
    {"flag": CELL_TYPE_ORDERING_CREATED_FLAG, "step": create_cell_type_ordering},
]


@log_func_runtime
def run_pipeline(corpus_path: Optional[str] = None):
    if corpus_path is None:
        corpus_path = str(int(time.time()))

    logger.info(f"Creating directory {corpus_path}")
    pathlib.Path(corpus_path).mkdir(parents=True, exist_ok=True)

    pipeline_state = load_pipeline_state(corpus_path)

    for pipeline_step in PIPELINE_STEPS:
        flag = pipeline_step["flag"]
        step = pipeline_step["step"]

        if not pipeline_state.get(flag):
            step(corpus_path)

    if os.environ.get(WMG_PIPELINE_TEST_RUN_KEY) != "true":
        try:
            is_valid = Validation(corpus_path).validate_cube()
        except Exception:
            is_valid = False

        snapshot_id = os.path.basename(os.path.normpath(corpus_path))
        cube_data_s3_path = upload_artifacts_to_s3(
            snapshot_source_path=corpus_path,
            snapshot_schema_version=WMG_DATA_SCHEMA_VERSION,
            snapshot_id=snapshot_id,
            is_snapshot_validation_successful=is_valid,
        )
        stats = _get_stats(corpus_path)

        if is_valid:
            logger.info(f"Updated latest_snapshot_identifier in s3. Current snapshot location: {cube_data_s3_path}")

        return cube_data_s3_path, stats, is_valid


def main():
    """
    To trigger the batch process in prod via the cli run
    AWS_PROFILE=single-cell-prod aws batch submit-job --job-name $JOB_NAME --job-queue dp-prod --job-definition dp-prod-prodstack-wmg-processing # noqa E501
    """
    # todo pass in validate_cubes as env arg
    try:
        corpus_path, stats, is_valid = run_pipeline()
        if not is_valid:
            pipeline_failure_message = gen_wmg_pipeline_failure_message(
                "Issue with cube validation, see logs for more detail"
            )
            data = format_failed_batch_issue_slack_alert(pipeline_failure_message)
            notify_slack(data)
        else:
            success_message = gen_wmg_pipeline_success_message(corpus_path, **stats)
            notify_slack(success_message)
    except Exception as e:
        logger.exception("Pipeline failed")
        failure_message = format_failed_batch_issue_slack_alert(
            gen_wmg_pipeline_failure_message(
                f"Issue with WMG snapshot generation pipeline: {e}. See logs for more detail."
            )
        )
        notify_slack(failure_message)


if __name__ == "__main__":
    main()
    sys.exit()


def _get_stats(corpus_path: str) -> dict[str, int]:
    pipeline_state = load_pipeline_state(corpus_path)
    if not pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG):
        raise PipelineStepMissing("cell_counts")

    # get dataset count
    with tiledb.open(os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)) as cc_cube:
        cell_counts_df = cc_cube.df[:]
    dataset_count = len(cell_counts_df["dataset_id"].unique())
    cell_count = int(cell_counts_df["n_cells"].sum())

    return {"dataset_count": dataset_count, "cell_count": cell_count}
