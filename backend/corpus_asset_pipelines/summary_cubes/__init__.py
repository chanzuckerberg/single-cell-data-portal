import sys

from backend.common.utils.result_notification import (
    format_failed_batch_issue_slack_alert,
    notify_slack,
    gen_wmg_pipeline_failure_message,
)

from backend.corpus_asset_pipelines.summary_cubes.expression_summary.job import create_expression_summary_cube
from backend.wmg.data.validation.validation import Validation
from backend.corpus_asset_pipelines.summary_cubes.cell_count import create_cell_count_cube


def run(corpus_path: str, validate_cube: bool) -> dict:
    """
    Build expression summary cube and cell count cube based
    on cell data stored in integrated corpus
    validate expression summary cube based on biological expectations
    if indicated by param
    """
    create_expression_summary_cube(corpus_path)
    create_cell_count_cube(corpus_path)
    if validate_cube:
        if Validation(corpus_path).validate_cube() is False:
            pipeline_failure_message = gen_wmg_pipeline_failure_message(
                "Issue with cube validation, see logs for more detail"
            )
            data = format_failed_batch_issue_slack_alert(pipeline_failure_message)
            notify_slack(data)
            sys.exit("Exiting due to cube validation failure")
