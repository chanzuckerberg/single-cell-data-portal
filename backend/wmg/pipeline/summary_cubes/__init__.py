from backend.common.utils.exceptions import CubeValidationException
from backend.common.utils.result_notification import (
    format_failed_batch_issue_slack_alert,
    gen_wmg_pipeline_failure_message,
    notify_slack,
)
from backend.wmg.data.validation.validation import Validation
from backend.wmg.pipeline.summary_cubes.cell_count import create_cell_count_cube
from backend.wmg.pipeline.summary_cubes.marker_genes import create_marker_genes_cube


def run(corpus_path: str, validate_cube: bool) -> dict:
    """
    Build expression summary cube and cell count cube based
    on cell data stored in integrated corpus
    validate expression summary cube based on biological expectations
    if indicated by param
    """
    create_cell_count_cube(corpus_path)
    create_marker_genes_cube(corpus_path)

    if validate_cube:
        try:
            is_valid = Validation(corpus_path).validate_cube()
        except Exception as e:
            raise CubeValidationException from e
        if is_valid is False:
            pipeline_failure_message = gen_wmg_pipeline_failure_message(
                "Issue with cube validation, see logs for more detail"
            )
            data = format_failed_batch_issue_slack_alert(pipeline_failure_message)
            notify_slack(data)
            raise CubeValidationException
