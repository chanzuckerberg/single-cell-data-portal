import logging

from backend.cellguide.pipeline import run_cellguide_pipeline
from backend.common.utils.result_notification import (
    format_failed_batch_issue_slack_alert,
    gen_cg_pipeline_failure_message,
    gen_cg_pipeline_success_message,
    notify_slack,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    try:
        output_path, description_output_path = run_cellguide_pipeline()
        success_message = gen_cg_pipeline_success_message(output_path, description_output_path)
        notify_slack(success_message)
    except Exception as e:
        logger.exception("Cell Guide Pipeline failed")
        failure_message = format_failed_batch_issue_slack_alert(
            gen_cg_pipeline_failure_message(f"Issue with Cell Guide pipeline run: {e}. See logs for more detail.")
        )
        notify_slack(failure_message)
