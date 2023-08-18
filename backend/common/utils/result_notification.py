import json
import logging
import os
from typing import Optional

import requests

from backend.common.corpora_config import CorporaConfig

logger = logging.getLogger(__name__)


aws_batch_job_url_fmt_str = (
    "https://{aws_region}.console.aws.amazon.com/batch/v2/home?region={aws_region}#jobs/detail/{job_id}"
)
aws_sfn_url_fmt_str = (
    "https://{aws_region}.console.aws.amazon.com/states/home?region={"
    "aws_region}#/v2/executions/details/{execution_arn}"
)


def notify_slack(data: dict, slack_webhook: Optional[str] = None):
    """
    This will only create a slack notification if called in the production env
    In all other envs (and in prod) it will simply log alert data
    """
    msg = json.dumps(data, indent=2)
    logger.info(f"Slack notification function called with message: {msg}")
    if os.getenv("DEPLOYMENT_STAGE") == "prod":
        if slack_webhook is None:
            slack_webhook = CorporaConfig().slack_webhook
        requests.post(slack_webhook, headers={"Content-type": "application/json"}, data=msg)


def upload_to_slack(file_name: str, content: str, initial_comment: str, channel: str, token: str) -> None:
    """
    Uploads a file to slack
    """
    slack_upload_url = "https://slack.com/api/files.upload"

    payload = {"channels": channel, "content": content, "filename": file_name, "initial_comment": initial_comment}
    headers = {"Authorization": f"Bearer {token}"}
    response = requests.post(slack_upload_url, data=payload, headers=headers)
    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError:
        logger.exception("Error uploading file to slack")


def format_failed_batch_issue_slack_alert(data: dict) -> dict:
    aws_region = os.getenv("AWS_DEFAULT_REGION")
    job_id = os.getenv("AWS_BATCH_JOB_ID")
    job_url = aws_batch_job_url_fmt_str.format(aws_region=aws_region, job_id=job_id)
    batch_data = {
        "type": "section",
        "text": {
            "type": "mrkdwn",
            "text": f"Batch processing job failed! @sc-oncall-eng\n" f"*Batch Job ID*:<{job_url}|{job_id}>\n",
        },
    }
    data["blocks"].append(batch_data)

    return data


def gen_wmg_pipeline_failure_message(failure_info: str) -> dict:
    return {
        "blocks": [
            {
                "type": "header",
                "text": {
                    "type": "plain_text",
                    "text": "WMG Snapshot Generation Pipeline FAILED:fire:",
                    "emoji": True,
                },
            },
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"WMG Snapshot Generation Pipeline failure @sc-oncall-eng \n{failure_info}",
                },
            },
        ]
    }


def gen_wmg_pipeline_success_message(snapshot_path: str, dataset_count: int, cell_count: int, gene_count: int) -> dict:
    return {
        "blocks": [
            {
                "type": "header",
                "text": {
                    "type": "plain_text",
                    "text": "WMG Snapshot Generation Pipeline Succeeded:tada: ",
                    "emoji": True,
                },
            },
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"\n* WMG snapshot stored in {snapshot_path}"
                    f"\n* The cube contains {cell_count} cells from {dataset_count} "
                    f"\n  datasets, with expression scores across {gene_count} genes.",
                },
            },
        ]
    }
