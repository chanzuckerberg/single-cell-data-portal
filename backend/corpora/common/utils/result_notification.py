import json
import logging
import os

import requests

from backend.corpora.common.corpora_config import CorporaConfig
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.json import CustomJSONEncoder

logger = logging.getLogger(__name__)


def dataset_processing_slack_notification(dataset_id):
    data = format_dataset_processing_failure_slack_message(dataset_id)
    notify_slack(data)


def notify_slack(data: dict):
    """
    This will only create a slack notification if called in the production env
    In all other envs (and in prod) it will simply log alert data
    """
    msg = json.dumps(data, indent=2)
    logger.info(f"Slack notification function called with message: {msg}")
    if os.getenv("DEPLOYMENT_STAGE") == "prod":
        slack_webhook = CorporaConfig().slack_webhook
        requests.post(slack_webhook, headers={"Content-type": "application/json"}, data=msg)


def format_dataset_processing_failure_slack_message(dataset_id):
    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_id, include_tombstones=True)
        collection = dataset.collection
        collection_id, collection_owner = collection.id, collection.owner
        processing_status = dataset.processing_status.to_dict(remove_relationships=True)

    data = {
        "blocks": [
            {
                "type": "header",
                "text": {
                    "type": "plain_text",
                    "text": "Dataset failed to process:fire:",
                    "emoji": True,
                },
            },
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"Dataset processing job failed! @sc-oncall-eng\n"
                    f"*Owner*: {collection_owner}\n"
                    f"*Collection*: https://cellxgene.cziscience.com/collections/{collection_id}\n"
                    f"*Processing Status*:\n",
                },
            },
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"```{json.dumps(processing_status, cls=CustomJSONEncoder, indent=2, sort_keys=True)}```",
                },
            },
        ]
    }
    batch_alert_data = format_failed_batch_issue_slack_alert(data)
    logger.info(batch_alert_data)
    return batch_alert_data


def format_failed_batch_issue_slack_alert(data: dict) -> dict:
    aws_region = os.getenv("AWS_DEFAULT_REGION")
    job_id = os.getenv("AWS_BATCH_JOB_ID")
    job_url = f"https://{aws_region}.console.aws.amazon.com/batch/v2/home?region={aws_region}#jobs/detail/{job_id}"
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
                    "text": f"WMG Snapshot Generation Pipeline FAILED:fire: @sc-oncall-eng \n{failure_info}",
                    "emoji": True,
                },
            }
        ]
    }


def gen_wmg_pipeline_success_message(snapshot_path: str, dataset_count: int, cell_count: int, gene_count: int) -> dict:
    return {
        "blocks": [
            {
                "type": "header",
                "text": {
                    "type": "plain_text",
                    "text": f"WMG Snapshot Generation Pipeline Succeeded:tada: "
                    f"\n* WMG snapshot stored in {snapshot_path}"
                    f"\n* The cube contains {cell_count} cells from {dataset_count} "
                    f"\n  datasets, with expression scores across {gene_count} genes.",
                    "emoji": True,
                },
            }
        ]
    }
