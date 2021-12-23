import json
import logging

import os

import requests

from backend.corpora.common.corpora_config import CorporaConfig
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.json import CustomJSONEncoder

logger = logging.getLogger(__name__)


def format_slack_message(dataset_id):
    # TODO: pass this info in so that the Slack message formatter doesn't have to know about the db
    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_id, include_tombstones=True)
        collection = dataset.collection
        collection_id, collection_owner = collection.id, collection.owner
        processing_status = dataset.processing_status.to_dict(remove_relationships=True)
    aws_region = os.getenv("AWS_DEFAULT_REGION")
    job_id = os.getenv("AWS_BATCH_JOB_ID")
    job_url = f"https://{aws_region}.console.aws.amazon.com/batch/v2/home?region={aws_region}#jobs/detail/{job_id}"
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
                    "text": f"Dataset processing job failed!\n"
                    f"*Batch Job ID*:<{job_url}|{job_id}>\n"
                    f"*Owner*: {collection_owner}\n"
                    f"*Collection*: https://cellxgene.cziscience.com/collections/{collection_id}/private\n"
                    f"*Processing Status*:\n",
                },
            },
            {
                "type": "section",
                "text": {
                    "type": "plain_text",
                    "text": json.dumps(processing_status, cls=CustomJSONEncoder, indent=2, sort_keys=True),
                    "emoji": False,
                },
            },
        ]
    }
    return json.dumps(data, indent=2)


def notify_slack_failure(dataset_id: str) -> None:
    data = format_slack_message(dataset_id)
    logger.info(data)

    # Slack notifications in non-prod environments (rdev, dev, staging) can be disabled by removing the
    # slack_webhook URL from corpora/backend/{env}/config in AWS Secrets Manager
    slack_webhook = CorporaConfig().slack_webhook
    if slack_webhook:
        response = requests.post(slack_webhook, headers={"Content-type": "application/json"}, data=data)

        if not response:
            logger.error(f"Failed to notify ingestion failure via Slack webhook! Response: {response.text}")
