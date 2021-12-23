import os

from backend.corpora.common.utils.aws import delete_many_from_s3
from backend.corpora.upload_failures.upload import update_dataset_processing_status_to_failed

from backend.corpora.dataset_processing.slack import notify_slack_failure


def handle_failure(event, context):
    dataset_uuid = event["dataset_uuid"]
    object_key = os.path.join(os.environ.get("REMOTE_DEV_PREFIX", ""), dataset_uuid).strip("/")
    delete_many_from_s3(os.environ["ARTIFACT_BUCKET"], object_key)

    deployment_stage = os.environ["DEPLOYMENT_STAGE"]

    cellxgene_bucket = f"hosted-cellxgene-{deployment_stage}"
    if os.getenv("CELLXGENE_BUCKET"):
        cellxgene_bucket = os.getenv("CELLXGENE_BUCKET")
    delete_many_from_s3(cellxgene_bucket, object_key)
    update_dataset_processing_status_to_failed(dataset_uuid, event["error"]["Cause"])

    notify_slack_failure(dataset_uuid)
