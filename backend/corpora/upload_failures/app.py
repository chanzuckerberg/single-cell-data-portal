import os
import sys

from backend.corpora.common.utils.aws import delete_many_from_s3
from backend.corpora.upload_failures.upload import update_dataset_processing_status_to_failed

def handle_failure(event, context):
    dataset_uuid = event["dataset_uuid"]
    object_key = os.path.join(os.environ.get("REMOTE_DEV_PREFIX", ""), dataset_uuid).strip("/")
    delete_many_from_s3(os.environ["ARTIFACT_BUCKET"], object_key)
    cellxgene_bucket = f"hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}"
    if os.getenv("CELLXGENE_BUCKET"):
        cellxgene_bucket = os.getenv("CELLXGENE_BUCKET")
    delete_many_from_s3(cellxgene_bucket, object_key)
    update_dataset_processing_status_to_failed(dataset_uuid, event["error"]["Cause"])
