from chalice import Chalice
import os
import sys

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "chalicelib"))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from upload_failures.upload import delete_many_from_s3, update_dataset_processing_status_to_failed


app = Chalice(app_name="upload_failures")


@app.lambda_function()
def handle_failure(event, context):
    dataset_uuid = event["dataset_uuid"]
    object_key = os.path.join(os.environ.get("REMOTE_DEV_PREFIX", ""), dataset_uuid).strip("/")
    delete_many_from_s3(os.environ["ARTIFACT_BUCKET"], object_key)
    delete_many_from_s3(f"hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}", object_key)
    update_dataset_processing_status_to_failed(dataset_uuid, event["error"]["Cause"])
