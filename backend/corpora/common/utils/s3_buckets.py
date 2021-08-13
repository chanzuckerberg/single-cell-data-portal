import os

import boto3

# s3_resource = None
# s3_client = None
# cxg_bucket = None

s3_resource = boto3.resource(
    "s3", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"), config=boto3.session.Config(signature_version="s3v4")
)
s3_client = boto3.client(
    "s3",
    endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None,
    config=boto3.session.Config(signature_version="s3v4"),
)
cxg_bucket = s3_resource.Bucket(os.getenv("CELLXGENE_BUCKET", f"hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}"))
